###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import sys
import logging
import operator
from collections import defaultdict
import multiprocessing as mp

import biolib.seq_io as seq_io

from genometreetk.exceptions import GenomeTreeTkError
from genometreetk.common import (read_gtdb_taxonomy,
                                 read_gtdb_ncbi_taxonomy,
                                 read_gtdb_ncbi_organism_name,
                                 species_label)
from genometreetk.aai import mismatches


class Cluster(object):
    """Cluster genomes based on AAI of concatenated alignment.

    Genomes are preferentially assigned to representatives from
    RefSeq, followed by GenBank, and final User genomes. This
    is done to ensure the majority of genomes are assigned to
    publicly available representatives. However, this means
    a genome will not necessarily be assigned to the representative
    with the highest AAI. Furthermore, genomes are not clustered
    if they have different valid species names. NCBI genomes are
    never assigned to User representatives to ensure they always
    appear in trees without User genomes.
    """

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Maximum number of cpus/threads to use.
        """

        self.logger = logging.getLogger()

        self.cpus = cpus

        self.source_order = {'R': 0,  # RefSeq
                                'G': 1,  # GenBank
                                'U': 2}  # User

    def _reassign_representative(self,
                                    cur_representative_id,
                                    cur_max_mismatches,
                                    new_representative_id,
                                    new_max_mismatches,
                                    trusted_user_genomes):
        """Determines best representative genome.

        Genomes are preferentially assigned to representatives based on
        source repository (RefSeq => GenBank => Trusted User => User) and AAI.
        """

        if not cur_representative_id:
            # no currently assigned representative
            return new_representative_id, new_max_mismatches

        cur_source = self.source_order[cur_representative_id[0]]
        new_source = self.source_order[new_representative_id[0]]
        if new_representative_id in trusted_user_genomes:
            new_source = self.source_order['G']
        
        if new_source < cur_source:
            # give preference to genome source
            return new_representative_id, new_max_mismatches
        elif new_source == cur_source:
            # for the same source, find the representative with the highest AAI
            if new_max_mismatches < cur_max_mismatches:
                return new_representative_id, new_max_mismatches

        return cur_representative_id, cur_max_mismatches

    def __worker(self,
                 representatives,
                 bac_seqs,
                 ar_seqs,
                 aai_threshold,
                 metadata_file,
                 trusted_user_genomes,
                 queue_in,
                 queue_out):
        """Process genomes in parallel."""

        # determine genus of each genome and representatives belonging to a genus
        gtdb_taxonomy = read_gtdb_taxonomy(metadata_file)
        ncbi_taxonomy = read_gtdb_ncbi_taxonomy(metadata_file)
        ncbi_organism_names = read_gtdb_ncbi_organism_name(metadata_file)
        species = species_label(gtdb_taxonomy, ncbi_taxonomy, ncbi_organism_names)

        # determine genus of all genomes and representatives
        genus = {}
        reps_from_genus = defaultdict(set)
        for genome_id, t in gtdb_taxonomy.iteritems():
            g = None
            if len(t) >= 6 and t[5] != 'g__':
                g = t[5]
            elif genome_id in species:
                g = species[genome_id].split(' ')[0][3:]

            if g:
                genus[genome_id] = g

                if genome_id in representatives:
                    reps_from_genus[g].add(genome_id)

        while True:
            genome_id = queue_in.get(block=True, timeout=None)
            if genome_id == None:
                break

            genome_genus = genus.get(genome_id, None)
            genome_species = species.get(genome_id, None)
            genome_bac_seq = bac_seqs[genome_id]
            genome_ar_seq = ar_seqs[genome_id]

            cur_aai_threshold = aai_threshold
            assigned_representative = None

            bac_max_mismatches = (1.0 - aai_threshold) * (len(genome_bac_seq) - genome_bac_seq.count('-'))
            ar_max_mismatches = (1.0 - aai_threshold) * (len(genome_ar_seq) - genome_ar_seq.count('-'))

            # speed up computation by first comparing genome
            # to representatives of the same genus
            cur_reps_from_genus = reps_from_genus.get(genome_genus, set())
            for rep_id in cur_reps_from_genus:
                # do not cluster genomes from different named species
                rep_species = species.get(rep_id, None)
                if rep_species and genome_species and rep_species != genome_species:
                    continue

                # do not cluster NCBI or trusted User genomes with User representatives
                if rep_id.startswith('U_') and (not genome_id.startswith('U__') or genome_id in trusted_user_genomes):
                    continue

                rep_bac_seq = bac_seqs[rep_id]
                rep_ar_seq = ar_seqs[rep_id]

                m = mismatches(rep_bac_seq, genome_bac_seq, bac_max_mismatches)
                if m is not None:  # necessary to distinguish None and 0
                    assigned_representative, bac_max_mismatches = self._reassign_representative(assigned_representative,
                                                                                                bac_max_mismatches,
                                                                                                rep_id,
                                                                                                m,
                                                                                                trusted_user_genomes)
                else:
                    m = mismatches(rep_ar_seq, genome_ar_seq, ar_max_mismatches)
                    if m is not None:  # necessary to distinguish None and 0
                        assigned_representative, ar_max_mismatches = self._reassign_representative(assigned_representative,
                                                                                                    ar_max_mismatches,
                                                                                                    rep_id,
                                                                                                    m,
                                                                                                    trusted_user_genomes)

            # compare genome to remaining representatives
            remaining_reps = representatives.difference(cur_reps_from_genus)
            for rep_id in remaining_reps:
                # do not cluster genomes from different named species
                rep_species = species.get(rep_id, None)
                if rep_species and genome_species and rep_species != genome_species:
                    continue

                # do not cluster genomes from different named groups
                skip = False
                for rep_taxon, taxon in zip(gtdb_taxonomy[rep_id], gtdb_taxonomy[genome_id]):
                    rep_taxon = rep_taxon[3:]
                    taxon = taxon[3:]
                    if rep_taxon and taxon and rep_taxon != taxon:
                        skip = True
                        break
                if skip:
                    continue

                # do not cluster NCBI or trusted User genomes with User representatives
                if rep_id.startswith('U_') and (not genome_id.startswith('U__') or genome_id in trusted_user_genomes):
                    continue

                rep_bac_seq = bac_seqs[rep_id]
                rep_ar_seq = ar_seqs[rep_id]

                m = mismatches(rep_bac_seq, genome_bac_seq, bac_max_mismatches)
                if m is not None:  # necessary to distinguish None and 0
                    assigned_representative, bac_max_mismatches = self._reassign_representative(assigned_representative,
                                                                                                bac_max_mismatches,
                                                                                                rep_id,
                                                                                                m,
                                                                                                trusted_user_genomes)
                else:
                    m = mismatches(rep_ar_seq, genome_ar_seq, ar_max_mismatches)
                    if m is not None:  # necessary to distinguish None and 0
                        assigned_representative, ar_max_mismatches = self._reassign_representative(assigned_representative,
                                                                                                    ar_max_mismatches,
                                                                                                    rep_id,
                                                                                                    m,
                                                                                                    trusted_user_genomes)

            queue_out.put((genome_id, assigned_representative))

    def __writer(self, representatives, num_genomes, output_file, writer_queue):
        """Process representative assignments from each process.

        Parameters
        ----------
        representatives : set
            Initial set of representative genomes.
        num_genomes : int
            Number of genomes being processed.
        output_file : str
            Output file specifying genome clustering.
        """

        # initialize clusters
        clusters = {}
        for rep_id in representatives:
            clusters[rep_id] = []

        # gather results for each genome
        processed_genomes = 0
        while True:
            genome_id, assigned_representative = writer_queue.get(block=True, timeout=None)
            if genome_id == None:
              break

            processed_genomes += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) genomes.' % (processed_genomes,
                                                                            num_genomes,
                                                                            float(processed_genomes) * 100 / num_genomes)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            if assigned_representative:
                clusters[assigned_representative].append(genome_id)

        sys.stdout.write('\n')

        # write out cluster
        fout = open(output_file, 'w')
        clustered_genomes = 0
        for c, cluster_rep in enumerate(sorted(clusters, key=lambda x: len(clusters[x]), reverse=True)):
            cluster_str = 'cluster_%d' % (c + 1)
            cluster = clusters[cluster_rep]
            clustered_genomes += len(cluster)
            fout.write('%s\t%s\t%d\t%s\n' % (cluster_rep, cluster_str, len(cluster) + 1, ','.join(cluster)))

        fout.close()

        self.logger.info('Assigned %d genomes to representatives.' % clustered_genomes)

    def _cluster(self,
                    representatives,
                    genomes_to_process,
                    aai_threshold,
                    ar_seqs,
                    bac_seqs,
                    metadata_file,
                    trusted_user_genomes,
                    output_file):
        """Assign genomes to representatives.


        Parameters
        ----------
        representatives : set
            Initial set of representative genomes.
        genomes_to_process : list
            Genomes to process for identification of new representatives.
        aai_threshold : float
              AAI threshold for assigning a genome to a representative.
        ar_seqs : d[genome_id] -> alignment
            Alignment of archaeal marker genes.
        bac_seqs : d[genome_id] -> alignment
            Alignment of bacterial marker genes.
        metadata_file : str
            Metadata, including CheckM estimates, for all genomes.
        trusted_user_genomes : set
            Trusted User genomes to treat as if they were in GenBank.
        output_file : str
            Output file specifying genome clustering.
        """

        # populate worker queue with data to process
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()

        for genome_id in genomes_to_process:
          worker_queue.put(genome_id)

        for _ in range(self.cpus):
          worker_queue.put(None)

        try:
          worker_proc = [mp.Process(target=self.__worker, args=(representatives,
                                                                    bac_seqs,
                                                                    ar_seqs,
                                                                    aai_threshold,
                                                                    metadata_file,
                                                                    trusted_user_genomes,
                                                                    worker_queue,
                                                                    writer_queue)) for _ in range(self.cpus)]
          write_proc = mp.Process(target=self.__writer, args=(representatives,
                                                              len(genomes_to_process),
                                                              output_file,
                                                              writer_queue))

          write_proc.start()

          for p in worker_proc:
              p.start()

          for p in worker_proc:
              p.join()

          writer_queue.put((None, None))
          write_proc.join()
        except:
          for p in worker_proc:
            p.terminate()

          write_proc.terminate()

    def run(self,
            representatives_file,
            trusted_user_file,
            ar_msa_file,
            bac_msa_file,
            aai_threshold,
            metadata_file,
            output_file):
        """Identify additional representatives based on AAI between aligned sequences.

        Parameters
        ----------
        representatives_file : str
            File listing genome identifiers as initial representatives.
        trusted_user_file : str
            File listing trusted User genomes that should be treated as if they are in GenBank.
        ar_msa_file : str
            Name of file containing canonical archaeal multiple sequence alignment.
        bac_msa_file : str
            Name of file containing canonical bacterial multiple sequence alignment.
        aai_threshold : float
              AAI threshold for clustering genomes to a representative.
        metadata_file : str
            Metadata, including CheckM estimates, for all genomes.
        output_file : str
            Output file specifying genome clustering.
        """

        # read sequences
        ar_seqs = seq_io.read_fasta(ar_msa_file)
        bac_seqs = seq_io.read_fasta(bac_msa_file)
        self.logger.info('Identified %d archaeal sequences in MSA.' % len(ar_seqs))
        self.logger.info('Identified %d bacterial sequences in MSA.' % len(bac_seqs))

        if len(ar_seqs) != len(bac_seqs):
            self.logger.error('Archaeal and bacterial MSA files do not contain the same number of sequences.')
            raise GenomeTreeTkError('Error with MSA input files.')

        genome_to_consider = set(ar_seqs.keys())

        # read initial representatives
        rep_genomes = set()
        for line in open(representatives_file):
            if line[0] == '#':
                continue

            genome_id = line.rstrip().split('\t')[0]
            if genome_id not in genome_to_consider:
                self.logger.error('Representative genome %s has no sequence data.' % genome_id)
            rep_genomes.add(genome_id)

        self.logger.info('Identified %d representatives.' % len(rep_genomes))
        
        # read trusted User genomes
        trusted_user_genomes = set()
        if trusted_user_file:
            for line in open(trusted_user_file):
                line_split = line.split('\t')
                trusted_user_genomes.add(line_split[0])
                
        self.logger.info('Identified %d trusted User genomes.' % len(trusted_user_genomes)) 

        # cluster genomes to representatives
        genomes_to_cluster = genome_to_consider - rep_genomes
        self.logger.info('Comparing %d genomes to %d representatives with threshold = %.3f.' % (len(genomes_to_cluster),
                                                                                                len(rep_genomes),
                                                                                                aai_threshold))
        self._cluster(rep_genomes,
                                    genomes_to_cluster,
                                    aai_threshold,
                                    ar_seqs,
                                    bac_seqs,
                                    metadata_file,
                                    trusted_user_genomes,
                                    output_file)
