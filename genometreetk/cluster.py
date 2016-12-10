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
from biolib.taxonomy import Taxonomy

from genometreetk.exceptions import GenomeTreeTkError
from genometreetk.common import (read_gtdb_taxonomy,
                                 read_gtdb_ncbi_taxonomy,
                                 read_gtdb_ncbi_organism_name,
                                 species_label,
                                 predict_bacteria,
                                 check_domain_assignment,
                                 assign_rep)


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
        
        gtdb_taxonomy = read_gtdb_taxonomy(metadata_file)
        ncbi_taxonomy = read_gtdb_ncbi_taxonomy(metadata_file)

        # read taxonomy information and determine 'best' species label for each genome
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
                    
        # predict domain of each representative
        rep_is_bacteria = {}
        for rep_id in representatives:
            rep_is_bacteria[rep_id], per_bac_aa, per_ar_aa = predict_bacteria(rep_id, bac_seqs, ar_seqs)
            if not check_domain_assignment(rep_id, gtdb_taxonomy, ncbi_taxonomy, rep_is_bacteria[rep_id]):
                print 'Problem with representative.'
                print 'Bac vs Ar:', per_bac_aa, per_ar_aa

        # assign representatives           
        while True:
            genome_id = queue_in.get(block=True, timeout=None)
            if genome_id == None:
                break
                
            cur_aai = aai_threshold
            assigned_representative = None

            genome_is_bacteria, per_bac_aa, per_ar_aa = predict_bacteria(genome_id, bac_seqs, ar_seqs)
            if not check_domain_assignment(genome_id, gtdb_taxonomy, ncbi_taxonomy, genome_is_bacteria):
                print 'Bac vs Ar:', per_bac_aa, per_ar_aa
            else:
                if genome_is_bacteria:
                    genome_seq = bac_seqs[genome_id]
                else:
                    genome_seq = ar_seqs[genome_id]

                genome_aa_count = len(genome_seq) - genome_seq.count('-')
                min_matches = max(0.5*genome_aa_count, 0.1*len(genome_seq))

                # speed up computation by first comparing genome 
                # to representatives of the same genus
                genome_genus = genus.get(genome_id, None)
                cur_reps_from_genus = reps_from_genus.get(genome_genus, set())
                remaining_reps = representatives.difference(cur_reps_from_genus)
                for rep_set in [cur_reps_from_genus, remaining_reps]:
                    for rep_id in rep_set:
                        assigned_representative, cur_aai = assign_rep(rep_id, 
                                                                        genome_id,
                                                                        rep_is_bacteria[rep_id], 
                                                                        genome_is_bacteria,
                                                                        bac_seqs, 
                                                                        ar_seqs,
                                                                        species,
                                                                        gtdb_taxonomy,
                                                                        genome_aa_count,
                                                                        trusted_user_genomes,
                                                                        aai_threshold,
                                                                        min_matches,
                                                                        assigned_representative,
                                                                        cur_aai)

            queue_out.put((genome_id, assigned_representative, cur_aai))

    def __writer(self, representatives, metadata_file, num_genomes, output_file, writer_queue):
        """Process representative assignments from each process.

        Parameters
        ----------
        representatives : set
            Initial set of representative genomes.
        metadata_file : str
            GTDB metadata file.
        num_genomes : int
            Number of genomes being processed.
        output_file : str
            Output file specifying genome clustering.
        """
        
        ncbi_taxonomy = read_gtdb_ncbi_taxonomy(metadata_file)

        # initialize clusters
        clusters = {}
        for rep_id in representatives:
            clusters[rep_id] = []

        # gather results for each genome
        fout_details = open(output_file + '.details', 'w')
        fout_details.write('Representative Id\tGenome Id\tAAI\tRep. NCBI taxonomy\tGenome NCBI Taxonomy\n')
        processed_genomes = 0
        while True:
            genome_id, assigned_representative, aai = writer_queue.get(block=True, timeout=None)
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
                
                fout_details.write('%s\t%s\t%.2f\t%s\t%s\n' % (assigned_representative,
                                                                genome_id,
                                                                aai * 100,
                                                                ';'.join(ncbi_taxonomy.get(assigned_representative, Taxonomy.rank_prefixes)),
                                                                ';'.join(ncbi_taxonomy.get(genome_id, Taxonomy.rank_prefixes))))

        sys.stdout.write('\n')
        fout_details.close()

        # write out clusters
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
                                                                metadata_file,
                                                                  len(genomes_to_process),
                                                                  output_file,
                                                                  writer_queue))

          write_proc.start()

          for p in worker_proc:
              p.start()

          for p in worker_proc:
              p.join()

          writer_queue.put((None, None, None))
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
