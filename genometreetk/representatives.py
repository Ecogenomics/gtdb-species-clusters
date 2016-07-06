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

import biolib.seq_io as seq_io

from genometreetk.exceptions import GenomeTreeTkError
from genometreetk.common import (read_gtdb_genome_quality,
                                    read_gtdb_metadata,
                                    read_gtdb_taxonomy,
                                    read_gtdb_ncbi_taxonomy,
                                    read_gtdb_ncbi_organism_name,
                                    species_label)
from genometreetk.aai import aai_test, mismatches


class Representatives(object):
    """Identify representative genomes.

    Representative genomes are identified in a
    greedy manner using an amino acid identity (AAI)
    criteria.

    To ensure good representatives are selected, genomes
    are order before processing. Genomes are order first
    based on their source: RefSeq, GenBank, User. Within
    each source, genomes are order by genome quality
    (completeness - contamination). A threshold is used
    to limit representative to genomes of sufficient
    quality. Furthermore, a genome is not clustered
    to an existing representative if they have different
    valid species names. NCBI genomes are also never assigned
    to User representatives.
    """

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger()

    def _order_genomes(self, genomes_to_consider, genome_quality, trusted_user_genomes):
        """Order genomes by source and genome quality.

        Parameters
        ----------
        genomes_to_consider : list
            Genomes to order.
        genome_quality : d[genome_id] -> genome quality
            Estimate quality (completeness - contamination) of each genome.
        trusted_user_genomes : set
            Trusted User genomes to treat as if they were in GenBank.

        Returns
        -------
        list
            Genomes order by source and quality.
        """
        # sort genomes by source repository followed by genome quality
        sorted_refseq_rep_genomes = []
        sorted_genbank_rep_genomes = []
        sorted_trusted_user_rep_genomes = []
        sorted_user_rep_genomes = []
        sorted_by_quality = sorted(genome_quality.items(), key=operator.itemgetter(1), reverse=True)
        for genome_id, _quality in sorted_by_quality:
            if genome_id not in genomes_to_consider:
                continue

            if genome_id.startswith('RS_'):
                sorted_refseq_rep_genomes.append(genome_id)
            elif genome_id.startswith('GB_'):
                sorted_genbank_rep_genomes.append(genome_id)
            elif genome_id in trusted_user_genomes:
                sorted_trusted_user_rep_genomes.append(genome_id)
            elif genome_id.startswith('U_'):
                sorted_user_rep_genomes.append(genome_id)
            else:
                self.logger.error('Unrecognized genome prefix: %s' % genome_id)
                sys.exit(-1)

        return (sorted_refseq_rep_genomes
                    + sorted_genbank_rep_genomes
                    + sorted_trusted_user_rep_genomes
                    + sorted_user_rep_genomes)

    def _greedy_representatives(self,
                                    representatives,
                                    genomes_to_process,
                                    aai_threshold,
                                    ar_seqs,
                                    bac_seqs,
                                    metadata_file,
                                    trusted_user_genomes):
        """Identify additional representative genomes in a greedy fashion.

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
            Metadata for all genomes, including NCBI taxonomy information.
        trusted_user_genomes : set
            Trusted User genomes to treat as if they were in GenBank.

        Returns
        -------
        set
            Representative genomes.
        """

        # read taxonomy information and determine 'best' species label for each genome
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

        total_genomes = len(genomes_to_process)
        processed_genomes = 0
        while len(genomes_to_process):
            processed_genomes += 1
            sys.stdout.write('==> Processed %d of %d genomes.\r' % (processed_genomes, total_genomes))
            sys.stdout.flush()

            genome_id = genomes_to_process.pop(0)

            genome_genus = genus.get(genome_id, None)
            genome_species = species.get(genome_id, None)
            genome_bac_seq = bac_seqs[genome_id]
            genome_ar_seq = ar_seqs[genome_id]

            # speed up computation by comparing genome to representatives from the same genus
            cur_reps_from_genus = reps_from_genus.get(genome_genus, set())
            bCluster = False
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

                bCluster = (aai_test(rep_bac_seq, genome_bac_seq, aai_threshold)
                            or aai_test(rep_ar_seq, genome_ar_seq, aai_threshold))
                if bCluster:
                    break

            if not bCluster:
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

                    bCluster = (aai_test(rep_bac_seq, genome_bac_seq, aai_threshold)
                                or aai_test(rep_ar_seq, genome_ar_seq, aai_threshold))
                    if bCluster:
                        break

            if not bCluster:
                # genome was not assigned to an existing representative,
                # so make it a new representative genome
                representatives.add(genome_id)
                if genome_genus:
                    reps_from_genus[genome_genus].add(genome_id)

        sys.stdout.write('\n')

        return representatives

    def run(self,
            refseq_representatives,
            trusted_user_file,
            ar_msa_file,
            bac_msa_file,
            aai_threshold,
            min_rep_comp,
            max_rep_cont,
            min_quality,
            max_contigs,
            min_N50,
            max_ambiguous,
            metadata_file,
            output_file):
        """Identify additional representatives based on AAI between aligned sequences.

        Parameters
        ----------
        refseq_representatives : str
            File listing RefSeq genome identifiers as initial representatives.
        trusted_user_file : str
            File listing trusted User genomes that should be treated as if they are in GenBank.
        ar_msa_file : str
            Name of file containing canonical archaeal multiple sequence alignment.
        bac_msa_file : str
            Name of file containing canonical bacterial multiple sequence alignment.
        aai_threshold : float
              AAI threshold for clustering genomes to a representative.
        min_rep_comp : float [0, 100]
            Minimum completeness for a genome to be a representative.
        max_rep_cont : float [0, 100]
            Maximum contamination for a genome to be a representative.
        min_quality : float [0, 100]
            Minimum quality (comp - 4*cont) for a genome to be a representative.
        max_contigs : int
            Maximum number of contigs for a genome to be a representative.
        min_N50 : int
            Minimum N50 for a genome to be a representative.
        max_ambiguous : int
            Maximum number of ambiguous bases for a genome to be a representative.
        metadata_file : str
            Metadata, including CheckM estimates, for all genomes.
        output_file : str
            Output file containing all genomes identified as representatives.
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
        refseq_rep_genomes = set()
        for line in open(refseq_representatives):
            if line[0] == '#':
                continue

            genome_id = line.rstrip().split('\t')[0]
            if genome_id not in genome_to_consider:
                self.logger.error('Representative genome %s has no sequence data.' % genome_id)
                sys.exit(-1)
            refseq_rep_genomes.add(genome_id)

        self.logger.info('Identified %d initial representatives.' % len(refseq_rep_genomes))
        
        # read trusted User genomes
        trusted_user_genomes = set()
        if trusted_user_file:
            for line in open(trusted_user_file):
                line_split = line.split('\t')
                trusted_user_genomes.add(line_split[0])
                
        self.logger.info('Identified %d trusted User genomes.' % len(trusted_user_genomes))        

        # read genome quality
        genome_quality = read_gtdb_genome_quality(metadata_file)

        missing_quality = genome_to_consider - set(genome_quality.keys())
        if missing_quality:
            self.logger.error('There are %d genomes with sequence data, but no genome quality information.' % len(missing_quality))
            sys.exit(-1)

        # remove existing representative genomes and genomes
        # of insufficient quality to be a representative
        genome_stats = read_gtdb_metadata(metadata_file, ['contig_count',
                                                            'n50_contigs',
                                                            'ambiguous_bases'])
        potential_reps = set()
        for genome_id in (genome_to_consider - refseq_rep_genomes):
            comp, cont, qual = genome_quality.get(genome_id, [-1, -1, -1])
            stats = genome_stats[genome_id]
            if (comp >= min_rep_comp 
                and cont <= max_rep_cont 
                and (comp - 4*cont) >= min_quality
                and stats.contig_count <= max_contigs
                and stats.n50_contigs >= min_N50
                and stats.ambiguous_bases <= max_ambiguous):
                    potential_reps.add(genome_id)

        # perform greedy identification of new representatives
        ordered_genomes = self._order_genomes(potential_reps, genome_quality, trusted_user_genomes)
        self.logger.info('Comparing %d genomes to %d representatives with threshold = %.3f.' % (len(ordered_genomes),
                                                                                                len(refseq_rep_genomes),
                                                                                                aai_threshold))
        representatives = self._greedy_representatives(refseq_rep_genomes,
                                                            ordered_genomes,
                                                            aai_threshold,
                                                            ar_seqs,
                                                            bac_seqs,
                                                            metadata_file,
                                                            trusted_user_genomes)

        self.logger.info('Identified %d representatives.' % len(representatives))

        # write out representative genomes
        fout = open(output_file, 'w')
        for rep in representatives:
            fout.write(rep + '\n')

        fout.close()
