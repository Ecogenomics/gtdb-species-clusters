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
import itertools

import biolib.seq_io as seq_io
from biolib.taxonomy import Taxonomy

import numpy as np

import genometreetk.ncbi as ncbi


class Cluster(object):
    """Cluster genomes based on AAI of concatenated alignment.

    Clustering is performed in a greedy manner similar to uClust,
    and the input order of sequences in important. Sequences are
    organized to process complete type strains, followed by
    partial type strains, complete genomes, and then all
    remaining genomes in a random order. Type strain and completeness
    information is taken from NCBI's type material and assembly level
    metadata.

    To form a cluster, the first sequence on the list of remaining genomes
    is taken and all genomes within a define AAI threshold and belonging
    to the same species as the genome are assigned to the same cluster.
    This process is repeated until all genomes have been assigned to a cluster.

    The constraint of belonging to the same species is used to ensure that similar
    genomes representing different species are not lost. In some (many?) cases
    this may represent an annotation error, but it is useful to identify when
    species annotation are incorrect.
    """

    def __init__(self, assembly_metadata_file, taxonomy_file, type_strain_file):
        """Initialization.

        Parameters
        ----------
        assembly_metadata_file : str
            File specifying assembly statistics of genomes.
        taxonomy_file : str
            File specifying 7 rank taxonomy of genomes.
        type_strain_file : str
            File specifying NCBI taxonomy ids representing type strains.
        """

        self.logger = logging.getLogger()

        self.assembly_metadata_file = assembly_metadata_file
        self.taxonomy_file = taxonomy_file
        self.type_strain_file = type_strain_file

    def _aai(self, seq1, seq2):
        """Calculate AAI between sequences.

        Parameters
        ----------
        seq1 : str
            First sequence.
        seq2 : float
            Second sequence.

        Returns
        -------
        float
            AAI between sequences.
        """

        mismatches = 0
        matches = 0
        for c1, c2 in itertools.izip(seq1, seq2):
            if c1 == '-' or c2 == '-':
                continue
            elif c1 != c2:
                mismatches += 1
            else:
                matches += 1

        return float(matches) / (matches + mismatches)

    def _order_genomes(self, genome_ids):
        """Determine order to process genomes based on type strain and completeness information.

        Parameters
        ----------
        genome_ids : set
            Unique identifier of genomes to order.

        Returns
        -------
        list
            Genomes placed in order for clustering.
        """

        self.logger.info('Sorting genomes into suitable order for clustering.')
        accession_to_taxid, complete_genomes = ncbi.read_metadata(self.assembly_metadata_file)
        type_strain_taxids = ncbi.read_type_strains(self.type_strain_file)
        type_strains = ncbi.get_type_strains(genome_ids, accession_to_taxid, type_strain_taxids)

        # limit genomes to those in the MSA file
        complete_genomes = complete_genomes.intersection(genome_ids)

        # identify genomes in different categories of interest
        complete_type_strains = complete_genomes.intersection(type_strains)
        remaining_type_strains = type_strains - complete_type_strains
        remaining_complete_genomes = complete_genomes - complete_type_strains
        remaining_genomes = genome_ids - complete_type_strains - remaining_type_strains - remaining_complete_genomes

        self.logger.info('Genomes to cluster = %d' % len(genome_ids))
        self.logger.info('Complete type strains = %d' % len(complete_type_strains))
        self.logger.info('Partial type strains = %d' % len(remaining_type_strains))
        self.logger.info('Remaining complete genomes= %d' % len(remaining_complete_genomes))
        self.logger.info('Remaining genomes = %d' % len(remaining_genomes))

        genomes_to_process = (list(complete_type_strains)
                                + list(remaining_type_strains)
                                + list(remaining_complete_genomes)
                                + list(remaining_genomes))

        return genomes_to_process

    def _species(self, genome_ids, taxonomy):
        """Determine species of each genome.

        Parameters
        ----------
        genome_ids : set
            Unique identifier of genomes of interest.
        taxonomy : d[genome_id] -> list of taxa

        Returns
        -------
        d[genome_id] -> species
            Species of each genome.
        """

        # determine genomes belonging to each named species
        self.logger.info('Determining species of each genome.')
        species = {}
        for genome_id, t in taxonomy.iteritems():
            if genome_id not in genome_ids:
                continue

            try:
                sp = t[Taxonomy.rank_index['s__']]
            except:
                self.logger.error('Genome does not have a complete taxonomy string: %s' % genome_id)
                sys.exit(0)

            if sp == 's__' or sp.lower() == 's__unclassified' or len(sp.split(' ')) < 2:
                # genomes is assigned to an unknown or incomplete species
                # so assign genome its own unique genome identifier to ensure
                # it is treated as its own species
                species[genome_id] = genome_id
            else:
                species[genome_id] = sp

        return species

    def run(self, msa_file, threshold, output_dir):
        """Perform clustering based on AAI between aligned sequences.

        Parameters
        ----------
        msa_file : str
            Name of file containing multiple sequence alignment in FASTA format.
        threshold : float
              AAI threshold for forming clusters.
        output_dir : str
            Output directory to store results.
        """

        # read sequences removing database prefix
        seqs = {}
        for seq_id, seq in seq_io.read_seq(msa_file):
            seqs[seq_id.replace('NCBI_', '')] = seq

        # determine species of each genome
        taxonomy = Taxonomy().read(self.taxonomy_file)
        species = self._species(set(seqs.keys()), taxonomy)

        # determine order to process genomes based on type strain and completeness information
        genomes_to_process = self._order_genomes(set(seqs.keys()))

        # perform greedy clustering of genomes
        self.logger.info('Performing greedy clustering on %d genomes with threshold = %.2f.' % (len(genomes_to_process), threshold))
        clusters = {}
        while len(genomes_to_process):
            genomeI = genomes_to_process.pop(0)
            seqI = seqs[genomeI]

            # get species if it exists, otherwise set to genome identifier
            # to ensure genome is processed as belonging to a unique species
            speciesI = species.get(genomeI, genomeI)

            cur_cluster = set()
            remaining_genome_to_process = []
            for genomeJ in genomes_to_process:
                if speciesI == species.get(genomeJ, genomeJ):
                    seqJ = seqs[genomeJ]

                    if self._aai(seqI, seqJ) >= threshold:
                        cur_cluster.add(genomeJ)
                    else:
                        remaining_genome_to_process.append(genomeJ)
                else:
                    remaining_genome_to_process.append(genomeJ)

            clusters[genomeI] = cur_cluster
            genomes_to_process = remaining_genome_to_process

            # report progress of clustering (extra spaces are to ensure line in completely overwritten)
            statusStr = '==> Cluster %d contains %d genomes. There are %d of %d genomes remaining.       ' % (len(clusters),
                                                                                                               len(cur_cluster) + 1,
                                                                                                               len(genomes_to_process),
                                                                                                               len(seqs))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

        self.logger.info('Writing out clusters.')
        cluster_file = os.path.join(output_dir, 'clusters.tsv')
        fout = open(cluster_file, 'w')
        for cluster_rep in sorted(clusters, key=lambda x: len(clusters[x]), reverse=True):
            cluster = clusters[cluster_rep]
            fout.write('%s\t%d\t%s\n' % (cluster_rep, len(cluster) + 1, ','.join(cluster)))
        fout.close()

        self.logger.info('Done.')
