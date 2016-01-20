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
import pickle
import logging
import operator
import itertools
from collections import defaultdict

import biolib.seq_io as seq_io
from biolib.parallel import Parallel

from genometreetk.exceptions import GenomeTreeTkError
from genometreetk.common import (read_gtdb_genome_quality,
                                    read_gtdb_ncbi_taxonomy,
                                    read_gtdb_phylum)


class Cluster(object):
    """Cluster genomes based on AAI of concatenated alignment.

    Clustering is performed in a greedy manner similar to uClust,
    and the input order of sequences is important. An initial set
    of representatives can be specified and clustering if first
    performed against these. Additional clusters are then formed
    from the remaining genomes.

    To form a cluster, the first sequence on the list of remaining genomes
    is taken and all genomes within a define AAI threshold are assigned to
    the cluster. This process is repeated until all genomes have been assigned
    to a cluster.

    It is useful to retain two genomes per cluster so the cluster can be named
    in ARB. A separate file is produced that keeps both the representative genome
    for a cluster along with another genome, preferable the one with the best
    genome quality.
    """

    def __init__(self, assembly_metadata_file, taxonomy_file, type_strain_file, cpus):
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

        self.cpus = cpus

    def _aai(self, seq1, seq2, threshold):
        """Calculate AAI between sequences.

        Parameters
        ----------
        seq1 : str
            First sequence.
        seq2 : float
            Second sequence.
        threshold : float
            Required AAI to cluster sequences.

        Returns
        -------
        bool
            True if AAI is greater than or equal to threshold, else False.
        """

        max_mismatches = (1.0 - threshold) * len(seq1)

        mismatches = 0
        matches = 0
        for c1, c2 in itertools.izip(seq1, seq2):
            if c1 == '-' or c2 == '-':
                continue
            elif c1 != c2:
                mismatches += 1
                if mismatches >= max_mismatches:
                    return False
            else:
                matches += 1

        aai = float(matches) / max(1, (matches + mismatches))

        return aai >= threshold

    def _order_genomes(self, genome_set, rep_genomes, genome_quality, ncbi_species):
        """Determine order to process genomes based on existing representatives and genome quality.

        Parameters
        ----------
        genome_set : set
            Unique identifier of genomes to order.
        rep_genomes : set
            Initial set of representative genomes.
        genome_quality : d[genome_id] -> genome quality
            Genome quality.
        ncbi_species : d[genome_id] -> species
            Species of NCBI genomes.

        Returns
        -------
        list
            Genomes placed in order for clustering.
        """

        self.logger.info('Sorting genomes into suitable order for clustering.')

        # determine most common species
        species_count = defaultdict(int)
        for species in ncbi_species.values():
            species_count[species] += 1

        # prioritize clustering on representatives from abundant species
        sorted_rep_genomes = []
        sorted_species_count = sorted(species_count.items(), key=operator.itemgetter(1), reverse=True)
        for species, count in sorted_species_count:
            if count < 100:
                break

            for genome_id in rep_genomes:
                if ncbi_species.get(genome_id, None) == species:
                    sorted_rep_genomes.append(genome_id)

        # sort genomes by initial representatives and then genome quality
        sorted_genomes_refeq = []
        sorted_genomes_genbank = []
        sorted_genomes_user = []
        sorted_by_quality = sorted(genome_quality.items(), key=operator.itemgetter(1), reverse=True)
        for genome_id, _quality in sorted_by_quality:
            if genome_id not in genome_set:
                continue

            if genome_id in rep_genomes:
                if genome_id not in sorted_rep_genomes:  # make sure genome hasn't already been added
                    sorted_rep_genomes.append(genome_id)
            else:
                if genome_id.startswith('RS_'):
                    sorted_genomes_refeq.append(genome_id)
                elif genome_id.startswith('GB_'):
                    sorted_genomes_genbank.append(genome_id)
                elif genome_id.startswith('U_'):
                    sorted_genomes_user.append(genome_id)
                else:
                    self.logger.error('Unrecognized genome prefix: %s' % genome_id)
                    sys.exit(-1)

        return sorted_rep_genomes + sorted_genomes_refeq + sorted_genomes_genbank + sorted_genomes_user

    def _producer(self, data):
        """Producer thread for parallel processing of AAIs."""

        index, genomeI, genomeJ = data

        ar_seqI = self.ar_seqs[genomeI]
        bac_seqI = self.bac_seqs[genomeI]

        ar_seqJ = self.ar_seqs[genomeJ]
        bac_seqJ = self.bac_seqs[genomeJ]

        bCluster = False
        if self._aai(bac_seqI, bac_seqJ, self.threshold) or self._aai(ar_seqI, ar_seqJ, self.threshold):
            bCluster = True

        return (index, bCluster, genomeJ)

    def _consumer(self, produced_data, consumer_data):
        """Consumer thread for parallel processing.

        Returns
        -------
        dict : d[genome_id] -> bool
            Indicates if a genome should be clustered.
        """

        if consumer_data == None:
            # setup structure for consumed data
            consumer_data = {}

        index, bCluster, genomeJ = produced_data
        consumer_data[index] = (genomeJ, bCluster)

        return consumer_data

    def _greedy_clustering(self,
                               genomes_to_process,
                               threshold,
                               min_rep_quality,
                               ar_seqs,
                               bac_seqs,
                               genome_quality,
                               genome_taxon,
                               strict_taxon_matching,
                               clusters):
        """Cluster genomes.

        Only genomes from the same taxon as defined by
        genome_taxon are compared. This drastically reduces
        processing time. If the strict_taxon_matching
        flag is set to True, only genomes from the same
        defined taxon are compared. When this flag is False,
        genomes are compared when they are from the same taxon
        or where one or both genomes have an undefined taxon.

        Parameters
        ----------
        genomes_to_process : list
            Genomes to cluster placed in desired clustering order.
        threshold : float
              AAI threshold for forming clusters.
        min_rep_quality : float
            Minimum genome quality for a genome to be a representative.
        ar_seqs : d[genome_id] -> alignment
            Alignment of archaeal marker genes.
        bac_seqs : d[genome_id] -> alignment
            Alignment of bacterial marker genes.
        genome_quality : d[genome_id] - genome quality
            Quality of each genome.
        genome_taxon : d[genome_id] -> taxon
            Taxon for genomes used to restrict comparisons.
        strict_taxon_matching : bool
            Determines how genomes without an assigned taxon are processed.
        clusters: d[genome_id] -> genomes in cluster
            Clustering of genomes to representatives.
        """

        remaining_genomes_to_cluster = list(genomes_to_process)

        self.threshold = threshold
        self.ar_seqs = ar_seqs
        self.bac_seqs = bac_seqs

        # perform greedy clustering of genomes
        parallel = Parallel(cpus=self.cpus)
        while len(genomes_to_process):
            genomeI = genomes_to_process.pop(0)

            if genome_quality.get(genomeI, 0) < min_rep_quality:
                # remaining genomes do not have sufficient quality
                # to become representatives
                break

            taxonI = genome_taxon.get(genomeI, None)
            if strict_taxon_matching and not taxonI:
                continue

            # determine genomes to compare with current genome
            genomes_to_compare = []
            results = {}
            for index, genomeJ in enumerate(genomes_to_process):
                # check taxon of both genomes
                taxonJ = genome_taxon.get(genomeJ, None)
                if strict_taxon_matching:
                    if not taxonJ or taxonI != taxonJ:
                        results[index] = (genomeJ, False)
                        continue
                else:
                    if taxonI and taxonJ and taxonI != taxonJ:
                        results[index] = (genomeJ, False)
                        continue

                if genomeJ in clusters:
                    # do not cluster representatives together
                    results[index] = (genomeJ, False)
                    continue

                genomes_to_compare.append((index, genomeI, genomeJ))

            # compare genomes in parallel
            if len(genomes_to_compare) > 10 * self.cpus:
                aai_results = parallel.run(self._producer, self._consumer, data_items=genomes_to_compare)
            else:
                aai_results = {}

                ar_seqI = self.ar_seqs[genomeI]
                bac_seqI = self.bac_seqs[genomeI]
                for index, genomeI, genomeJ in genomes_to_compare:
                    ar_seqJ = self.ar_seqs[genomeJ]
                    bac_seqJ = self.bac_seqs[genomeJ]

                    bCluster = self._aai(bac_seqI, bac_seqJ, self.threshold) or self._aai(ar_seqI, ar_seqJ, self.threshold)
                    aai_results[index] = (genomeJ, bCluster)

            if aai_results:
                results.update(aai_results)

            cur_cluster = []
            remaining_genomes_to_process = []
            for index in sorted(results.keys()):
                genomeJ, bCluster = results[index]
                if bCluster:
                    cur_cluster.append(genomeJ)
                    remaining_genomes_to_cluster.remove(genomeJ)
                else:
                    remaining_genomes_to_process.append(genomeJ)

            clusters[genomeI].extend(cur_cluster)
            genomes_to_process = remaining_genomes_to_process

            # report progress of clustering (extra spaces are to ensure line in completely overwritten)
            #***self.logger.info('%s\t%d\t%s\t%d' % (genomeI, len(cur_cluster), genome_taxon.get(genomeI, None), len(genomes_to_compare)))
            statusStr = '==> Cluster %s contains %d genomes. There are %d of %d genomes remaining.       ' % (genomeI,
                                                                                                               len(cur_cluster) + 1,
                                                                                                               len(genomes_to_process),
                                                                                                               len(ar_seqs))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

        return remaining_genomes_to_cluster, clusters

    def run(self,
            ar_msa_file,
            bac_msa_file,
            representative_genomes,
            metadata_file,
            threshold,
            min_rep_quality,
            output_dir):
        """Perform clustering based on AAI between aligned sequences.

        Parameters
        ----------
        ar_msa_file : str
            Name of file containing canonical archaeal multiple sequence alignment.
        bac_msa_file : str
            Name of file containing canonical bacterial multiple sequence alignment.
        representative_genomes : str
            File listing genome identifiers for initial representatives.
        metadata_file : str
            Metadata, including CheckM estimates, for all genomes.
        threshold : float
              AAI threshold for forming clusters.
        min_rep_quality : float
            Minimum genome quality for a genome to be a representative.
        output_dir : str
            Output directory to store results.
        """

        # read sequences
        ar_seqs = seq_io.read_fasta(ar_msa_file)
        bac_seqs = seq_io.read_fasta(bac_msa_file)
        self.logger.info('Identified %d archaeal sequences in MSA.' % len(ar_seqs))
        self.logger.info('Identified %d bacterial sequences in MSA.' % len(bac_seqs))

        if len(ar_seqs) != len(bac_seqs):
            self.logger.error('Archaeal and bacterial MSA files do not contain the same number of sequences.')
            raise GenomeTreeTkError('Error with MSA input files.')

        genome_set = set(ar_seqs.keys())

        # read initial representatives
        rep_genomes = set()
        for line in open(representative_genomes):
            if line[0] == '#':
                continue

            genome_id = line.rstrip().split('\t')[0]
            if genome_id not in genome_set:
                self.logger.error('Representative genomes %s has no sequence data.' % genome_id)
            rep_genomes.add(genome_id)

        self.logger.info('Identified %d initial representatives.' % len(rep_genomes))

        # read genome quality
        genome_quality = read_gtdb_genome_quality(metadata_file, keep_db_prefix=True)

        missing_quality = genome_set - set(genome_quality.keys())
        if missing_quality:
            self.logger.warning('There are %d genomes with sequence data, but no metadata information.' % len(missing_quality))

        # read NCBI taxonomy
        ncbi_taxonomy = read_gtdb_ncbi_taxonomy(metadata_file, keep_db_prefix=True)
        ncbi_species = {}
        for genome_id, t in ncbi_taxonomy.iteritems():
            if len(t) == 7 and t[6] != 's__':
                ncbi_species[genome_id] = t[6]

        phyla = read_gtdb_phylum(metadata_file, keep_db_prefix=True)
        gtdb_phyla = {}
        for genome_id, phyla in phyla.iteritems():
            if phyla and phyla != 'p__':
                gtdb_phyla[genome_id] = phyla
            else:
                t = ncbi_taxonomy.get(genome_id, None)
                if t and len(t) >= 2:
                    gtdb_phyla[genome_id] = t[1]

        self.logger.info('Identified %d genomes with defined species, and %d genomes with defined phyla.' % (len(ncbi_species), len(gtdb_phyla)))

        # order genomes
        genomes_to_process = self._order_genomes(genome_set, rep_genomes, genome_quality, ncbi_species)

        # ensure specified representative genomes remain representatives
        clusters = defaultdict(list)
        for genome_id in rep_genomes:
            clusters[genome_id] = []

        # perform greedy clustering on genomes from same species
        self.logger.info('Performing species-specific greedy clustering of %d genomes with threshold = %.2f.' % (len(genomes_to_process), threshold))
        genomes_to_process, clusters = self._greedy_clustering(genomes_to_process,
                                                                       threshold,
                                                                       min_rep_quality,
                                                                       ar_seqs,
                                                                       bac_seqs,
                                                                       genome_quality,
                                                                       ncbi_species,
                                                                       True,
                                                                       clusters)

        with open('clusters.pkl', 'wb') as f:
            pickle.dump(clusters, f, pickle.HIGHEST_PROTOCOL)

        with open('genomes_to_process.pkl', 'wb') as f:
            pickle.dump(genomes_to_process, f, pickle.HIGHEST_PROTOCOL)

        with open('clusters.pkl', 'rb') as f:
            clusters = pickle.load(f)

        with open('genomes_to_process.pkl', 'rb') as f:
            genomes_to_process = pickle.load(f)

        sorted_clusters = sorted(clusters, key=lambda k: len(clusters[k]), reverse=True)
        for cluster_rep in sorted_clusters[0:10]:
            print cluster_rep, len(clusters[cluster_rep]), ncbi_species.get(cluster_rep, 'none')

        # perform greedy clustering on genomes from same species
        self.logger.info('Performing phylum-specific greedy clustering on %d genomes with threshold = %.2f.' % (len(genomes_to_process), threshold))
        genomes_to_process, clusters = self._greedy_clustering(genomes_to_process,
                                                                       threshold,
                                                                       min_rep_quality,
                                                                       ar_seqs,
                                                                       bac_seqs,
                                                                       genome_quality,
                                                                       gtdb_phyla,
                                                                       False,
                                                                       clusters)

        # write out representative genomes, genomes clustered to each representative, and
        # genomes to retain (which may, or may not, have a representative)
        self.logger.info('Writing out clusters.')
        cluster_file = os.path.join(output_dir, 'representatives.tsv')
        genome_to_retain_file = os.path.join(output_dir, 'genomes_to_retain.tsv')

        fout = open(cluster_file, 'w')
        fout_retain = open(genome_to_retain_file, 'w')
        for c, cluster_rep in enumerate(sorted(clusters, key=lambda x: len(clusters[x]), reverse=True)):
            cluster_str = 'cluster_%d' % (c + 1)
            cluster = clusters[cluster_rep]
            fout.write('%s\t%s\t%d\t%s\n' % (cluster_rep, cluster_str, len(cluster) + 1, ','.join(cluster)))

            fout_retain.write('%s\t%d\n' % (cluster_rep, c))

            reps = set(cluster_rep).intersection(rep_genomes)
            for r in reps:
                fout_retain.write('%s\t%d\n' % (r, c))

            if not len(reps) and len(cluster) >= 1:
                # two genomes should be retained for every cluster in order to aid annotating
                # in ARB. The first genome in cluster list  will be the highest quality genome.
                fout_retain.write('%s\t%d\n' % (cluster[0], c))

        # the remaining genomes to be processed are of insufficient quality to form representatives
        # but should be retained as they are novelty relative to the representatives
        for genome_id in set(genomes_to_process).difference(clusters.keys()):
            fout_retain.write('%s\t%d\n' % (genome_id, -1))

        fout.close()
        fout_retain.close()

        self.logger.info('Done.')
