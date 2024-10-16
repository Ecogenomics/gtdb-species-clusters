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
import random
import pickle
from itertools import combinations
from collections import defaultdict, namedtuple

from numpy import (mean as np_mean,
                   median as np_median,
                   std as np_std,
                   zeros as np_zeros,
                   argmin as np_argmin)

from gtdblib.util.shell.execute import check_dependencies
from gtdblib.util.bio.accession import canonical_gid

from gtdb_species_clusters.genome_utils import read_genome_path
from gtdb_species_clusters.taxon_utils import read_gtdb_taxonomy
from gtdb_species_clusters.type_genome_utils import GenomeRadius
from gtdb_species_clusters.skani import Skani
from gtdb_species_clusters import defaults as Defaults


class PMC_ClusterStats(object):
    """Calculate statistics for species cluster."""

    def __init__(self, af_sp, max_genomes, ani_cache_file, cpus, output_dir):
        """Initialization."""

        check_dependencies(['skani'])

        self.cpus = cpus
        self.output_dir = output_dir

        self.log = logging.getLogger('rich')

        self.af_sp = af_sp

        self.skani = Skani(ani_cache_file, cpus)

        # maximum number of randomly selected genomes to
        self.max_genomes_for_stats = max_genomes
        # consider when calculating pairwise statistics

        self.RepStats = namedtuple(
            'RepStats', 'min_ani mean_ani std_ani median_ani')
        self.PairwiseStats = namedtuple('PairwiseStats', ('min_ani',
                                                          'mean_ani',
                                                          'std_ani',
                                                          'median_ani',
                                                          'ani_to_medoid',
                                                          'mean_ani_to_medoid',
                                                          'mean_ani_to_rep',
                                                          'ani_below_sp_threshold'))

    def find_multiple_reps(self, clusters, cluster_radius):
        """Determine number of non-rep genomes within ANI radius of multiple rep genomes.

        This method assumes the ANI cache contains all relevant ANI calculations between
        representative and non-representative genomes. This is the case once the de novo
        clustering has been performed.
        """

        self.log.info(
            'Determine number of non-rep genomes within ANI radius of multiple rep genomes.')

        # get clustered genomes IDs
        clustered_gids = []
        for rid in clusters:
            clustered_gids += clusters[rid]

        self.log.info('Considering {:,} representatives and {:,} non-representative genomes.'.format(
            len(clusters),
            len(clustered_gids)))

        nonrep_rep_count = defaultdict(set)
        for idx, gid in enumerate(clustered_gids):
            cur_ani_cache = self.skani.ani_cache[gid]
            for rid in clusters:
                if rid not in cur_ani_cache:
                    continue

                ani, af = Skani.symmetric_ani(
                    self.skani.ani_cache, gid, rid)
                if af >= self.af_sp and ani >= cluster_radius[rid].ani:
                    nonrep_rep_count[gid].add((rid, ani))

            if (idx+1) % 100 == 0 or (idx+1) == len(clustered_gids):
                statusStr = '-> Processing %d of %d (%.2f%%) clusters genomes.'.ljust(86) % (
                    idx+1,
                    len(clustered_gids),
                    float((idx+1)*100)/len(clustered_gids))
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()

        sys.stdout.write('\n')

        return nonrep_rep_count

    def intragenus_pairwise_ani(self, clusters, species, genome_files, gtdb_taxonomy):
        """Determine pairwise intra-genus ANI between representative genomes."""

        self.log.info(
            'Calculating pairwise intra-genus ANI values between GTDB representatives.')

        # get genus for each representative
        genus = {}
        for rid, sp in species.items():
            genus[rid] = sp.split()[0].replace('s__', '')
            assert genus[rid] == gtdb_taxonomy[rid][5].replace('g__', '')

        # determine all intra-genus pairs between representative genomes
        self.log.info('Determining intra-genus genome pairs.')
        ani_pairs = []
        all_ani_pairs = []
        for qid in clusters:
            for rid in clusters:
                if qid == rid:
                    continue

                genusA = genus[qid]
                genusB = genus[rid]
                if genusA != genusB:
                    continue

                # skani ANI is symmetric and the AF is calculated in both
                # directions so only need to run a given pair once
                if qid < rid:
                    ani_pairs.append((qid, rid))
                else:
                    ani_pairs.append((rid, qid))

                all_ani_pairs.append((qid, rid))
                all_ani_pairs.append((rid, qid))

        self.log.info(
            'Identified {:,} intra-genus genome pairs.'.format(len(ani_pairs)))

        # calculate ANI between pairs
        self.log.info(
            'Calculating ANI between {:,} genome pairs:'.format(len(ani_pairs)))
        if True:  # ***DEBUGGING
            ani_af = self.skani.pairs(ani_pairs, genome_files, preset = Defaults.SKANI_PRESET,)
            pickle.dump(ani_af, open(os.path.join(
                self.output_dir, 'type_genomes_ani_af.pkl'), 'wb'))
        else:
            ani_af = pickle.load(
                open(os.path.join(self.output_dir, 'type_genomes_ani_af.pkl'), 'rb'))

        # find closest intra-genus pair for each rep
        fout = open(os.path.join(self.output_dir,
                                 'intra_genus_pairwise_ani.tsv'), 'w')
        fout.write(
            'Genus\tSpecies 1\tGenome ID 1\tSpecies 2\tGenome ID2\tANI\tAF\n')
        closest_intragenus_rep = {}
        for qid, rid in all_ani_pairs:
            genusA = genus[qid]
            genusB = genus[rid]
            assert genusA == genusB

            cur_ani = 0
            ani, af = ('n/a', 'n/a')
            if qid in ani_af and rid in ani_af[qid]:
                ani, af = Skani.symmetric_ani(ani_af, qid, rid)
                cur_ani = ani
                cur_af = af

            fout.write('{}\t{}\t{}\t{}\t{}\t{:.2f}\t{:.3f}\n'.format(
                genusA,
                species[qid],
                qid,
                species[rid],
                rid,
                ani,
                af))

            if cur_ani > closest_intragenus_rep.get(qid, ['dummy', 0, 0])[1]:
                closest_intragenus_rep[qid] = (
                    rid, 
                    cur_ani, 
                    cur_af)

        fout.close()

        # write out closest intra-genus species to each representative
        fout = open(os.path.join(self.output_dir,
                                 'closest_intragenus_rep.tsv'), 'w')
        fout.write(
            'Genome ID\tSpecies\tIntra-genus neighbour\tIntra-genus species\tANI\tAF\n')
        for qid in closest_intragenus_rep:
            rid, ani, af = closest_intragenus_rep[qid]

            fout.write('{}\t{}\t{}\t{}\t{:.2f}\t{:.3f}\n'.format(
                qid,
                species[qid],
                rid,
                species[rid],
                ani,
                af))
        fout.close()

    def parse_sp_clusters(self, cluster_file):
        """Parse species clustering information."""

        species = {}
        clusters = {}
        cluster_radius = {}
        with open(cluster_file) as f:
            headers = f.readline().strip().split('\t')

            rid_index = headers.index('Representative genome')
            sp_index = headers.index('GTDB species')
            num_clustered_index = headers.index('No. clustered genomes')
            clustered_genomes_index = headers.index('Clustered genomes')
            ani_radius_index = headers.index('ANI circumscription radius')

            for line in f:
                line_split = line.strip().split('\t')

                rid = line_split[rid_index]
                rid = canonical_gid(rid)

                species[rid] = line_split[sp_index]

                clusters[rid] = set()
                num_clustered = int(line_split[num_clustered_index])
                if num_clustered > 0:
                    for gid in [g.strip() for g in line_split[clustered_genomes_index].split(',')]:
                        gid = canonical_gid(gid)
                        clusters[rid].add(gid)

                cluster_radius[rid] = GenomeRadius(ani=float(line_split[ani_radius_index]),
                                                   af=self.af_sp,
                                                   neighbour_gid=None)

        return clusters, species, cluster_radius

    def rep_genome_stats(self, clusters, genome_files):
        """Calculate statistics relative to representative genome."""

        self.log.info('Calculating statistics to cluster representatives:')
        stats = {}
        for idx, (rid, cids) in enumerate(clusters.items()):
            if len(cids) == 0:
                stats[rid] = self.RepStats(min_ani=-1,
                                           mean_ani=-1,
                                           std_ani=-1,
                                           median_ani=-1)
            else:
                # calculate ANI to representative genome
                gid_pairs = []
                for cid in cids:
                    gid_pairs.append((cid, rid))
                    gid_pairs.append((rid, cid))

                if True:  # *** DEBUGGING
                    ani_af = self.skani.pairs(gid_pairs,
                                                genome_files,
                                                preset = Defaults.SKANI_PRESET,
                                                report_progress=False)
                else:
                    ani_af = self.skani.ani_cache

                # calculate statistics
                anis = [Skani.symmetric_ani(ani_af, cid, rid)[
                    0] for cid in cids]

                stats[rid] = self.RepStats(min_ani=min(anis),
                                           mean_ani=np_mean(anis),
                                           std_ani=np_std(anis),
                                           median_ani=np_median(anis))

            statusStr = '-> Processing %d of %d (%.2f%%) clusters.'.ljust(86) % (
                idx+1,
                len(clusters),
                float((idx+1)*100)/len(clusters))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

        return stats

    def pairwise_stats(self, clusters, genome_files):
        """Calculate statistics for all pairwise comparisons in a species cluster."""

        self.log.info(
            f'Restricting pairwise comparisons to {self.max_genomes_for_stats:,} randomly selected genomes.')
        self.log.info(
            'Calculating statistics for all pairwise comparisons in a species cluster:')

        stats = {}
        for idx, (rid, cids) in enumerate(clusters.items()):
            statusStr = '-> Processing {:,} of {:,} ({:2f}%) clusters (size = {:,}).'.ljust(86).format(
                idx+1,
                len(clusters),
                float((idx+1)*100)/len(clusters),
                len(cids))
            sys.stdout.write('{}\r'.format(statusStr))
            sys.stdout.flush()

            if len(cids) == 0:
                stats[rid] = self.PairwiseStats(min_ani=-1,
                                                mean_ani=-1,
                                                std_ani=-1,
                                                median_ani=-1,
                                                ani_to_medoid=-1,
                                                mean_ani_to_medoid=-1,
                                                mean_ani_to_rep=-1,
                                                ani_below_sp_threshold=-1)
            else:
                if len(cids) > self.max_genomes_for_stats:
                    cids = set(random.sample(cids, self.max_genomes_for_stats))

                # calculate ANI to representative genome
                gid_pairs = []
                gids = list(cids.union([rid]))
                for gid1, gid2 in combinations(gids, 2):
                    gid_pairs.append((gid1, gid2))
                    gid_pairs.append((gid2, gid1))

                if True:  # ***DEBUGGING
                    ani_af = self.skani.pairs(gid_pairs,
                                                genome_files,
                                                preset = Defaults.SKANI_PRESET,
                                                report_progress=False)
                else:
                    ani_af = self.skani.ani_cache

                # calculate medoid point
                if len(gids) > 2:
                    dist_mat = np_zeros((len(gids), len(gids)))
                    for i, gid1 in enumerate(gids):
                        for j, gid2 in enumerate(gids):
                            if i < j:
                                ani, _af = Skani.symmetric_ani(
                                    ani_af, gid1, gid2)
                                dist_mat[i, j] = 100 - ani
                                dist_mat[j, i] = 100 - ani

                    medoid_idx = np_argmin(dist_mat.sum(axis=0))
                    medoid_gid = gids[medoid_idx]
                else:
                    # with only 2 genomes in a cluster, the representative is the
                    # natural medoid at least for reporting statistics for the
                    # individual species cluster
                    medoid_gid = rid

                mean_ani_to_medoid = np_mean([Skani.symmetric_ani(ani_af, gid, medoid_gid)[0]
                                              for gid in gids if gid != medoid_gid])

                mean_ani_to_rep = np_mean([Skani.symmetric_ani(ani_af, gid, rid)[0]
                                           for gid in gids if gid != rid])

                if mean_ani_to_medoid < mean_ani_to_rep:
                    self.log.error('mean_ani_to_medoid < mean_ani_to_rep')
                    sys.exit(-1)

                # calculate statistics
                anis = []
                for gid1, gid2 in combinations(gids, 2):
                    ani, _af = Skani.symmetric_ani(ani_af, gid1, gid2)
                    anis.append(ani)

                stats[rid] = self.PairwiseStats(
                    min_ani=min(anis),
                    mean_ani=np_mean(anis),
                    std_ani=np_std(anis),
                    median_ani=np_median(anis),
                    ani_to_medoid=Skani.symmetric_ani(
                        ani_af, rid, medoid_gid)[0],
                    mean_ani_to_medoid=mean_ani_to_medoid,
                    mean_ani_to_rep=mean_ani_to_rep,
                    ani_below_sp_threshold=sum([1 for ani in anis if ani < Defaults.ANI_SP]))

        sys.stdout.write('\n')

        return stats

    def write_cluster_stats(self,
                            stats_file,
                            clusters,
                            species,
                            cluster_radius,
                            rep_stats,
                            pairwise_stats):
        """Write file with cluster statistics."""

        fout = open(stats_file, 'w')
        fout.write('Species\tRep genome\tNo. clustered genomes')
        fout.write(
            '\tMin ANI to rep\tMean ANI to rep\tStd ANI to rep\tMedian ANI to rep')
        fout.write(
            '\tMin pairwise ANI\tMean pairwise ANI\tStd pairwise ANI\tMedian pairwise ANI')
        fout.write(
            f'\tANI to medoid\tMean ANI to medoid\tMean ANI to rep (w/ subsampling)\tANI pairs <{Defaults.ANI_SP}%')
        fout.write(
            '\tClosest species\tClosest rep genome\tANI radius\tAF closest')
        fout.write('\tClustered genomes\n')

        for rid in clusters:
            fout.write('%s\t%s\t%d' % (species[rid], rid, len(clusters[rid])))
            fout.write('\t%.2f\t%.2f\t%.3f\t%.2f' % (
                rep_stats[rid].min_ani,
                rep_stats[rid].mean_ani,
                rep_stats[rid].std_ani,
                rep_stats[rid].median_ani))

            fout.write('\t%.2f\t%.2f\t%.3f\t%.2f\t%.2f\t%.2f\t%.2f\t%d' % (
                pairwise_stats[rid].min_ani,
                pairwise_stats[rid].mean_ani,
                pairwise_stats[rid].std_ani,
                pairwise_stats[rid].median_ani,
                pairwise_stats[rid].ani_to_medoid,
                pairwise_stats[rid].mean_ani_to_medoid,
                pairwise_stats[rid].mean_ani_to_rep,
                pairwise_stats[rid].ani_below_sp_threshold))

            if cluster_radius[rid].neighbour_gid != 'N/A':
                fout.write('\t%s\t%s\t%.2f\t%.2f' % (
                    species[cluster_radius[rid].neighbour_gid],
                    cluster_radius[rid].neighbour_gid,
                    cluster_radius[rid].ani,
                    cluster_radius[rid].af))
            else:
                fout.write('\t%s\t%s\t%.2f\t%.2f' %
                           ('N/A', 'N/A', Defaults.ANI_SP, 0))

            fout.write('\t%s\n' % ','.join(clusters[rid]))

        fout.close()

    def run(self, sp_cluster_file, genome_path_file, metadata_file):
        """Calculate statistics for species cluster."""

        # read the GTDB taxonomy
        self.log.info('Reading GTDB taxonomy from metadata file.')
        gtdb_taxonomy = read_gtdb_taxonomy(metadata_file)

        # get path to genome FASTA files
        self.log.info('Reading path to genome FASTA files.')
        genome_files = read_genome_path(genome_path_file)
        self.log.info(f'Read path for {len(genome_files):,} genomes.')

        # determine type genomes and genomes clustered to type genomes
        self.log.info('Reading species clusters.')
        clusters, species, cluster_radius = self.parse_sp_clusters(
            sp_cluster_file)
        self.log.info(f'Identified {len(clusters):,} species clusters.')

        # determine species assignment for clustered genomes
        clustered_species = {}
        for rid, cids in clusters.items():
            for cid in cids:
                clustered_species[cid] = species[rid]

        # find closest representative genome to each representative genome
        self.intragenus_pairwise_ani(
            clusters,
            species,
            genome_files,
            gtdb_taxonomy)

        sys.exit() #***

        # determine number of non-rep genomes within ANI radius of multiple rep genomes
        nonrep_rep_count = self.find_multiple_reps(clusters, cluster_radius)

        fout = open(os.path.join(self.output_dir,
                                 'nonrep_rep_ani_radius_count.tsv'), 'w')
        fout.write('Genome ID\tSpecies\tNo. rep radii\tMean radii')
        fout.write('\t<0.25%\t<0.5%\t<0.75%\t<1%\t<1.5%\t<2%')
        fout.write('\tRep genomes IDs\n')
        for gid, rid_info in nonrep_rep_count.items():
            rids = [rid for rid, ani in rid_info]
            anis = [ani for rid, ani in rid_info]

            fout.write('%s\t%s\t%d\t%.2f' % (
                gid,
                clustered_species[gid],
                len(rids),
                np_mean([cluster_radius[rid].ani for rid in rids])))

            if len(anis) >= 2:
                max_ani = max(anis)
                ani_2nd = sorted(anis, reverse=True)[1]
                diff = max_ani - ani_2nd
                fout.write('\t%s' % (diff < 0.25))
                fout.write('\t%s' % (diff < 0.5))
                fout.write('\t%s' % (diff < 0.75))
                fout.write('\t%s' % (diff < 1.0))
                fout.write('\t%s' % (diff < 1.5))
                fout.write('\t%s' % (diff < 2.0))
            else:
                fout.write('\tFalse\tFalse\tFalse\tFalse\tFalse\tFalse')

            fout.write('\t%s\n' % ','.join(rids))
        fout.close()

        # identify statistics relative to representative genome
        rep_stats = self.rep_genome_stats(clusters, genome_files)

        # identify pairwise statistics
        pairwise_stats = self.pairwise_stats(clusters, genome_files)

        # report statistics
        stats_file = os.path.join(self.output_dir, 'cluster_stats.tsv')
        self.write_cluster_stats(stats_file,
                                 clusters,
                                 species,
                                 cluster_radius,
                                 rep_stats,
                                 pairwise_stats)
