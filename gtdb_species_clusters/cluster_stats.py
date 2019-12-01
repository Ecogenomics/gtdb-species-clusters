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
import re
import random
import shutil
import tempfile
import ntpath
import pickle
from itertools import combinations
from collections import defaultdict, namedtuple, Counter

from biolib.external.execute import check_dependencies

from numpy import (mean as np_mean,
                    median as np_median,
                    std as np_std,
                    zeros as np_zeros,
                    argmin as np_argmin)

from gtdb_species_clusters.genome_utils import read_genome_path
from gtdb_species_clusters.taxon_utils import read_gtdb_ncbi_taxonomy
                                    
from gtdb_species_clusters.type_genome_utils import (GenomeRadius,
                                                        symmetric_ani)
                                            
from gtdb_species_clusters.mash import Mash
from gtdb_species_clusters.fastani import FastANI

class ClusterStats(object):
    """Calculate statistics for species cluster."""

    def __init__(self, af_sp, max_genomes, ani_cache_file, cpus, output_dir):
        """Initialization."""
        
        check_dependencies(['fastANI', 'mash'])
        
        self.cpus = cpus
        self.output_dir = output_dir

        self.logger = logging.getLogger('timestamp')
        
        self.af_sp = af_sp
        
        self.fastani = FastANI(ani_cache_file, cpus)
        
        self.max_genomes_for_stats = max_genomes    # maximum number of randomly selected genomes to
                                                    # consider when calculating pairwise statistics
        
        self.RepStats = namedtuple('RepStats', 'min_ani mean_ani std_ani median_ani')
        self.PairwiseStats = namedtuple('PairwiseStats', ('min_ani',
                                                           'mean_ani', 
                                                           'std_ani', 
                                                           'median_ani', 
                                                           'ani_to_medoid',
                                                           'mean_ani_to_medoid',
                                                           'mean_ani_to_rep',
                                                           'ani_below_95'))
                                                           
    def _find_multiple_reps(self, clusters, cluster_radius):
        """Determine number of non-rep genomes within ANI radius of multiple rep genomes.
        
        This method assumes the ANI cache contains all relevant ANI calculations between
        representative and non-representative genomes. This is the case once the de novo
        clustering has been performed.
        """
        
        self.logger.info('Determine number of non-rep genomes within ANI radius of multiple rep genomes.')
        
        # get clustered genomes IDs
        clustered_gids = []
        for rid in clusters:
            clustered_gids += clusters[rid]
            
        self.logger.info('Considering %d representatives and %d non-representative genomes.' % (len(clusters), len(clustered_gids)))
            
        nonrep_rep_count = defaultdict(set)
        for idx, gid in enumerate(clustered_gids):
            cur_ani_cache = self.fastani.ani_cache[gid]
            for rid in clusters:
                if rid not in cur_ani_cache:
                    continue
                    
                ani, af = symmetric_ani(self.fastani.ani_cache, gid, rid)
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

    def _find_closest_intragenus_rep(self, clusters, genome_files, ncbi_taxonomy):
        """Find closest intra-genus representative genomes."""
        
        self.logger.info('Identifying closest intra-genus neighbours for representative genomes.')
        
        min_mash_ani = 85.0
        
        mash = Mash(self.cpus)
        
        # create Mash sketch for potential representative genomes
        genome_list_file = os.path.join(self.output_dir, 'gtdb_type_genomes.lst')
        sketch = os.path.join(self.output_dir, 'gtdb_type_genomes.msh')
        mash.sketch(clusters.keys(), genome_files, genome_list_file, sketch)

        # get Mash distances
        mash_dist_file = os.path.join(self.output_dir, 'gtdb_type_genomes.dst')
        mash.dist_pairwise((100-min_mash_ani)/100, sketch, mash_dist_file)

        # read Mash distances
        mash_ani = mash.read_ani(mash_dist_file)

        # get pairs above Mash threshold
        mash_ani_pairs = []
        for qid in mash_ani:
            for rid in mash_ani[qid]:
                if qid == rid:
                    continue
                
                ncbi_genusA = ncbi_taxonomy[qid][5]
                ncbi_genusB = ncbi_taxonomy[rid][5]
                if ncbi_genusA != ncbi_genusB:
                    continue
            
                if mash_ani[qid][rid] >= min_mash_ani:
                    if qid != rid:
                        mash_ani_pairs.append((qid, rid))
                        mash_ani_pairs.append((rid, qid))
                
        self.logger.info('Identified %d intra-genus genome pairs with a Mash ANI >= %.1f%%.' % (len(mash_ani_pairs), min_mash_ani))
        
        # calculate ANI between pairs
        self.logger.info('Calculating ANI between %d genome pairs:' % len(mash_ani_pairs))
        if True: #***
            ani_af = self.fastani.pairs(mash_ani_pairs, genome_files)
            pickle.dump(ani_af, open(os.path.join(self.output_dir, 'type_genomes_ani_af.pkl'), 'wb'))
        else:
            ani_af = pickle.load(open(os.path.join(self.output_dir, 'type_genomes_ani_af.pkl'), 'rb'))
            
        # find closest intra-genus pair for each rep
        closest_intragenus_rep = {}
        for rid in clusters:
            if rid in ani_af:
                closest_ani = 0
                closest_gid = None
                for n_gid in ani_af[rid]:
                    ani, af = symmetric_ani(ani_af, rid, n_gid)
                    if ani > closest_ani:
                        closest_ani = ani
                        closest_gid = n_gid
                        
                if n_gid:
                    closest_intragenus_rep[rid] = (closest_gid, closest_ani)
            
        return closest_intragenus_rep
    
    def _parse_clusters(self, cluster_file):
        """Parse species clustering information."""
        
        species = {}
        clusters = {}
        cluster_radius = {}
        with open(cluster_file) as f:
            headers = f.readline().strip().split('\t')
            
            type_sp_index = headers.index('NCBI species')
            type_genome_index = headers.index('Type genome')
            num_clustered_index = headers.index('No. clustered genomes')
            clustered_genomes_index = headers.index('Clustered genomes')
            closest_type_index = headers.index('Closest type genome')
            ani_radius_index = headers.index('ANI radius')
            af_index = headers.index('AF closest')

            for line in f:
                line_split = line.strip().split('\t')
                
                rid = line_split[type_genome_index]
                species[rid] = line_split[type_sp_index]
                
                clusters[rid] = set()
                num_clustered = int(line_split[num_clustered_index])
                if num_clustered > 0:
                    for gid in [g.strip() for g in line_split[clustered_genomes_index].split(',')]:
                        clusters[rid].add(gid)
                        
                cluster_radius[rid] = GenomeRadius(ani = float(line_split[ani_radius_index]), 
                                                     af = float(line_split[af_index]),
                                                     neighbour_gid = line_split[closest_type_index])
                        
        return clusters, species, cluster_radius
        
    def _rep_genome_stats(self, clusters, genome_files):
        """Calculate statistics relative to representative genome."""
        
        self.logger.info('Calculating statistics to cluster representatives:')
        stats = {}
        for idx, (rid, cids) in enumerate(clusters.items()):
            if len(cids) == 0:
                stats[rid] = self.RepStats(min_ani = -1,
                                            mean_ani = -1,
                                            std_ani = -1,
                                            median_ani = -1)
            else:
                # calculate ANI to representative genome
                gid_pairs = []
                for cid in cids:
                    gid_pairs.append((cid, rid))
                    gid_pairs.append((rid, cid))
                
                if False: #***
                    ani_af = self.fastani.pairs(gid_pairs, 
                                                            genome_files, 
                                                            report_progress=False)
                else:
                    ani_af = self.fastani.ani_cache
                
                # calculate statistics
                anis = [symmetric_ani(ani_af, cid, rid)[0] for cid in cids]

                stats[rid] = self.RepStats(min_ani = min(anis),
                                            mean_ani = np_mean(anis),
                                            std_ani = np_std(anis),
                                            median_ani = np_median(anis))
                                            
            statusStr = '-> Processing %d of %d (%.2f%%) clusters.'.ljust(86) % (
                                idx+1, 
                                len(clusters), 
                                float((idx+1)*100)/len(clusters))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
                
        sys.stdout.write('\n')
            
        return stats
        
    def _pairwise_stats(self, clusters, genome_files, rep_stats):
        """Calculate statistics for all pairwise comparisons in a species cluster."""
        
        self.logger.info('Restricting pairwise comparisons to %d randomly selected genomes.' % self.max_genomes_for_stats)
        self.logger.info('Calculating statistics for all pairwise comparisons in a species cluster:')
        
        stats = {}
        for idx, (rid, cids) in enumerate(clusters.items()):
            statusStr = '-> Processing %d of %d (%.2f%%) clusters (size = %d).'.ljust(86) % (
                                idx+1, 
                                len(clusters), 
                                float((idx+1)*100)/len(clusters),
                                len(cids))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
                                
            if len(cids) == 0:
                stats[rid] = self.PairwiseStats(min_ani = -1,
                                                mean_ani = -1,
                                                std_ani = -1,
                                                median_ani = -1,
                                                ani_to_medoid = -1,
                                                mean_ani_to_medoid = -1,
                                                mean_ani_to_rep = -1,
                                                ani_below_95 = -1)
            else:
                if len(cids) > self.max_genomes_for_stats:
                    cids = set(random.sample(cids, self.max_genomes_for_stats))
                
                # calculate ANI to representative genome
                gid_pairs = []
                gids = list(cids.union([rid]))
                for gid1, gid2 in combinations(gids, 2):
                    gid_pairs.append((gid1, gid2))
                    gid_pairs.append((gid2, gid1))
                    
                if False: #***
                    ani_af = self.fastani.pairs(gid_pairs, 
                                                        genome_files, 
                                                        report_progress=False)
                else:
                    ani_af = self.fastani.ani_cache
                                                        
                # calculate medoid point
                if len(gids) > 2:
                    dist_mat = np_zeros((len(gids), len(gids)))
                    for i, gid1 in enumerate(gids):
                        for j, gid2 in enumerate(gids):
                            if i < j:
                                ani, af = symmetric_ani(ani_af, gid1, gid2)
                                dist_mat[i, j] = 100 - ani
                                dist_mat[j, i] = 100 - ani

                    medoid_idx = np_argmin(dist_mat.sum(axis=0))
                    medoid_gid = gids[medoid_idx]
                else:
                    # with only 2 genomes in a cluster, the representative is the
                    # natural medoid at least for reporting statistics for the
                    # individual species cluster
                    medoid_gid = rid
                    
                mean_ani_to_medoid = np_mean([symmetric_ani(ani_af, gid, medoid_gid)[0] 
                                                for gid in gids if gid != medoid_gid])
                                                
                mean_ani_to_rep = np_mean([symmetric_ani(ani_af, gid, rid)[0] 
                                                for gid in gids if gid != rid])
                                                
                if mean_ani_to_medoid < mean_ani_to_rep:
                    self.logger.error('mean_ani_to_medoid < mean_ani_to_rep')
                    sys.exit(-1)

                # calculate statistics
                anis = []
                for gid1, gid2 in combinations(gids, 2):
                    ani, af = symmetric_ani(ani_af, gid1, gid2)
                    anis.append(ani)
                    
                stats[rid] = self.PairwiseStats(min_ani = min(anis),
                                                mean_ani = np_mean(anis),
                                                std_ani = np_std(anis),
                                                median_ani = np_median(anis),
                                                ani_to_medoid = symmetric_ani(ani_af, rid, medoid_gid)[0],
                                                mean_ani_to_medoid = mean_ani_to_medoid,
                                                mean_ani_to_rep = mean_ani_to_rep,
                                                ani_below_95 = sum([1 for ani in anis if ani < 95]))

        sys.stdout.write('\n')
            
        return stats
        
    def _write_cluster_stats(self, 
                                stats_file,
                                clusters,
                                species,
                                cluster_radius,
                                rep_stats,
                                pairwise_stats):
        """Write file with cluster statistics."""
        
        fout = open(stats_file, 'w')
        fout.write('Species\tRep genome\tNo. clustered genomes')
        fout.write('\tMin ANI to rep\tMean ANI to rep\tStd ANI to rep\tMedian ANI to rep')
        fout.write('\tMin pairwise ANI\tMean pairwise ANI\tStd pairwise ANI\tMedian pairwise ANI')
        fout.write('\tANI to medoid\tMean ANI to medoid\tMean ANI to rep (w/ subsampling)\tANI pairs <95%')
        fout.write('\tClosest species\tClosest rep genome\tANI radius\tAF closest')
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
                        pairwise_stats[rid].ani_below_95))
                        
            if cluster_radius[rid].neighbour_gid != 'N/A':
                fout.write('\t%s\t%s\t%.2f\t%.2f' % (
                            species[cluster_radius[rid].neighbour_gid],
                            cluster_radius[rid].neighbour_gid,
                            cluster_radius[rid].ani,
                            cluster_radius[rid].af))
            else:
                fout.write('\t%s\t%s\t%.2f\t%.2f' % ('N/A', 'N/A', 95, 0))
                
            fout.write('\t%s\n' % ','.join(clusters[rid]))

        fout.close()
 
    def run(self, cluster_file, genome_path_file, metadata_file):
        """Calculate statistics for species cluster."""
        
        # read the NCBI taxonomy
        self.logger.info('Reading NCBI taxonomy from GTDB metadata file.')
        ncbi_taxonomy, ncbi_update_count = read_gtdb_ncbi_taxonomy(metadata_file, None)

        # get path to genome FASTA files
        self.logger.info('Reading path to genome FASTA files.')
        genome_files = read_genome_path(genome_path_file)
        self.logger.info('Read path for %d genomes.' % len(genome_files))
        
        # determine type genomes and genomes clustered to type genomes
        self.logger.info('Reading species clusters.')
        clusters, species, cluster_radius = self._parse_clusters(cluster_file)
        self.logger.info('Identified %d species clusters.' % len(clusters))
        
        # determine species assignment for clustered genomes
        clustered_species = {}
        for rid, cids in clusters.items():
            for cid in cids:
                clustered_species[cid] = species[rid]
        
        # determine number of non-rep genomes with ANI radius of multiple rep genomes
        if False:
            nonrep_rep_count = self._find_multiple_reps(clusters, cluster_radius)
            
            fout = open(os.path.join(self.output_dir, 'nonrep_rep_ani_radius_count.tsv'), 'w')
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

            sys.exit(0)
            
            # find closest representative genome to each representative genome
            closest_intragenus_rep = self._find_closest_intragenus_rep(clusters, genome_files, ncbi_taxonomy)
            
            fout = open(os.path.join(self.output_dir, 'closest_intragenus_rep.tsv'), 'w')
            fout.write('Genome ID\tSpecies\tIntra-genus neighbour\tIntra-genus species\tANI\n')
            for rid in closest_intragenus_rep:
                nid, ani = closest_intragenus_rep[rid]
                
                fout.write('%s\t%s\t%s\t%s\t%.2f\n' % (
                                rid,
                                species[rid],
                                nid,
                                species[nid],
                                ani))
            fout.close()
            
            sys.exit(0)

        # identify statistics relative to representative genome
        rep_stats = self._rep_genome_stats(clusters, genome_files)
        
        # identify pairwise statistics
        pairwise_stats = self._pairwise_stats(clusters, genome_files, rep_stats)
        
        # report statistics
        stats_file = os.path.join(self.output_dir, 'cluster_stats.tsv')
        self._write_cluster_stats(stats_file,
                                    clusters,
                                    species,
                                    cluster_radius,
                                    rep_stats,
                                    pairwise_stats)
        