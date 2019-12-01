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
import shutil
import tempfile
import ntpath
import pickle
from itertools import combinations
from collections import defaultdict, namedtuple, Counter

from biolib.external.execute import check_dependencies

from numpy import (mean as np_mean,
                    std as np_std)

from gtdb_species_clusters.genome_utils import read_genome_path
                                    
from gtdb_species_clusters.type_genome_utils import (GenomeRadius,
                                                        read_quality_metadata,
                                                        read_clusters,
                                                        symmetric_ani)
                                    
from gtdb_species_clusters.fastani import FastANI
from gtdb_species_clusters.mash import Mash

class ClusterUser(object):
    """Cluster User genomes to GTDB species clusters."""

    def __init__(self, ani_cache_file, cpus, output_dir):
        """Initialization."""
        
        check_dependencies(['fastANI', 'mash'])
        
        self.cpus = cpus
        self.output_dir = output_dir

        self.logger = logging.getLogger('timestamp')

        self.min_mash_ani = 90.0
        
        self.af_sp = 0.65

        self.fastani = FastANI(ani_cache_file, cpus)

    def _mash_ani(self, genome_files, user_genomes, sp_clusters):
        """Calculate Mash ANI estimates between User genomes and species clusters."""
        
        mash = Mash(self.cpus)
        
        # create Mash sketch for User genomes
        mash_user_sketch_file = os.path.join(self.output_dir, 'gtdb_user_genomes.msh')
        genome_list_file = os.path.join(self.output_dir, 'gtdb_user_genomes.lst')
        mash.sketch(user_genomes, genome_files, genome_list_file, mash_user_sketch_file)
        
        # create Mash sketch for species clusters
        mash_sp_sketch_file = os.path.join(self.output_dir, 'gtdb_sp_genomes.msh')
        genome_list_file = os.path.join(self.output_dir, 'gtdb_sp_genomes.lst')
        mash.sketch(sp_clusters, genome_files, genome_list_file, mash_sp_sketch_file)

        # get Mash distances
        mash_dist_file = os.path.join(self.output_dir, 'gtdb_user_vs_sp.dst')
        mash.dist(float(100 - self.min_mash_ani)/100, mash_sp_sketch_file, mash_user_sketch_file, mash_dist_file)

        # read Mash distances
        mash_ani = mash.read_ani(mash_dist_file)
        
        # report pairs above Mash threshold
        mash_ani_pairs = []
        for qid in mash_ani:
            for rid in mash_ani[qid]:
                if mash_ani[qid][rid] >= self.min_mash_ani:
                    if qid != rid:
                        mash_ani_pairs.append((qid, rid))
                        mash_ani_pairs.append((rid, qid))
                
        self.logger.info('Identified %d genome pairs with a Mash ANI >= %.1f%%.' % (len(mash_ani_pairs), self.min_mash_ani))

        return mash_ani
        
    def _cluster(self,
                    genome_files,
                    sp_clusters,
                    rep_radius, 
                    user_genomes, 
                    mash_anis):
        """Cluster User genomes to existing species clusters."""
        
        # assign User genomes to closest species cluster
        
        for idx, cur_gid in enumerate(user_genomes):
            # determine species cluster to calculate ANI between
            ani_pairs = []
            if cur_gid in mash_anis:
                for rep_gid in sp_clusters:
                    if mash_anis[cur_gid].get(rep_gid, 0) >= self.min_mash_ani:
                        ani_pairs.append((cur_gid, rep_gid))
                        ani_pairs.append((rep_gid, cur_gid))

                # determine if genome clusters with representative
                clustered = False
                if ani_pairs:
                    ani_af = self.fastani.pairs(ani_pairs, genome_files, report_progress=False)

                    closest_rep_gid = None
                    closest_rep_ani = 0
                    closest_rep_af = 0
                    for rep_gid in sp_clusters:
                        ani, af = symmetric_ani(ani_af, cur_gid, rep_gid)
                        
                        if af >= self.af_sp:
                            if ani > closest_rep_ani or (ani == closest_rep_ani and af > closest_rep_af):
                                closest_rep_gid = rep_gid
                                closest_rep_ani = ani
                                closest_rep_af = af
                        
                    if closest_rep_gid and closest_rep_ani > rep_radius[closest_rep_gid].ani:
                        sp_clusters[closest_rep_gid].append(cur_gid)
                    else:
                        self.logger.warning('Failed to assign genome %s to representative.' % cur_gid)
                     
            statusStr = '-> Assigned %d of %d (%.2f%%) genomes.'.ljust(86) % (idx+1, 
                                                                                len(user_genomes), 
                                                                                float(idx+1)*100/len(user_genomes))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
        sys.stdout.write('\n')

    def run(self, 
                metadata_file,
                genome_path_file,
                final_cluster_file):
        """Cluster User genomes to GTDB species clusters."""

        # get path to genome FASTA files
        self.logger.info('Reading path to genome FASTA files.')
        genome_files = read_genome_path(genome_path_file)
        
        # read existing cluster information
        self.logger.info('Reading already established species clusters.')
        sp_clusters, species, rep_radius = read_clusters(final_cluster_file)
        
        clustered_genomes = set()
        for rep_id in sp_clusters:
            clustered_genomes.add(rep_id)
            clustered_genomes.update(sp_clusters[rep_id])
        
        self.logger.info('Identified %d species clusters spanning %d genomes.' % (len(sp_clusters), len(clustered_genomes)))

        # get User genomes to cluster
        self.logger.info('Parse quality statistics for all genomes.')
        quality_metadata = read_quality_metadata(metadata_file)
        
        user_genomes = set()
        for gid in quality_metadata:
            if gid in clustered_genomes:
                continue
                
            if (quality_metadata[gid].checkm_completeness > 50
                and quality_metadata[gid].checkm_contamination < 10):
                    user_genomes.add(gid)
                    
        self.logger.info('Identified %d User genomes to cluster.' % len(user_genomes))

        # calculate Mash ANI estimates between unclustered genomes
        self.logger.info('Calculating Mash ANI estimates between User genomes and species clusters.')
        mash_anis = self._mash_ani(genome_files, user_genomes, sp_clusters)

        # cluster User genomes to species clusters
        self.logger.info('Assigning User genomes to closest species cluster.')
        self._cluster(genome_files,
                        sp_clusters,
                        rep_radius, 
                        user_genomes, 
                        mash_anis)
                        
        clustered_genomes = 0
        for rep_id in sp_clusters:
            clustered_genomes += 1
            clustered_genomes += len(sp_clusters[rep_id])
            
        self.logger.info('The %d species clusters span %d genomes, including User genomes.' % (len(sp_clusters), clustered_genomes))

        # report clustering
        user_cluster_file = os.path.join(self.output_dir, 'gtdb_user_clusters.tsv')
        fout = open(user_cluster_file, 'w')
        fout.write('Type genome\tNo. clustered genomes\tClustered genomes\n')
        for rep_id in sp_clusters:
            fout.write('%s\t%d\t%s\n' % (rep_id, len(sp_clusters[rep_id]), ','.join(sp_clusters[rep_id])))
        fout.close()
            