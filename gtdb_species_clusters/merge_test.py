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
import shutil
import tempfile
import ntpath
import pickle
from itertools import combinations
from collections import defaultdict, namedtuple

from biolib.taxonomy import Taxonomy
from biolib.external.execute import check_dependencies

from numpy import (mean as np_mean)

from gtdb_species_clusters.genome_utils import (read_gtdb_sp_clusters,
                                                read_genome_path)
from gtdb_species_clusters.taxon_utils import read_gtdb_taxonomy
                                    
from gtdb_species_clusters.type_genome_utils import symmetric_ani
                                    
from gtdb_species_clusters.ani_cache import ANI_Cache


class MergeTest(object):
    """Produce information relevant to merging two sister species."""

    def __init__(self, ani_cache_file, cpus, output_dir):
        """Initialization."""
        
        check_dependencies(['fastANI'])
        
        self.cpus = cpus
        self.output_dir = output_dir
        
        self.logger = logging.getLogger('timestamp')

        self.ani_cache = ANI_Cache(ani_cache_file, cpus)
        
    def top_hits(self, species, rid, ani_af, gtdb_taxonomy):
        """Report top 5 hits to species."""
        
        results = {}
        for qid in ani_af[rid]:
            ani, af = symmetric_ani(ani_af, rid, qid)
            results[qid] = (ani, af)
            
        self.logger.info(f'Closest 5 species to {species} ({rid}):')
        idx = 0
        for qid, (ani, af) in sorted(results.items(), key=lambda x: x[1], reverse=True):
            q_species = gtdb_taxonomy[qid][6]
            self.logger.info(f'{q_species} ({qid}): ANI={ani:.1f}%, AF={af:.2f}')
            if idx == 5:
                break
                
            idx += 1

    def run(self, sp_cluster_file,
                    gtdb_metadata_file,
                    genome_path_file,
                    species1,
                    species2):
        """Produce information relevant to merging two sister species."""
 
        # read GTDB species clusters
        self.logger.info('Reading GTDB species clusters.')
        sp_clusters, species = read_gtdb_sp_clusters(sp_cluster_file)
        self.logger.info(' ... identified {:,} species clusters spanning {:,} genomes.'.format(
                            len(sp_clusters),
                            sum([len(cids) for cids in sp_clusters.values()])))
                            
        # find species of interest
        gid1 = None
        gid2 = None
        for gid, species in species.items():
            if species == species1:
                gid1 = gid
            elif species == species2:
                gid2 = gid
                
        if gid1 is None:
            self.logger.error(f'Unable to find representative genome for {species1}.')
            sys.exit(-1)
            
        if gid2 is None:
            self.logger.error(f'Unable to find representative genome for {species2}.')
            sys.exit(-1)
            
        # read GTDB taxonomy
        gtdb_taxonomy = read_gtdb_taxonomy(gtdb_metadata_file)
            
        # calculate ANI between all genome in genus
        genus1 = gtdb_taxonomy[gid1][5]
        genus2 = gtdb_taxonomy[gid2][5]
        if genus1 != genus2:
            self.logger.error(f'Genomes must be from same genus: {genus1} {genus2}')
            sys.exit(-1)
        
        self.logger.info(f'Identifying {genus1} species representatives.')
        reps_in_genera = set()
        for gid, taxa in gtdb_taxonomy.items():
            if gid not in sp_clusters:
                continue
                
            if gtdb_taxonomy[gid][5] == genus1:
                reps_in_genera.add(gid)
        
        self.logger.info(f' ... identified {len(reps_in_genera):,} representatives.')
        
        # read path to genomic FASTA files
        genomic_files = read_genome_path(genome_path_file)
        
        # calculate ANI between genomes
        self.logger.info(f'Calculating ANI to {species1}.')
        gid_pairs = []
        for gid in reps_in_genera:
            if gid != gid1:
                gid_pairs.append((gid1, gid))
                gid_pairs.append((gid, gid1))
        ani_af1 = self.ani_cache.fastani_pairs(gid_pairs, genomic_files)
        
        self.logger.info(f'Calculating ANI to {species2}.')
        gid_pairs = []
        for gid in reps_in_genera:
            if gid != gid2:
                gid_pairs.append((gid2, gid))
                gid_pairs.append((gid, gid2))
        ani_af2 = self.ani_cache.fastani_pairs(gid_pairs, genomic_files)
        
        # report results
        ani12, af12 = ani_af1[gid1][gid2]
        ani21, af21 = ani_af2[gid2][gid1]
        ani, af = symmetric_ani(ani_af1, gid1, gid2)
        
        self.logger.info(f'{species1} ({gid1}) -> {species2} ({gid2}): ANI={ani12:.1f}%, AF={af12:.2f}')
        self.logger.info(f'{species2} ({gid2}) -> {species1} ({gid1}): ANI={ani21:.1f}%, AF={af21:.2f}')
        self.logger.info(f'Max. ANI={ani:.1f}%, Max. AF={af:.2f}')
        
        # report top hits
        self.top_hits(species1, gid1, ani_af1, gtdb_taxonomy)
        self.top_hits(species2, gid2, ani_af2, gtdb_taxonomy)
        