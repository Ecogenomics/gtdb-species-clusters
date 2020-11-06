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
import ntpath
import pickle
import operator
from itertools import combinations, permutations
from collections import defaultdict, namedtuple

from biolib.external.execute import check_dependencies

from numpy import (mean as np_mean,
                    std as np_std)

from gtdb_species_clusters.mash import Mash
from gtdb_species_clusters.fastani import FastANI
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.type_genome_utils import (ClusteredGenome,
                                                        GenomeRadius,
                                                        write_rep_radius,
                                                        write_clusters)
                                                        
from gtdb_species_clusters.genome_utils import canonical_gid


class IntraGenusANI(object):
    """Calculate intra-genus ANI/AF values between GTDB representative genomes."""

    def __init__(self,
                    ani_cache_file, 
                    cpus, 
                    output_dir):
        """Initialization."""
        
        check_dependencies(['fastANI'])
        
        self.cpus = cpus
        self.output_dir = output_dir

        self.logger = logging.getLogger('timestamp')

        self.fastani = FastANI(ani_cache_file, cpus)
        
        self.user_id_map = {}

    def run(self, target_genus,
                    gtdb_metadata_file,
                    genomic_path_file,
                    uba_gid_table):
        """Dereplicate GTDB species clusters using ANI/AF criteria."""
        
        # map user IDs to UBA IDs
        with open(uba_gid_table) as f:
            for line in f:
                tokens = line.strip().split('\t')
                
                if len(tokens) == 3:
                    self.user_id_map[tokens[0]] = tokens[2]
                else:
                    self.user_id_map[tokens[0]] = tokens[1]
        
        # create GTDB genome sets
        self.logger.info('Creating GTDB genome set.')
        genomes = Genomes()
        genomes.load_from_metadata_file(gtdb_metadata_file,
                                        uba_genome_file=uba_gid_table)
        genomes.load_genomic_file_paths(genomic_path_file)
        self.logger.info(' - genome set has {:,} species clusters spanning {:,} genomes.'.format(
                            len(genomes.sp_clusters),
                            genomes.sp_clusters.total_num_genomes()))
                            
        # identify GTDB representatives from target genus
        self.logger.info('Identifying GTDB representatives from target genus.')
        target_gids = set()
        for gid in genomes:
            if genomes[gid].is_gtdb_sp_rep() and genomes[gid].gtdb_taxa.genus == target_genus:
                target_gids.add(gid)
        self.logger.info(' - identified {:,} genomes.'.format(len(target_gids)))

        # calculate FastANI ANI/AF between target genomes
        self.logger.info('Calculating pairwise ANI between target genomes.')
        ani_af = self.fastani.pairwise(target_gids, 
                                        genomes.genomic_files, 
                                        check_cache=True)
        self.fastani.write_cache(silence=True)
    
        # write out results
        genus_label = target_genus.replace('g__', '').lower()
        fout = open(os.path.join(self.output_dir, '{}_rep_ani.tsv'.format(genus_label)), 'w')
        fout.write('Query ID\tQuery species\tTarget ID\tTarget species\tANI\tAF\n')
        for qid in target_gids:
            for rid in target_gids:
                ani, af = FastANI.symmetric_ani(ani_af, qid, rid)
                
                fout.write('{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\n'.format(
                            qid,
                            genomes[qid].gtdb_taxa.species,
                            rid,
                            genomes[rid].gtdb_taxa.species,
                            ani,
                            af))
        fout.close()
                            
                