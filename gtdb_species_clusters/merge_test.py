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

from gtdb_species_clusters.taxon_utils import read_gtdb_taxonomy
from gtdb_species_clusters.type_genome_utils import symmetric_ani
from gtdb_species_clusters.fastani import FastANI
from gtdb_species_clusters.genomes import Genomes


class MergeTest(object):
    """Produce information relevant to merging two sister species."""

    def __init__(self, ani_cache_file, cpus, output_dir):
        """Initialization."""
        
        check_dependencies(['fastANI'])
        
        self.cpus = cpus
        self.output_dir = output_dir
        
        self.logger = logging.getLogger('timestamp')

        self.fastani = FastANI(ani_cache_file, cpus)
        
    def top_hits(self, species, rid, ani_af, genomes):
        """Report top 5 hits to species."""
        
        results = {}
        for qid in ani_af[rid]:
            ani, af = symmetric_ani(ani_af, rid, qid)
            results[qid] = (ani, af)
            
        self.logger.info(f'Closest 5 species to {species} ({rid}):')
        idx = 0
        for qid, (ani, af) in sorted(results.items(), key=lambda x: x[1], reverse=True):
            q_species = genomes[qid].gtdb_species
            self.logger.info(f'{q_species} ({qid}): ANI={ani:.1f}%, AF={af:.2f}')
            if idx == 5:
                break
                
            idx += 1
            
    def merge_ani_radius(self, species, rid, merged_sp_cluster, genomic_files):
        """Determine ANI radius if species were merged."""
        
        self.logger.info(f'Calculating ANI from {species} to all genomes in merged species cluster.')
        
        gid_pairs = []
        for gid in merged_sp_cluster:
            gid_pairs.append((rid, gid))
            gid_pairs.append((gid, rid))
        merged_ani_af1 = self.fastani.pairs(gid_pairs, genomic_files)
        
        ani_radius = 100
        for gid in merged_sp_cluster:
            ani, af = symmetric_ani(merged_ani_af1, rid, gid)
            if ani < ani_radius:
                ani_radius = ani
                af_radius = af
        self.logger.info(f'Merged cluster with {species} rep: ANI radius={ani_radius:.1f}%, AF={af_radius:.2f}')

    def run(self, gtdb_metadata_file,
                    genome_path_file,
                    species1,
                    species2):
        """Produce information relevant to merging two sister species."""
 
        # read GTDB species clusters
        self.logger.info('Reading GTDB species clusters.')
        genomes = Genomes()
        genomes.load_from_metadata_file(gtdb_metadata_file)
        genomes.load_genomic_file_paths(genome_path_file)
        self.logger.info(' ... identified {:,} species clusters spanning {:,} genomes.'.format(
                            len(genomes.sp_clusters),
                            genomes.sp_clusters.total_num_genomes()))
                            
        # find species of interest
        gid1 = None
        gid2 = None
        for gid, species in genomes.sp_clusters.species():
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
            
        self.logger.info(' ... identified {:,} genomes in {}.'.format(
                            len(genomes.sp_clusters[gid1]),
                            species1))
        self.logger.info(' ... identified {:,} genomes in {}.'.format(
                            len(genomes.sp_clusters[gid2]),
                            species2))

        # calculate ANI between all genome in genus
        genus1 = genomes[gid1].gtdb_genus
        genus2 = genomes[gid2].gtdb_genus
        if genus1 != genus2:
            self.logger.error(f'Genomes must be from same genus: {genus1} {genus2}')
            sys.exit(-1)
        
        self.logger.info(f'Identifying {genus1} species representatives.')
        reps_in_genera = set()
        for rid in genomes.sp_clusters:
            if genomes[rid].gtdb_genus == genus1:
                reps_in_genera.add(rid)
        
        self.logger.info(f' ... identified {len(reps_in_genera):,} representatives.')
        
        # calculate ANI between genomes
        self.logger.info(f'Calculating ANI to {species1}.')
        gid_pairs = []
        for gid in reps_in_genera:
            if gid != gid1:
                gid_pairs.append((gid1, gid))
                gid_pairs.append((gid, gid1))
        ani_af1 = self.fastani.pairs(gid_pairs, genomes.genomic_files)
        
        self.logger.info(f'Calculating ANI to {species2}.')
        gid_pairs = []
        for gid in reps_in_genera:
            if gid != gid2:
                gid_pairs.append((gid2, gid))
                gid_pairs.append((gid, gid2))
        ani_af2 = self.fastani.pairs(gid_pairs, genomes.genomic_files)

        # report results
        ani12, af12 = ani_af1[gid1][gid2]
        ani21, af21 = ani_af2[gid2][gid1]
        ani, af = symmetric_ani(ani_af1, gid1, gid2)
        
        self.logger.info(f'{species1} ({gid1}) -> {species2} ({gid2}): ANI={ani12:.1f}%, AF={af12:.2f}')
        self.logger.info(f'{species2} ({gid2}) -> {species1} ({gid1}): ANI={ani21:.1f}%, AF={af21:.2f}')
        self.logger.info(f'Max. ANI={ani:.1f}%, Max. AF={af:.2f}')
        
        # report top hits
        self.top_hits(species1, gid1, ani_af1, genomes)
        self.top_hits(species2, gid2, ani_af2, genomes)
        
        # calculate ANI from species to all genomes in merged species cluster
        merged_sp_cluster = genomes.sp_clusters[gid1].union(genomes.sp_clusters[gid2])
        self.merge_ani_radius(species1, gid1, merged_sp_cluster, genomes.genomic_files)
        self.merge_ani_radius(species2, gid2, merged_sp_cluster, genomes.genomic_files)
        