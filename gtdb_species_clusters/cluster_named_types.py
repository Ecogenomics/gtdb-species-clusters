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

from gtdb_species_clusters.common import read_gtdb_metadata
                                            
from gtdb_species_clusters.genome_utils import (read_genome_path,
                                                read_qc_file)
                                            
from gtdb_species_clusters.taxon_utils import (read_gtdb_taxonomy,
                                                read_gtdb_ncbi_taxonomy)
                                    
from gtdb_species_clusters.type_genome_utils import (GenomeRadius,
                                            symmetric_ani,
                                            write_clusters,
                                            write_rep_radius)
                                    
from gtdb_species_clusters.fastani import FastANI
from gtdb_species_clusters.mash import Mash

class ClusterNamedTypes(object):
    """Cluster genomes to selected GTDB type genomes."""

    def __init__(self, ani_sp, af_sp, ani_cache_file, cpus, output_dir):
        """Initialization."""
        
        check_dependencies(['fastANI', 'mash'])
        
        self.cpus = cpus
        self.output_dir = output_dir

        self.logger = logging.getLogger('timestamp')

        self.ani_sp = ani_sp
        self.af_sp = af_sp

        self.max_ani_neighbour = 97.0
        self.min_mash_ani = 90.0
        
        self.ClusteredGenome = namedtuple('ClusteredGenome', 'ani af gid')
        
        self.fastani = FastANI(ani_cache_file, cpus)
        
    def _type_genome_radius(self, type_gids, type_genome_ani_file):
        """Calculate circumscription radius for type genomes."""
        
        # set type radius for all type genomes to default values
        type_radius = {}
        for gid in type_gids:
            type_radius[gid] = GenomeRadius(ani = self.ani_sp, 
                                                 af = None,
                                                 neighbour_gid = None)
        
        # determine closest ANI neighbour and restrict ANI radius as necessary
        with open(type_genome_ani_file) as f:
            header = f.readline().strip().split('\t')
            
            type_gid1_index = header.index('Type genome 1')
            type_gid2_index = header.index('Type genome 2')
            ani_index = header.index('ANI')
            af_index = header.index('AF')

            for line in f:
                line_split = line.strip().split('\t')
                
                type_gid1 = line_split[type_gid1_index]
                type_gid2 = line_split[type_gid2_index]

                if type_gid1 not in type_gids or type_gid2 not in type_gids:
                    continue

                ani = float(line_split[ani_index])
                af = float(line_split[af_index])

                if ani > type_radius[type_gid1].ani:
                    if af < self.af_sp:
                        if ani >= self.ani_sp:
                            self.logger.warning('ANI for %s and %s is >%.2f, but AF <%.2f [pair skipped].' % (
                                                    type_gid1,
                                                    type_gid2,
                                                    ani, af))
                        continue
                    
                    if ani > self.max_ani_neighbour:
                        self.logger.error('ANI neighbour %s is >%.2f for %s.' % (type_gid2, ani, type_gid1))
 
                    type_radius[type_gid1] = GenomeRadius(ani = ani, 
                                                                 af = af,
                                                                 neighbour_gid = type_gid2)
                    
        self.logger.info('ANI circumscription radius: min=%.2f, mean=%.2f, max=%.2f' % (
                                min([d.ani for d in type_radius.values()]), 
                                np_mean([d.ani for d in type_radius.values()]), 
                                max([d.ani for d in type_radius.values()])))
                        
        return type_radius
        
    def _calculate_ani(self, type_gids, genome_files, ncbi_taxonomy, type_genome_sketch_file):
        """Calculate ANI between type and non-type genomes."""
        
        mash = Mash(self.cpus)
        
        # create Mash sketch for type genomes
        if not type_genome_sketch_file or not os.path.exists(type_genome_sketch_file):
            type_genome_list_file = os.path.join(self.output_dir, 'gtdb_type_genomes.lst')
            type_genome_sketch_file = os.path.join(self.output_dir, 'gtdb_type_genomes.msh')
            mash.sketch(type_gids, genome_files, type_genome_list_file, type_genome_sketch_file)
            
        # create Mash sketch for non-type genomes
        nontype_gids = set()
        for gid in genome_files:
            if gid not in type_gids:
                nontype_gids.add(gid)
                
        nontype_genome_list_file = os.path.join(self.output_dir, 'gtdb_nontype_genomes.lst')
        nontype_genome_sketch_file = os.path.join(self.output_dir, 'gtdb_nontype_genomes.msh')
        mash.sketch(nontype_gids, genome_files, nontype_genome_list_file, nontype_genome_sketch_file)

        # get Mash distances
        mash_dist_file = os.path.join(self.output_dir, 'gtdb_type_vs_nontype_genomes.dst')
        mash.dist(float(100 - self.min_mash_ani)/100, 
                                type_genome_sketch_file, 
                                nontype_genome_sketch_file, 
                                mash_dist_file)

        # read Mash distances
        mash_ani = mash.read_ani(mash_dist_file)

        # get pairs above Mash threshold
        mash_ani_pairs = []
        for qid in mash_ani:
            for rid in mash_ani[qid]:
                if mash_ani[qid][rid] >= self.min_mash_ani:
                    if qid != rid:
                        mash_ani_pairs.append((qid, rid))
                        mash_ani_pairs.append((rid, qid))
                
        self.logger.info('Identified %d genome pairs with a Mash ANI >= %.1f%%.' % (len(mash_ani_pairs), self.min_mash_ani))
        
        # calculate ANI between pairs
        self.logger.info('Calculating ANI between %d genome pairs:' % len(mash_ani_pairs))
        if True: #***
            ani_af = self.fastani.pairs(mash_ani_pairs, genome_files)
            pickle.dump(ani_af, open(os.path.join(self.output_dir, 'ani_af_type_vs_nontype.pkl'), 'wb'))
        else:
            ani_af = pickle.load(open(os.path.join(self.output_dir, 'ani_af_type_vs_nontype.pkl'), 'rb'))

        return ani_af

    def _cluster(self, ani_af, nontype_gids, type_radius):
        """Cluster non-type genomes to type genomes using species specific ANI thresholds."""
        
        clusters = {}
        for rep_id in type_radius:
            clusters[rep_id] = []
            
        for idx, nontype_gid in enumerate(nontype_gids):
            if idx % 100 == 0:
                sys.stdout.write('==> Processed %d of %d genomes.\r' % (idx+1, len(nontype_gids)))
                sys.stdout.flush()
                
            if nontype_gid not in ani_af:
                continue

            closest_type_gid = None
            closest_ani = 0
            closest_af = 0
            for type_gid in type_radius:
                if type_gid not in ani_af[nontype_gid]:
                    continue

                ani, af = symmetric_ani(ani_af, type_gid, nontype_gid)
                
                if af >= self.af_sp:
                    if ani > closest_ani or (ani == closest_ani and af > closest_af):
                        closest_type_gid = type_gid
                        closest_ani = ani
                        closest_af = af
                
            if closest_type_gid:
                if closest_ani > type_radius[closest_type_gid].ani:
                    clusters[closest_type_gid].append(self.ClusteredGenome(gid=nontype_gid, 
                                                                            ani=closest_ani, 
                                                                            af=closest_af))
                
        sys.stdout.write('==> Processed %d of %d genomes.\r' % (idx, 
                                                                len(nontype_gids)))
        sys.stdout.flush()
        sys.stdout.write('\n')

        self.logger.info('Assigned %d genomes to representatives.' % sum([len(clusters[type_gid]) for type_gid in clusters]))
        
        return clusters

    def run(self, qc_file,
                    metadata_file,
                    genome_path_file,
                    named_type_genome_file,
                    type_genome_ani_file,
                    mash_sketch_file,
                    species_exception_file):
        """Cluster genomes to selected GTDB type genomes."""
        
        # identify genomes failing quality criteria
        self.logger.info('Reading QC file.')
        passed_qc = read_qc_file(qc_file)
        self.logger.info('Identified %d genomes passing QC.' % len(passed_qc))

        # get type genomes
        type_gids = set()
        species_type_gid = {}
        with open(named_type_genome_file) as f:
            header = f.readline().strip().split('\t')
            type_gid_index = header.index('Type genome')
            sp_index = header.index('NCBI species')
            
            for line in f:
                line_split = line.strip().split('\t')
                type_gids.add(line_split[type_gid_index])
                species_type_gid[line_split[type_gid_index]] = line_split[sp_index]
        self.logger.info('Identified type genomes for %d species.' % len(species_type_gid))

        # calculate circumscription radius for type genomes
        self.logger.info('Determining ANI species circumscription for %d type genomes.' % len(type_gids))
        type_radius = self._type_genome_radius(type_gids, type_genome_ani_file)
        assert(len(type_radius) == len(species_type_gid))
        
        write_rep_radius(type_radius, species_type_gid, os.path.join(self.output_dir, 'gtdb_type_genome_ani_radius.tsv'))
        
        # get path to genome FASTA files
        self.logger.info('Reading path to genome FASTA files.')
        genome_files = read_genome_path(genome_path_file)
        self.logger.info('Read path for %d genomes.' % len(genome_files))
        for gid in set(genome_files):
            if gid not in passed_qc:
                genome_files.pop(gid)
        self.logger.info('Considering %d genomes after removing unwanted User genomes.' % len(genome_files))
        assert(len(genome_files) == len(passed_qc))
        
        # get GTDB and NCBI taxonomy strings for each genome
        self.logger.info('Reading NCBI taxonomy from GTDB metadata file.')
        ncbi_taxonomy, ncbi_update_count = read_gtdb_ncbi_taxonomy(metadata_file, species_exception_file)
        self.logger.info('Read NCBI taxonomy for %d genomes with %d manually defined updates.' % (len(ncbi_taxonomy), ncbi_update_count))
        
        # calculate ANI between type and non-type genomes
        self.logger.info('Calculating ANI between type and non-type genomes.')
        ani_af = self._calculate_ani(type_gids, genome_files, ncbi_taxonomy, mash_sketch_file)

        # cluster remaining genomes to type genomes
        nontype_gids = set(genome_files) - set(type_radius)
        self.logger.info('Clustering %d non-type genomes to type genomes using species specific ANI radii.' % len(nontype_gids))
        clusters = self._cluster(ani_af, nontype_gids, type_radius)
        
        # write out clusters
        write_clusters(clusters, type_radius, species_type_gid, os.path.join(self.output_dir, 'gtdb_type_genome_clusters.tsv'))
