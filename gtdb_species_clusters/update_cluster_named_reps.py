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

from gtdb_species_clusters.mash import Mash
from gtdb_species_clusters.fastani import FastANI

from gtdb_species_clusters.genome import Genome
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_clusters import SpeciesClusters
                                    
from gtdb_species_clusters.genome_utils import (canonical_gid,
                                                read_qc_file,
                                                exclude_from_refseq)

from gtdb_species_clusters.type_genome_utils import (GenomeRadius,
                                                        symmetric_ani,
                                                        write_rep_radius,
                                                        write_clusters)
                                                        


class UpdateClusterNamedReps(object):
    """Cluster genomes to selected GTDB representatives."""

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
        
    def _rep_radius(self, rep_gids, rep_ani_file):
        """Calculate circumscription radius for representative genomes."""
        
        # set radius for all representative genomes to default values
        rep_radius = {}
        for gid in rep_gids:
            rep_radius[gid] = GenomeRadius(ani = self.ani_sp, 
                                                 af = None,
                                                 neighbour_gid = None)
        
        # determine closest ANI neighbour and restrict ANI radius as necessary
        with open(rep_ani_file) as f:
            header = f.readline().strip().split('\t')
            
            rep_gid1_index = header.index('Representative 1')
            rep_gid2_index = header.index('Representative 2')
            ani_index = header.index('ANI')
            af_index = header.index('AF')

            for line in f:
                line_split = line.strip().split('\t')
                
                rep_gid1 = line_split[rep_gid1_index]
                rep_gid2 = line_split[rep_gid2_index]

                if rep_gid1 not in rep_gids or rep_gid2 not in rep_gids:
                    continue

                ani = float(line_split[ani_index])
                af = float(line_split[af_index])

                if ani > rep_radius[rep_gid1].ani:
                    if af < self.af_sp:
                        if ani >= self.ani_sp:
                            self.logger.warning('ANI for {} and {} is >{:.2f}, but AF <{:.2f} [pair skipped].'.format(
                                                    rep_gid1,
                                                    rep_gid2,
                                                    ani, af))
                        continue
                    
                    if ani > self.max_ani_neighbour:
                        self.logger.error('ANI neighbour {} is >{:.2f} for {}.'.format(rep_gid2, ani, rep_gid1))
 
                    rep_radius[rep_gid1] = GenomeRadius(ani = ani, 
                                                         af = af,
                                                         neighbour_gid = rep_gid2)
                    
        self.logger.info('ANI circumscription radius: min={:.2f}, mean={:.2f}, max={:.2f}'.format(
                                min([d.ani for d in rep_radius.values()]), 
                                np_mean([d.ani for d in rep_radius.values()]), 
                                max([d.ani for d in rep_radius.values()])))
                        
        return rep_radius
        
    def _calculate_ani(self, cur_genomes, rep_gids, rep_mash_sketch_file):
        """Calculate ANI between representative and non-representative genomes."""
        
        if False: #***
            mash = Mash(self.cpus)
            
            # create Mash sketch for representative genomes
            if not rep_mash_sketch_file or not os.path.exists(rep_mash_sketch_file):
                rep_genome_list_file = os.path.join(self.output_dir, 'gtdb_reps.lst')
                rep_mash_sketch_file = os.path.join(self.output_dir, 'gtdb_reps.msh')
                mash.sketch(rep_gids, cur_genomes.genomic_files, rep_genome_list_file, rep_mash_sketch_file)
                
            # create Mash sketch for non-representative genomes
            nonrep_gids = set()
            for gid in cur_genomes:
                if gid not in rep_gids:
                    nonrep_gids.add(gid)
                    
            nonrep_genome_list_file = os.path.join(self.output_dir, 'gtdb_nonreps.lst')
            nonrep_genome_sketch_file = os.path.join(self.output_dir, 'gtdb_nonreps.msh')
            mash.sketch(nonrep_gids, cur_genomes.genomic_files, nonrep_genome_list_file, nonrep_genome_sketch_file)

            # get Mash distances
            mash_dist_file = os.path.join(self.output_dir, 'gtdb_reps_vs_nonreps.dst')
            mash.dist(float(100 - self.min_mash_ani)/100, 
                                    rep_mash_sketch_file, 
                                    nonrep_genome_sketch_file, 
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
                    
            self.logger.info('Identified {:,} genome pairs with a Mash ANI >= {:.1f}%.'.format(len(mash_ani_pairs), self.min_mash_ani))
            
            # calculate ANI between pairs
            self.logger.info('Calculating ANI between {:,} genome pairs:'.format(len(mash_ani_pairs)))
        
            ani_af = self.fastani.pairs(mash_ani_pairs, cur_genomes.genomic_files)
            pickle.dump(ani_af, open(os.path.join(self.output_dir, 'ani_af_rep_vs_nonrep.pkl'), 'wb'))
        else:
            self.logger.warning('Using previously calculated results in: {}'.format('ani_af_rep_vs_nonrep.pkl'))
            ani_af = pickle.load(open(os.path.join(self.output_dir, 'ani_af_rep_vs_nonrep.pkl'), 'rb'))

        return ani_af

    def _cluster(self, ani_af, non_reps, rep_radius):
        """Cluster non-representative to representative genomes using species specific ANI thresholds."""
        
        clusters = {}
        for rep_id in rep_radius:
            clusters[rep_id] = []
            
        for idx, non_rep in enumerate(non_reps):
            if idx % 100 == 0:
                sys.stdout.write('==> Processed {:,} of {:,} genomes.\r'.format(idx+1, len(non_reps)))
                sys.stdout.flush()
                
            if non_rep not in ani_af:
                continue

            closest_rep = None
            closest_ani = 0
            closest_af = 0
            for rep in rep_radius:
                if rep not in ani_af[non_rep]:
                    continue

                ani, af = symmetric_ani(ani_af, rep, non_rep)
                
                if af >= self.af_sp:
                    if ani > closest_ani or (ani == closest_ani and af > closest_af):
                        closest_rep = rep
                        closest_ani = ani
                        closest_af = af
                
            if closest_rep:
                if closest_ani > rep_radius[closest_rep].ani:
                    clusters[closest_rep].append(self.ClusteredGenome(gid=non_rep, 
                                                                        ani=closest_ani, 
                                                                        af=closest_af))
                
        sys.stdout.write('==> Processed {:,} of {:,} genomes.\r'.format(idx, len(non_reps)))
        sys.stdout.flush()
        sys.stdout.write('\n')

        self.logger.info('Assigned {:,} genomes to representatives.'.format(sum([len(clusters[rep]) for reps in clusters])))
        
        return clusters

    def run(self, qc_passed_file,
                    cur_gtdb_metadata_file,
                    cur_genomic_path_file,
                    uba_genome_paths,
                    named_rep_file,
                    rep_mash_sketch_file,
                    rep_ani_file,
                    species_exception_file,
                    genus_exception_file,
                    gtdb_type_strains_ledger):
        """Cluster genomes to selected GTDB representatives."""
        
        # create current GTDB genome sets
        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                                species_exception_file,
                                                genus_exception_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                create_sp_clusters=False,
                                                uba_genome_file=uba_genome_paths,
                                                qc_passed_file=qc_passed_file)
        self.logger.info(f' ... current genome set contains {len(cur_genomes):,} genomes.')
        
        # get path to previous and current genomic FASTA files
        self.logger.info('Reading path to current genomic FASTA files.')
        cur_genomes.load_genomic_file_paths(cur_genomic_path_file)
        cur_genomes.load_genomic_file_paths(uba_genome_paths)

        # get representative genomes
        rep_gids = set()
        proposed_sp = {}
        with open(named_rep_file) as f:
            header = f.readline().strip().split('\t')
            rep_index = header.index('Representative')
            sp_index = header.index('Proposed species')

            for line in f:
                line_split = line.strip().split('\t')
                gid = line_split[rep_index]
                rep_gids.add(gid)
                proposed_sp[gid] = line_split[sp_index]
        self.logger.info('Identified representative genomes for {:,} species.'.format(len(rep_gids)))

        # calculate circumscription radius for representative genomes
        self.logger.info('Determining ANI species circumscription for {:,} representative genomes.'.format(len(rep_gids)))
        rep_radius = self._rep_radius(rep_gids, rep_ani_file)      
        write_rep_radius(rep_radius, proposed_sp, cur_genomes, os.path.join(self.output_dir, 'gtdb_rep_ani_radius.tsv'))

        # calculate ANI between representative and non-representative genomes
        self.logger.info('Calculating ANI between representative and non-representative genomes.')
        ani_af = self._calculate_ani(cur_genomes, rep_gids, rep_mash_sketch_file)

        # cluster remaining genomes to representatives
        non_reps = set(cur_genomes.genomic_files) - set(rep_radius)
        self.logger.info('Clustering {:,} non-representatives to {:,} representatives using species specific ANI radii.'.format(len(non_reps), len(rep_radius)))
        if True: #***
            clusters = self._cluster(ani_af, non_reps, rep_radius)
        else:
            self.logger.warning('Using previously calculated results in: {}'.format('clusters.pkl'))
            clusters = pickle.load(open(os.path.join(self.output_dir, 'clusters.pkl'), 'rb'))
        
        # write out clusters
        write_clusters(clusters, 
                        rep_radius,
                        proposed_sp,
                        cur_genomes,
                        os.path.join(self.output_dir, 'gtdb_rep_clusters.tsv'))
