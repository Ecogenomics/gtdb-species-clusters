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
        self.max_af_neighbour = 0.65
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
        af_warning_count = 0
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
                
                if ani >= self.max_ani_neighbour and af >= self.max_af_neighbour:
                    # typically, representative genomes should not exceed this ANI and AF
                    # criteria as they should have been declared synonyms in 
                    # the u_sel_reps step if they are this similar to each other.
                    # However, a 'fudge factor' is used to allow previous GTDB clusters
                    # to remain as seperate clusters if they exceed these thresholds by
                    # a small margin as this can simply be due to differences in the 
                    # version of FastANI used to calculate ANI and AF.
                    self.logger.warning('ANI neighbours {} and {} have ANI={:.2f} and AF={:.2f}.'.format(
                                        rep_gid1, rep_gid2,
                                        ani, af))

                if ani > rep_radius[rep_gid1].ani:
                    if af < self.af_sp:
                        af_warning_count += 1
                        #self.logger.warning('ANI for {} and {} is >{:.2f}, but AF <{:.2f} [pair skipped].'.format(
                        #                        rep_gid1,
                        #                        rep_gid2,
                        #                        ani, af))
                        continue

                    rep_radius[rep_gid1] = GenomeRadius(ani = ani, 
                                                         af = af,
                                                         neighbour_gid = rep_gid2)
                    
        self.logger.info('ANI circumscription radius: min={:.2f}, mean={:.2f}, max={:.2f}'.format(
                                min([d.ani for d in rep_radius.values()]), 
                                np_mean([d.ani for d in rep_radius.values()]), 
                                max([d.ani for d in rep_radius.values()])))
        
        self.logger.warning('Identified {:,} genome pairs meeting ANI radius criteria, but with an AF <{:.2f}'.format(
                                af_warning_count,
                                self.af_sp))
                                
        
        return rep_radius
        
    def _calculate_ani(self, cur_genomes, rep_gids, rep_mash_sketch_file):
        """Calculate ANI between representative and non-representative genomes."""
        
        if True: #***
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
                        n_qid = cur_genomes.user_uba_id_map.get(qid, qid)
                        n_rid = cur_genomes.user_uba_id_map.get(rid, rid)
                        if n_qid != n_rid:
                            mash_ani_pairs.append((n_qid, n_rid))
                            mash_ani_pairs.append((n_rid, n_qid))
                    
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
            
        num_clustered = 0
        for idx, non_rid in enumerate(non_reps):
            if idx % 100 == 0:
                sys.stdout.write('==> Processed {:,} of {:,} genomes [no. clustered = {:,}].\r'.format(
                                    idx+1, 
                                    len(non_reps),
                                    num_clustered))
                sys.stdout.flush()

            if non_rid not in ani_af:
                continue

            closest_rid = None
            closest_ani = 0
            closest_af = 0
            for rid in rep_radius:
                if rid not in ani_af[non_rid]:
                    continue

                ani, af = symmetric_ani(ani_af, rid, non_rid)
                
                if af >= self.af_sp:
                    if ani > closest_ani or (ani == closest_ani and af > closest_af):
                        closest_rid = rid
                        closest_ani = ani
                        closest_af = af
                
            if closest_rid:
                if closest_ani > rep_radius[closest_rid].ani:
                    num_clustered += 1
                    clusters[closest_rid].append(self.ClusteredGenome(gid=non_rid, 
                                                                        ani=closest_ani, 
                                                                        af=closest_af))
                
        sys.stdout.write('==> Processed {:,} of {:,} genomes [no. clustered = {:,}].\r'.format(
                                    idx+1, 
                                    len(non_reps),
                                    num_clustered))
        sys.stdout.flush()
        sys.stdout.write('\n')

        num_unclustered = len(non_reps) - num_clustered
        self.logger.info('Assigned {:,} genomes to {:,} representatives; {:,} genomes remain unclustered.'.format(
                            sum([len(clusters[rid]) for rid in clusters]),
                            len(clusters),
                            num_unclustered))
        
        return clusters

    def run(self, named_rep_file,
                    cur_gtdb_metadata_file,
                    cur_genomic_path_file,
                    uba_genome_paths,
                    qc_passed_file,
                    gtdbtk_classify_file,
                    ncbi_genbank_assembly_file,
                    untrustworthy_type_file,
                    rep_mash_sketch_file,
                    rep_ani_file,
                    gtdb_type_strains_ledger):
        """Cluster genomes to selected GTDB representatives."""
        
        # create current GTDB genome sets
        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                create_sp_clusters=False,
                                                uba_genome_file=uba_genome_paths,
                                                qc_passed_file=qc_passed_file,
                                                gtdbtk_classify_file=gtdbtk_classify_file,
                                                ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                                untrustworthy_type_ledger=untrustworthy_type_file)
        self.logger.info(f' ... current genome set contains {len(cur_genomes):,} genomes.')
        
        # get path to previous and current genomic FASTA files
        self.logger.info('Reading path to current genomic FASTA files.')
        cur_genomes.load_genomic_file_paths(cur_genomic_path_file)
        cur_genomes.load_genomic_file_paths(uba_genome_paths)

        # get representative genomes
        rep_gids = set()
        with open(named_rep_file) as f:
            header = f.readline().strip().split('\t')
            rep_index = header.index('Representative')
            sp_index = header.index('Proposed species')

            for line in f:
                line_split = line.strip().split('\t')
                gid = line_split[rep_index]
                assert gid in cur_genomes
                rep_gids.add(gid)
                
        self.logger.info('Identified representative genomes for {:,} species.'.format(len(rep_gids)))

        # calculate circumscription radius for representative genomes
        self.logger.info('Determining ANI species circumscription for {:,} representative genomes.'.format(len(rep_gids)))
        rep_radius = self._rep_radius(rep_gids, rep_ani_file)
        write_rep_radius(rep_radius, cur_genomes, os.path.join(self.output_dir, 'gtdb_rep_ani_radius.tsv'))

        # calculate ANI between representative and non-representative genomes
        self.logger.info('Calculating ANI between representative and non-representative genomes.')
        ani_af = self._calculate_ani(cur_genomes, rep_gids, rep_mash_sketch_file)
        self.logger.info(' ... ANI values determined for {:,} query genomes.'.format(len(ani_af)))
        self.logger.info(' ... ANI values determined for {:,} genome pairs.'.format(
                            sum([len(ani_af[qid]) for qid in ani_af])))

        # cluster remaining genomes to representatives
        non_reps = set(cur_genomes.genomes) - set(rep_radius)
        self.logger.info('Clustering {:,} non-representatives to {:,} representatives using species-specific ANI radii.'.format(len(non_reps), len(rep_radius)))
        clusters = self._cluster(ani_af, non_reps, rep_radius)
        
        # write out clusters
        write_clusters(clusters, 
                        rep_radius,
                        cur_genomes,
                        os.path.join(self.output_dir, 'gtdb_named_rep_clusters.tsv'))
