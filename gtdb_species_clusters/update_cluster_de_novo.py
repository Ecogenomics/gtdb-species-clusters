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
from itertools import combinations
from collections import defaultdict, namedtuple

from biolib.external.execute import check_dependencies

from numpy import (mean as np_mean,
                    std as np_std)

from gtdb_species_clusters.mash import Mash
from gtdb_species_clusters.fastani import FastANI
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.type_genome_utils import (ClusteredGenome,
                                                        GenomeRadius,
                                                        symmetric_ani,
                                                        write_rep_radius,
                                                        write_clusters)
                                                        


class UpdateClusterDeNovo(object):
    """Infer de novo species clusters and representatives for remaining genomes."""

    def __init__(self, ani_sp, af_sp, ani_cache_file, cpus, output_dir):
        """Initialization."""
        
        check_dependencies(['fastANI', 'mash'])
        
        self.cpus = cpus
        self.output_dir = output_dir

        self.logger = logging.getLogger('timestamp')
        
        self.true_str = ['t', 'T', 'true', 'True']
        
        self.ani_sp = ani_sp
        self.af_sp = af_sp

        self.min_mash_ani = 90.0

        self.fastani = FastANI(ani_cache_file, cpus)
        
    def _parse_named_clusters(self, named_cluster_file):
        """Parse named GTDB species clusters."""
        
        rep_gids = set()
        rep_clustered_gids = set()
        rep_radius = {}
        with open(named_cluster_file) as f:
            headers = f.readline().strip().split('\t')
            
            rep_index = headers.index('Representative')
            num_clustered_index = headers.index('No. clustered genomes')
            clustered_genomes_index = headers.index('Clustered genomes')
            closest_type_index = headers.index('Closest representative')
            ani_radius_index = headers.index('ANI radius')
            af_index = headers.index('AF closest')

            for line in f:
                line_split = line.strip().split('\t')

                rep_gid = line_split[rep_index]
                rep_gids.add(rep_gid)
                
                num_clustered = int(line_split[num_clustered_index])
                if num_clustered > 0:
                    for gid in [g.strip() for g in line_split[clustered_genomes_index].split(',')]:
                        rep_clustered_gids.add(gid)
                        
                rep_radius[rep_gid] = GenomeRadius(ani = float(line_split[ani_radius_index]), 
                                                     af = float(line_split[af_index]),
                                                     neighbour_gid = line_split[closest_type_index])
                        
        return rep_gids, rep_clustered_gids, rep_radius

    def _nonrep_radius(self, unclustered_gids, rep_gids, ani_af_rep_vs_nonrep):
        """Calculate circumscription radius for unclustered, nontype genomes."""
        
        # set radius for genomes to default values
        nonrep_radius = {}
        for gid in unclustered_gids:
            nonrep_radius[gid] = GenomeRadius(ani = self.ani_sp, 
                                                     af = None,
                                                     neighbour_gid = None)

        # determine closest type ANI neighbour and restrict ANI radius as necessary
        ani_af = pickle.load(open(ani_af_rep_vs_nonrep, 'rb'))
        for nonrep_gid in unclustered_gids:
            if nonrep_gid not in ani_af:
                continue
                    
            for rep_gid in rep_gids:
                if rep_gid not in ani_af[nonrep_gid]:
                    continue
                    
                ani, af = symmetric_ani(ani_af, nonrep_gid, rep_gid)

                if ani > nonrep_radius[nonrep_gid].ani and af >= self.af_sp:
                    nonrep_radius[nonrep_gid] = GenomeRadius(ani = ani, 
                                                             af = af,
                                                             neighbour_gid = rep_gid)
                    
        self.logger.info('ANI circumscription radius: min={:.2f}, mean={:.2f}, max={:.2f}'.format(
                                min([d.ani for d in nonrep_radius.values()]), 
                                np_mean([d.ani for d in nonrep_radius.values()]), 
                                max([d.ani for d in nonrep_radius.values()])))
                        
        return nonrep_radius
        
    def _mash_ani_unclustered(self, cur_genomes, gids):
        """Calculate pairwise Mash ANI estimates between genomes."""
        
        mash = Mash(self.cpus)
        
        # create Mash sketch for potential representative genomes
        mash_nontype_sketch_file = os.path.join(self.output_dir, 'gtdb_unclustered_genomes.msh')
        genome_list_file = os.path.join(self.output_dir, 'gtdb_unclustered_genomes.lst')
        mash.sketch(gids, cur_genomes.genomic_files, genome_list_file, mash_nontype_sketch_file)

        # get Mash distances
        mash_dist_file = os.path.join(self.output_dir, 'gtdb_unclustered_genomes.dst')
        mash.dist_pairwise( float(100 - self.min_mash_ani)/100, mash_nontype_sketch_file, mash_dist_file)

        # read Mash distances
        mash_ani = mash.read_ani(mash_dist_file)
        
        # report pairs above Mash threshold
        mash_ani_pairs = []
        for qid in mash_ani:
            for rid in mash_ani[qid]:
                if mash_ani[qid][rid] >= self.min_mash_ani:
                    n_qid = cur_genomes.user_uba_id_map.get(qid, qid)
                    n_rid = cur_genomes.user_uba_id_map.get(rid, rid)
                    if n_qid != n_rid:
                        mash_ani_pairs.append((n_qid, n_rid))
                        mash_ani_pairs.append((n_rid, n_qid))
                
        self.logger.info('Identified {:,} genome pairs with a Mash ANI >= {:.1f}%.'.format(
                            len(mash_ani_pairs), 
                            self.min_mash_ani))

        return mash_ani
        
    def _selected_rep_genomes(self,
                                cur_genomes,
                                nonrep_radius, 
                                unclustered_qc_gids, 
                                mash_ani):
        """Select de novo representatives for species clusters in a greedy fashion using species-specific ANI thresholds."""

        # sort genomes by quality score
        self.logger.info('Selecting de novo representatives in a greedy manner based on quality.')
        q = {gid:cur_genomes[gid].score_type_strain() for gid in unclustered_qc_gids}
        q_sorted = sorted(q.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)

        # greedily determine representatives for new species clusters
        cluster_rep_file = os.path.join(self.output_dir, 'cluster_reps.tsv')
        clusters = set()
        if not os.path.exists(cluster_rep_file):
            clustered_genomes = 0
            max_ani_pairs = 0
            for idx, (cur_gid, _score) in enumerate(q_sorted):

                # determine reference genomes to calculate ANI between
                ani_pairs = []
                if cur_gid in mash_ani:
                    for rep_gid in clusters:
                        if mash_ani[cur_gid].get(rep_gid, 0) >= self.min_mash_ani:
                            ani_pairs.append((cur_gid, rep_gid))
                            ani_pairs.append((rep_gid, cur_gid))

                # determine if genome clusters with representative
                clustered = False
                if ani_pairs:
                    if len(ani_pairs) > max_ani_pairs:
                        max_ani_pairs = len(ani_pairs)
                    
                    ani_af = self.fastani.pairs(ani_pairs, cur_genomes.genomic_files, report_progress=False)

                    closest_rep_gid = None
                    closest_rep_ani = 0
                    closest_rep_af = 0
                    for rep_gid in clusters:
                        ani, af = symmetric_ani(ani_af, cur_gid, rep_gid)

                        if af >= self.af_sp:
                            if ani > closest_rep_ani or (ani == closest_rep_ani and af > closest_rep_af):
                                closest_rep_gid = rep_gid
                                closest_rep_ani = ani
                                closest_rep_af = af

                        if ani > nonrep_radius[cur_gid].ani and af >= self.af_sp:
                            nonrep_radius[cur_gid] = GenomeRadius(ani = ani, 
                                                                         af = af,
                                                                         neighbour_gid = rep_gid)
                                                                         
                    if closest_rep_gid and closest_rep_ani > nonrep_radius[closest_rep_gid].ani:
                        clustered = True
                    
                if not clustered:
                    # genome is a new species cluster representative
                    clusters.add(cur_gid)
                else:
                    clustered_genomes += 1
                
                if (idx+1) % 10 == 0 or idx+1 == len(q_sorted):
                    statusStr = '-> Clustered {:,} of {:,} ({:.2f}%) genomes [ANI pairs: {:,}; clustered genomes: {:,}; clusters: {:,}].'.format(
                                    idx+1, 
                                    len(q_sorted), 
                                    float(idx+1)*100/len(q_sorted),
                                    max_ani_pairs,
                                    clustered_genomes,
                                    len(clusters)).ljust(96)
                    sys.stdout.write('{}\r'.format(statusStr))
                    sys.stdout.flush()
                    max_ani_pairs = 0
            sys.stdout.write('\n')
            
            # write out selected cluster representative
            fout = open(cluster_rep_file, 'w')
            for gid in clusters:
                fout.write('{}\n'.format(gid))
            fout.close()
        else:
            # read cluster reps from file
            self.logger.warning('Using previously determined cluster representatives.')
            for line in open(cluster_rep_file):
                gid = line.strip()
                clusters.add(gid)
                
        self.logger.info('Selected {:,} representative genomes for de novo species clusters.'.format(len(clusters)))
        
        return clusters
        
    def _cluster_genomes(self,
                            cur_genomes,
                            de_novo_rep_gids,
                            named_rep_gids, 
                            final_cluster_radius):
        """Cluster new representatives to representatives of named GTDB species clusters."""
        
        all_reps = de_novo_rep_gids.union(named_rep_gids)
        nonrep_gids = set(cur_genomes.genomes.keys()) - all_reps
        self.logger.info('Clustering {:,} genomes to {:,} named and de novo representatives.'.format(
                            len(nonrep_gids), len(all_reps)))

        if True: #***
            # calculate MASH distance between non-representatives and representatives genomes
            mash = Mash(self.cpus)
            
            mash_rep_sketch_file = os.path.join(self.output_dir, 'gtdb_rep_genomes.msh')
            rep_genome_list_file = os.path.join(self.output_dir, 'gtdb_rep_genomes.lst')
            mash.sketch(all_reps, cur_genomes.genomic_files, rep_genome_list_file, mash_rep_sketch_file)

            mash_none_rep_sketch_file = os.path.join(self.output_dir, 'gtdb_nonrep_genomes.msh')
            non_rep_file = os.path.join(self.output_dir, 'gtdb_nonrep_genomes.lst')
            mash.sketch(nonrep_gids, cur_genomes.genomic_files, non_rep_file, mash_none_rep_sketch_file)

            # get Mash distances
            mash_dist_file = os.path.join(self.output_dir, 'gtdb_rep_vs_nonrep_genomes.dst')
            mash.dist(float(100 - self.min_mash_ani)/100, 
                        mash_rep_sketch_file, 
                        mash_none_rep_sketch_file, 
                        mash_dist_file)

            # read Mash distances
            mash_ani = mash.read_ani(mash_dist_file)
            
            # calculate ANI between non-representatives and representatives genomes
            clusters = {}
            for gid in all_reps:
                clusters[gid] = []

            if False: #***
                mash_ani_pairs = []
                for gid in nonrep_gids:
                    if gid in mash_ani:
                        for rid in clusters:
                            if mash_ani[gid].get(rid, 0) >= self.min_mash_ani:
                                n_gid = cur_genomes.user_uba_id_map.get(gid, gid)
                                n_rid = cur_genomes.user_uba_id_map.get(rid, rid)
                                if n_gid != n_rid:
                                    mash_ani_pairs.append((n_gid, n_rid))
                                    mash_ani_pairs.append((n_rid, n_gid))
                                    
            mash_ani_pairs = []
            for qid in mash_ani:
                n_qid = cur_genomes.user_uba_id_map.get(qid, qid)
                assert n_qid in nonrep_gids
                
                for rid in mash_ani[qid]:
                    n_rid = cur_genomes.user_uba_id_map.get(rid, rid)
                    assert n_rid in all_reps
                    
                    if (mash_ani[qid][rid] >= self.min_mash_ani
                        and n_qid != n_rid):
                        mash_ani_pairs.append((n_qid, n_rid))
                        mash_ani_pairs.append((n_rid, n_qid))
                            
            self.logger.info('Calculating ANI between {:,} species clusters and {:,} unclustered genomes ({:,} pairs):'.format(
                                len(clusters), 
                                len(nonrep_gids),
                                len(mash_ani_pairs)))
            ani_af = self.fastani.pairs(mash_ani_pairs, cur_genomes.genomic_files)

            # assign genomes to closest representatives 
            # that is within the representatives ANI radius
            self.logger.info('Assigning genomes to closest representative.')
            for idx, cur_gid in enumerate(nonrep_gids):
                closest_rep_gid = None
                closest_rep_ani = 0
                closest_rep_af = 0
                for rep_gid in clusters:
                    ani, af = symmetric_ani(ani_af, cur_gid, rep_gid)
                    
                    if ani >= final_cluster_radius[rep_gid].ani and af >= self.af_sp:
                        if ani > closest_rep_ani or (ani == closest_rep_ani and af > closest_rep_af):
                            closest_rep_gid = rep_gid
                            closest_rep_ani = ani
                            closest_rep_af = af
                    
                if closest_rep_gid:
                    clusters[closest_rep_gid].append(ClusteredGenome(gid=cur_gid, 
                                                                            ani=closest_rep_ani, 
                                                                            af=closest_rep_af))
                else:
                    self.logger.warning('Failed to assign genome {} to representative.'.format(cur_gid))
                    if closest_rep_gid:
                        self.logger.warning(' ...closest_rep_gid = {}'.format(closest_rep_gid))
                        self.logger.warning(' ...closest_rep_ani = {:.2f}'.format(closest_rep_ani))
                        self.logger.warning(' ...closest_rep_af = {:.2f}'.format(closest_rep_af))
                        self.logger.warning(' ...closest rep radius = {:.2f}'.format(final_cluster_radius[closest_rep_gid].ani))
                    else:
                        self.logger.warning(' ...no representative with an AF >{:.2f} identified.'.format(self.af_sp))
                 
                statusStr = '-> Assigned {:,} of {:,} ({:.2f}%) genomes.'.format(idx+1, 
                                                                                    len(nonrep_gids), 
                                                                                    float(idx+1)*100/len(nonrep_gids)).ljust(86)
                sys.stdout.write('{}\r'.format(statusStr))
                sys.stdout.flush()
            sys.stdout.write('\n')
            
            pickle.dump(clusters, open(os.path.join(self.output_dir, 'clusters.pkl'), 'wb'))
            pickle.dump(ani_af, open(os.path.join(self.output_dir, 'ani_af_rep_vs_nonrep.de_novo.pkl'), 'wb'))
        else:
            self.logger.warning('Using previously calculated results in: {}'.format('clusters.pkl'))
            clusters = pickle.load(open(os.path.join(self.output_dir, 'clusters.pkl'), 'rb'))
            
            self.logger.warning('Using previously calculated results in: {}'.format('ani_af_rep_vs_nonrep.de_novo.pkl'))
            ani_af = pickle.load(open(os.path.join(self.output_dir, 'ani_af_rep_vs_nonrep.de_novo.pkl'), 'rb'))

        return clusters, ani_af

    def run(self, named_cluster_file,
                    cur_gtdb_metadata_file,
                    cur_genomic_path_file,
                    uba_genome_paths,
                    qc_passed_file,
                    gtdbtk_classify_file,
                    ncbi_genbank_assembly_file,
                    untrustworthy_type_file,
                    ani_af_rep_vs_nonrep,
                    gtdb_type_strains_ledger):
        """Infer de novo species clusters and representatives for remaining genomes."""
        
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

        # determine representatives and genomes clustered to each representative
        self.logger.info('Reading named GTDB species clusters.')
        named_rep_gids, rep_clustered_gids, rep_radius = self._parse_named_clusters(named_cluster_file)
        self.logger.info(' ... identified {:,} representative genomes.'.format(len(named_rep_gids)))
        self.logger.info(' ... identified {:,} clustered genomes.'.format(len(rep_clustered_gids)))
        
        # determine genomes left to be clustered
        unclustered_gids = set(cur_genomes.genomes.keys()) - named_rep_gids - rep_clustered_gids
        self.logger.info('Identified {:,} unclustered genomes passing QC.'.format(len(unclustered_gids)))

        # establish closest representative for each unclustered genome
        self.logger.info('Determining ANI circumscription for {:,} unclustered genomes.'.format(len(unclustered_gids)))
        nonrep_radius = self._nonrep_radius(unclustered_gids, named_rep_gids, ani_af_rep_vs_nonrep)

        # calculate Mash ANI estimates between unclustered genomes
        self.logger.info('Calculating Mash ANI estimates between unclustered genomes.')
        mash_anis = self._mash_ani_unclustered(cur_genomes, unclustered_gids)

        # select de novo species representatives in a greedy fashion based on genome quality
        de_novo_rep_gids = self._selected_rep_genomes(cur_genomes,
                                                        nonrep_radius, 
                                                        unclustered_gids, 
                                                        mash_anis)

        # cluster all non-representative genomes to representative genomes
        final_cluster_radius = rep_radius.copy()
        final_cluster_radius.update(nonrep_radius)
        
        final_clusters, ani_af = self._cluster_genomes(cur_genomes,
                                                        de_novo_rep_gids,
                                                        named_rep_gids, 
                                                        final_cluster_radius)

        # remove genomes that are not representatives of a species cluster and then write out representative ANI radius
        for gid in set(final_cluster_radius) - set(final_clusters):
            del final_cluster_radius[gid]

        self.logger.info('Writing {:,} species clusters to file.'.format(len(final_clusters)))
        self.logger.info('Writing {:,} cluster radius information to file.'.format(len(final_cluster_radius)))
        
        write_clusters(final_clusters, 
                        final_cluster_radius, 
                        cur_genomes,
                        os.path.join(self.output_dir, 'gtdb_clusters_de_novo.tsv'))

        write_rep_radius(final_cluster_radius, 
                            cur_genomes,
                            os.path.join(self.output_dir, 'gtdb_ani_radius_de_novo.tsv'))
