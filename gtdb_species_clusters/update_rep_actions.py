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
import argparse
import logging
from collections import defaultdict

from gtdb_species_clusters.fastani import FastANI
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_clusters import SpeciesClusters
from gtdb_species_clusters.type_genome_utils import symmetric_ani


class RepActions(object):
    """Perform initial actions required for changed representatives."""

    def __init__(self, ani_cache_file, cpus, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
        self.fastani = FastANI(ani_cache_file, cpus)
        
        # action parameters
        self.genomic_update_ani = 99.0
        self.genomic_update_af = 0.80
        
        self.new_rep_ani = 99.0
        self.new_rep_af = 0.80
        self.new_rep_qs_threshold = 10      # increase in ANI score require to select 
                                            # new representative
        
        self.action_log = open(os.path.join(self.output_dir, 'action_log.tsv'), 'w')
        self.action_log.write('Genome ID\tPrevious GTDB species\tAction\tParameters\n')
        
    def rep_change_gids(self, rep_change_summary_file, field, value):
        """Get genomes with a specific change."""
        
        gids = {}
        with open(rep_change_summary_file) as f:
            header = f.readline().strip().split('\t')
            
            field_index = header.index(field)
            prev_sp_index = header.index('Previous GTDB species')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                v = line_split[field_index]
                if v == value:
                    prev_sp = line_split[prev_sp_index]
                    gids[line_split[0]] = prev_sp
                    
        return gids
        
    def top_ani_score_prev_rep(self,
                                prev_rid,
                                sp_cids, 
                                prev_genomes,
                                cur_genomes):
        """Identify genome in cluster with highest balanced ANI score to genomic file of representative in previous GTDB release."""

        max_score = -1e6
        max_rid = None
        max_ani = None
        max_af = None
        for cid in sp_cids:
            ani, af = self.fastani.symmetric_ani_cached(f'{prev_rid}-P', f'{cid}-C', 
                                                        prev_genomes[prev_rid].genomic_file, 
                                                        cur_genomes[cid].genomic_file)

            cur_score = cur_genomes[cid].score_ani(ani)
            if (cur_score > max_score
                or (cur_score == max_score and ani > max_ani)):
                max_score = cur_score
                max_rid = cid
                max_ani = ani
                max_af = af
                        
        return max_rid, max_score, max_ani, max_af
        
    def top_ani_score(self,
                        prev_rid,
                        sp_cids, 
                        cur_genomes):
        """Identify genome in cluster with highest balanced ANI score to representative genome."""
        
        # calculate ANI between representative and genomes in species cluster
        gid_pairs = []
        for cid in sp_cids:
            gid_pairs.append((cid, prev_rid))
            gid_pairs.append((prev_rid, cid))

        ani_af = self.fastani.pairs(gid_pairs, cur_genomes.genomic_files, report_progress=False)

        # find genome with top ANI score
        max_score = -1e6
        max_rid = None
        max_ani = None
        max_af = None
        for cid in sp_cids:
            ani, af = symmetric_ani(ani_af, prev_rid, cid)

            if cur_genomes[cid].score_ani(ani) > max_score:
                max_score = cur_genomes[cid].score_ani(ani)
                max_rid = cid
                max_ani = ani
                max_af = af
                        
        return max_rid, max_score, max_ani, max_af

    def genomes_in_current_sp_cluster(self, 
                                        prev_rid, 
                                        prev_genomes, 
                                        new_updated_sp_clusters, 
                                        cur_genomes):
        """Get genomes in current species cluster."""

        sp_cids = prev_genomes.sp_clusters[prev_rid]
        if prev_rid in new_updated_sp_clusters:
            sp_cids = sp_cids.union(new_updated_sp_clusters[prev_rid])
        sp_cids = sp_cids.intersection(cur_genomes)

        return sp_cids
        
    def action_genomic_lost(self, 
                                rep_change_summary_file,
                                prev_genomes, 
                                cur_genomes,
                                new_updated_sp_clusters,
                                new_reps):
        """Handle species with lost representative genome."""
        
        # get genomes with specific changes
        self.logger.info('Identifying species with lost representative genome.')
        genomic_lost_rids = self.rep_change_gids(rep_change_summary_file, 
                                                    'GENOMIC_CHANGE', 
                                                    'LOST')
        self.logger.info(f' ...identified {len(genomic_lost_rids):,} genomes.')
        
        # calculate ANI between previous and current genomes
        for prev_rid, prev_gtdb_sp in genomic_lost_rids.items():
            sp_cids = self.genomes_in_current_sp_cluster(prev_rid, 
                                                            prev_genomes, 
                                                            new_updated_sp_clusters, 
                                                            cur_genomes)

            params = {}
            if sp_cids:
                action = 'GENOMIC_CHANGE:LOST:REPLACED'
                
                new_rid, top_score, ani, af = self.top_ani_score_prev_rep(prev_rid,
                                                                            sp_cids,
                                                                            prev_genomes,
                                                                            cur_genomes)
                assert(new_rid != prev_rid)

                params['new_rid'] = new_rid
                params['ani'] = ani
                params['af'] = af
                params['new_rid_assembly_quality'] = cur_genomes[new_rid].score_assembly
                params['prev_rid_assembly_quality'] = prev_genomes[prev_rid].score_assembly
                
                new_reps[prev_gtdb_sp] = new_rid
            else:
                action = 'GENOMIC_CHANGE:LOST:SPECIES_RETIRED'
                new_reps[prev_gtdb_sp] = None
                
            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                        prev_rid, 
                                        prev_gtdb_sp, 
                                        action, 
                                        params))
        
    def action_genomic_update(self, 
                                rep_change_summary_file,
                                prev_genomes, 
                                cur_genomes,
                                new_updated_sp_clusters,
                                new_reps):
        """Handle representatives with updated genomes."""
        
        # get genomes with specific changes
        self.logger.info('Identifying representatives with updated genomic files.')
        genomic_update_gids = self.rep_change_gids(rep_change_summary_file, 
                                                    'GENOMIC_CHANGE', 
                                                    'UPDATED')
        self.logger.info(f' ...identified {len(genomic_update_gids):,} genomes.')
        
        # calculate ANI between previous and current genomes
        for prev_rid, prev_gtdb_sp in genomic_update_gids.items():
            if prev_rid not in cur_genomes:
                # indicates genome has been lost
                continue
                
            ani, af = self.fastani.symmetric_ani_cached(f'{prev_rid}-P', f'{prev_rid}-C', 
                                                        prev_genomes[prev_rid].genomic_file, 
                                                        cur_genomes[prev_rid].genomic_file)

            params = {}
            params['ani'] = ani
            params['af'] = af
            params['prev_ncbi_accession'] = prev_genomes[prev_rid].ncbi_accn
            params['cur_ncbi_accession'] = cur_genomes[prev_rid].ncbi_accn
            
            if ani >= self.genomic_update_ani and af >= self.genomic_update_af:
                action = 'GENOMIC_CHANGE:UPDATED:MINOR_CHANGE'
            else:
                sp_cids = self.genomes_in_current_sp_cluster(prev_rid, 
                                                            prev_genomes, 
                                                            new_updated_sp_clusters, 
                                                            cur_genomes)

                if sp_cids:
                    new_rid, top_score, ani, af = self.top_ani_score_prev_rep(prev_rid,
                                                                                sp_cids,
                                                                                prev_genomes,
                                                                                cur_genomes)
                    
                    if new_rid == prev_rid:
                        action = 'GENOMIC_CHANGE:UPDATED:RETAINED'
                    else:
                        action = 'GENOMIC_CHANGE:UPDATED:REPLACED'
                        params['new_rid'] = new_rid
                        params['ani'] = ani
                        params['af'] = af
                        params['new_rid_assembly_quality'] = cur_genomes[new_rid].score_assembly
                        params['prev_rid_assembly_quality'] = prev_genomes[prev_rid].score_assembly
                        
                        new_reps[prev_gtdb_sp] = new_rid
                else:
                    action = 'GENOMIC_CHANGE:UPDATED:SPECIES_RETIRED'
                    new_reps[prev_gtdb_sp] = None
                
            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                        prev_rid, 
                                        prev_gtdb_sp, 
                                        action, 
                                        params))
                                        
    def action_species_reassigned(self,
                                    rep_change_summary_file,
                                    prev_genomes, 
                                    cur_genomes):
        """Handle representatives with new NCBI species assignments."""
        
        # get genomes with new NCBI species assignments
        self.logger.info('Identifying representatives with new NCBI species assignments.')
        ncbi_sp_change_gids = self.rep_change_gids(rep_change_summary_file, 
                                                    'NCBI_SPECIES_CHANGE', 
                                                    'REASSIGNED')
        self.logger.info(f' ...identified {len(ncbi_sp_change_gids):,} genomes.')
        
        # get species assignment of previous GTDB species clusters
        prev_sp_clusters = prev_genomes.sp_clusters
        gtdb_specific_epithets = defaultdict(set)
        specific_epithets_rid = defaultdict(lambda: {})
        for rid, sp in prev_sp_clusters.species():
            gtdb_genus = prev_genomes[rid].gtdb_genus
            gtdb_specific_epithet = prev_genomes[rid].gtdb_specific_epithet
            gtdb_specific_epithets[gtdb_genus].add(gtdb_specific_epithet)
            specific_epithets_rid[gtdb_genus][gtdb_specific_epithet] = rid

        # determine if change can be made without causing a conflict with 
        # other GTDB species clusters
        for prev_rid, prev_gtdb_sp in ncbi_sp_change_gids.items():
            prev_ncbi_sp = prev_genomes[prev_rid].ncbi_species
            cur_ncbi_sp = cur_genomes[prev_rid].ncbi_species
            prev_gtdb_genus = prev_genomes[prev_rid].gtdb_genus
            assert(prev_ncbi_sp != cur_ncbi_sp)
            assert(prev_gtdb_sp == prev_genomes[prev_rid].gtdb_species)
            
            params = {}
            params['prev_ncbi_sp'] = prev_ncbi_sp
            params['cur_ncbi_sp'] = cur_ncbi_sp
            params['prev_gtdb_sp'] = prev_gtdb_sp
            
            cur_ncbi_specific_epithet = cur_genomes[prev_rid].ncbi_specific_epithet
            prev_gtdb_specific_epithet = prev_genomes[prev_rid].gtdb_specific_epithet
            if cur_ncbi_specific_epithet == prev_gtdb_specific_epithet:
                action = 'NCBI_SPECIES_CHANGE:REASSIGNED:UNCHANGED'
            elif cur_ncbi_specific_epithet not in gtdb_specific_epithets[prev_gtdb_genus]:
                action = 'NCBI_SPECIES_CHANGE:REASSIGNED:UPDATED'
            else:
                # check if species cluster should be merged
                ncbi_proposed_rid = specific_epithets_rid[prev_gtdb_genus][cur_ncbi_specific_epithet]
                ani, af = self.fastani.symmetric_ani_cached(prev_rid, ncbi_proposed_rid,
                                                            cur_genomes[prev_rid].genomic_file, 
                                                            cur_genomes[ncbi_proposed_rid].genomic_file)
                                                        
                params['ani'] = ani
                params['af'] = af
                params['ncbi_subspecies'] = cur_genomes[prev_rid].ncbi_subspecies
                                                        
                if ani >= 95 and af >= 0.65:
                    action = 'NCBI_SPECIES_CHANGE:REASSIGNED:MERGE'
                else:
                    action = 'NCBI_SPECIES_CHANGE:REASSIGNED:CONFLICT'
                
            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                        prev_rid, 
                                        prev_gtdb_sp, 
                                        action, 
                                        params))
                                        
    def action_type_strain_lost(self, 
                                    rep_change_summary_file,
                                    prev_genomes, 
                                    cur_genomes,
                                    new_updated_sp_clusters,
                                    new_reps):
        """Handle representatives which have lost type strain genome status."""
        
        # get genomes with new NCBI species assignments
        self.logger.info('Identifying representative that lost type strain genome status.')
        ncbi_type_species_lost = self.rep_change_gids(rep_change_summary_file, 
                                                    'TYPE_STRAIN_CHANGE', 
                                                    'LOST')
        self.logger.info(f' ...identified {len(ncbi_type_species_lost):,} genomes.')

        for prev_rid, prev_gtdb_sp in ncbi_type_species_lost.items():
            if prev_rid not in cur_genomes:
                # indicates genome has been lost
                continue
                
            sp_cids = self.genomes_in_current_sp_cluster(prev_rid, 
                                                            prev_genomes, 
                                                            new_updated_sp_clusters, 
                                                            cur_genomes)
                
            prev_rep_score = cur_genomes[prev_rid].score_ani(100)
            new_rid, top_score, ani, af = self.top_ani_score(prev_rid,
                                                                sp_cids,
                                                                cur_genomes)

            params = {}
            params['prev_rid_prev_strain_ids'] = prev_genomes[prev_rid].ncbi_strain_identifiers
            params['prev_rid_cur_strain_ids'] = cur_genomes[prev_rid].ncbi_strain_identifiers
            params['prev_rid_prev_gtdb_type_designation'] = prev_genomes[prev_rid].gtdb_type_designation
            params['prev_rid_cur_gtdb_type_designation'] = cur_genomes[prev_rid].gtdb_type_designation
            params['prev_rid_prev_gtdb_type_designation_sources'] = prev_genomes[prev_rid].gtdb_type_designation_sources
            params['prev_rid_cur_gtdb_type_designation_sources'] = cur_genomes[prev_rid].gtdb_type_designation_sources

            if top_score > prev_rep_score:
                action = 'TYPE_STRAIN_CHANGE:LOST:REPLACED'
                assert(prev_rid != new_rid)
                
                params['new_rid'] = new_rid
                params['ani'] = ani
                params['af'] = af
                params['new_rid_assembly_quality'] = cur_genomes[new_rid].score_assembly
                params['prev_rid_assembly_quality'] = prev_genomes[prev_rid].score_assembly

                params['new_rid_strain_ids'] = prev_genomes[new_rid].ncbi_strain_identifiers
                params['new_rid_gtdb_type_designation'] = prev_genomes[new_rid].gtdb_type_designation
                params['new_rid_gtdb_type_designation_sources'] = prev_genomes[new_rid].gtdb_type_designation_sources
                
                new_reps[prev_gtdb_sp] = new_rid
            else:
                action = 'TYPE_STRAIN_CHANGE:LOST:RETAINED'
                
            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                            prev_rid, 
                                            prev_gtdb_sp, 
                                            action, 
                                            params))
                                            
    def action_improved_rep(self, 
                            prev_genomes, 
                            cur_genomes, 
                            new_updated_sp_clusters,
                            new_reps):
        """Check if representative should be replace with higher quality genome."""

        self.logger.info('Identifying improved representatives for GTDB species clusters.')
        num_improved = 0
        num_gtdb_type_sp = 0
        num_ncbi_type_sp = 0
        num_complete = 0
        num_isolate = 0
        for idx, (prev_rid, cids) in enumerate(new_updated_sp_clusters.clusters()):
            if prev_rid not in cur_genomes:
                # indicates genome has been lost
                continue
                
            prev_gtdb_sp = new_updated_sp_clusters.get_species(prev_rid)
            statusStr = '-> Processing %d of %d (%.2f%%) species [%s: %d new/updated genomes].'.ljust(86) % (
                                idx+1, 
                                len(new_updated_sp_clusters), 
                                float(idx+1)*100/len(new_updated_sp_clusters),
                                prev_gtdb_sp,
                                len(cids))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            prev_rep_score = cur_genomes[prev_rid].score_ani(100)
            new_rid, top_score, ani, af = self.top_ani_score(prev_rid,
                                                                cids,
                                                                cur_genomes)

            params = {}
            action = None
            if top_score > prev_rep_score + self.new_rep_qs_threshold:
                action = 'IMPROVED_REP:REPLACED:HIGHER_QS'
                assert(prev_rid != new_rid)
                
                if cur_genomes[prev_rid].genomic_file is None:
                    self.logger.error(f'No genomic file: {prev_rid}')
                    sys.exit(-1)
                    
                if cur_genomes[new_rid].genomic_file is None:
                    self.logger.error(f'No genomic file for: {new_rid}')
                    sys.exit(-1)

                new_reps[prev_gtdb_sp] = new_rid
                num_improved += 1

                params['new_rid'] = new_rid
                params['ani'] = ani
                params['af'] = af
                params['new_rid_assembly_quality'] = cur_genomes[new_rid].score_assembly
                params['prev_rid_assembly_quality'] = prev_genomes[prev_rid].score_assembly

                improvement_list = []
                if cur_genomes[new_rid].is_gtdb_type_strain and not cur_genomes[prev_rid].is_gtdb_type_strain:
                    num_gtdb_type_sp += 1
                    improvement_list.append('replaced with genome assembled from type strain according to GTDB')
                elif cur_genomes[new_rid].is_ncbi_type_strain and not cur_genomes[prev_rid].is_ncbi_type_strain:
                    num_ncbi_type_sp += 1
                    improvement_list.append('replaced with genome assembled from type strain according to NCBI')

                if cur_genomes[new_rid].is_isolate and not cur_genomes[prev_rid].is_isolate:
                    num_isolate += 1
                    improvement_list.append('MAG/SAG replaced with isolate')

                if cur_genomes[new_rid].is_complete_genome and not cur_genomes[prev_rid].is_complete_genome:
                    num_complete += 1
                    improvement_list.append('replaced with complete genome')
                    
                params['improvements'] = ':'.join(improvement_list)
                    
                self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                                prev_rid, 
                                                prev_gtdb_sp, 
                                                action, 
                                                params))
            
        sys.stdout.write('\n')
        self.logger.info(f' ...identified {num_improved:,} species with improved representatives.')
        self.logger.info(f'   ...{num_gtdb_type_sp:,} replaced with GTDB genome from type strain.')
        self.logger.info(f'   ...{num_ncbi_type_sp:,} replaced with NCBI genome from type strain.')
        self.logger.info(f'   ...{num_isolate:,} replaced MAG/SAG with isolate.')
        self.logger.info(f'   ...{num_complete:,} replaced with complete genome assembly.')
            
    def run(self, 
            rep_change_summary_file,
            prev_gtdb_metadata_file,
            prev_genomic_path_file,
            cur_gtdb_metadata_file,
            cur_genomic_path_file,
            uba_genome_paths,
            genomes_new_updated_file,
            qc_passed_file,
            gtdbtk_classify_file,
            species_exception_file,
            genus_exception_file,
            gtdb_type_strains_ledger):
        """Perform initial actions required for changed representatives."""
        
        # create previous and current GTDB genome sets
        self.logger.info('Creating previous GTDB genome set.')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                                species_exception_file,
                                                genus_exception_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                uba_genome_file=uba_genome_paths)
        self.logger.info(f' ...previous genome set contains {len(prev_genomes):,} genomes.')
        self.logger.info(' ...previous genome set has {:,} species clusters spanning {:,} genomes.'.format(
                            len(prev_genomes.sp_clusters),
                            prev_genomes.sp_clusters.total_num_genomes()))

        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                                species_exception_file,
                                                genus_exception_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                create_sp_clusters=False,
                                                uba_genome_file=uba_genome_paths,
                                                qc_passed_file=qc_passed_file)
        self.logger.info(f' ...current genome set contains {len(cur_genomes):,} genomes.')
        
        # get path to previous and current genomic FASTA files
        self.logger.info('Reading path to previous and current genomic FASTA files.')
        prev_genomes.load_genomic_file_paths(prev_genomic_path_file)
        prev_genomes.load_genomic_file_paths(uba_genome_paths)
        cur_genomes.load_genomic_file_paths(cur_genomic_path_file)
        cur_genomes.load_genomic_file_paths(uba_genome_paths)

        # created expanded previous GTDB species clusters
        new_updated_sp_clusters = SpeciesClusters()
        new_reps = {}

        self.logger.info('Creating species clusters of new and updated genomes based on GTDB-Tk classifications.')
        new_updated_sp_clusters.create_expanded_clusters(prev_genomes.sp_clusters,
                                                            genomes_new_updated_file,
                                                            qc_passed_file,
                                                            gtdbtk_classify_file)

        self.logger.info('Identified {:,} expanded species clusters spanning {:,} genomes.'.format(
                            len(new_updated_sp_clusters),
                            new_updated_sp_clusters.total_num_genomes()))
        
        # take required action for each changed representatives
        self.action_genomic_lost(rep_change_summary_file,
                                    prev_genomes, 
                                    cur_genomes,
                                    new_updated_sp_clusters,
                                    new_reps)

        self.action_genomic_update(rep_change_summary_file,
                                    prev_genomes, 
                                    cur_genomes,
                                    new_updated_sp_clusters,
                                    new_reps)

        self.action_type_strain_lost(rep_change_summary_file,
                                            prev_genomes, 
                                            cur_genomes,
                                            new_updated_sp_clusters,
                                            new_reps)

        self.action_improved_rep(prev_genomes, 
                                    cur_genomes,
                                    new_updated_sp_clusters,
                                    new_reps)

        if False:   
            self.action_species_reassigned(rep_change_summary_file,
                                            prev_genomes, 
                                            cur_genomes)

        num_retired_sp = sum([1 for v in new_reps.values() if v is None])  
        num_replaced_rids = sum([1 for v in new_reps.values() if v is not None])                         
        self.logger.info(f'Identified {num_retired_sp:,} retired species.')
        self.logger.info(f'Identified {num_replaced_rids:,} species with a modified representative genome.')
        
        self.action_log.close()
        self.fastani.write_cache()

        # write out representatives for existing species clusters
        fout = open(os.path.join(self.output_dir, 'species_reps.tsv'), 'w')
        for rid, sp_name in prev_genomes.sp_clusters.species():
            if sp_name in new_reps:
                if new_reps[sp_name] is not None:
                    fout.write(f'{sp_name}\t{new_reps[sp_name]}\tREPLACED\n')
            else:
                fout.write(f'{sp_name}\t{rid}\tUNCHANGED\n')
            
        fout.close()