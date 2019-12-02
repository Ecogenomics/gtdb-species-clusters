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


class RepActions(object):
    """Perform initial actions required for changed representatives."""

    def __init__(self, ani_cache_file, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
        self.fastani = FastANI(ani_cache_file, 1)
        
        # action parameters
        self.genomic_update_ani = 99.0
        self.genomic_update_af = 0.80
        
        self.new_rep_ani = 99.0
        self.new_rep_af = 0.80
        self.new_close_rep_qs_threshold = 5  # increase in quality score to replace with similar genome
        self.new_rep_qs_threshold = 20      # increase in quality require to select 
                                            # new representative regardless of similarity
        
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
        
    def action_genomic_update(self, 
                                rep_change_summary_file,
                                prev_genomes, 
                                cur_genomes):
        """Handle representatives with updated genomes."""
        
        # get genomes with specific changes
        self.logger.info('Identifying representatives with updated genomic files.')
        genomic_update_gids = self.rep_change_gids(rep_change_summary_file, 
                                                    'GENOMIC_CHANGE', 
                                                    'UPDATED')
        self.logger.info(f' ...identified {len(genomic_update_gids):,} genomes.')
        
        # calculate ANI between previous and current genomes
        for gid, prev_gtdb_sp in genomic_update_gids.items():
            ani, af = self.fastani.symmetric_ani_cached(f'{gid}-P', f'{gid}-C', 
                                                        prev_genomes[gid].genomic_file, 
                                                        cur_genomes[gid].genomic_file)
            
            params = {}
            params['ani'] = ani
            params['af'] = af
            
            if ani >= self.genomic_update_ani and af >= self.genomic_update_af:
                action = 'GENOMIC_CHANGE:UPDATED:RETAINED'
            else:
                action = 'GENOMIC_CHANGE:UPDATED:REPLACED'
                
            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                        gid, 
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
        for gid, prev_gtdb_sp in ncbi_sp_change_gids.items():
            prev_ncbi_sp = prev_genomes[gid].ncbi_species
            cur_ncbi_sp = cur_genomes[gid].ncbi_species
            prev_gtdb_genus = prev_genomes[gid].gtdb_genus
            assert(prev_ncbi_sp != cur_ncbi_sp)
            assert(prev_gtdb_sp == prev_genomes[gid].gtdb_species)
            
            params = {}
            params['prev_ncbi_sp'] = prev_ncbi_sp
            params['cur_ncbi_sp'] = cur_ncbi_sp
            params['prev_gtdb_sp'] = prev_gtdb_sp
            
            cur_ncbi_specific_epithet = cur_genomes[gid].ncbi_specific_epithet
            prev_gtdb_specific_epithet = prev_genomes[gid].gtdb_specific_epithet
            if cur_ncbi_specific_epithet == prev_gtdb_specific_epithet:
                action = 'NCBI_SPECIES_CHANGE:REASSIGNED:UNCHANGED'
            elif cur_ncbi_specific_epithet not in gtdb_specific_epithets[prev_gtdb_genus]:
                action = 'NCBI_SPECIES_CHANGE:REASSIGNED:UPDATED'
            else:
                # check if species cluster should be merged
                ncbi_proposed_rid = specific_epithets_rid[prev_gtdb_genus][cur_ncbi_specific_epithet]
                ani, af = self.fastani.symmetric_ani_cached(gid, ncbi_proposed_rid,
                                                            cur_genomes[gid].genomic_file, 
                                                            cur_genomes[ncbi_proposed_rid].genomic_file)
                                                        
                params['ani'] = ani
                params['af'] = af
                params['ncbi_subspecies'] = cur_genomes[gid].ncbi_subspecies
                                                        
                if ani >= 95 and af >= 0.65:
                    action = 'NCBI_SPECIES_CHANGE:REASSIGNED:MERGE'
                else:
                    action = 'NCBI_SPECIES_CHANGE:REASSIGNED:CONFLICT'
                
            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                        gid, 
                                        prev_gtdb_sp, 
                                        action, 
                                        params))
                                        
    def action_type_strain_lost(self, 
                                    rep_change_summary_file,
                                    prev_genomes, 
                                    cur_genomes):
        """..."""
        # get genomes with new NCBI species assignments
        self.logger.info('Identifying representative that lost type strain genome status.')
        ncbi_type_species_lost = self.rep_change_gids(rep_change_summary_file, 
                                                    'TYPE_STRAIN_CHANGE', 
                                                    'LOST')
        self.logger.info(f' ...identified {len(ncbi_type_species_lost):,} genomes.')
        
        prev_sp_clusters = prev_genomes.sp_clusters
        for rid, prev_gtdb_sp in ncbi_type_species_lost.items():
            cur_rep_quality = cur_genomes[rid].score_update
            
            highest_genome_quality = cur_rep_quality
            highest_quality_gid = rid
            for cid in prev_sp_clusters[rid]:
                if cid not in cur_genomes:
                    continue
                    
                cid_q = cur_genomes[cid].score_update
                if cid_q > highest_genome_quality:
                    highest_genome_quality = cid_q
                    highest_quality_gid = cid

            params = {}
            params['prev_gtdb_sp'] = prev_gtdb_sp
            params['prev_rep_quality'] = cur_rep_quality
            params['highest_genome_quality'] = highest_genome_quality
            params['highest_quality_gid'] = highest_quality_gid
            
            if highest_genome_quality > cur_rep_quality + self.new_rep_qs_threshold:
                action = 'TYPE_STRAIN_CHANGE:LOST:REPLACED'
                
                ani, af = self.fastani.symmetric_ani_cached(rid, highest_quality_gid,
                                                            cur_genomes[rid].genomic_file, 
                                                            cur_genomes[highest_quality_gid].genomic_file)
                                                        
                params['ani'] = ani
                params['af'] = af
            else:
                action = 'TYPE_STRAIN_CHANGE:LOST:RETAINED'
                
            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                            rid, 
                                            prev_gtdb_sp, 
                                            action, 
                                            params))
                                            
    def action_improved_rep(self, prev_genomes, cur_genomes):
        """Check if representative should be replace with higher quality genome."""
        
        self.logger.info('Identifying improved representatives for GTDB species clusters.')
        prev_sp_clusters = prev_genomes.sp_clusters
        new_updated_gids = prev_sp_clusters.new_gids.union(prev_sp_clusters.updated_gids)
        for idx, (rid, cids) in enumerate(prev_sp_clusters.clusters()):
            if rid not in cur_genomes:
                # indicates genome has been lost
                continue
                
            prev_gtdb_sp = prev_sp_clusters.get_species(rid)
            new_updated_cids = set(cids).intersection(new_updated_gids)
            prev_rep_quality = cur_genomes[rid].score_update
            
            statusStr = '-> Processing %d of %d (%.2f%%) species [%s: %d new/updated genomes].'.ljust(86) % (
                                idx+1, 
                                len(prev_sp_clusters), 
                                float(idx+1)*100/len(prev_sp_clusters),
                                prev_gtdb_sp,
                                len(new_updated_cids))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            highest_genome_quality = prev_rep_quality
            highest_quality_gid = rid
            improved_cids = {}
            for cid in new_updated_cids:
                cid_qs = cur_genomes[cid].score_update
                if cid_qs > highest_genome_quality:
                    highest_genome_quality = cid_qs
                    highest_quality_gid = cid
                    
                if cid_qs > prev_rep_quality + self.new_close_rep_qs_threshold:
                    improved_cids[cid] = cid_qs

            params = {}
            params['prev_gtdb_sp'] = prev_gtdb_sp
            params['prev_rep_quality'] = prev_rep_quality

            action = None
            if highest_genome_quality > prev_rep_quality + self.new_rep_qs_threshold:
                action = 'IMPROVED_REP:REPLACED:HIGHER_QS'
                
                if cur_genomes[rid].genomic_file is None:
                    print('No genomic file', rid)
                    
                if cur_genomes[highest_quality_gid].genomic_file is None:
                    print('No genomic file', highest_quality_gid)
                
                ani, af = self.fastani.symmetric_ani_cached(rid, highest_quality_gid,
                                                            cur_genomes[rid].genomic_file, 
                                                            cur_genomes[highest_quality_gid].genomic_file)
                params['new_rep_quality'] = highest_genome_quality
                params['new_rid'] = highest_quality_gid
                params['ani'] = ani
                params['af'] = af
            elif improved_cids:
                sorted_improved_cids = sorted(improved_cids.items(), 
                                                key = lambda kv: kv[1], 
                                                reverse=True)
                for cid, qs in sorted_improved_cids:
                
                    if cur_genomes[rid].genomic_file is None:
                        print('No genomic file', rid)
                    
                    if cur_genomes[cid].genomic_file is None:
                        print('No genomic file', cid)
                
                    ani, af = self.fastani.symmetric_ani_cached(rid, cid,
                                                                cur_genomes[rid].genomic_file, 
                                                                cur_genomes[cid].genomic_file)

                    if ani >= self.new_rep_ani and af >= self.new_rep_af:
                        action = 'IMPROVED_REP:REPLACED:HIGHER_QS_HIGH_SIMILARITY'
                        params['new_rep_quality'] = qs
                        params['new_rid'] = cid
                        params['ani'] = ani
                        params['af'] = af
                        break
                        
            if action:
                new_rid = params['new_rid']
                improvement_list = []
                if cur_genomes[new_rid].is_gtdb_type_strain and not cur_genomes[rid].is_gtdb_type_strain:
                    improvement_list.append('replaced with genome assembled from type strain according to GTDB')
                elif cur_genomes[new_rid].is_ncbi_type_strain and not cur_genomes[new_rid].is_ncbi_type_strain:
                    improvement_list.append('replaced with genome assembled from type strain according to NCBI')

                if cur_genomes[new_rid].is_isolate and not cur_genomes[rid].is_isolate:
                    improvement_list.append('MAG/SAG replaced with isolate')

                if cur_genomes[new_rid].is_complete_genome and not cur_genomes[rid].is_complete_genome:
                    improvement_list.append('replaced with complete genome')
                    
                params['improvements'] = ':'.join(improvement_list)
                    
                self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                                rid, 
                                                prev_gtdb_sp, 
                                                action, 
                                                params))
            
        sys.stdout.write('\n')
            
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
            genus_exception_file):
        """Perform initial actions required for changed representatives."""
        
        # create previous and current GTDB genome sets
        self.logger.info('Creating previous GTDB genome set.')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                                species_exception_file,
                                                genus_exception_file)
        self.logger.info(f' ...previous genome set contains {len(prev_genomes):,} genomes.')
        self.logger.info(' ...previous genome set has {:,} species clusters spanning {:,} genomes.'.format(
                            len(prev_genomes.sp_clusters),
                            prev_genomes.sp_clusters.total_num_genomes()))

        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                                species_exception_file,
                                                genus_exception_file,
                                                create_sp_clusters=False)
        self.logger.info(f' ...current genome set contains {len(cur_genomes):,} genomes.')
        
        # get path to previous and current genomic FASTA files
        self.logger.info('Reading path to previous and current genomic FASTA files.')
        prev_genomes.load_genomic_file_paths(prev_genomic_path_file)
        prev_genomes.load_genomic_file_paths(uba_genome_paths)
        cur_genomes.load_genomic_file_paths(cur_genomic_path_file)
        cur_genomes.load_genomic_file_paths(uba_genome_paths)

        # expand previous GTDB species clusters to contain new genomes, 
        # and verify assignment of updated genomes
        self.logger.info('Expanding previous species clusters to contain new genomes based on GTDB-Tk classifications.')
        prev_genomes.sp_clusters.expand_sp_clusters(genomes_new_updated_file,
                                                    qc_passed_file,
                                                    gtdbtk_classify_file)
        self.logger.info(' ...expanded genome set has {:,} species clusters spanning {:,} genomes.'.format(
                            len(prev_genomes.sp_clusters),
                            prev_genomes.sp_clusters.total_num_genomes()))
        
        # take required action for each changed representatives
        #if False:
        self.action_genomic_update(rep_change_summary_file,
                                    prev_genomes, 
                                    cur_genomes)
                                
        self.action_species_reassigned(rep_change_summary_file,
                                        prev_genomes, 
                                        cur_genomes)
                                        
        self.action_type_strain_lost(rep_change_summary_file,
                                        prev_genomes, 
                                        cur_genomes)
                                            
        self.action_improved_rep(prev_genomes, cur_genomes)
                                        
        self.action_log.close()
        self.fastani.write_cache()
