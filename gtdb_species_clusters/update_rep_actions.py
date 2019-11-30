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
from copy import deepcopy
from collections import defaultdict

from gtdb_species_clusters.common import (generic_name,
                                            specific_epithet,
                                            canonical_generic_name,
                                            read_genome_path,
                                            read_gtdb_sp_clusters,
                                            read_gtdb_ncbi_taxonomy,
                                            read_ncbi_subsp,
                                            read_gtdb_metadata,
                                            read_gtdb_taxonomy,
                                            rep_change_gids,
                                            read_cur_new_updated,
                                            read_qc_file,
                                            read_gtdbtk_classifications,
                                            expand_sp_clusters)
                                            
from gtdb_species_clusters.type_genome_utils import (gtdb_type_strain_of_species,
                                                        quality_score_update,
                                                        read_quality_metadata,
                                                        is_isolate,
                                                        is_type_strain,
                                                        is_ncbi_type_strain,
                                                        is_complete_genome)

from gtdb_species_clusters.ani_cache import ANI_Cache

                                                        
class RepActions(object):
    """Perform initial actions required for changed representatives."""

    def __init__(self, ani_cache_file, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
        self.fastani = ANI_Cache(ani_cache_file, 1)
        
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
        
    def action_genomic_update(self, 
                                rep_change_summary_file,
                                prev_genomic_files,
                                cur_genomic_files):
        """Handle representatives with updated genomes."""
        
        # get genomes with specific changes
        self.logger.info('Identifying representatives with updated genomic files.')
        genomic_update_gids = rep_change_gids(rep_change_summary_file, 
                                                'GENOMIC_CHANGE', 
                                                'UPDATED')
        self.logger.info(f' ...identified {len(genomic_update_gids):,} genomes.')
        
        # calculate ANI between previous and current genomes
        for gid, prev_gtdb_sp in genomic_update_gids.items():
            ani, af = self.fastani.symmetric_ani(prev_genomic_files[gid], cur_genomic_files[gid])
            
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
                                    prev_gtdb_metadata_file,
                                    prev_sp_cluster_file,
                                    cur_gtdb_metadata_file,
                                    species_exception_file,
                                    genus_exception_file,
                                    cur_genomic_files):
        """Handle representatives with new NCBI species assignments."""
        
        # get genomes with new NCBI species assignments
        self.logger.info('Identifying representatives with new NCBI species assignments.')
        ncbi_sp_change_gids = rep_change_gids(rep_change_summary_file, 
                                                'NCBI_SPECIES_CHANGE', 
                                                'REASSIGNED')
        self.logger.info(f' ...identified {len(ncbi_sp_change_gids):,} genomes.')
        
        # get species assignment of previous GTDB species clusters
        self.logger.info('Reading previous GTDB species clusters.')
        prev_sp_clusters, gtdb_rep_sp = read_gtdb_sp_clusters(prev_sp_cluster_file)
        gtdb_rep_species = set(gtdb_rep_sp.values())
        
        gtdb_specific_epithets = defaultdict(set)
        specific_epithets_rid = defaultdict(lambda: {})
        for rid, sp in gtdb_rep_sp.items():
            gtdb_specific_epithets[generic_name(sp)].add(specific_epithet(sp))
            specific_epithets_rid[generic_name(sp)][specific_epithet(sp)] = rid
            
        assert(len(gtdb_rep_species) == len(prev_sp_clusters))
        self.logger.info(' ... identified {:,} species clusters spanning {:,} genomes.'.format(
                            len(prev_sp_clusters),
                            sum([len(cids) for cids in prev_sp_clusters.values()])))
        
        # get previous NCBI taxonomy
        self.logger.info('Reading previous NCBI and GTDB taxonomies from GTDB metadata file.')
        prev_ncbi_taxonomy, _ = read_gtdb_ncbi_taxonomy(prev_gtdb_metadata_file, 
                                                        species_exception_file,
                                                        genus_exception_file)
                                                        
        prev_gtdb_taxonomy = read_gtdb_taxonomy(prev_gtdb_metadata_file)
        
        # get current NCBI taxonomy
        self.logger.info('Reading current NCBI taxonomy from GTDB metadata file.')
        cur_ncbi_taxonomy, _ = read_gtdb_ncbi_taxonomy(cur_gtdb_metadata_file, 
                                                        species_exception_file,
                                                        genus_exception_file)
        cur_ncbi_subsp = read_ncbi_subsp(cur_gtdb_metadata_file)
                                                        
        # determine if change can be made without causing a conflict with 
        # other GTDB species clusters
        for gid, prev_gtdb_sp in ncbi_sp_change_gids.items():
            prev_ncbi_sp = prev_ncbi_taxonomy[gid][6]
            cur_ncbi_sp = cur_ncbi_taxonomy[gid][6]
            prev_gtdb_genus = generic_name(prev_gtdb_sp)
            assert(prev_ncbi_sp != cur_ncbi_sp)
            assert(prev_gtdb_sp == prev_gtdb_taxonomy[gid][6])
            
            params = {}
            params['prev_ncbi_sp'] = prev_ncbi_sp
            params['cur_ncbi_sp'] = cur_ncbi_sp
            params['prev_gtdb_sp'] = prev_gtdb_sp
            
            if specific_epithet(cur_ncbi_sp) == specific_epithet(prev_gtdb_sp):
                action = 'NCBI_SPECIES_CHANGE:REASSIGNED:UNCHANGED'
            elif specific_epithet(cur_ncbi_sp) not in gtdb_specific_epithets[prev_gtdb_genus]:
                action = 'NCBI_SPECIES_CHANGE:REASSIGNED:UPDATED'
            else:
                # check if species cluster should be merged
                ncbi_proposed_rid = specific_epithets_rid[prev_gtdb_genus][specific_epithet(cur_ncbi_sp)]
                ani, af = self.fastani.symmetric_ani_cached(gid, rid,
                                                            cur_genomic_files[gid], 
                                                            cur_genomic_files[ncbi_proposed_rid])
                                                        
                params['ani'] = ani
                params['af'] = af
                params['ncbi_subspecies'] = cur_ncbi_subsp.get(gid, 'none')
                                                        
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
                                    expanded_sp_clusters,
                                    cur_genome_quality,
                                    cur_genomic_files):
        """..."""
        # get genomes with new NCBI species assignments
        self.logger.info('Identifying representative that lost type strain genome status.')
        ncbi_type_species_lost = rep_change_gids(rep_change_summary_file, 
                                                'TYPE_STRAIN_CHANGE', 
                                                'LOST')
        self.logger.info(f' ...identified {len(ncbi_type_species_lost):,} genomes.')
        
        for rid, prev_gtdb_sp in ncbi_type_species_lost.items():
            cur_rep_quality = cur_genome_quality.get(rid, -1e6)
            
            highest_genome_quality = cur_rep_quality
            highest_quality_gid = rid
            for cid in expanded_sp_clusters[rid]:
                cid_q = cur_genome_quality.get(cid, -1e6)
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
                                                            cur_genomic_files[rid], 
                                                            cur_genomic_files[highest_quality_gid])
                                                        
                params['ani'] = ani
                params['af'] = af
            else:
                action = 'TYPE_STRAIN_CHANGE:LOST:RETAINED'
                
            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                            rid, 
                                            prev_gtdb_sp, 
                                            action, 
                                            params))
                                            
    def action_improved_rep(self,
                            expanded_sp_clusters,
                            prev_gtdb_species,
                            cur_genome_quality,
                            cur_genomic_files,
                            cur_new, 
                            cur_updated,
                            metadata):
        """Check if representative should be replace with higher quality genome."""
        
        self.logger.info('Identifying improved representatives for GTDB species clusters.')
        new_updated_gids = cur_new.union(cur_updated)
        for idx, (rid, cids) in enumerate(expanded_sp_clusters.items()):
            if rid.startswith('U'): #***
                continue
                
            if rid not in cur_genomic_files:
                # indicates genome has been lost
                continue
                
            prev_gtdb_sp = prev_gtdb_species[rid]
            new_updated_cids = set(cids).intersection(new_updated_gids)
            prev_rep_quality = cur_genome_quality.get(rid, -1e6)
            
            statusStr = '-> Processing %d of %d (%.2f%%) species [%s: %d new/updated genomes].'.ljust(86) % (
                                idx+1, 
                                len(expanded_sp_clusters), 
                                float(idx+1)*100/len(expanded_sp_clusters),
                                prev_gtdb_sp,
                                len(new_updated_cids))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            highest_genome_quality = prev_rep_quality
            highest_quality_gid = rid
            improved_cids = {}
            for cid in new_updated_cids:
                cid_qs = cur_genome_quality.get(cid, -1e6)
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
                ani, af = self.fastani.symmetric_ani(cur_genomic_files[rid], 
                                                        cur_genomic_files[highest_quality_gid])
                params['new_rep_quality'] = highest_genome_quality
                params['new_rid'] = highest_quality_gid
                params['ani'] = ani
                params['af'] = af
            elif improved_cids:
                sorted_improved_cids = sorted(improved_cids.items(), 
                                                key = lambda kv: kv[1], 
                                                reverse=True)
                for cid, qs in sorted_improved_cids:
                    ani, af = self.fastani.symmetric_ani_cached(rid, cid,
                                                                cur_genomic_files[rid], 
                                                                cur_genomic_files[cid])

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
                if is_type_strain(new_rid, metadata) and not is_type_strain(rid, metadata):
                    improvement_list.append('replaced with genome assembled from type strain according to GTDB')
                elif is_ncbi_type_strain(new_rid, metadata) and not is_ncbi_type_strain(rid, metadata):
                    improvement_list.append('replaced with genome assembled from type strain according to NCBI')

                if is_isolate(new_rid, metadata) and not is_isolate(rid, metadata):
                    improvement_list.append('MAG/SAG replaced with isolate')

                if is_complete_genome(new_rid, metadata) and not is_complete_genome(rid, metadata):
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
            prev_genomic_path_file,
            cur_genomic_path_file,
            prev_gtdb_metadata_file,
            prev_sp_cluster_file,
            cur_gtdb_metadata_file,
            genomes_new_updated_file,
            qc_passed_file,
            gtdbtk_classify_file,
            species_exception_file,
            genus_exception_file):
        """Perform initial actions required for changed representatives."""
        
        # get path to previous and current genomic FASTA files
        self.logger.info('Reading path to previous and current genomic FASTA files.')
        prev_genomic_files = read_genome_path(prev_genomic_path_file)
        cur_genomic_files = read_genome_path(cur_genomic_path_file)
        self.logger.info(' ...read path for {:,} previous and {:,} current genomes.'.format(
                            len(prev_genomic_files),
                            len(cur_genomic_files)))
                            
        # get previous GTDB species clusters
        self.logger.info('Reading previous GTDB species clusters.')
        prev_sp_clusters, prev_gtdb_species = read_gtdb_sp_clusters(prev_sp_cluster_file)
        self.logger.info(' ... identified {:,} species clusters spanning {:,} genomes.'.format(
                            len(prev_sp_clusters),
                            sum([len(cids) for cids in prev_sp_clusters.values()])))
                            
        # read GTDB-Tk classifications for new and updated genomes
        self.logger.info('Reading GTDB-Tk classifications.')
        gtdbtk_classifications = read_gtdbtk_classifications(gtdbtk_classify_file)
        self.logger.info(f' ... identified {len(gtdbtk_classifications):,} classifications.')
        
        # get new and updated genomes in current GTDB release
        self.logger.info('Reading new and updated genomes in current GTDB release.')
        cur_new, cur_updated = read_cur_new_updated(genomes_new_updated_file)
        self.logger.info(f' ... identified {len(cur_new):,} new and {len(cur_updated):,} updated genomes.')
        
        # get list of genomes passing QC
        self.logger.info('Reading genomes passing QC.')
        gids_pass_qc = read_qc_file(qc_passed_file)
        self.logger.info(f' ... identified {len(gids_pass_qc):,} genomes.')
                            
        # expand previous GTDB species clusters to contain new genomes, 
        # and verify assignment of updated genomes
        self.logger.info('Expanding previous species clusters to contain new genomes based on GTDB-Tk classifications.')
        expanded_sp_clusters, sp_assignments = expand_sp_clusters(prev_sp_clusters, 
                                                                    prev_gtdb_species, 
                                                                    gtdbtk_classifications,
                                                                    cur_new,
                                                                    cur_updated,
                                                                    gids_pass_qc)
                                                                        
        # calculate quality score for genomes
        self.logger.info('Calculating genome quality score.')
        quality_metadata = read_quality_metadata(cur_gtdb_metadata_file)
        cur_genome_quality = quality_score_update(quality_metadata.keys(), 
                                                    quality_metadata)
        
        # take required action for each changed representatives
        if False:
            self.action_genomic_update(rep_change_summary_file,
                                        prev_genomic_files,
                                        cur_genomic_files)
                                    
            self.action_species_reassigned(rep_change_summary_file,
                                            prev_gtdb_metadata_file,
                                            prev_sp_cluster_file,
                                            cur_gtdb_metadata_file,
                                            species_exception_file,
                                            genus_exception_file,
                                            cur_genomic_files)
                                            
            self.action_type_strain_lost(rep_change_summary_file,
                                            expanded_sp_clusters,
                                            cur_genome_quality,
                                            cur_genomic_files)
                                            
        self.action_improved_rep(expanded_sp_clusters,
                                    prev_gtdb_species,
                                    cur_genome_quality,
                                    cur_genomic_files,
                                    cur_new, 
                                    cur_updated,
                                    quality_metadata)
                                        

        self.action_log.close()
        self.fastani.write_cache()
