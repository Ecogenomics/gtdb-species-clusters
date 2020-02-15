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
import pickle
from collections import defaultdict

from numpy import (mean as np_mean, std as np_std)

from gtdb_species_clusters.fastani import FastANI
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_clusters import SpeciesClusters
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters.type_genome_utils import symmetric_ani
from gtdb_species_clusters.taxon_utils import is_placeholder_taxon
from gtdb_species_clusters.genome_utils import select_highest_quality


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

        self.new_reps = {}
        
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

        ani_af = self.fastani.pairs(gid_pairs,
                                    cur_genomes.genomic_files, 
                                    report_progress=False,
                                    check_cache=True)

        # find genome with top ANI score
        max_score = -1e6
        max_rid = None
        max_ani = None
        max_af = None
        for cid in sp_cids:
            ani, af = symmetric_ani(ani_af, prev_rid, cid)

            cur_score = cur_genomes[cid].score_ani(ani)
            if cur_score > max_score:
                max_score = cur_score
                max_rid = cid
                max_ani = ani
                max_af = af
                        
        return max_rid, max_score, max_ani, max_af
        
    def get_updated_rid(self, prev_rid):
        """Get updated representative."""
        
        if prev_rid in self.new_reps:
            gid, action = self.new_reps[prev_rid]
            return gid
        
        return prev_rid

    def update_rep(self, prev_rid, new_rid, action):
        """Update representative genome for GTDB species cluster."""

        if prev_rid in self.new_reps and self.new_reps[prev_rid][0] != new_rid:
            self.logger.warning('Representative {} was reassigned multiple times: {} {}.'.format(prev_rid, self.new_reps[prev_rid], (new_rid, action)))
            self.logger.warning('Assuming last reassignment of {}: {} has priority.'.format(new_rid, action))
        
        self.new_reps[prev_rid] = (new_rid, action)

    def genomes_in_current_sp_cluster(self, 
                                        prev_rid, 
                                        prev_genomes, 
                                        new_updated_sp_clusters, 
                                        cur_genomes):
        """Get genomes in current species cluster."""
        
        assert prev_rid in prev_genomes.sp_clusters

        sp_cids = prev_genomes.sp_clusters[prev_rid]
        if prev_rid in new_updated_sp_clusters:
            sp_cids = sp_cids.union(new_updated_sp_clusters[prev_rid])
        sp_cids = sp_cids.intersection(cur_genomes)

        return sp_cids
        
    def action_genomic_lost(self, 
                            rep_change_summary_file,
                            prev_genomes, 
                            cur_genomes,
                            new_updated_sp_clusters):
        """Handle species with lost representative genome."""
        
        # get genomes with specific changes
        self.logger.info('Identifying species with lost representative genome.')
        genomic_lost_rids = self.rep_change_gids(rep_change_summary_file, 
                                                    'GENOMIC_CHANGE', 
                                                    'LOST')
        self.logger.info(f' ... identified {len(genomic_lost_rids):,} genomes.')
        
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
                params['new_assembly_quality'] = cur_genomes[new_rid].score_assembly()
                params['prev_assembly_quality'] = prev_genomes[prev_rid].score_assembly()
                
                self.update_rep(prev_rid, new_rid, action)
            else:
                action = 'GENOMIC_CHANGE:LOST:SPECIES_RETIRED'
                self.update_rep(prev_rid, None, action)
                
            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                        prev_rid, 
                                        prev_gtdb_sp, 
                                        action, 
                                        params))
        
    def action_genomic_update(self, 
                                rep_change_summary_file,
                                prev_genomes, 
                                cur_genomes,
                                new_updated_sp_clusters):
        """Handle representatives with updated genomes."""
        
        # get genomes with specific changes
        self.logger.info('Identifying representatives with updated genomic files.')
        genomic_update_gids = self.rep_change_gids(rep_change_summary_file, 
                                                    'GENOMIC_CHANGE', 
                                                    'UPDATED')
        self.logger.info(f' ... identified {len(genomic_update_gids):,} genomes.')
        
        # calculate ANI between previous and current genomes
        assembly_score_change = []
        for prev_rid, prev_gtdb_sp in genomic_update_gids.items():
            # check that genome hasn't been lost which should
            # be handled differently
            assert prev_rid in cur_genomes
                
            ani, af = self.fastani.symmetric_ani_cached(f'{prev_rid}-P', f'{prev_rid}-C', 
                                                        prev_genomes[prev_rid].genomic_file, 
                                                        cur_genomes[prev_rid].genomic_file)

            params = {}
            params['ani'] = ani
            params['af'] = af
            params['prev_ncbi_accession'] = prev_genomes[prev_rid].ncbi_accn
            params['cur_ncbi_accession'] = cur_genomes[prev_rid].ncbi_accn
            assert prev_genomes[prev_rid].ncbi_accn != cur_genomes[prev_rid].ncbi_accn
            
            if ani >= self.genomic_update_ani and af >= self.genomic_update_af:
                params['prev_assembly_quality'] = prev_genomes[prev_rid].score_assembly()
                params['new_assembly_quality'] = cur_genomes[prev_rid].score_assembly()
                action = 'GENOMIC_CHANGE:UPDATED:MINOR_CHANGE'
                
                d = cur_genomes[prev_rid].score_assembly() - prev_genomes[prev_rid].score_assembly()
                assembly_score_change.append(d)
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
                        params['prev_assembly_quality'] = prev_genomes[prev_rid].score_assembly()
                        params['new_assembly_quality'] = cur_genomes[prev_rid].score_assembly()
                        action = 'GENOMIC_CHANGE:UPDATED:RETAINED'
                    else:
                        action = 'GENOMIC_CHANGE:UPDATED:REPLACED'
                        params['new_rid'] = new_rid
                        params['ani'] = ani
                        params['af'] = af
                        params['new_assembly_quality'] = cur_genomes[new_rid].score_assembly()
                        params['prev_assembly_quality'] = prev_genomes[prev_rid].score_assembly()
                        
                        self.update_rep(prev_rid, new_rid, action)
                else:
                    action = 'GENOMIC_CHANGE:UPDATED:SPECIES_RETIRED'
                    self.update_rep(prev_rid, None, action)
                
            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                        prev_rid, 
                                        prev_gtdb_sp, 
                                        action, 
                                        params))
                                        
        self.logger.info(' ... change in assembly score for updated genomes: {:.2f} +/- {:.2f}'.format(
                            np_mean(assembly_score_change),
                            np_std(assembly_score_change)))

    def action_type_strain_lost(self, 
                                    rep_change_summary_file,
                                    prev_genomes, 
                                    cur_genomes,
                                    new_updated_sp_clusters):
        """Handle representatives which have lost type strain genome status."""
        
        # get genomes with new NCBI species assignments
        self.logger.info('Identifying representative that lost type strain genome status.')
        ncbi_type_species_lost = self.rep_change_gids(rep_change_summary_file, 
                                                    'TYPE_STRAIN_CHANGE', 
                                                    'LOST')
        self.logger.info(f' ... identified {len(ncbi_type_species_lost):,} genomes.')

        for prev_rid, prev_gtdb_sp in ncbi_type_species_lost.items():
            # check that genome hasn't been lost which should
            # be handled differently
            assert prev_rid in cur_genomes
                
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
                params['new_assembly_quality'] = cur_genomes[new_rid].score_assembly()
                params['prev_assembly_quality'] = prev_genomes[prev_rid].score_assembly()

                params['new_rid_strain_ids'] = prev_genomes[new_rid].ncbi_strain_identifiers
                params['new_rid_gtdb_type_designation'] = prev_genomes[new_rid].gtdb_type_designation
                params['new_rid_gtdb_type_designation_sources'] = prev_genomes[new_rid].gtdb_type_designation_sources

                self.update_rep(prev_rid, new_rid, action)
            else:
                action = 'TYPE_STRAIN_CHANGE:LOST:RETAINED'
                
            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                            prev_rid, 
                                            prev_gtdb_sp, 
                                            action, 
                                            params))

    def action_domain_change(self, 
                                rep_change_summary_file,
                                prev_genomes, 
                                cur_genomes):
        """Handle representatives which have new domain assignments."""
        
        # get genomes with new NCBI species assignments
        self.logger.info('Identifying representative with new domain assignments.')
        domain_changed = self.rep_change_gids(rep_change_summary_file, 
                                                    'DOMAIN_CHECK', 
                                                    'REASSIGNED')
        self.logger.info(f' ... identified {len(domain_changed):,} genomes.')

        for prev_rid, prev_gtdb_sp in domain_changed.items():
            action = 'DOMAIN_CHECK:REASSIGNED'
            params = {}
            params['prev_gtdb_domain'] = prev_genomes[prev_rid].gtdb_taxa.domain
            params['cur_gtdb_domain'] = cur_genomes[prev_rid].gtdb_taxa.domain

            self.update_rep(prev_rid, None, action)
            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                            prev_rid, 
                                            prev_gtdb_sp, 
                                            action, 
                                            params))
                                            
    def action_improved_rep(self, 
                            prev_genomes, 
                            cur_genomes, 
                            new_updated_sp_clusters):
        """Check if representative should be replace with higher quality genome."""

        self.logger.info('Identifying improved representatives for GTDB species clusters.')
        num_gtdb_ncbi_type_sp = 0
        num_gtdb_type_sp = 0
        num_ncbi_type_sp = 0
        num_complete = 0
        num_isolate = 0
        anis = []
        afs = []
        improved_reps = {}
        for idx, (prev_rid, cids) in enumerate(new_updated_sp_clusters.clusters()):
            if prev_rid not in cur_genomes:
                # indicates genome has been lost
                continue

            prev_gtdb_sp = new_updated_sp_clusters.get_species(prev_rid)
            statusStr = '-> Processing {:,} of {:,} ({:.2f}%) species [{}: {:,} new/updated genomes].'.format(
                                idx+1, 
                                len(new_updated_sp_clusters), 
                                float(idx+1)*100/len(new_updated_sp_clusters),
                                prev_gtdb_sp,
                                len(cids)).ljust(86)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            # get latest representative of GTDB species clusters as it may
            # have been updated by a previous update rule
            prev_updated_rid = self.get_updated_rid(prev_rid)
            
            prev_rep_score = cur_genomes[prev_updated_rid].score_ani(100)
            new_rid, top_score, ani, af = self.top_ani_score(prev_updated_rid,
                                                                cids,
                                                                cur_genomes)

            params = {}
            action = None

            if top_score > prev_rep_score + self.new_rep_qs_threshold:
                assert(prev_updated_rid != new_rid)
                
                if (cur_genomes[prev_updated_rid].is_gtdb_type_strain()
                    and cur_genomes[prev_updated_rid].ncbi_taxa.specific_epithet != cur_genomes[new_rid].ncbi_taxa.specific_epithet
                    and self.sp_priority_mngr.has_priority(cur_genomes, prev_updated_rid, new_rid)):
                    # GTDB species cluster should not be moved to a different type strain genome 
                    # that has lower naming priority
                    self.logger.warning('Reassignments to type strain genome with lower naming priority is not allowed: {}/{}/{}, {}/{}/{}'.format(
                                            prev_updated_rid,
                                            cur_genomes[prev_updated_rid].ncbi_taxa.species,
                                            cur_genomes[prev_updated_rid].year_of_priority(),
                                            new_rid,
                                            cur_genomes[new_rid].ncbi_taxa.species, 
                                            cur_genomes[new_rid].year_of_priority()))
                    continue
            
                action = 'IMPROVED_REP:REPLACED:HIGHER_QS'

                params['new_rid'] = new_rid
                params['ani'] = ani
                params['af'] = af
                params['new_assembly_quality'] = cur_genomes[new_rid].score_assembly()
                params['prev_assembly_quality'] = cur_genomes[prev_updated_rid].score_assembly()
                params['new_gtdb_type_strain'] = cur_genomes[new_rid].is_gtdb_type_strain()
                params['prev_gtdb_type_strain'] = cur_genomes[prev_updated_rid].is_gtdb_type_strain()
                params['new_ncbi_type_strain'] = cur_genomes[new_rid].is_ncbi_type_strain()
                params['prev_ncbi_type_strain'] = cur_genomes[prev_updated_rid].is_ncbi_type_strain()

                anis.append(ani)
                afs.append(af)

                improvement_list = []
                gtdb_type_improv = cur_genomes[new_rid].is_gtdb_type_strain() and not cur_genomes[prev_updated_rid].is_gtdb_type_strain()
                ncbi_type_improv = cur_genomes[new_rid].is_ncbi_type_strain() and not cur_genomes[prev_updated_rid].is_ncbi_type_strain()

                if gtdb_type_improv and ncbi_type_improv:
                    num_gtdb_ncbi_type_sp += 1
                    improvement_list.append('replaced with genome from type strain according to GTDB and NCBI')
                elif gtdb_type_improv:
                    num_gtdb_type_sp += 1
                    improvement_list.append('replaced with genome from type strain according to GTDB')
                elif ncbi_type_improv:
                    num_ncbi_type_sp += 1
                    improvement_list.append('replaced with genome from type strain according to NCBI')

                if cur_genomes[new_rid].is_isolate() and not cur_genomes[prev_updated_rid].is_isolate():
                    num_isolate += 1
                    improvement_list.append('MAG/SAG replaced with isolate')

                if cur_genomes[new_rid].is_complete_genome() and not cur_genomes[prev_updated_rid].is_complete_genome():
                    num_complete += 1
                    improvement_list.append('replaced with complete genome')
                    
                if len(improvement_list) == 0:
                    improvement_list.append('replaced with higher quality genome')
                    
                params['improvements'] = '; '.join(improvement_list)
                    
                self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                                prev_rid, 
                                                prev_gtdb_sp, 
                                                action, 
                                                params))
                                                
                improved_reps[prev_rid] = (new_rid, action)

        sys.stdout.write('\n')
        self.logger.info(f' ... identified {len(improved_reps):,} species with improved representatives.')
        self.logger.info(f'   ... {num_gtdb_ncbi_type_sp:,} replaced with GTDB/NCBI genome from type strain.')
        self.logger.info(f'   ... {num_gtdb_type_sp:,} replaced with GTDB genome from type strain.')
        self.logger.info(f'   ... {num_ncbi_type_sp:,} replaced with NCBI genome from type strain.')
        self.logger.info(f'   ... {num_isolate:,} replaced MAG/SAG with isolate.')
        self.logger.info(f'   ... {num_complete:,} replaced with complete genome assembly.')
        self.logger.info(f' ... ANI = {np_mean(anis):.2f} +/- {np_std(anis):.2f}%; AF = {np_mean(afs)*100:.2f} +/- {np_std(afs)*100:.2f}%.')
        
        return improved_reps
        
    def action_naming_priority(self,
                                prev_genomes,
                                cur_genomes,
                                new_updated_sp_clusters):
        """Check if representative should be replace with genome with higher nomenclatural priority."""
        
        self.logger.info('Identifying genomes with naming priority in GTDB species clusters.')
        
        out_file = os.path.join(self.output_dir, 'update_priority.tsv')
        fout = open(out_file, 'w')
        fout.write('NCBI species\tGTDB species\tRepresentative\tStrain IDs\tRepresentative type sources\tPriority year\tGTDB type species\tGTDB type strain\tNCBI assembly type')
        fout.write('\tNCBI synonym\tGTDB synonym\tSynonym genome\tSynonym strain IDs\tSynonym type sources\tPriority year\tGTDB type species\tGTDB type strain\tSynonym NCBI assembly type')
        fout.write('\tANI\tAF\tPriority note\n')

        num_higher_priority = 0
        assembly_score_change = []
        anis = []
        afs = []
        for idx, prev_rid in enumerate(prev_genomes.sp_clusters):
            # get type strain genomes in GTDB species cluster, including genomes new to this release
            type_strain_gids = [gid for gid in prev_genomes.sp_clusters[prev_rid] 
                                    if gid in cur_genomes and cur_genomes[gid].is_effective_type_strain()]
            if prev_rid in new_updated_sp_clusters:
                new_type_strain_gids = [gid for gid in new_updated_sp_clusters[prev_rid] if cur_genomes[gid].is_effective_type_strain()]
                type_strain_gids.extend(new_type_strain_gids)
                
            if len(type_strain_gids) == 0:
                continue
                
            # check if representative has already been updated
            updated_rid = self.get_updated_rid(prev_rid)

            type_strain_sp = set([cur_genomes[gid].ncbi_taxa.species for gid in type_strain_gids])
            if len(type_strain_sp) == 1 and updated_rid in type_strain_gids:
                continue
                
            updated_sp = cur_genomes[updated_rid].ncbi_taxa.species
            highest_priority_gid = updated_rid
            
            if updated_rid not in type_strain_gids:
                highest_priority_gid = None
                if updated_sp in type_strain_sp:
                    sp_gids = [gid for gid in type_strain_gids 
                                if cur_genomes[gid].ncbi_taxa.species == updated_sp]
                    hq_gid = select_highest_quality(sp_gids, cur_genomes)
                    highest_priority_gid = hq_gid

                #self.logger.warning('Representative is a non-type strain genome even though type strain genomes exist in species cluster: {}: {}, {}: {}'.format(
                #                    prev_rid, cur_genomes[prev_rid].is_effective_type_strain(), updated_rid, cur_genomes[updated_rid].is_effective_type_strain()))
                #self.logger.warning('Type strain genomes: {}'.format(','.join(type_strain_gids)))

            # find highest priority genome
            for sp in type_strain_sp:
                if sp == updated_sp:
                    continue
                    
                # get highest quality genome from species
                sp_gids = [gid for gid in type_strain_gids 
                            if cur_genomes[gid].ncbi_taxa.species == sp]
                hq_gid = select_highest_quality(sp_gids, cur_genomes)
                
                if highest_priority_gid is None:
                    highest_priority_gid = hq_gid
                else:
                    highest_priority_gid, note = self.sp_priority_mngr.priority(cur_genomes, 
                                                                                highest_priority_gid, 
                                                                                hq_gid)
                    
            # check if representative should be updated
            if highest_priority_gid != updated_rid:
                num_higher_priority += 1
                
                ani, af = self.fastani.symmetric_ani_cached(updated_rid, 
                                                            highest_priority_gid, 
                                                            cur_genomes[updated_rid].genomic_file, 
                                                            cur_genomes[highest_priority_gid].genomic_file)

                anis.append(ani)
                afs.append(af)
                
                d = cur_genomes[highest_priority_gid].score_assembly() - cur_genomes[updated_rid].score_assembly()
                assembly_score_change.append(d)
                
                action = 'NOMENCLATURE_PRIORITY:REPLACED'
                params = {}
                params['prev_ncbi_species'] = cur_genomes[updated_rid].ncbi_taxa.species
                params['prev_year_of_priority'] = cur_genomes[updated_rid].year_of_priority()
                params['new_ncbi_species'] = cur_genomes[highest_priority_gid].ncbi_taxa.species
                params['new_year_of_priority'] = cur_genomes[highest_priority_gid].year_of_priority()
                params['new_rid'] = highest_priority_gid
                params['ani'] = ani
                params['af'] = af
                params['priority_note'] = note

                self.update_rep(prev_rid, highest_priority_gid, action)
                self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                                prev_rid, 
                                                cur_genomes[updated_rid].gtdb_taxa.species, 
                                                action, 
                                                params))
                                                
                fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                            cur_genomes[highest_priority_gid].ncbi_taxa.species,
                            cur_genomes[highest_priority_gid].gtdb_taxa.species,
                            highest_priority_gid,
                            ','.join(sorted(cur_genomes[highest_priority_gid].strain_ids())),
                            ','.join(sorted(cur_genomes[highest_priority_gid].gtdb_type_sources())).upper().replace('STRAININFO', 'StrainInfo'),
                            cur_genomes[highest_priority_gid].year_of_priority(),
                            cur_genomes[highest_priority_gid].is_gtdb_type_species(),
                            cur_genomes[highest_priority_gid].is_gtdb_type_strain(),
                            cur_genomes[highest_priority_gid].ncbi_type_material))
                fout.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                            cur_genomes[updated_rid].ncbi_taxa.species,
                            cur_genomes[updated_rid].gtdb_taxa.species,
                            updated_rid,
                            ','.join(sorted(cur_genomes[updated_rid].strain_ids())),
                            ','.join(sorted(cur_genomes[updated_rid].gtdb_type_sources())).upper().replace('STRAININFO', 'StrainInfo'),
                            cur_genomes[updated_rid].year_of_priority(),
                            cur_genomes[updated_rid].is_gtdb_type_species(),
                            cur_genomes[updated_rid].is_gtdb_type_strain(),
                            cur_genomes[updated_rid].ncbi_type_material))
                fout.write('\t{:.3f}\t{:.4f}\t{}\n'.format(ani, af, note))
                                                
        fout.close()

        self.logger.info(f' ... identified {num_higher_priority:,} species with representative changed to genome with higher nomenclatural priority.')
        self.logger.info(' ... change in assembly score for new representatives: {:.2f} +/- {:.2f}'.format(
                            np_mean(assembly_score_change),
                            np_std(assembly_score_change)))
        self.logger.info(' ... ANI: {:.2f} +/- {:.2f}'.format(
                            np_mean(anis),
                            np_std(anis)))
        self.logger.info(' ... AF: {:.2f} +/- {:.2f}'.format(
                            np_mean(afs),
                            np_std(afs)))

    def write_updated_clusters(self, 
                                prev_genomes,
                                cur_genomes,
                                new_reps, 
                                new_updated_sp_clusters, 
                                out_file):
        """Write out updated GTDB species clusters."""
        
        self.logger.info('Writing updated GTDB species clusters to file: {}'.format(out_file))
        
        fout = open(out_file, 'w')
        fout.write('Representative genome\tGTDB species\tNo. clustered genomes\tClustered genomes\n')
        
        cur_genome_set = set(cur_genomes)

        num_clusters = 0
        for idx, prev_rid in enumerate(prev_genomes.sp_clusters):

            new_rid, action = new_reps.get(prev_rid, [prev_rid, None])
            if new_rid is None:
                continue

            sp_cids = self.genomes_in_current_sp_cluster(prev_rid, 
                                                            prev_genomes, 
                                                            new_updated_sp_clusters, 
                                                            cur_genome_set)

            fout.write('{}\t{}\t{}\t{}\n'.format(
                        new_rid,
                        prev_genomes.sp_clusters.get_species(prev_rid),
                        len(sp_cids),
                        ','.join(sp_cids)))
            num_clusters += 1
            
        fout.close()

        self.logger.info(f' ... wrote {num_clusters:,} clusters.')

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
            ncbi_genbank_assembly_file,
            untrustworthy_type_file,
            gtdb_type_strains_ledger,
            sp_priority_ledger):
        """Perform initial actions required for changed representatives."""
        
        # create previous and current GTDB genome sets
        self.logger.info('Creating previous GTDB genome set.')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                uba_genome_file=uba_genome_paths,
                                                ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                                untrustworthy_type_ledger=untrustworthy_type_file)
        self.logger.info(' ... previous genome set has {:,} species clusters spanning {:,} genomes.'.format(
                            len(prev_genomes.sp_clusters),
                            prev_genomes.sp_clusters.total_num_genomes()))

        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                create_sp_clusters=False,
                                                uba_genome_file=uba_genome_paths,
                                                qc_passed_file=qc_passed_file,
                                                ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                                untrustworthy_type_ledger=untrustworthy_type_file)
        self.logger.info(f' ... current genome set contains {len(cur_genomes):,} genomes.')

        # get path to previous and current genomic FASTA files
        self.logger.info('Reading path to previous and current genomic FASTA files.')
        prev_genomes.load_genomic_file_paths(prev_genomic_path_file)
        prev_genomes.load_genomic_file_paths(uba_genome_paths)
        cur_genomes.load_genomic_file_paths(cur_genomic_path_file)
        cur_genomes.load_genomic_file_paths(uba_genome_paths)

        # created expanded previous GTDB species clusters
        new_updated_sp_clusters = SpeciesClusters()

        self.logger.info('Creating species clusters of new and updated genomes based on GTDB-Tk classifications.')
        new_updated_sp_clusters.create_expanded_clusters(prev_genomes.sp_clusters,
                                                            genomes_new_updated_file,
                                                            qc_passed_file,
                                                            gtdbtk_classify_file)

        self.logger.info('Identified {:,} expanded species clusters spanning {:,} genomes.'.format(
                            len(new_updated_sp_clusters),
                            new_updated_sp_clusters.total_num_genomes()))
                            
        # initialize species priority manager
        self.sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger)

        # take required action for each changed representatives
        self.action_genomic_lost(rep_change_summary_file,
                                    prev_genomes, 
                                    cur_genomes,
                                    new_updated_sp_clusters)

        self.action_genomic_update(rep_change_summary_file,
                                    prev_genomes, 
                                    cur_genomes,
                                    new_updated_sp_clusters)

        self.action_type_strain_lost(rep_change_summary_file,
                                            prev_genomes, 
                                            cur_genomes,
                                            new_updated_sp_clusters)

        self.action_domain_change(rep_change_summary_file,
                                    prev_genomes, 
                                    cur_genomes)

        if True: #***
            improved_reps = self.action_improved_rep(prev_genomes, 
                                                    cur_genomes,
                                                    new_updated_sp_clusters)
                                                    
            pickle.dump(improved_reps, open(os.path.join(self.output_dir, 'improved_reps.pkl'), 'wb'))
        else:
            self.logger.warning('Reading improved_reps for pre-cached file. Generally used only for debugging.')
            improved_reps = pickle.load(open(os.path.join(self.output_dir, 'improved_reps.pkl'), 'rb'))
            
        for prev_rid, (new_rid, action) in improved_reps.items():
            self.update_rep(prev_rid, new_rid, action)
            
        self.action_naming_priority(prev_genomes,
                                    cur_genomes,
                                    new_updated_sp_clusters)

        # report basic statistics
        num_retired_sp = sum([1 for v in self.new_reps.values() if v[0] is None])  
        num_replaced_rids = sum([1 for v in self.new_reps.values() if v[0] is not None])
        self.logger.info(f'Identified {num_retired_sp:,} retired species.')
        self.logger.info(f'Identified {num_replaced_rids:,} species with a modified representative genome.')
        
        self.action_log.close()

        # write out representatives for existing species clusters
        fout = open(os.path.join(self.output_dir, 'updated_species_reps.tsv'), 'w')
        fout.write('Previous representative ID\tNew representative ID\tAction\tRepresentative status\n')
        for rid in prev_genomes.sp_clusters:
            if rid in self.new_reps:
                new_rid, action = self.new_reps[rid]
                if new_rid is not None:
                    fout.write(f'{rid}\t{new_rid}\t{action}\tREPLACED\n')
                else:
                    fout.write(f'{rid}\t{new_rid}\t{action}\tLOST\n') 
            else:
                fout.write(f'{rid}\t{rid}\tNONE\tUNCHANGED\n')
            
        fout.close()
        
        # write out updated species clusters
        out_file = os.path.join(self.output_dir, 'updated_sp_clusters.tsv')
        self.write_updated_clusters(prev_genomes, 
                                    cur_genomes,
                                    self.new_reps, 
                                    new_updated_sp_clusters,
                                    out_file)