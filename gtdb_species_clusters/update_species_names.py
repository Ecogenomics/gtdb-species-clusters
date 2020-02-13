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
from gtdb_species_clusters.species_name_manager import SpeciesNameManager
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters.type_genome_utils import (symmetric_ani, 
                                                        read_clusters, 
                                                        write_clusters, 
                                                        write_rep_radius)
from gtdb_species_clusters.taxon_utils import is_placeholder_taxon


class UpdateSpeciesNames(object):
    """Update names of GTDB species clusters."""

    def __init__(self, ani_cache_file, cpus, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
        self.fastani = FastANI(ani_cache_file, cpus)
        
        self.sp_name_log = open(os.path.join(self.output_dir, 'sp_name_log.tsv'), 'w')
        self.sp_name_log.write('Genome ID\tPrevious GTDB species\tNew GTDB species\tAction\tParameters\n')
        
        self.new_sp_name = {}
        
    def parse_species_reps(self, species_reps_file):
        """Get GTDB species clusters with a new representative."""
        
        new_rids = {}
        with open(species_reps_file) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                
                prev_rid = tokens[0]
                new_rid = tokens[1]
                if prev_rid == new_rid:
                    continue
                
                if new_rid.lower() != 'none':
                    new_rids[prev_rid] = new_rid
                
        return new_rids
                
    def update_sp(self, gid, species):
        """Update name assigned to GTDB species cluster."""
    
        if gid in self.new_sp_name and self.new_sp_name[gid] != species:
            self.logger.warning('Representative {} assigned multiple distinct species names: {} {}.'.format(gid, self.new_sp_name[gid], species))
            self.logger.warning('Assuming first reassignment of {} has priority.'.format(self.new_sp_name[gid]))
            return
            
        if species in self.new_sp_name.values():
            self.logger.warning(f'Species name {species} assigned multiple times.')
    
        self.new_sp_name[gid] = species
        
    def action_improved_rep_sp_assignments(self, prev_genomes, cur_genomes, new_reps):
        """Check if GTDB species name needs to be updated due to changed representative."""

        self.logger.info('Identifying improved representatives which have new NCBI specific epithet assignments.')

        action_count = defaultdict(int)
        fout = open(os.path.join(self.output_dir, 'species_changes_new_reps.tsv'), 'w')
        fout.write('Action')
        fout.write('\tOriginal ID\tType strain\tPriority year\tGTDB species\tPrevious NCBI species\tCurrent NCBI species')
        fout.write('\tNew ID\tType strain\tPriority year\tPrevious NCBI species\tCurrent NCBI species\tUpdated GTDB species')
        fout.write('\tConflicting ID\tANI to new ID\tAF to new ID\n')
        new_rids = set(new_reps.values())
        for new_rid in cur_genomes.sort_by_assembly_score():
            if new_rid not in new_rids:
                continue
                
            orig_rids = [rid for rid, nrid in new_reps.items() if nrid == new_rid]
            if len(orig_rids) > 1:
                self.logger.error('Multiple representatives assigned to the same genome {}: {}'.format(new_rid, orig_rids))
                sys.exit(-1)
                
            orig_rid = orig_rids[0]
            actions = self.sp_name_manager.update(orig_rid, new_rid)
            for action, prev_rid, cur_rid, cur_gtdb_sp, conflicting_rid, ani, af in actions:
                action_count[action] += 1

                if action != 'UNCHANGED':
                    self.update_sp(prev_rid, cur_gtdb_sp)

                orig_prev_gtdb_sp = prev_genomes[prev_rid].gtdb_taxa.species
                orig_prev_ncbi_sp = prev_genomes[prev_rid].ncbi_taxa.species
                orig_cur_ncbi_sp = 'n/a'
                prev_is_type_strain = 'n/a'
                prev_year = 'n/a'
                if prev_rid in cur_genomes:
                    orig_cur_ncbi_sp = cur_genomes[prev_rid].ncbi_taxa.species
                    prev_is_type_strain = cur_genomes[prev_rid].is_effective_type_strain()
                    prev_year = cur_genomes[prev_rid].year_of_priority()

                new_cur_ncbi_sp = cur_genomes[cur_rid].ncbi_taxa.species
                new_prev_ncbi_sp = 'n/a'
                if cur_rid in prev_genomes:
                    new_prev_ncbi_sp = prev_genomes[cur_rid].ncbi_taxa.species

                cur_is_type_strain = cur_genomes[cur_rid].is_effective_type_strain()
                    
                fout.write(f'{action}')
                fout.write(f'\t{prev_rid}\t{prev_is_type_strain}\t{prev_year}\t{orig_prev_gtdb_sp}\t{orig_prev_ncbi_sp}\t{orig_cur_ncbi_sp}')
                fout.write(f'\t{cur_rid}\t{cur_is_type_strain}\t{cur_genomes[cur_rid].year_of_priority()}\t{new_prev_ncbi_sp}\t{new_cur_ncbi_sp}\t{cur_gtdb_sp}')
                
                if conflicting_rid:
                    fout.write('\t{}\t{:.2f}\t{:.3f}\n'.format(conflicting_rid, ani, af))
                else:
                    fout.write('\tn/a\tn/a\tn/a\n')
                
                params = {}
                params['new_rid'] = cur_rid
                params['new_ncbi_subspecies'] = cur_genomes[cur_rid].ncbi_unfiltered_taxa.subspecies
                params['orig_type_strain'] = prev_is_type_strain
                params['new_type_strain'] = cur_is_type_strain

                action = 'IMPROVED_REP:SPECIES_CHANGE:{}'.format(action)
                self.sp_name_log.write('{}\t{}\t{}\t{}\t{}\n'.format(
                                            prev_rid, 
                                            orig_prev_gtdb_sp, 
                                            cur_gtdb_sp,
                                            action, 
                                            params))

        fout.close()

        self.logger.info(f"  ... {action_count['UNCHANGED']:,} unchanged.")
        self.logger.info(f"  ... {action_count['UPDATED']:,} updated.")
        self.logger.info(f"  ... {action_count['MANUAL_CURATION']:,} require manual curation.")
        self.logger.info(f"  ... {action_count['CONFLICT']:,} conflict.")

    def action_species_reassigned(self,
                                    rep_change_summary_file,
                                    prev_genomes, 
                                    cur_genomes,
                                    improved_reps):
        """Handle representatives with new NCBI species assignments."""
        
        # get genomes with new NCBI species assignments
        self.logger.info('Identifying representatives with new NCBI specific epithet assignments.')
        ncbi_sp_change_gids = self.rep_change_gids(rep_change_summary_file, 
                                                    'NCBI_SPECIES_CHANGE', 
                                                    'REASSIGNED')
        self.logger.info(f' ... identified {len(ncbi_sp_change_gids):,} genomes.')
        
        # remove cases resolved by selecting an improved genome
        self.logger.info('Disregarding representatives resolved via an improved representative.')
        for gid in improved_reps:
            if gid in ncbi_sp_change_gids:
                del ncbi_sp_change_gids[gid]
        self.logger.info(f' ... {len(ncbi_sp_change_gids):,} genomes after filtering.')
        
        # determine if change can be made without causing a conflict with 
        # other GTDB species clusters
        action_count = defaultdict(int)
        fout = open(os.path.join(self.output_dir, 'species_changes.ncbi_sp_update.tsv'), 'w')
        fout.write('Action')
        fout.write('\tGenome ID\tGTDB species\tPrevious NCBI species\tCurrent NCBI species\tUpdated GTDB species')
        fout.write('\tPrevious type strain\tCurrent type strain')
        fout.write('\tConflicting ID\tANI to new ID\tAF to new ID\n')
        for rid in cur_genomes.sort_by_assembly_score():
            if rid not in ncbi_sp_change_gids:
                continue
                
            actions = self.sp_name_manager.update(rid, rid)
            for action, prev_rid, cur_rid, cur_gtdb_sp, conflicting_rid, ani, af in actions:
                assert prev_rid == cur_rid
                
                action_count[action] += 1
                self.update_rep(prev_rid, cur_rid, action)
                
                if action != 'UNCHANGED':
                    self.update_sp(prev_rid, cur_gtdb_sp)

                prev_gtdb_sp = prev_genomes[prev_rid].gtdb_taxa.species
                prev_ncbi_sp = prev_genomes[prev_rid].ncbi_taxa.species
                cur_ncbi_sp = cur_genomes[prev_rid].ncbi_taxa.species

                prev_is_type_strain = prev_genomes[prev_rid].is_effective_type_strain()
                cur_is_type_strain = cur_genomes[prev_rid].is_effective_type_strain()

                fout.write(f'{action}')
                fout.write(f'\t{prev_rid}\t{prev_gtdb_sp}\t{prev_ncbi_sp}\t{cur_ncbi_sp}\t{cur_gtdb_sp}')
                fout.write(f'\t{prev_is_type_strain}\t{cur_is_type_strain}')
                
                if conflicting_rid:
                    fout.write('\t{}\t{:.2f}\t{:.3f}\n'.format(conflicting_rid, ani, af))
                else:
                    fout.write('\tn/a\tn/a\tn/a\n')

                params = {}
                params['orig_gtdb_sp'] = prev_gtdb_sp
                params['new_gtdb_sp'] = cur_gtdb_sp
                params['new_ncbi_subspecies'] = cur_genomes[cur_rid].ncbi_unfiltered_taxa.subspecies
                params['orig_type_strain'] = cur_genomes[prev_rid].is_effective_type_strain()
                params['new_type_strain'] = cur_genomes[cur_rid].is_effective_type_strain()
                    
                action = 'NCBI_SPECIES_CHANGE:REASSIGNED::{}'.format(action)
                self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                                            prev_rid, 
                                            prev_gtdb_sp, 
                                            action, 
                                            params))
                                        
        fout.close()

        self.logger.info(f"  ... {action_count['UNCHANGED']:,} unchanged.")
        self.logger.info(f"  ... {action_count['UPDATED']:,} updated.")
        self.logger.info(f"  ... {action_count['MANUAL_CURATION']:,} require manual curation.")
        self.logger.info(f"  ... {action_count['CONFLICT']:,} conflict.")
        
    def _parse_explicit_sp_updates(self, gtdb_species_updates_ledger):
        """Get explicit updates to previous GTDB species names."""
        
        explicit_sp_updates = {}
        with open(gtdb_species_updates_ledger) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                
                prev_sp = tokens[0].strip()
                if not prev_sp.startswith('s__'):
                    prev_sp = 's__' + prev_sp
                    
                new_sp = tokens[1].strip()
                if not new_sp.startswith('s__'):
                    new_sp = 's__' + new_sp
                    
                explicit_sp_updates[prev_sp] = new_sp
                
        return explicit_sp_updates

    def run(self, 
                gtdb_clusters_file,
                prev_gtdb_metadata_file,
                prev_genomic_path_file,
                cur_gtdb_metadata_file,
                cur_genomic_path_file,
                uba_genome_paths,
                qc_passed_file,
                gtdbtk_classify_file,
                ncbi_genbank_assembly_file,
                untrustworthy_type_file,
                gtdb_type_strains_ledger,
                sp_priority_ledger,
                gtdb_species_updates_ledger):
        """Perform initial actions required for changed representatives."""
        
        # read named GTDB species clusters
        self.logger.info('Reading GTDB species clusters.')
        clusters, rep_radius = read_clusters(gtdb_clusters_file)
        self.logger.info(' ... identified {:,} clusters spanning {:,} genomes.'.format(
                            len(clusters),
                            sum([len(gids) + 1 for gids in clusters.values()])))
                            
        # get explicit updates to previous GTDB species names
        self.logger.info('Reading explicit species updates.')
        explicit_sp_updates = self._parse_explicit_sp_updates(gtdb_species_updates_ledger)
        self.logger.info(f' ... identified {len(explicit_sp_updates):,} updates.')
        
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
                                                untrustworthy_type_ledger=untrustworthy_type_file,
                                                gtdbtk_classify_file=gtdbtk_classify_file)
        self.logger.info(f' ... current genome set contains {len(cur_genomes):,} genomes.')
        
        # sanity check that current genomes have same GTDB assignments as the
        # current genomes
        
        # *** TBD

        # get path to previous and current genomic FASTA files
        self.logger.info('Reading path to previous and current genomic FASTA files.')
        prev_genomes.load_genomic_file_paths(prev_genomic_path_file)
        prev_genomes.load_genomic_file_paths(uba_genome_paths)
        cur_genomes.load_genomic_file_paths(cur_genomic_path_file)
        cur_genomes.load_genomic_file_paths(uba_genome_paths)
                            
        # create species name manager
        self.logger.info('Initializing species name manager.')
        self.sp_name_manager = SpeciesNameManager(prev_genomes, 
                                                    cur_genomes,
                                                    self.fastani)
                                                    
        # initialize species priority manager
        self.sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger)
 
        # establish appropriate species names for GTDB clusters with new representatives
        self.action_improved_rep_sp_assignments(prev_genomes,
                                                cur_genomes,
                                                new_reps)

        if False:
            self.action_species_reassigned(rep_change_summary_file,
                                            prev_genomes, 
                                            cur_genomes,
                                            improved_reps)

        self.sp_name_log.close()

        # write out cluster information with finalized GTDB cluster names
        self.logger.info('Writing {:,} species clusters to file.'.format(len(all_species)))
        self.logger.info('Writing {:,} cluster radius information to file.'.format(len(final_cluster_radius)))
        
        write_clusters(clusters, 
                        rep_radius, 
                        cur_genomes,
                        os.path.join(self.output_dir, 'gtdb_clusters_de_novo.tsv'))

        write_rep_radius(rep_radius, 
                            cur_genomes,
                            os.path.join(self.output_dir, 'gtdb_ani_radius_de_novo.tsv'))