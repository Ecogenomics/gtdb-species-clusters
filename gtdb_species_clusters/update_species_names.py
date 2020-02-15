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
from collections import defaultdict, Counter

from numpy import (mean as np_mean, std as np_std)

from gtdb_species_clusters.mash import Mash
from gtdb_species_clusters.fastani import FastANI
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_clusters import SpeciesClusters
from gtdb_species_clusters.species_name_manager import SpeciesNameManager
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters.type_genome_utils import (symmetric_ani, 
                                                        read_clusters, 
                                                        write_clusters, 
                                                        write_rep_radius)
from gtdb_species_clusters.taxon_utils import parse_synonyms, is_placeholder_taxon


class UpdateSpeciesNames(object):
    """Update names of GTDB species clusters."""

    def __init__(self, ani_cache_file, cpus, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
        self.fastani = FastANI(ani_cache_file, cpus)
        mash = Mash(cpus)
        
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
        
    def assign_species_names(self, cur_genomes, clusters, names_in_use):
        """Assign a species name to each species cluster."""
        
        orig_names_in_use = set(names_in_use)

        fout = open(os.path.join(self.output_dir, 'gtdb_assigned_sp.tsv'), 'w')
        fout.write('Representative genome\tAssigned species\tGTDB taxonomy\tNo. clustered genomes\tClustered GTDB genera\tClustered GTDB species\tSpecies name in use\tMost common name in use\tClustered genomes\n')
        cluster_sp_names = {}
        for rid in sorted(clusters, key=lambda x: len(clusters[x]), reverse=True):
            clustered_gids = [c.gid for c in clusters[rid]]
            
            # find most common genus name in cluster
            gtdb_genera = [cur_genomes[gid].gtdb_taxa.genus for gid in clustered_gids] + [cur_genomes[rid].gtdb_taxa.genus]
            gtdb_genus_counter = Counter(gtdb_genera)
            gtdb_common_genus = None 
            gtdb_common_genus_count = 0
            for genus, count in gtdb_genus_counter.most_common(): 
                if genus != 'g__':
                    gtdb_common_genus = genus
                    gtdb_common_genus_count = count
                    break
                    
            # in case of ties involving genus of representative genome, 
            # defer to classification of representative
            rep_genus = cur_genomes[rid].gtdb_taxa.genus
            if gtdb_genus_counter[rep_genus] == gtdb_common_genus_count and rep_genus != 'g__':
                gtdb_common_genus = rep_genus
            
            # get most common GTDB species name 
            gtdb_sp = [cur_genomes[gid].gtdb_taxa.species for gid in clustered_gids] + [cur_genomes[rid].gtdb_taxa.species]
            gtdb_sp_counter = Counter(gtdb_sp)
            gtdb_common_sp = None
            gtdb_common_sp_count = 0
            for sp, count in gtdb_sp_counter.most_common(): 
                if sp != 's__':
                    gtdb_common_sp = sp
                    gtdb_common_sp_count = count
                    break
                    
            most_common_in_use = gtdb_common_sp in names_in_use

            min_req_genomes = 0.5*(sum(gtdb_sp_counter.values()) - gtdb_sp_counter.get('s__', 0))
            if gtdb_common_sp_count >= min_req_genomes and not most_common_in_use:
                # assign common species if it occurs in >=50% of the clustered genomes,
                # excluding genomes with no species assignment
                names_in_use.add(gtdb_common_sp)
                cluster_sp_names[rid] = gtdb_common_sp
            else:
                # derive new species name from genus, if possible, 
                # and accession number of representative genome
                genus = '{unresolved}'
                if gtdb_common_genus and gtdb_common_genus != 'g__':
                    genus = gtdb_common_genus[3:]
                
                acc = rid
                if rid.startswith('U_'):
                    self.logger.error('User genome was selected as the representative for a new species cluster.')
                    sys.exit(-1)

                derived_sp = 's__' + '{} sp{}'.format(genus, acc[acc.rfind('_')+1:acc.rfind('.')])
                if derived_sp in names_in_use:
                    self.logger.error('Derived species name already in use: {}, {}'.format(derived_sp, acc))
                    sys.exit(-1)

                names_in_use.add(derived_sp)
                cluster_sp_names[rid] = derived_sp
                
            fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        rid, 
                        cluster_sp_names[rid],
                        cur_genomes[rid].gtdb_taxa,
                        len(clustered_gids),
                        ', '.join("%s=%r" % (genus, count) for (genus, count) in gtdb_genus_counter.most_common()),
                        ', '.join("%s=%r" % (sp, count) for (sp, count) in gtdb_sp_counter.most_common()),
                        ', '.join("%s=%s" % (sp, sp in names_in_use) for sp, _count in gtdb_sp_counter.most_common()),
                        '%s=%d' % (gtdb_common_sp, gtdb_common_sp_count) if most_common_in_use else 'n/a',
                        ', '.join(clustered_gids)))
                
        fout.close()
        
        return cluster_sp_names
        
    def name_type_strain_clusters(self,
                                        clusters,
                                        prev_genomes, 
                                        cur_genomes, 
                                        explicit_sp_updates,
                                        synonyms):
        """Assign names to clusters represented by a validally or effectively published type strain genome."""
        
        # identify NCBI species with multiple genomes assembled from type strain of species
        self.logger.info('Determining effective type strain genomes in each NCBI species.')
        ncbi_sp_type_strain_genomes = cur_genomes.ncbi_sp_effective_type_genomes()
        self.logger.info(' ... identified effective type strain genomes for {:,} NCBI species.'.format(
                            len(ncbi_sp_type_strain_genomes)))
        
        # verify that type genomes for a species are contained in a 
        # single GTDB species cluster
        rid_map = {}
        for rid, gids in clusters.items():
            rid_map[rid] = rid
            for gid in gids:
                rid_map[gid] = rid
                
        for ncbi_sp, type_gids in ncbi_sp_type_strain_genomes.items():
            gtdb_rids = set([rid_map[gid] for gid in type_gids])
            if len(gtdb_rids) > 1:
                self.logger.warning('Type strain genomes from NCBI species {} were assigned to {:,} GTDB species clusters: {}.'.format(
                                    ncbi_sp,
                                    len(gtdb_rids),
                                    [(rid, cur_genomes[rid].gtdb_taxa.species) for rid in gtdb_rids]))
                                    
                                    
        # assign names
        cluster_sp_names = {}
        used_sp_names = set()
        for rid, cids in clusters.items():
            if not cur_genomes[rid].is_effective_type_strain():
                continue
                
            ncbi_sp = cur_genomes[rid].ncbi_taxa.species
            
            if ncbi_sp in used_sp_names:
                dup_gid = [gid for gid, sp in cluster_sp_names.values() if sp == ncbi_sp][0]
                self.logger.error('NCBI species name of type strain genome already in use: {}: {}, {}'.format(
                                    ncbi_sp, rid, dup_gid))
            
            if ncbi_sp in synonyms:
                self.logger.error('Assigning synonym to type strain species cluster: {}: {}'.format(ncbi_sp, rid))
                                    
            cluster_sp_names[rid] = ncbi_sp
            used_sp_names.add(ncbi_sp)
        
        return cluster_sp_names, used_sp_names
        
    def name_binomial_clusters(self,
                                cluster_sp_names,
                                used_sp_names,
                                clusters,
                                prev_genomes, 
                                cur_genomes, 
                                explicit_sp_updates,
                                synonyms):
        """Assign names to clusters with suitable binomial Latin names which lack type strain genomes."""
        
        for rid, cids in clusters.items():
            if rid in cluster_sp_names:
                continue

            ncbi_sp = cur_genomes[rid].ncbi_taxa.species
            if is_placeholder_taxon(ncbi_sp):
                continue
            
            if ncbi_sp in used_sp_names:
                dup_gid = [gid for gid, sp in cluster_sp_names.values() if sp == ncbi_sp][0]
                self.logger.error('NCBI species name of non-type strain genome already in use: {}: {}, {}'.format(
                                    ncbi_sp, rid, dup_gid))
                                    
            if ncbi_sp in synonyms:
                self.logger.error('Assigning synonym to non-type strain species cluster: {}: {}'.format(ncbi_sp, rid))
                                    
            cluster_sp_names[rid] = ncbi_sp
            used_sp_names.add(ncbi_sp)
        
        return cluster_sp_names, used_sp_names
        
    def name_placeholder_clusters(self,
                                    cluster_sp_names,
                                    used_sp_names,
                                    clusters,
                                    prev_genomes, 
                                    cur_genomes, 
                                    explicit_sp_updates,
                                    synonyms):
        """Assign placeholder names to clusters that could not be assigned a binomial Latin name."""
        
        for rid, cids in clusters.items():
            if rid in cluster_sp_names:
                continue

            gtdb_sp = cur_genomes[rid].gtdb_taxa.species
            if gtdb_sp != 's__':
                placeholder_name = gtdb_sp
            else:
                genus = '{unresolved}'
                gtdb_genus = cur_genomes[rid].gtdb_taxa.genus
                if gtdb_genus != 'g__':
                    genus = gtdb_genus
                    
                placeholder_name = 's__' + '{} sp{}'.format(genus, rid[1:])
            
            if placeholder_name in used_sp_names:
                dup_gid = [gid for gid, sp in cluster_sp_names.values() if sp == placeholder_name][0]
                self.logger.error('Placeholder name already in use: {}: {}, {}'.format(
                                    placeholder_name, rid, dup_gid))
                                    
            if ncbi_sp in synonyms:
                self.logger.error('Assigning synonym to placeholder species cluster: {}: {}'.format(ncbi_sp, rid))
                                    
            cluster_sp_names[rid] = placeholder_name
            used_sp_names.add(placeholder_name)
        
        return cluster_sp_names, used_sp_names
        
    def update_species_names(self, 
                                clusters,
                                prev_genomes, 
                                cur_genomes, 
                                explicit_sp_updates,
                                synonyms):
        """Determine most suitable name for each GTDB species cluster."""
        
        # determine names for species clusters represented by a validally
        # or effectively published species name with a designated type strain
        cluster_sp_names, used_sp_names = self.name_type_strain_clusters(clusters,
                                                                            prev_genomes, 
                                                                            cur_genomes, 
                                                                            explicit_sp_updates,
                                                                            synonyms)
            
        # determine name for species clusters where a NCBI-define 
        # binomial latin name exists, but there is no designated type strain
        # (e.g., Candidatus species)
        cluster_sp_names, used_sp_names = self.name_binomial_clusters(cluster_sp_names,
                                                                            used_sp_names,
                                                                            clusters,
                                                                            prev_genomes, 
                                                                            cur_genomes, 
                                                                            explicit_sp_updates,
                                                                            synonyms)
        
        
        # determine names for remaining clusters which should all have a
        # placeholder name 
        cluster_sp_names, used_sp_names = self.name_placeholder_clusters(cluster_sp_names,
                                                                            used_sp_names,
                                                                            clusters,
                                                                            prev_genomes, 
                                                                            cur_genomes, 
                                                                            explicit_sp_updates,
                                                                            synonyms)
                                                                            
        return cluster_sp_names

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
        
    def write_rep_info(self, 
                        cur_genomes,
                        clusters, 
                        cluster_sp_names, 
                        excluded_from_refseq_note,
                        ani_af,
                        output_file):
        """Write out information about selected representative genomes."""
                                            
        fout = open(output_file, 'w')
        fout.write('Species\tRepresentative\tNCBI assembly level\tNCBI genome category')
        fout.write('\tGenome size (bp)\tQuality score\tCompleteness (%)\tContamination (%)\tNo. scaffolds\tNo. contigs\tN50 contigs\tAmbiguous bases\tSSU count\tSSU length (bp)')
        fout.write('\tNo. genomes in cluster\tMean ANI\tMean AF\tMin ANI\tMin AF\tNCBI exclude from RefSeq\n')
        
        for gid in clusters:
            fout.write('{}\t{}\t{}\t{}'.format(
                        cluster_sp_names[gid], 
                        gid, 
                        cur_genomes[gid].ncbi_assembly_level,
                        cur_genomes[gid].ncbi_genome_category))

            fout.write('\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{}\t{:.1f}\t{}\t{}\t{}'.format(
                            cur_genomes[gid].length,
                            cur_genomes[gid].score_assembly(), 
                            cur_genomes[gid].comp,
                            cur_genomes[gid].cont,
                            cur_genomes[gid].scaffold_count,
                            cur_genomes[gid].contig_count,
                            cur_genomes[gid].contig_n50,
                            cur_genomes[gid].ambiguous_bases,
                            cur_genomes[gid].ssu_count,
                            cur_genomes[gid].ssu_length))
                            
            anis = []
            afs = []
            for cluster_id in clusters[gid]:
                ani, af = symmetric_ani(ani_af, gid, cluster_id)
                anis.append(ani)
                afs.append(af)
            
            if anis:
                fout.write('\t{}\t{:.1f}\t{:.2f}\t{:.1f}\t{:.2f}\t{}\n'.format(len(clusters[gid]),
                                                                    np_mean(anis), np_mean(afs),
                                                                    min(anis), min(afs),
                                                                    excluded_from_refseq_note.get(gid, '')))
            else:
                fout.write('\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(len(clusters[gid]),
                                                            'n/a', 'n/a', 'n/a', 'n/a',
                                                            excluded_from_refseq_note.get(gid, '')))
        fout.close()

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
                synonym_file,
                gtdb_type_strains_ledger,
                sp_priority_ledger,
                gtdb_species_updates_ledger):
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
                                                untrustworthy_type_ledger=untrustworthy_type_file,
                                                gtdbtk_classify_file=gtdbtk_classify_file)
        self.logger.info(f' ... current genome set contains {len(cur_genomes):,} genomes.')
        
        # get path to previous and current genomic FASTA files
        self.logger.info('Reading path to previous and current genomic FASTA files.')
        prev_genomes.load_genomic_file_paths(prev_genomic_path_file)
        prev_genomes.load_genomic_file_paths(uba_genome_paths)
        cur_genomes.load_genomic_file_paths(cur_genomic_path_file)
        cur_genomes.load_genomic_file_paths(uba_genome_paths)
        
        # read named GTDB species clusters
        clusters = {}
        if False: #***
            self.logger.info('Reading GTDB species clusters.')
            clusters, rep_radius = read_clusters(gtdb_clusters_file)
            self.logger.info(' ... identified {:,} clusters spanning {:,} genomes.'.format(
                                len(clusters),
                                sum([len(gids) + 1 for gids in clusters.values()])))
                            
        # get explicit updates to previous GTDB species names
        self.logger.info('Reading explicit species updates.')
        explicit_sp_updates = self._parse_explicit_sp_updates(gtdb_species_updates_ledger)
        self.logger.info(f' ... identified {len(explicit_sp_updates):,} updates.')
        
        # get list of synonyms in order to restrict usage of species names
        self.logger.info('Reading GTDB synonyms.')
        synonyms = parse_synonyms(synonym_file)
        self.logger.info(' ... identified {:,} synonyms from {:,} distinct species.'.format(
                            len(synonyms),
                            len(set(synonyms.values()))))
        
        # set current genomes to have same GTDB assignments as in previous
        # GTDB release. This is necessary, since genomes may have different
        # NCBI accession numbers between releases and thus the previous GTDB
        # taxonomy will not be reflected in the latest GTDB database. The 
        # exception is if a genome has changed domains, in which case the
        # previous assignment is invalid.
        self.logger.info('Setting GTDB taxonomy of genomes in current genome set.')
        update_count = 0
        conflicting_domain_count = 0
        for prev_gid in prev_genomes:
            if prev_gid in cur_genomes:
                if prev_genomes[prev_gid].gtdb_taxa != cur_genomes[prev_gid].gtdb_taxa:
                    if prev_genomes[prev_gid].gtdb_taxa.domain == cur_genomes[prev_gid].gtdb_taxa.domain:
                        update_count += 1
                        cur_genomes[prev_gid].gtdb_taxa.update_taxa(prev_genomes[prev_gid].gtdb_taxa)
                    else:
                        conflicting_domain_count += 1
        self.logger.info(f' ... updated {update_count:,} genomes.')
        self.logger.info(f' ... identified {conflicting_domain_count:,} genomes with conflicting domain assignments.')

        # create species name manager
        self.logger.info('Initializing species name manager.')
        self.sp_name_manager = SpeciesNameManager(prev_genomes, 
                                                    cur_genomes,
                                                    self.fastani)
                                                    
        # initialize species priority manager
        self.sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger)
 
        # establish appropriate species names for GTDB clusters with new representatives
        self.update_species_names(clusters,
                                    prev_genomes, 
                                    cur_genomes, 
                                    explicit_sp_updates,
                                    synonyms)
                                    
                                    
                                    
        if False:
                                        
                                        
            # assign species names to de novo species clusters
            names_in_use = synonyms.union(rep_species)
            self.logger.info('Identified {:,} species names already in use.'.format(len(names_in_use)))
            self.logger.info('Assigning species name to each de novo species cluster.')
            cluster_sp_names = self.assign_species_names(cur_genomes,
                                                            rep_clusters, 
                                                            names_in_use)
                                        
                                        
                                        
                                        
                                        
            # write out file with details about selected representative genomes
            self._write_rep_info(cur_genomes,
                                    rep_clusters, 
                                    cluster_sp_names,
                                    excluded_from_refseq_note,
                                    ani_af,
                                    os.path.join(self.output_dir, 'gtdb_rep_genome_info.tsv'))
                                        
                                    

        # write out cluster information with finalized GTDB cluster names
        self.logger.info('Writing {:,} species clusters to file.'.format(len(clusters)))
        self.logger.info('Writing {:,} cluster radius information to file.'.format(len(rep_radius)))
        
        write_clusters(clusters, 
                        rep_radius, 
                        cur_genomes,
                        os.path.join(self.output_dir, 'gtdb_clusters_de_novo.tsv'))

        write_rep_radius(rep_radius, 
                            cur_genomes,
                            os.path.join(self.output_dir, 'gtdb_ani_radius_de_novo.tsv'))