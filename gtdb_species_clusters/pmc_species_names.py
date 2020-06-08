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
import copy
import argparse
import logging
import pickle
from collections import defaultdict, Counter

from numpy import (mean as np_mean, std as np_std)

from biolib.taxonomy import Taxonomy

from gtdb_species_clusters.fastani import FastANI
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_clusters import SpeciesClusters

from gtdb_species_clusters.species_name_manager import SpeciesNameManager
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters.specific_epithet_manager import SpecificEpithetManager
from gtdb_species_clusters.genome_utils import canonical_gid
from gtdb_species_clusters.type_genome_utils import (symmetric_ani, 
                                                        read_clusters, 
                                                        write_clusters, 
                                                        write_rep_radius)
from gtdb_species_clusters.taxon_utils import (generic_name,
                                                specific_epithet,
                                                canonical_taxon,
                                                canonical_species,
                                                taxon_suffix,
                                                parse_synonyms,
                                                gtdb_merged_genera,
                                                sort_by_naming_priority,
                                                longest_common_prefix,
                                                is_placeholder_taxon,
                                                is_placeholder_sp_epithet,
                                                is_alphanumeric_taxon,
                                                is_suffixed_taxon,
                                                test_same_epithet)


class PMC_SpeciesNames(object):
    """Establish final species names based on manual curation."""

    def __init__(self, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
        self.final_name_log = open(os.path.join(output_dir, 'final_sp_name.log'), 'w')
        self.final_name_log.write('Genome ID\tGTDB species\tCase\tNote\n')

    def parse_species_rep_status(self, updated_species_reps_file):
        """Parse file indicating status of previous GTDB species representatives."""
        
        prev_gtdb_rep_status = {}
        new_to_prev_rid = {}
        with open(updated_species_reps_file) as f:
            header = f.readline().strip().split('\t')
            
            prev_rid_index = header.index('Previous representative ID')
            new_rid_index = header.index('New representative ID')
            rep_status_index = header.index('Representative status')
            
            for line in f:
                tokens = line.strip().split('\t')
                
                prev_rid = tokens[prev_rid_index]
                new_rid = tokens[new_rid_index]
                if new_rid.lower() != 'none':
                    if new_rid in prev_gtdb_rep_status:
                        self.logger.error(f'New representative indicated multiple times: {new_rid}')
                        sys.exit(-1)
                        
                    prev_gtdb_rep_status[new_rid] = tokens[rep_status_index]
                    new_to_prev_rid[new_rid] = prev_rid
                else:
                    prev_gtdb_rep_status[prev_rid] = tokens[rep_status_index]
        
        return prev_gtdb_rep_status, new_to_prev_rid

    def key_taxon(self, gid, taxonomy):
        """Get key taxon for determine species assignments."""
        
        genus = taxonomy[gid][Taxonomy.GENUS_INDEX]
        species = taxonomy[gid][Taxonomy.SPECIES_INDEX]
        generic = generic_name(species)
        specific = specific_epithet(species)

        return genus, species, generic, specific
        
    def create_placeholder_species_name(self, gid, generic, specific):
        """Create most appropriate placeholder species name.
        
        Species names with a Latin specific name are suffixed, other names
        are turned into sp<accn> IDs.
        """
        
        sp = 's__{} {}'.format(generic, specific)
        
        if is_placeholder_sp_epithet(specific):
            specific = self.sp_name_mngr.numeric_placeholder_sp_epithet(gid)
        else:
            specific = self.sp_name_mngr.suffixed_placeholder_sp_epithet(generic, specific)
            
        return 's__{} {}'.format(generic, specific)
        
    def create_numeric_sp_placeholder(self, gid, generic):
        """Explicitly create a numeric sp<accn> placeholder name."""
        
        specific = self.sp_name_mngr.numeric_placeholder_sp_epithet(gid)
        
        return 's__{} {}'.format(generic, specific)

    def set_type_species(self, 
                            rid, 
                            final_taxonomy, 
                            cur_genomes, 
                            case_count):
        """Establish final species name for type species of genus genome."""
        
        note = ''
        
        gtdb_genus, gtdb_species, gtdb_generic, gtdb_specific = self.key_taxon(rid, 
                                                                                final_taxonomy)
        ncbi_genus = cur_genomes[rid].ncbi_taxa.genus
        ncbi_sp = cur_genomes[rid].ncbi_taxa.species
        
        if gtdb_species != ncbi_sp:
            gtdb_priority = self.sp_priority_mngr.genus_priority(gtdb_genus, ncbi_genus)

            if gtdb_genus != ncbi_genus and (gtdb_priority is None or gtdb_priority == ncbi_genus):
                case_count['INCONGRUENT_TYPE_SPECIES_GENERIC'] += 1
                note = 'incongruent type species'
                print('INCONGRUENT_TYPE_SPECIES_GENERIC', rid, gtdb_genus, ncbi_genus)
            
            if not test_same_epithet(gtdb_specific, specific_epithet(ncbi_sp)):
                case_count['INCONGRUENT_TYPE_SPECIES_SPECIFIC'] += 1
                note = 'incongruent type species'
                print('INCONGRUENT_TYPE_SPECIES_SPECIFIC', rid, gtdb_species, ncbi_sp)

        return gtdb_species, note
        
    def set_type_strain(self, rid, final_taxonomy, cur_genomes, case_count):
        """Establish final species name for type strain of species genome."""
        
        note = ''
        
        gtdb_genus, gtdb_species, gtdb_generic, gtdb_specific = self.key_taxon(rid, 
                                                                                final_taxonomy)
                                                                                
        ncbi_sp = cur_genomes[rid].ncbi_taxa.species
        ncbi_specific = specific_epithet(ncbi_sp)

        final_sp = gtdb_species
        if not test_same_epithet(gtdb_specific, ncbi_specific):
            if test_same_epithet(canonical_taxon(gtdb_specific), ncbi_specific):
                # epithet are a match if suffix is removed, and the suffix
                # is by definition in error since this is a type strain genome
                case_count['ERRONEOUS_SUFFIX_TYPE_STRAIN'] += 1
                note = 'type strain genome had erroneous suffix'
                print('ERRONEOUS_SUFFIX_TYPE_STRAIN', rid, gtdb_species, ncbi_sp)
                final_sp = 's__{} {}'.format(gtdb_generic, canonical_taxon(gtdb_specific))
            else:
                # epithets disagree and require manual curation
                case_count['INCONGRUENT_TYPE_STRAIN'] += 1
                note = 'incongruent type strain'
                print('INCONGRUENT_TYPE_STRAIN', rid, gtdb_species, ncbi_sp)

        return final_sp, note
        
    def set_nontype_cluster(self,
                            rid,
                            final_taxonomy, 
                            cur_genomes, 
                            prev_genomes,
                            cur_clusters,
                            new_to_prev_rid,
                            misclassified_gids,
                            ncbi_synonyms,
                            case_count):
        """Establish final species name for non-type material genome."""
        
        forbidden_names = set(['cyanobacterium'])
        
        # get taxon names
        gtdb_genus, gtdb_species, gtdb_generic, gtdb_specific = self.key_taxon(rid, 
                                                                                final_taxonomy)
                                                                                
        ncbi_species = cur_genomes[rid].ncbi_taxa.species
        ncbi_species = ncbi_synonyms.get(ncbi_species, ncbi_species)
        ncbi_specific = specific_epithet(ncbi_species)
        
        final_species = gtdb_species
        note = ''
        
        # check if genome has a forbidden specific epithet so should be treated
        # as having no species assignment at NCBI
        if gtdb_specific in forbidden_names:
            note = 'forbidden specific epithet'
            
            prev_rid = new_to_prev_rid.get(rid, rid)
            
            if prev_rid in prev_genomes:
                prev_sp = prev_genomes[prev_rid].gtdb_taxa.species
                prev_gtdb_specific = specific_epithet(prev_sp)
                if prev_gtdb_specific.startswith('sp'):
                    final_species = 's__{} {}'.format(gtdb_generic, prev_gtdb_specific)
                else:
                    self.logger.warning('Unusual case for {} with forbidden specific epithet: {}'.format(rid, prev_sp))
                    final_species = self.create_numeric_sp_placeholder(rid, gtdb_generic)
            else:
                final_species = self.create_numeric_sp_placeholder(rid, gtdb_generic)
                
        # Latin specific names are not allowed with alphanumeric generic/genus names
        if is_alphanumeric_taxon(gtdb_genus) and not is_placeholder_sp_epithet(gtdb_specific):
            note = 'Latin specific name not allowed with placeholder generic name'
            
            prev_rid = new_to_prev_rid.get(rid, rid)
            if prev_rid in prev_genomes:
                prev_sp = prev_genomes[prev_rid].gtdb_taxa.species
                prev_gtdb_specific = specific_epithet(prev_sp)
                if prev_gtdb_specific.startswith('sp'):
                    final_species = 's__{} {}'.format(gtdb_generic, prev_gtdb_specific)
                else:
                    self.logger.warning('Unusual case for {} Latin specific name with non-suffixed, placeholder generic name: {}'.format(rid, prev_sp))
                    final_species = self.create_numeric_sp_placeholder(rid, gtdb_generic)
            else:
                final_species = self.create_numeric_sp_placeholder(rid, gtdb_generic)
                
        # check if NCBI species of genome is misclassified and thus specific name needs to be a placeholder
        if rid in misclassified_gids:
            note = 'misclassified NCBI genome'

            using_prev_specific_placeholder = False
            prev_rid = new_to_prev_rid.get(rid, rid)
            if prev_rid in prev_genomes:
                prev_sp = prev_genomes[prev_rid].gtdb_taxa.species
                prev_gtdb_specific = specific_epithet(prev_sp)
                if is_placeholder_sp_epithet(prev_gtdb_specific):
                    final_species = 's__{} {}'.format(gtdb_generic, prev_gtdb_specific)
                    using_prev_specific_placeholder = True
                    
            # need to generate a specific placeholder name which should be a suffixed
            # placeholder if the misclassified genome is still in the expected family,
            # and a alphanumeric placeholder is the misclassified genome is in a different
            # family
            if not using_prev_specific_placeholder:
                type_gid = misclassified_gids[rid]
                is_alphanumeric_genus = is_alphanumeric_taxon(gtdb_genus)
                is_same_family = final_taxonomy[rid][Taxonomy.FAMILY_INDEX] == final_taxonomy[type_gid][Taxonomy.FAMILY_INDEX]
                if is_same_family:
                    specific = self.sp_name_mngr.suffixed_placeholder_sp_epithet(gtdb_generic, ncbi_specific)
                    final_species = 's__{} {}'.format(gtdb_generic, specific)
                else:
                    final_species = self.create_numeric_sp_placeholder(rid, gtdb_generic)

        return final_species, note
        
    def resolve_specific_epithet_suffixes(self, final_taxonomy, cur_genomes):
        """Resolve cases where specific epithet needs to be modified to account for genus transfer."""
        
        for rid in final_taxonomy:
            gtdb_genus, gtdb_species, gtdb_generic, gtdb_specific = self.key_taxon(rid, final_taxonomy)
            
            canonical_gtdb_species = canonical_taxon(gtdb_species)
            gtdb_species_corr = self.sp_epithet_mngr.translate_species(canonical_gtdb_species)
            gtdb_specific_corr = specific_epithet(gtdb_species_corr)

            if canonical_taxon(gtdb_specific) != canonical_taxon(gtdb_specific_corr):
                suffix = taxon_suffix(gtdb_specific)
                if suffix is None:
                    final_taxonomy[rid][Taxonomy.SPECIES_INDEX] = gtdb_species_corr
                else:
                    final_taxonomy[rid][Taxonomy.SPECIES_INDEX] = '{}_{}'.format(gtdb_species_corr, suffix)
                    
    def resolve_synonyms(self, final_taxonomy, cur_genomes, ncbi_synonyms):
        """Resolve cases where specific suffix should be changed as NCBI species are considered synonyms under the GTDB."""
        
        for rid in final_taxonomy:
            gtdb_genus, gtdb_species, gtdb_generic, gtdb_specific = self.key_taxon(rid, final_taxonomy)

            if is_placeholder_sp_epithet(gtdb_specific):
                # should not replace a placeholder name with a Latin name, since placeholder
                # names should be preserved
                continue
                
            ncbi_species = cur_genomes[rid].ncbi_taxa.species
            if ncbi_species in ncbi_synonyms:
                synonym = ncbi_synonyms[ncbi_species]

                synonym_specific = specific_epithet(synonym)
                
                if gtdb_specific != synonym_specific:
                    final_taxonomy[rid][Taxonomy.SPECIES_INDEX] = 's__{} {}'.format(gtdb_generic, synonym_specific)

    def resolve_duplicate_type_strain_species(self, mc_taxonomy, cur_genomes):
        """Resolve cases where 2 or more type strain genomes have the same species name."""
        
        sp_gids = defaultdict(list)
        for rid, taxa in mc_taxonomy.items():
            if cur_genomes[rid].is_gtdb_type_strain():
                sp = taxa[Taxonomy.SPECIES_INDEX]
                sp_gids[sp].append(rid)
                
        resolved_gids = set()
        for sp, rids in sp_gids.items():
            if len(rids) > 1:
                # if only 1 genome has a congruent NCBI and GTDB genus
                # genus assignment, the situation can be automatically resolved
                common_genera = []
                for rid in rids:
                    genus = mc_taxonomy[rid][Taxonomy.GENUS_INDEX]
                    if cur_genomes[rid].ncbi_taxa.genus == canonical_taxon(genus):
                        common_genera.append(rid)
                        
                if len(common_genera) == 1:
                    for rid in rids:
                        if rid not in common_genera:
                            resolved_gids.add(rid)
                            new_sp = self.create_placeholder_species_name(rid, generic_name(sp), specific_epithet(sp))
                            mc_taxonomy[rid][Taxonomy.SPECIES_INDEX] = new_sp
                            
                            print('Modified {} from {} to {} as multiple type strain species had the proposed name and this genome was transferred from a different NCBI genera, {}.'.format( 
                                    rid,
                                    sp,
                                    mc_taxonomy[rid][Taxonomy.SPECIES_INDEX],
                                    cur_genomes[rid].ncbi_taxa.species))
                                    
        return resolved_gids
        
    def amend_multiple_nontype_species(self, 
                                        sorted_rids, 
                                        final_taxonomy, 
                                        cur_genomes, 
                                        mc_species, 
                                        reassigned_rid_mapping):
        """Amend species assignments for genera with multiple, identical species names.
        
        Amends species assignments where there are multiple species clusters in
        a genus with the same specific name, . Ideally, this should be handled during
        initial species assignment, but this wasn't done in R95 and this acts as a good
        safety check.
        """
        
        #test_genera = ['g__Crenothrix'] #'g__Crenothrix', 'g__Pauljensenia', 'g__Acinetobacter', 'g__Thermosynechococcus'] #***
        
        starting_taxonomy = copy.deepcopy(final_taxonomy)

        # get all representatives in a genus
        genus_reps = defaultdict(list)
        for rid in sorted_rids:
            genus = final_taxonomy[rid][Taxonomy.GENUS_INDEX]

            #if genus in test_genera: #***
            genus_reps[genus].append(rid)
            
        # check if canonical, Latin specific names occurs multiple times
        for genus in genus_reps:
            generic = genus.replace('g__', '')
            
            canonical_specific_names = defaultdict(list)
            for rid in genus_reps[genus]:
                species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
                specific_name = specific_epithet(species)
                canonical_specific_name = canonical_taxon(specific_name)
                
                if not is_placeholder_sp_epithet(canonical_specific_name):
                    canonical_specific_names[canonical_specific_name].append(rid)
                
            for canonical_specific_name, rids in canonical_specific_names.items():
                if len(rids) > 1:
                    # check if there is a type strain genome
                    type_strain_genomes = [rid for rid in rids if cur_genomes[rid].is_effective_type_strain()]
                    if len(type_strain_genomes) > 1:
                        # make sure only 1 of these genomes doesn't already have a suffix
                        non_suffixed_count = 0
                        for rid in type_strain_genomes:
                            specific_name = specific_epithet(final_taxonomy[rid][Taxonomy.SPECIES_INDEX])
                            if not is_placeholder_sp_epithet(specific_name):
                                non_suffixed_count += 1
                    
                        if non_suffixed_count > 1:
                            self.logger.error('Identified multiple type strain genomes for {} {}: {}.'.format(
                                                    generic, 
                                                    canonical_specific_name, 
                                                    type_strain_genomes))
                            sys.exit(-1)
                    elif len(type_strain_genomes) == 1:
                        if final_taxonomy[type_strain_genomes[0]][Taxonomy.SPECIES_INDEX] != 's__{} {}'.format(generic, canonical_specific_name):
                            self.logger.error('Type strain genome has incorrect assignment? (rid = {}, final_taxonomy = {}, proposed = {})'.format( 
                                                    type_strain_genomes[0],
                                                    final_taxonomy[type_strain_genomes[0]][Taxonomy.SPECIES_INDEX], 
                                                    's__{} {}'.format(generic, canonical_specific_name)))
                            #***sys.exit(-1)
                    
                    # update species assignment if required
                    updated = False
                    for rid in rids:
                        if rid in type_strain_genomes:
                            continue

                        # check if species name defined by manual curation,
                        # or in previous GTDB release
                        if rid in mc_species:
                            cur_species = mc_species[rid]
                        else:
                            prev_rid = reassigned_rid_mapping.get(rid, rid)
                            cur_species = cur_genomes[prev_rid].gtdb_taxa.species
                            
                        cur_specific = specific_epithet(cur_species)
                        congruent_canonical = canonical_taxon(cur_specific) == canonical_specific_name
                        
                        if rid not in mc_species and (not congruent_canonical and cur_species != 's__'):
                            self.logger.warning('Incongruent canonical names identified when amending names requiring polyphyletic suffixes: {} {} {} {}'.format(
                                                    rid, genus, canonical_specific_name, cur_species))
                        
                        if rid in mc_species or (congruent_canonical and is_placeholder_sp_epithet(cur_specific)):
                            # use previously assigned suffix
                            updated_specific = cur_specific
                        else:
                            # need to generate a new suffix for the species
                            updated_specific = self.sp_name_mngr.suffixed_placeholder_sp_epithet(generic, canonical_specific_name)
                        
                        final_taxonomy[rid][Taxonomy.SPECIES_INDEX] = 's__{} {}'.format(generic, updated_specific)
                        updated = True
                            
                    if updated:
                        rows = ['{}_{}'.format(generic, canonical_specific_name)]
                        print_table = False
                        for rid in rids:
                            if rid in reassigned_rid_mapping:
                                prev_rid = reassigned_rid_mapping[rid]
                                rows.append('  {} -> {} = {} -> {} (effective type strain = {})'.format(
                                                prev_rid,
                                                rid, 
                                                cur_genomes[prev_rid].gtdb_taxa.species, 
                                                final_taxonomy[rid][Taxonomy.SPECIES_INDEX],
                                                cur_genomes[rid].is_effective_type_strain()))
                            else:
                                prev_rid = rid
                                rows.append('  {} = {} -> {} (effective type strain = {})'.format(
                                                rid, 
                                                cur_genomes[rid].gtdb_taxa.species, 
                                                final_taxonomy[rid][Taxonomy.SPECIES_INDEX],
                                                cur_genomes[rid].is_effective_type_strain()))
                                                
                            final_specific = specific_epithet(final_taxonomy[rid][Taxonomy.SPECIES_INDEX])
                            if cur_genomes[prev_rid].gtdb_taxa.species != final_taxonomy[rid][Taxonomy.SPECIES_INDEX]:
                                print_table = True
                            elif is_placeholder_sp_epithet(final_specific) and cur_genomes[rid].is_gtdb_type_strain():
                                print_table = True
                                                                                
                        if False: #print_table or ('g__' + generic in test_genera): #***
                            for row in rows:
                                print(row)
                                
    def identified_ncbi_misclassification(self, final_taxonomy, cur_genomes, prev_genomes, ncbi_synonyms):
        """Identify genomes that are misclassified at NCBI."""

        # get NCBI species defined by a type strain genome
        self.logger.info('Identifying misclassified genomes.')
        ncbi_anchored_species = {}
        for gid in final_taxonomy:
            ncbi_species = cur_genomes[gid].ncbi_taxa.species
            ncbi_species = ncbi_synonyms.get(ncbi_species, ncbi_species)
            
            if ncbi_species != 's__':
                gtdb_genus = final_taxonomy[gid][Taxonomy.GENUS_INDEX]
                ncbi_generic = generic_name(ncbi_species)
                anchored_sp = (cur_genomes[gid].is_effective_type_strain()
                                or gtdb_genus == 'g__' + ncbi_generic)
                                
                if anchored_sp:
                    gtdb_species = final_taxonomy[gid][Taxonomy.SPECIES_INDEX]
                    
                    if ncbi_species not in ncbi_anchored_species:
                        ncbi_anchored_species[ncbi_species] = (gtdb_species, gid)
                    elif cur_genomes[gid].is_effective_type_strain():
                        # favor having the GTDB species which is a type strain genome
                        ncbi_anchored_species[ncbi_species] = (gtdb_species, gid)
                    
        self.logger.info(' - identified {:,} NCBI anchored species.'.format(len(ncbi_anchored_species)))

        fout = open(os.path.join(self.output_dir, 'misclassified_genomes.tsv'), 'w')
        fout.write('Genome ID\tType strain genome\tProposed GTDB species\tNCBI species\tNCBI later synonym')
        fout.write('\tType strain genome ID\tGTDB species of type strain genome')
        fout.write('\tIs alphanumeric genus\tIn same family\tIn R89\tGenus requires update\n')
        misclassified_gids = {}
        for cur_gid in final_taxonomy:
            ncbi_sp = cur_genomes[cur_gid].ncbi_taxa.species
            ncbi_sp = ncbi_synonyms.get(ncbi_sp, ncbi_sp)
            
            gtdb_genus = final_taxonomy[cur_gid][Taxonomy.GENUS_INDEX]
            gtdb_sp = final_taxonomy[cur_gid][Taxonomy.SPECIES_INDEX]

            if ncbi_sp in ncbi_anchored_species:
                gtdb_type_sp, type_gid = ncbi_anchored_species[ncbi_sp]
                
                if generic_name(gtdb_sp) != generic_name(gtdb_type_sp):
                    is_alphanumeric = is_alphanumeric_taxon(gtdb_genus)
                    is_same_family = final_taxonomy[cur_gid][Taxonomy.FAMILY_INDEX] == final_taxonomy[type_gid][Taxonomy.FAMILY_INDEX]
                    in_r89 = cur_gid in prev_genomes
                    requires_update = False
                    if not in_r89:
                        if is_alphanumeric and is_same_family:
                            requires_update = True
                        if is_suffixed_taxon(gtdb_genus) and not is_same_family:
                            requires_update = True
                    
                    fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                cur_gid,
                                cur_genomes[cur_gid].is_effective_type_strain(),
                                gtdb_sp,
                                cur_genomes[cur_gid].ncbi_taxa.species,
                                ncbi_sp,
                                type_gid,
                                gtdb_type_sp,
                                is_alphanumeric,
                                is_same_family,
                                in_r89,
                                requires_update))
                    misclassified_gids[cur_gid] = type_gid
        
        self.logger.info(' - identified {:,} misclassified genomes.'.format(len(misclassified_gids)))

        fout.close()
    
        return misclassified_gids
    
    def finalize_species_name(self, 
                                rid, 
                                final_sp, 
                                case,
                                note,
                                final_taxonomy, 
                                cur_gtdb_sp):
    
        if final_sp in cur_gtdb_sp:
            prev_rid, case = cur_gtdb_sp[final_sp]
            self.logger.error('Manually-curated GTDB species name already exists: {} {} {} {}'.format(
                                final_sp, rid, prev_rid, case))
                                    
        cur_gtdb_sp[final_sp] = (rid, case)
        
        final_taxonomy[rid][Taxonomy.SPECIES_INDEX] = final_sp
        
        self.final_name_log.write('{}\t{}\t{}\t{}\n'.format(
                                    rid,
                                    final_sp,
                                    case,
                                    note))

    def finalize_species_names(self,
                                mc_taxonomy, 
                                mc_species,
                                cur_clusters,
                                reassigned_rid_mapping,
                                prev_genomes, 
                                cur_genomes, 
                                gtdb_type_strain_ledger,
                                prev_gtdb_rep_status,
                                new_to_prev_rid,
                                ncbi_synonyms):
        """Establish final species names."""
        
        final_taxonomy = copy.deepcopy(mc_taxonomy)

        # resolve specific epithets requiring changes due to genus transfers
        self.resolve_specific_epithet_suffixes(final_taxonomy, cur_genomes)

        # resolve synonyms
        self.resolve_synonyms(final_taxonomy, cur_genomes, ncbi_synonyms)

        # resolve type strain genomes with the same name
        resolved_type_strain_gids = self.resolve_duplicate_type_strain_species(final_taxonomy, 
                                                                                cur_genomes)

        # get all previous GTDB species
        prev_gtdb_sp = set()
        for gid in prev_genomes:
            prev_gtdb_sp.add(prev_genomes[gid].gtdb_taxa.species)
            
        # order representatives by naming priority
        rtn = sort_by_naming_priority(final_taxonomy,
                                        cur_clusters,
                                        prev_genomes, 
                                        cur_genomes, 
                                        gtdb_type_strain_ledger,
                                        mc_species)
        (manual_curation_rids, type_species_rids, type_strains_rids, binomial_rids, placeholder_rids) = rtn
        sorted_rids = manual_curation_rids + type_species_rids + type_strains_rids + binomial_rids + placeholder_rids

        # establish final name for each GTDB species cluster
        case_count = defaultdict(int)
        cur_gtdb_sp = {}
        
        for rid in manual_curation_rids:
            case = 'MANUAL_CURATION'
            note = 'Species name set by manual curation'
            final_sp = mc_species[rid]
            
            self.finalize_species_name(rid,
                                        final_sp,
                                        case,
                                        note,
                                        final_taxonomy, 
                                        cur_gtdb_sp)

        for rid in type_species_rids:
            case = 'TYPE_SPECIES_OF_GENUS'
            final_sp, note = self.set_type_species(rid,
                                                    final_taxonomy, 
                                                    cur_genomes, 
                                                    case_count)
                                                    
            self.finalize_species_name(rid,
                                        final_sp,
                                        case,
                                        note,
                                        final_taxonomy, 
                                        cur_gtdb_sp)
                                        
        for rid in type_strains_rids:
            case = 'TYPE_STRAIN_OF_SPECIES'
            if rid in resolved_type_strain_gids:
                note = 'Resolved duplicate type strain species name'
                final_sp = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            else:
                final_sp, note = self.set_type_strain(rid,
                                                        final_taxonomy, 
                                                        cur_genomes, 
                                                        case_count)
                                                        
            self.finalize_species_name(rid,
                                        final_sp,
                                        case,
                                        note,
                                        final_taxonomy, 
                                        cur_gtdb_sp)
                                        
        # amend species assignments where there are multiple species clusters in
        # a genus with the same specific name. Ideally, this should be handled during 
        # initial species assignment, but this wasn't done in R95 and this acts as a 
        # good safety check.
        print('###5', final_taxonomy['G000974745'][Taxonomy.SPECIES_INDEX])
        self.amend_multiple_nontype_species(sorted_rids,
                                                final_taxonomy, 
                                                cur_genomes,
                                                mc_species,
                                                reassigned_rid_mapping)

        # identify genomes that are misclassified at NCBI
        misclassified_gids = self.identified_ncbi_misclassification(final_taxonomy, 
                                                                    cur_genomes,
                                                                    prev_genomes,
                                                                    ncbi_synonyms)

        # establish final names for non-type genomes
        for rid in binomial_rids + placeholder_rids:
            case = 'NONTYPE_SPECIES_CLUSTER'
            final_sp, note = self.set_nontype_cluster(rid,
                                                        final_taxonomy, 
                                                        cur_genomes, 
                                                        prev_genomes,
                                                        cur_clusters,
                                                        new_to_prev_rid,
                                                        misclassified_gids,
                                                        ncbi_synonyms,
                                                        case_count)
                
            self.finalize_species_name(rid,
                                        final_sp,
                                        case,
                                        note,
                                        final_taxonomy, 
                                        cur_gtdb_sp)

        # write out final taxonomy
        fout = open(os.path.join(self.output_dir, 'final_taxonomy.tsv'), 'w')
        for rid, taxa in final_taxonomy.items():
            fout.write('{}\t{}\t{}\t{}\t{}\n'.format(
                        rid,
                        ';'.join(taxa),
                        cur_genomes[rid].is_gtdb_type_strain(),
                        cur_genomes[rid].is_effective_type_strain(),
                        cur_genomes[rid].is_gtdb_type_species()))
        fout.close()
        
        for case, count in case_count.items():
            print('{}\t{}'.format(case, count))

        if False:
            # status string should indicate:
            #  status of GTDB species cluster: EXISTING or NEW
            #  status of representative genome: UNCHANGED, REPLACED, or NEW
            #  status of GTDB species name: UNCHANGED, CHANGED, or NEW
            status = 'NEW_SPECIES_CLUSTER:NEW_REPRESENTATIVE:NEW_SPECIES_NAME'
            if rid in prev_gtdb_rep_status:
                status = 'EXISTING_SPECIES_CLUSTER'
                status += ':{}_REPRESENTATIVE'.format(prev_gtdb_rep_status[rid])
                if prev_gtdb_rep_status[rid] not in ['REPLACED', 'UNCHANGED', 'LOST']:
                    self.logger.error(f'Unexpected status for previous GTDB representative: {rid}, {prev_gtdb_rep_status[rid]}')
                    sys.exit(-1)

                prev_sp_rid = new_to_prev_rid.get(rid, rid)
                if final_sp == prev_genomes[prev_sp_rid].gtdb_taxa.species:
                    status += ':UNCHANGED_SPECIES_NAME'
                else:
                    status += ':CHANGED_SPECIES_NAME'
            else:
                # new GTDB species cluster so sanity check species name is new
                if rid in prev_genomes.sp_clusters:
                    self.logger.error(f'New species cluster assigned previous representative genomes: {rid}')
                                    
                if final_sp in prev_gtdb_sp:
                    self.logger.error(f'New species cluster assigned previous GTDB species name: {rid}, {final_sp}')

    def parse_gtdb_type_strain_ledger(self, gtdb_type_strains_ledger, cur_genomes):
        """Read and validate GTDB type strain ledger."""
        
        gtdb_type_strain_ledger = {}
        with open(gtdb_type_strains_ledger) as f:
            header = f.readline().strip().split('\t')
            
            gid_index = header.index('Genome ID')
            sp_name_index = header.index('Proposed species name')
            
            for line in f:
                tokens = line.strip().split('\t')
                
                gid = canonical_gid(tokens[gid_index].strip())
                gtdb_sp_name = tokens[sp_name_index].strip()
                if not gtdb_sp_name.startswith('s__'):
                    gtdb_sp_name = 's__' + gtdb_sp_name
                gtdb_type_strain_ledger[gid] = gtdb_sp_name
                
                # validate assignment
                ncbi_sp = cur_genomes[gid].ncbi_taxa.species
                if ncbi_sp != 's__' and ncbi_sp != gtdb_sp_name:
                    self.logger.warning('GTDB type strain ledger disagrees with NCBI species assignment: {} {} {}'.format(
                                        gid, gtdb_sp_name, ncbi_sp))

        return gtdb_type_strain_ledger

    def run(self,
                manual_taxonomy,
                manual_sp_names,
                pmc_custom_species,
                gtdb_clusters_file,
                prev_gtdb_metadata_file,
                cur_gtdb_metadata_file,
                uba_genome_paths,
                qc_passed_file,
                updated_species_reps_file,
                ncbi_genbank_assembly_file,
                untrustworthy_type_file,
                synonym_file,
                gtdb_type_strains_ledger,
                sp_priority_ledger,
                genus_priority_ledger,
                specific_epithet_ledger,
                dsmz_bacnames_file):
        """Finalize species names based on results of manual curation."""

        # read manually-curated taxonomy
        self.logger.info('Parsing manually-curated taxonomy.')
        mc_taxonomy = Taxonomy().read(manual_taxonomy, use_canonical_gid=True)
        self.logger.info(' - identified taxonomy strings for {:,} genomes.'.format(
                            len(mc_taxonomy)))
                            
        # read species names explicitly set via manual curation
        self.logger.info('Parsing manually-curated species.')
        mc_species = {}
        with open(manual_sp_names) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                mc_species[tokens[0]] = tokens[2]
        self.logger.info(' - identified manually-curated species names for {:,} genomes.'.format(
                            len(mc_species)))
                            
        # read post-curation, manually defined species
        self.logger.info('Parsing post-curation, manually-curated species.')
        pmc_species = {}
        with open(pmc_custom_species) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                gid = tokens[0]
                species = tokens[1]
                if gid in mc_species:
                    self.logger.warning('Manually-curated genome {} reassigned from {} to {}.'.format(
                                            gid, mc_species[gid], species))
                pmc_species[gid] = species
                mc_species[gid] = species
        self.logger.info(' - identified post-curation, manually-curated species names for {:,} genomes.'.format(
                            len(pmc_species)))

        # create previous and current GTDB genome sets
        self.logger.info('Creating previous GTDB genome set.')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                uba_genome_file=uba_genome_paths,
                                                ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                                untrustworthy_type_ledger=untrustworthy_type_file)
        self.logger.info(' - previous genome set has {:,} species clusters spanning {:,} genomes.'.format(
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
        self.logger.info(f' - current genome set contains {len(cur_genomes):,} genomes.')
        
        cur_genomes.set_prev_gtdb_classifications(prev_genomes)

        # read ledger with explicit names for type strain genomes
        self.logger.info('Validating species names given in GTDB type strain ledger.')
        gtdb_type_strain_ledger = self.parse_gtdb_type_strain_ledger(gtdb_type_strains_ledger, cur_genomes)
        self.logger.info(' - identified {:,} genomes in ledger.'.format(len(gtdb_type_strain_ledger)))
        num_conflict = 0
        for gid in gtdb_type_strain_ledger:
            if gid not in mc_taxonomy:
                continue
                
            cur_gtdb_sp = mc_taxonomy[gid][Taxonomy.SPECIES_INDEX]
            ledger_sp = gtdb_type_strain_ledger[gid]
            if cur_gtdb_sp != ledger_sp:
                num_conflict += 1
                self.logger.warning('Assigned species name for {} conflicts with GTDB type strain ledger: {} {}'.format(
                                    gid, cur_gtdb_sp, ledger_sp))
        
        # read named GTDB species clusters
        self.logger.info('Reading GTDB species clusters.')
        cur_clusters, rep_radius = read_clusters(gtdb_clusters_file)
        self.logger.info(' - identified {:,} clusters spanning {:,} genomes.'.format(
                            len(cur_clusters),
                            sum([len(gids) + 1 for gids in cur_clusters.values()])))
        assert len(set(mc_taxonomy) - set(cur_clusters)) == 0

        # get list of synonyms in order to restrict usage of species names
        self.logger.info('Reading GTDB synonyms.')
        ncbi_synonyms = parse_synonyms(synonym_file)
        self.logger.info(' - identified {:,} synonyms from {:,} distinct species.'.format(
                            len(ncbi_synonyms),
                            len(set(ncbi_synonyms.values()))))

        # create specific epithet manager
        self.sp_epithet_mngr = SpecificEpithetManager()
        self.sp_epithet_mngr.parse_specific_epithet_ledger(specific_epithet_ledger)
        self.sp_epithet_mngr.infer_epithet_map(mc_taxonomy, 
                                                mc_species, 
                                                cur_genomes, 
                                                cur_clusters)
        self.sp_epithet_mngr.write_diff_epithet_map(os.path.join(self.output_dir, 'specific_epithet_diff_map.tsv'))
        self.sp_epithet_mngr.write_epithet_map(os.path.join(self.output_dir, 'specific_epithet_map.tsv'))
        
        # create species name manager
        self.logger.info('Initializing species name manager.')
        self.sp_name_mngr = SpeciesNameManager(prev_genomes, 
                                                    cur_genomes,
                                                    None)
                                                    
        # initialize species priority manager
        self.sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger,
                                                        genus_priority_ledger,
                                                        dsmz_bacnames_file)
                                                        
        # parse file indicating status of previous GTDB species representatives
        self.logger.info('Reading status of previous GTDB species representatives.')
        prev_gtdb_rep_status, new_to_prev_rid = self.parse_species_rep_status(updated_species_reps_file)
        self.logger.info(' - read status for {:,} representatives.'.format(
                            len(prev_gtdb_rep_status)))
        assert len(prev_gtdb_rep_status) == len(prev_genomes.sp_clusters)
        
        # get species clusters with reassigned representative
        self.logger.info('Identifying species clusters with reassigned representatives.')
        reassigned_rid_mapping, _ = prev_genomes.sp_clusters.reassigned_rids(cur_clusters)
        self.logger.info(' - identified {:,} species clusters with a reassigned representative.'.format(
                            len(reassigned_rid_mapping)))

        # establish appropriate species names for GTDB clusters with new representatives
        self.finalize_species_names(mc_taxonomy, 
                                        mc_species,
                                        cur_clusters,
                                        reassigned_rid_mapping,
                                        prev_genomes, 
                                        cur_genomes,
                                        gtdb_type_strain_ledger,
                                        prev_gtdb_rep_status,
                                        new_to_prev_rid,
                                        ncbi_synonyms)
                                        
        self.final_name_log.close()
                                        