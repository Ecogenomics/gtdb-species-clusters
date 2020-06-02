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

from gtdb_species_clusters.taxon_suffix_manager import TaxonSuffixManager
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
                                                taxon_suffix,
                                                parse_synonyms,
                                                gtdb_merged_genera,
                                                sort_by_naming_priority,
                                                longest_common_prefix,
                                                is_placeholder_taxon,
                                                is_placeholder_sp_epithet,
                                                test_same_epithet)

class PMC_Validation(object):
    """Validate final species names."""

    def __init__(self, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
    def report_validation(self, invalid_cases, validated_count, description):
        """Report validation results."""
        
        print('Validated {:,} genomes.'.format(validated_count))
        if len(invalid_cases) > 0:
            print(description)
            for gid in invalid_cases:
                if len(invalid_cases[gid]) == 3:
                    f1, f2, f3 = invalid_cases[gid]
                    print('{}\t{}\t{}\t{}'.format(gid, f1, f2, f3))
                elif len(invalid_cases[gid]) == 2:
                    f1, f2 = invalid_cases[gid]
                    print('{}\t{}\t{}'.format(gid, f1, f2))
                else:
                    f = invalid_cases[gid]
                    print('{}\t{}'.format(gid, f))
            print('')
    
    def validate_genus_generic_match(self, final_taxonomy):
        """Validate that genus and generic names match."""
        
        invalid_names = {}
        validated_count = 0
        for gid in final_taxonomy:
            gtdb_genus = final_taxonomy[gid][Taxonomy.GENUS_INDEX]
            gtdb_species = final_taxonomy[gid][Taxonomy.SPECIES_INDEX]
            gtdb_generic = generic_name(gtdb_species)

            if gtdb_genus != 'g__' + gtdb_generic: 
                invalid_names[gid] = (gtdb_genus, gtdb_species)
            else:
                validated_count += 1
                    
        self.report_validation(
            invalid_names, 
            validated_count, 
            'Identified {:,} genomes with unmatched genus and generic names (GTDB genus, GTDB species):'.format(len(invalid_names)))
        
    def validate_type_species(self, final_taxonomy, cur_genomes, sp_priority_mngr):
        """Validate generic name of type species."""
        
        invalid_type_sp = {}
        validated_count = 0
        for gid in final_taxonomy:
            if cur_genomes[gid].is_gtdb_type_species():
                gtdb_genus = final_taxonomy[gid][Taxonomy.GENUS_INDEX]
                gtdb_species = final_taxonomy[gid][Taxonomy.SPECIES_INDEX]
                gtdb_generic = generic_name(gtdb_species)
                assert 'g__' + gtdb_generic == gtdb_genus
                
                ncbi_genus = cur_genomes[gid].ncbi_taxa.genus

                if (gtdb_genus != ncbi_genus 
                    and sp_priority_mngr.genus_priority(gtdb_genus, ncbi_genus) != gtdb_genus):
                    invalid_type_sp[gid] = (gtdb_species, ncbi_genus)
                else:
                    validated_count += 1
                    
        self.report_validation(
            invalid_type_sp, 
            validated_count, 
            'Identified {:,} invalid type species genomes (GTDB species, NCBI genus):'.format(len(invalid_type_sp)))

    def validate_type_strain(self, final_taxonomy, cur_genomes):
        """Validate specific name of type strains."""
        
        invalid_type_strain = {}
        validated_count = 0
        for gid in final_taxonomy:
            if cur_genomes[gid].is_gtdb_type_strain():
                gtdb_species = final_taxonomy[gid][Taxonomy.SPECIES_INDEX]
                gtdb_specific = specific_epithet(gtdb_species)

                ncbi_species = cur_genomes[gid].ncbi_taxa.species
                ncbi_specific = specific_epithet(ncbi_species)
                
                if not test_same_epithet(gtdb_specific, ncbi_specific):
                    invalid_type_strain[gid] = (gtdb_species, ncbi_species)
                else:
                    validated_count += 1

        self.report_validation(
            invalid_type_strain, 
            validated_count, 
            'Identified {:,} invalid type strain genomes (GTDB species, NCBI species):'.format(len(invalid_type_strain)))
            
    def validate_manually_curated_species(self, final_taxonomy, manual_sp_names):
        """Validate final assignment of species that were set by manual curation."""
        
        # read species names explicitly set via manual curation
        invalid_mc = {}
        validated_count = 0
        with open(manual_sp_names) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                gid = tokens[0]
                mc_species = tokens[2]
                
                if final_taxonomy[gid][Taxonomy.SPECIES_INDEX] != mc_species:
                    invalid_mc[gid] = (final_taxonomy[gid][Taxonomy.SPECIES_INDEX], mc_species)
                else:
                    validated_count += 1

        self.report_validation(
            invalid_mc, 
            validated_count, 
            'Identified {:,} genomes that did not match manual curation (GTDB species, manually curated species):'.format(len(invalid_mc)))
            
    def validate_ground_truth(self, final_taxonomy, ground_truth_test_cases):
        """Validate final assignments against a set of ground truth assignments."""
        
        invalid_gt = {}
        validated_count = 0
        with open(ground_truth_test_cases) as f:
            for line in f:
                tokens = line.strip().split('\t')
                gid = tokens[0]
                gt_taxa_str = tokens[1]
                gt_taxa = [t.strip() for t in gt_taxa_str.split(';')]
                
                if final_taxonomy[gid] != gt_taxa:
                    invalid_gt[gid] = (final_taxonomy[gid][Taxonomy.SPECIES_INDEX], gt_taxa[Taxonomy.SPECIES_INDEX])
                else:
                    validated_count += 1

        self.report_validation(
            invalid_gt, 
            validated_count, 
            'Identified {:,} genomes that did not match ground truth assignments (GTDB species, ground truth species):'.format(len(invalid_gt)))

    def validate_placeholder_suffixes(self, final_taxonomy, cur_genomes, reassigned_rid_mapping):
        """Validate that species clusters have same placeholder suffixes as previously assigned."""

        invalid_suffix = {}
        validated_count = 0
        suffix_count = 0
        for rid in final_taxonomy:
            gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            gtdb_specific = specific_epithet(gtdb_species)
            
            suffix = taxon_suffix(gtdb_specific)
            if suffix is None:
                continue
                
            suffix_count += 1
            
            prev_rid = reassigned_rid_mapping.get(rid, rid)
            if prev_rid in cur_genomes:
                prev_species = cur_genomes[prev_rid].gtdb_taxa.species
                prev_specific = specific_epithet(prev_species)
                prev_suffix = taxon_suffix(prev_specific)
                
                if prev_species == 's__':
                    # since genome has no previous GTDB assignment it should have a
                    # suffix higher than observed in the previous taxonomy
                    highest_suffix = self.taxon_suffix_manager.highest_suffix(canonical_taxon(gtdb_species))
                    
                    if highest_suffix is None or self.taxon_suffix_manager.is_higher_suffix(suffix, highest_suffix):
                        validated_count += 1
                    else:
                        invalid_suffix[rid] = (gtdb_species, prev_species)
                elif prev_suffix is None:
                    # this should only occur if a species name was previously unique,
                    # and is now associated with multiple species clusters
                    num_canonical_sp = sum([1 for gid in final_taxonomy 
                                            if canonical_taxon(final_taxonomy[gid][Taxonomy.SPECIES_INDEX]) == canonical_taxon(gtdb_species)])
                    #***print(rid, gtdb_species, prev_species, num_canonical_sp)
                    
                    if num_canonical_sp > 1:
                        validated_count += 1
                    else:
                        invalid_suffix[rid] = (gtdb_species, prev_species)
                else:
                    # suffix should be the same between releases
                    if suffix != prev_suffix:
                        invalid_suffix[rid] = (gtdb_species, prev_species)
                    else:
                        validated_count += 1

        print('Identified {:,} genomes with a suffixed specific name.'.format(suffix_count))
        self.report_validation(
            invalid_suffix, 
            validated_count, 
            'Identified {:,} genomes with an invalid placeholder suffix (Current GTDB species, Previous GTDB species):'.format(len(invalid_suffix)))

    def validate_no_suffix(self, final_taxonomy, cur_genomes, reassigned_rid_mapping):
        """Validate GTDB species assignments without a polyphyletic suffixes."""
        
        # determine number of GTDB species with same canonical name
        canonical_sp_count = defaultdict(set)
        for rid in final_taxonomy:
            gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            canonical_species = canonical_taxon(gtdb_species)
            
            canonical_sp_count[canonical_species].add(gtdb_species)
        
        invalid_names = {}
        validated_count = 0
        case_count = 0
        for rid in final_taxonomy:
            gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            canonical_species = canonical_taxon(gtdb_species)
            gtdb_specific = specific_epithet(gtdb_species)
            
            suffix = taxon_suffix(gtdb_specific)
            if suffix is not None:
                continue
                
            case_count += 1
            
            if cur_genomes[rid].is_effective_type_strain():
                validated_count += 1
            elif len(canonical_sp_count[canonical_species]) == 1:
                validated_count += 1
            else:
                invalid_names[rid] = (gtdb_species, 
                                        cur_genomes[rid].is_effective_type_strain(), 
                                        ', '.join(canonical_sp_count[canonical_species]))

        print('Identified {:,} genomes without a suffixed specific name.'.format(case_count))
        self.report_validation(
            invalid_names, 
            validated_count, 
            'Identified {:,} genomes which should have a polyphyletic suffix (Current GTDB species, effective type material, GTDB canonical species):'.format(len(invalid_names)))
        
    def validate_forbidden_specific_names(self, final_taxonomy):
        """Validate that no species names contain a forbidden specific name."""
        
        forbidden_names = set(['cyanobacterium'])
        
        invalid_specific = {}
        validated_count = 0
        for gid, taxa in final_taxonomy.items():
            gtdb_species = taxa[Taxonomy.SPECIES_INDEX]
            specific = specific_epithet(gtdb_species)
            
            if specific in forbidden_names:
                invalid_specific[gid] = gtdb_species
            else:
                validated_count += 1
        
        self.report_validation(
            invalid_specific, 
            validated_count, 
            'Identified {:,} genomes with a forbidden specific name:'.format(len(invalid_specific)))
    
    def validate_latin_sp_placeholder_generic(self, final_taxonomy):
        """Validate that Latinized specific names are never used with placeholder generic names."""

        invalid_sp_name = {}
        validated_count = 0
        for gid, taxa in final_taxonomy.items():
            gtdb_species = taxa[Taxonomy.SPECIES_INDEX]
            canonical_generic = canonical_taxon(generic_name(gtdb_species))
            specific = specific_epithet(gtdb_species)

            if (is_placeholder_taxon('g__' + canonical_generic) and 
                not is_placeholder_sp_epithet(specific)):
                    invalid_sp_name[gid] = gtdb_species
            else:
                validated_count += 1

        self.report_validation(
            invalid_sp_name, 
            validated_count, 
            'Identified {:,} genomes with a non-suffixed, placeholder generic name, but Latin specific name:'.format(len(invalid_sp_name)))
        
    def report_modified_species_names(self, final_taxonomy, prev_genomes, prev_to_new_rids):
        """Report GTDB species clusters with modified names."""
        
        domains_of_interest = set([taxa[Taxonomy.DOMAIN_INDEX] for taxa in final_taxonomy.values()])
        
        fout = open(os.path.join(self.output_dir, 'modified_sp_names.tsv'), 'w')
        fout.write('Previous ID\tNew ID\tPrevious name\tNew name\tModified generic name\tModified specific name\tUpdated representative\n')
        
        modified_generic_names = 0
        modified_specific_names = 0
        same_species_names = 0
        lost_sp_cluster = 0
        total_names = 0
        for prev_rid in prev_genomes.sp_clusters:
            prev_domain = prev_genomes[prev_rid].gtdb_taxa.domain
            if prev_domain not in domains_of_interest:
                continue
                
            total_names += 1

            prev_sp = prev_genomes[prev_rid].gtdb_taxa.species
            prev_genus = generic_name(prev_sp)
            prev_specific = specific_epithet(prev_sp)
            
            new_rid = prev_to_new_rids.get(prev_rid, prev_rid)
            if new_rid in final_taxonomy:
                new_sp = final_taxonomy[new_rid][Taxonomy.SPECIES_INDEX]
                new_genus = generic_name(new_sp)
                new_specific = specific_epithet(new_sp)
            else:
                lost_sp_cluster += 1
                fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            prev_rid,
                            'n/a',
                            prev_sp,
                            'n/a',
                            'n/a',
                            'n/a',
                            'n/a'))
                continue
            
            if prev_sp == new_sp:
                same_species_names += 1
            else:
                fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            prev_rid,
                            new_rid,
                            prev_sp,
                            new_sp,
                            prev_genus != new_genus,
                            prev_specific != new_specific,
                            prev_rid != new_rid))
                
            if prev_genus != new_genus:
                modified_generic_names += 1
                
            if prev_specific != new_specific:
                modified_specific_names += 1
            
        self.logger.info(' - total previous species names: {:,}'.format(total_names))
        self.logger.info(' - unmodified names: {:,} ({:.2f}%)'.format(
                        same_species_names,
                        same_species_names*100.0/total_names))
        self.logger.info(' - modified generic names: {:,} ({:.2f}%)'.format(
                        modified_generic_names,
                        modified_generic_names*100.0/total_names))
        self.logger.info(' - modified specific names: {:,} ({:.2f}%)'.format(
                        modified_specific_names,
                        modified_specific_names*100.0/total_names))
        self.logger.info(' - lost species clusters: {:,} ({:.2f}%)'.format(
                        lost_sp_cluster,
                        lost_sp_cluster*100.0/total_names))

        fout.close()
        
    def run(self,
                final_taxonomy,
                manual_sp_names,
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
                dsmz_bacnames_file,
                ground_truth_test_cases):
        """Validate final species names."""
        
        # read manually-curated taxonomy
        self.logger.info('Parsing final taxonomy.')
        final_taxonomy = Taxonomy().read(final_taxonomy, use_canonical_gid=True)
        self.logger.info(' - identified taxonomy strings for {:,} genomes.'.format(
                            len(final_taxonomy)))

        # get priority of species and generic (genus) names
        sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger,
                                                    genus_priority_ledger,
                                                    dsmz_bacnames_file)
                            
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
        
        # sanity check
        self.logger.info('Checking that previous and current genomes have same GTDB assignments.')
        dup_genomes = 0
        for gid in set(cur_genomes).intersection(prev_genomes):
            if cur_genomes[gid].gtdb_taxa != prev_genomes[gid].gtdb_taxa:
                self.logger.warning('Current GTDB assignments differ from previous assignment: {}, {}, {}'.format(
                                        gid, cur_genomes[gid].gtdb_taxa, prev_genomes[gid].gtdb_taxa))
            dup_genomes += 1
        self.logger.info(' - validate GTDB assignments for {:,} genomes in common.'.format(dup_genomes))
        
        # read named GTDB species clusters
        self.logger.info('Reading GTDB species clusters.')
        cur_clusters, rep_radius = read_clusters(gtdb_clusters_file)
        self.logger.info(' - identified {:,} clusters spanning {:,} genomes.'.format(
                            len(cur_clusters),
                            sum([len(gids) + 1 for gids in cur_clusters.values()])))
        assert len(set(final_taxonomy) - set(cur_clusters)) == 0
        
        # get species clusters with reassigned representative
        self.logger.info('Identifying species clusters with reassigned representatives.')
        reassigned_rid_mapping, prev_to_new_rids = prev_genomes.sp_clusters.reassigned_rids(cur_clusters)
        self.logger.info(' - identified {:,} species clusters with a reassigned representative.'.format(
                            len(reassigned_rid_mapping)))
                            
        # get polyphyletic suffixes for GTDB taxa
        self.taxon_suffix_manager = TaxonSuffixManager()
        
        # validate properly formatted species names
        self.logger.info('Performing basic validation of taxonomy.')
        Taxonomy().validate(final_taxonomy,
                            check_prefixes=True, 
                            check_ranks=True, 
                            check_hierarchy=True, 
                            check_species=True, 
                            check_group_names=True,
                            check_duplicate_names=True,
                            check_capitalization=True,
                            report_errors=True)
                            
        # validate genus and generic names are the same
        self.logger.info('Validating that genus and generic names are the same.')
        self.validate_genus_generic_match(final_taxonomy)
        
        # validate type species
        self.logger.info('Validating type species.')
        self.validate_type_species(final_taxonomy, cur_genomes, sp_priority_mngr)
        
        # validate type strain
        self.logger.info('Validating type strains.')
        self.validate_type_strain(final_taxonomy, cur_genomes)
        
        # validate manually-curated species names
        self.logger.info('Validating manually-curated species names.')
        self.validate_manually_curated_species(final_taxonomy, manual_sp_names)
        
        # validate ground truth test cases
        self.logger.info('Validating ground truth test cases.')
        self.validate_ground_truth(final_taxonomy, ground_truth_test_cases)
        
        # validate placeholder suffixes
        self.logger.info('Validating placeholder suffixes.')
        self.validate_placeholder_suffixes(final_taxonomy, cur_genomes, reassigned_rid_mapping)
        
        # validate GTDB species assignments without a polyphyletic suffixes
        self.logger.info('Validating species assignments without a placeholder suffix .')
        self.validate_no_suffix(final_taxonomy, cur_genomes, reassigned_rid_mapping)
        
        # validate that no forbidden specific names are used
        self.logger.info('Validating that no forbidden specific names were used.')
        self.validate_forbidden_specific_names(final_taxonomy)
       
        # validate that Latinized specific names are never used with non-suffixed, placeholder generic names
        self.logger.info('Validating that Latinized specific names are not used with non-suffixed, placeholder generic names.')
        self.validate_latin_sp_placeholder_generic(final_taxonomy)
        
        # create table with modified GTDB species names for manual inspection
        self.logger.info('Creating table with modified GTDB species names.')
        self.report_modified_species_names(final_taxonomy, prev_genomes, prev_to_new_rids)
