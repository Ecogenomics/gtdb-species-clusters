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
from collections import defaultdict, Counter

import dendropy 

from numpy import (mean as np_mean, 
                    std as np_std,
                    median as np_median)

from Levenshtein import (distance as lvn_distance)

from biolib.taxonomy import Taxonomy
from biolib.newick import parse_label

from gtdb_species_clusters.genomes import Genomes

from gtdb_species_clusters.taxon_suffix_manager import TaxonSuffixManager
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters.ncbi_species_manager import NCBI_SpeciesManager
from gtdb_species_clusters.type_genome_utils import (read_clusters)
from gtdb_species_clusters.taxon_utils import (generic_name,
                                                specific_epithet,
                                                canonical_taxon,
                                                canonical_species,
                                                taxon_suffix,
                                                is_placeholder_taxon,
                                                is_placeholder_sp_epithet,
                                                is_alphanumeric_taxon,
                                                is_suffixed_taxon,
                                                test_same_epithet)

class PMC_Validation(object):
    """Validate final species names."""

    def __init__(self, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
    def ncbi_specific_name_support(self, rid, gtdb_specific, cur_genomes, cur_clusters):
        """Validate assignment of specific name based on NCBI assignment of genomes in species cluster."""
        
        assert rid in cur_clusters[rid]

        support = 0
        total_ncbi_sp = 0
        for gid in cur_clusters[rid]:
            ncbi_sp = cur_genomes[gid].ncbi_taxa.species
            ncbi_specific = specific_epithet(ncbi_sp)
            
            if ncbi_sp != 's__':
                total_ncbi_sp += 1
            
                if ncbi_specific == gtdb_specific:
                    support += 1
            
        if total_ncbi_sp == 0:
            return 0
            
        return float(support)/total_ncbi_sp
        
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
            'Identified {:,} genomes with unmatched genus and generic names (Proposed GTDB genus, GTDB species):'.format(len(invalid_names)))
        
    def validate_type_species(self, final_taxonomy, cur_genomes, sp_priority_mngr):
        """Validate generic name of type species."""
        
        # cases reviewed and accepted by MC
        # (https://uq-my.sharepoint.com/:x:/g/personal/uqpchaum_uq_edu_au/EU5Elv0IHkVBr2NK1iJgzTsBqEvlXJrYwrtsRJjss7WmeQ?e=NBoPia)
        accepted_violations = set(['G900167455',
                                    'G900109875'])

        invalid_type_sp = {}
        validated_count = 0
        for gid in final_taxonomy:
            if gid in accepted_violations:
                continue
                
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
            'Identified {:,} invalid type species genomes (Proposed GTDB species, NCBI genus):'.format(len(invalid_type_sp)))

    def validate_type_strain(self, final_taxonomy, cur_genomes):
        """Validate specific name of type strains."""
        
        # cases reviewed and accepted by MC or DHP
        accepted_violations = set(['G002964605', # Vallitalea okinawensis S15 is type strain
                                    'G003149745', # Flagellimonas aquimarina is basonym of M. koreensis
                                    ])

        # get all GTDB species represented by a type strain:
        gtdb_type_species = set()
        for rid in final_taxonomy:
            if cur_genomes[rid].is_effective_type_strain():
                gtdb_type_species.add(final_taxonomy[rid][Taxonomy.SPECIES_INDEX])

        # check that species clusters represented by type strain genomes
        # have the specific name of the type strain
        invalid_type_strain = {}
        validated_count = 0
        for rid in final_taxonomy:
            if rid in accepted_violations:
                continue
                
            if cur_genomes[rid].is_effective_type_strain():
                gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
                gtdb_generic = generic_name(gtdb_species)
                gtdb_specific = specific_epithet(gtdb_species)

                ncbi_species = cur_genomes[rid].ncbi_taxa.species
                ncbi_generic = generic_name(ncbi_species)
                ncbi_specific = specific_epithet(ncbi_species)
                
                # check if genome is a valid genus transfer into a genus
                # that already contains a species with the specific
                # name which results in a polyphyletic suffix being required
                # e.g. G002240355 is Prauserella marina at NCBI and is
                # transferred into Saccharomonospora under the GTDB. However,
                # Saccharomonospora marina already exists so this genome
                # needs to be S. marina_A.
                if (is_placeholder_taxon(gtdb_species) 
                    and gtdb_generic != ncbi_generic
                    and canonical_species(gtdb_species) in gtdb_type_species):
                    validated_count += 1
                elif not test_same_epithet(gtdb_specific, ncbi_specific):
                    invalid_type_strain[rid] = (gtdb_species, ncbi_species)
                else:
                    validated_count += 1

        self.report_validation(
            invalid_type_strain, 
            validated_count, 
            'Identified {:,} invalid type strain genomes (Proposed GTDB species, NCBI species):'.format(len(invalid_type_strain)))
            
    def validate_manually_curated_species(self, final_taxonomy, manual_sp_names, pmc_custom_species):
        """Validate final assignment of species that were set by manual curation."""
        
        # get expected species names with post-manual curation assignments 
        # taken precedence
        mc_species = {}
        with open(manual_sp_names) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                gid = tokens[0]
                mc_species[gid] = tokens[2]
                
        with open(pmc_custom_species) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                gid = tokens[0]
                mc_species[gid] = tokens[1]
        
        # read species names explicitly set via manual curation
        invalid_mc = {}
        validated_count = 0
        for gid, mc_sp in mc_species.items():
            if final_taxonomy[gid][Taxonomy.SPECIES_INDEX] != mc_sp:
                invalid_mc[gid] = (final_taxonomy[gid][Taxonomy.SPECIES_INDEX], mc_sp)
            else:
                validated_count += 1

        self.report_validation(
            invalid_mc, 
            validated_count, 
            'Identified {:,} genomes that did not match manual curation (Proposed GTDB species, manually curated species):'.format(len(invalid_mc)))
            
        return mc_species
            
    def validate_mannually_curated_specifix_suffix(self, final_taxonomy, cur_genomes, specific_epithet_ledger):
        """Validate manually-curated specific name suffix changes resulting from genus transfers."""
        
        # read manually curated suffix changes
        curated_epithet_changes = defaultdict(lambda: {})
        with open(specific_epithet_ledger) as f:
            header = f.readline().strip().split('\t')
            
            generic_idx = header.index('GTDB generic')
            sp_corr_idx = header.index('GTDB specific_corrected')
            sp_idx = header.index('GTDB specific')
            ncbi_sp_idx = header.index('NCBI specific')
            
            for line in f:
                tokens = [t.strip() for t in line.strip().split('\t')]
                
                gtdb_generic = tokens[generic_idx]
                specific_corr = tokens[sp_corr_idx]
                gtdb_specific = tokens[sp_idx]
                ncbi_specific = tokens[ncbi_sp_idx]
                
                if specific_corr:
                    curated_epithet_changes[gtdb_generic][ncbi_specific] = specific_corr
                else:
                    curated_epithet_changes[gtdb_generic][ncbi_specific] = gtdb_specific

        # validate names
        invalid_suffix = {}
        validated_count = 0
        for rid in final_taxonomy:
            gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            gtdb_generic = generic_name(gtdb_species)
            gtdb_specific = specific_epithet(gtdb_species)
            gtdb_canonical_specific = canonical_taxon(gtdb_specific)
            
            ncbi_species = cur_genomes[rid].ncbi_taxa.species
            ncbi_specific = specific_epithet(ncbi_species)
            
            if gtdb_generic in curated_epithet_changes:
                if ncbi_specific in curated_epithet_changes[gtdb_generic]:
                    mc_specific = curated_epithet_changes[gtdb_generic][ncbi_specific]
                    if (test_same_epithet(gtdb_canonical_specific, curated_epithet_changes[gtdb_generic][ncbi_specific])
                        and gtdb_canonical_specific != curated_epithet_changes[gtdb_generic][ncbi_specific]):
                        # epithets are highly similar, but not identical indicating a likely curation issue
                        invalid_suffix[rid] = (gtdb_species, curated_epithet_changes[gtdb_generic][ncbi_specific])
                    else:
                        validated_count += 1

        self.report_validation(
            invalid_suffix, 
            validated_count, 
            'Identified {:,} genomes that did not match manually curated changes to the suffix of specific name (Proposed GTDB species, manually curated specific name):'.format(len(invalid_suffix)))
        
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
            'Identified {:,} genomes that did not match ground truth assignments (Proposed GTDB species, ground truth species):'.format(len(invalid_gt)))

    def validate_synonyms(self, final_taxonomy, cur_genomes, ncbi_synonyms):
        """Validate that GTDB synonyms are respected."""
        
        invalid_synonym = {}
        validated_count = 0
        for rid in final_taxonomy:
            ncbi_species = cur_genomes[rid].ncbi_taxa.species

            if ncbi_species in ncbi_synonyms:
                ncbi_specific = specific_epithet(ncbi_species)
                
                gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
                gtdb_specific = specific_epithet(gtdb_species)
                
                if (not is_placeholder_sp_epithet(gtdb_specific)
                    and test_same_epithet(ncbi_specific, gtdb_specific)):
                    # GTDB specific name is the specific name of a GTDB synonym so is in error
                    invalid_synonym[rid] = (gtdb_species, 
                                                ncbi_species, 
                                                '{} is a synonym of {}'.format(ncbi_species, ncbi_synonyms[ncbi_species]))
                else:
                    validated_count += 1

        self.report_validation(
            invalid_synonym, 
            validated_count, 
            'Identified {:,} genomes that do not respect GTDB synonyms (Proposed GTDB species, NCBI species, NCBI synonym):'.format(len(invalid_synonym)))
            
    def validate_single_type_species(self, final_taxonomy, cur_genomes):
        """Validating that each GTDB genus has a single type species cluster."""
        
        # get type species in each genus
        genus_gtdb_type_species = defaultdict(set)
        for rid in final_taxonomy:
            if cur_genomes[rid].is_gtdb_type_species():
                gtdb_genus = final_taxonomy[rid][Taxonomy.GENUS_INDEX]
                gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
                ncbi_genus = cur_genomes[rid].ncbi_taxa.genus
                ncbi_species = cur_genomes[rid].ncbi_taxa.species
                
                if gtdb_genus == ncbi_genus:
                    genus_gtdb_type_species[gtdb_genus].add((gtdb_species, ncbi_species))
                
        # report genus with multiple type species
        invalid_num_type_species = {}
        validated_count = 0
        for genus, type_species in genus_gtdb_type_species.items():
            if len(type_species) > 1:
                invalid_num_type_species[genus] = ', '.join([str(gtdb_ncbi_sp) for gtdb_ncbi_sp in type_species])
            else:
                validated_count += 1

        self.report_validation(
            invalid_num_type_species, 
            validated_count, 
            'Identified {:,} genera with multiple type species (GTDB genera, (GTDB species, NCBI species)):'.format(len(invalid_num_type_species)))
            
    def validate_prev_specific_placeholder_suffixes(self, final_taxonomy, cur_genomes, misclassified_gids, reassigned_rid_mapping):
        """Validate that species clusters have same specific name placeholder suffixes as previously assigned."""

        invalid_suffix = {}
        validated_count = 0
        suffix_count = 0
        for rid in final_taxonomy:
            if rid in self.mc_gids:
                continue
                
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
                
                if rid in misclassified_gids:
                    # specific name should be a placeholder
                    if is_placeholder_sp_epithet(gtdb_specific):
                        validated_count += 1
                    else:
                        invalid_suffix[rid] = (gtdb_species, prev_species)
                elif prev_species == 's__':
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
            'Identified {:,} genomes with an suffixed specific name that disagrees with previous GTDB assignment (Proposed GTDB species, Previous GTDB species):'.format(len(invalid_suffix)))

    def validate_specific_placeholder_suffixes(self,
                                                final_taxonomy, 
                                                cur_genomes, 
                                                misclassified_gids,
                                                ncbi_synonyms):
        """Validate suffix assignment rules for specific names."""
        
        # identify cases where a specific name is used multiple times in a genus
        rid_synonyms = {}
        sp_name_map = defaultdict(lambda: defaultdict(list))
        for rid in final_taxonomy:
            gtdb_genus = final_taxonomy[rid][Taxonomy.GENUS_INDEX]
            gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            gtdb_specific = specific_epithet(gtdb_species)
            canonical_specific = canonical_taxon(gtdb_specific)
            
            ncbi_species = cur_genomes[rid].ncbi_taxa.species
            if ncbi_species in ncbi_synonyms:
                # change specific name for purposes of establishing suffix rules
                ncbi_species = ncbi_synonyms[ncbi_species]
                canonical_specific = specific_epithet(ncbi_species)
                rid_synonyms[rid] = ncbi_species

            if not is_placeholder_sp_epithet(canonical_specific):
                sp_name_map[gtdb_genus][canonical_specific].append(rid)
            
        # validate suffixes on specific names follow assignment rules:
        # 1) single occurrence of name should be unsuffixed unless it is a misclassification
        # 2) multiple occurrences should all be suffixed, except type strain genome
        invalid_suffix = {}
        validated_count = 0
        for gtdb_genus in sp_name_map:
            for canonical_specific, rids in sp_name_map[gtdb_genus].items():
                if len(rids) == 1:
                    rid = rids[0]
                    if rid in self.mc_gids:
                        continue
                        
                    ncbi_species = cur_genomes[rid].ncbi_taxa.species
                    gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
                    specific_name = specific_epithet(gtdb_species)
                    if (taxon_suffix(specific_name) is not None
                        and rid not in misclassified_gids):
                        # violation rule 1 as single instance of specific epithet should
                        # be unsuffixed unless it is a misclassification
                        invalid_suffix[rid] = (gtdb_species, 
                                                ncbi_species,
                                                'Suffixed single instance of specific name; genome is not misclassified')
                    else:
                        # in general, a single instance of a specific epithet should be unsuffixed
                        validated_count += 1
                else:
                    type_strain_genomes = [rid for rid in rids if cur_genomes[rid].is_effective_type_strain()]
                    for rid in rids:
                        if rid in self.mc_gids:
                            continue
                            
                        ncbi_species = cur_genomes[rid].ncbi_taxa.species
                        gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
                        specific_name = specific_epithet(gtdb_species)
                        if (taxon_suffix(specific_name) is None
                            and not cur_genomes[rid].is_effective_type_strain()):
                            # violates rule 2 as only type strain genomes should lack a suffix
                            note = 'Specific name used multiple times; non-type strain genome lacks suffix'
                            if rid in rid_synonyms:
                                note += '; NCBI species assignment is a synonym of {}'.format(rid_synonyms[rid])
                            invalid_suffix[rid] = (gtdb_species, 
                                                    ncbi_species,
                                                    note)
                        elif (taxon_suffix(specific_name) is not None
                                and cur_genomes[rid].is_effective_type_strain()
                                and len(type_strain_genomes) <= 1):
                            # violates rule 2 as type strain genomes should not have a suffix unless
                            # there are multiple type strains in species clusters 
                            note = 'Specific name used multiple times; type strain genome is suffixed'
                            if rid in rid_synonyms:
                                note += '; NCBI species assignment is a synonym of {}'.format(rid_synonyms[rid])
                            invalid_suffix[rid] = (gtdb_species, 
                                                    ncbi_species,
                                                    note)
                        else:
                            validated_count += 1

        self.report_validation(
            invalid_suffix, 
            validated_count, 
            'Identified {:,} genomes which violate the suffixing rules for specific names (Proposed GTDB species, NCBI species, Note):'.format(len(invalid_suffix)))
            
        sys.exit(-1) #***
            
    def validate_genera_placeholder_suffixes(self, final_taxonomy, cur_genomes, reassigned_rid_mapping):
        """Validate placeholder suffixes of genera."""
        
        invalid_suffix = {}
        validated_count = 0
        suffix_count = 0
        for rid in final_taxonomy:
            gtdb_genus = final_taxonomy[rid][Taxonomy.GENUS_INDEX]
            gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            gtdb_generic = generic_name(gtdb_species)
            
            suffix = taxon_suffix(gtdb_generic)
            if suffix is None:
                continue
                
            suffix_count += 1
            
            prev_rid = reassigned_rid_mapping.get(rid, rid)
            if prev_rid in cur_genomes:
                prev_species = cur_genomes[prev_rid].gtdb_taxa.species
                prev_generic = generic_name(prev_species)
                prev_suffix = taxon_suffix(prev_generic)
                
                if prev_species == 's__':
                    # might fall under any GTDB genus, so difficult to do any
                    # sort of validation
                    continue
                elif prev_suffix is None:
                    # this should only occur if a genus name was previously unique,
                    # and is now associated with multiple species clusters
                    num_canonical_sp = sum([1 for gid in final_taxonomy 
                                            if canonical_taxon(final_taxonomy[gid][Taxonomy.GENUS_INDEX]) == canonical_taxon(gtdb_genus)])

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
            'Identified {:,} genomes with a modified placeholder suffix which should be manually verified (Proposed GTDB species, Previous GTDB species):'.format(len(invalid_suffix)))

    def validate_genera_spelling(self, final_taxonomy, cur_genomes, reassigned_rid_mapping):
        """Validate spelling of genus names."""
        
        invalid_spelling = {}
        validated_count = 0
        suffix_count = 0
        for rid in final_taxonomy:
            gtdb_genus = final_taxonomy[rid][Taxonomy.GENUS_INDEX]
            canonical_genus = canonical_taxon(gtdb_genus)
            
            prev_gtdb_genus = cur_genomes[rid].gtdb_taxa.genus
            canonical_prev_genus  = canonical_taxon(prev_gtdb_genus)
            
            if (canonical_genus != canonical_prev_genus
                and lvn_distance(canonical_genus, canonical_prev_genus) <= 2
                and not (canonical_genus.startswith('g__UBA') and canonical_prev_genus.startswith('g__UBA'))):
                # Names are highly similar, but not identical so may be a typo. 
                # UBA genera are ignored since these appear similar to each other so are often flagged as false positives,
                # given that valid genus transfers do occur.
                invalid_spelling[gtdb_genus] = prev_gtdb_genus
            else:
                validated_count += 1
                
        self.report_validation(
            invalid_spelling, 
            validated_count, 
            'Identified {:,} genera with potential spelling error (Proposed GTDB genus, Previous GTDB genus):'.format(len(invalid_spelling)))
        
    def validate_specific_no_suffix(self, final_taxonomy, cur_genomes, reassigned_rid_mapping):
        """Validate GTDB specific name assignments without a polyphyletic suffixes."""
        
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
            'Identified {:,} genomes which should have a polyphyletic suffix on specific name (Proposed GTDB species, effective type material, GTDB canonical species):'.format(len(invalid_names)))
        
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
    
    def validate_latin_sp_placeholder_generic(self, final_taxonomy, cur_genomes):
        """Validate that Latinized specific names are never used with placeholder generic names."""

        invalid_sp_name = {}
        validated_count = 0
        for gid, taxa in final_taxonomy.items():
            gtdb_species = taxa[Taxonomy.SPECIES_INDEX]
            canonical_generic = canonical_taxon(generic_name(gtdb_species))
            specific = specific_epithet(gtdb_species)

            if (is_placeholder_taxon('g__' + canonical_generic) and 
                not is_placeholder_sp_epithet(specific)):
                    ncbi_species = cur_genomes[gid].ncbi_taxa.species
                    invalid_sp_name[gid] = (gtdb_species, ncbi_species)
            else:
                validated_count += 1

        self.report_validation(
            invalid_sp_name, 
            validated_count, 
            'Identified {:,} genomes with a non-suffixed, placeholder generic name, but Latin specific name (Proposed GTDB species, NCBI species):'.format(len(invalid_sp_name)))

    def validate_missing_ncbi_genera(self, final_taxonomy, cur_genomes, curation_tree):
        """Validate missing NCBI genera."""

        # get maps between genome IDs and genera to nodes in the curation tree
        leaf_map = {}
        for leaf in curation_tree.leaf_node_iter():
            leaf_map[leaf.taxon.label] = leaf
            
        genus_map = {}
        red_genus = {}
        for node in curation_tree.preorder_node_iter():
            _support, taxon_name, _auxiliary_info = parse_label(node.label)
            if taxon_name:
                taxa = [t.strip() for t in taxon_name.split(';')]
                for taxon in taxa:
                    if taxon.startswith('g__'):
                        genus_map[taxon] = node
                        
                if taxa[-1].startswith('g__'):
                    red_genus[taxa[-1]] = node.root_distance

        median_red_genera = np_median(list(red_genus.values()))
        self.logger.info(' - median RED of genera: {:.3f}'.format(median_red_genera))

        # get all genera present in the GTDB
        gtdb_genera = set([final_taxonomy[rid][Taxonomy.GENUS_INDEX] for rid in final_taxonomy])
        
        # get NCBI genera not present in the GTDB
        missing_ncbi_genera = defaultdict(list)
        for rid in final_taxonomy:
            ncbi_genus = cur_genomes[rid].ncbi_taxa.genus
            if ncbi_genus != 'g__' and ncbi_genus not in gtdb_genera:
                missing_ncbi_genera[ncbi_genus].append(rid)
        self.logger.info(' - found {:,} NCBI genera not represented in the GTDB.'.format(len(missing_ncbi_genera)))
                
        invalid_missing_ncbi_genera = {}
        validated_count = 0
        for ncbi_genus, rids in missing_ncbi_genera.items():
            found = False
            for rid in rids:
                gtdb_genus = final_taxonomy[rid][Taxonomy.GENUS_INDEX]
                if canonical_taxon(gtdb_genus) == ncbi_genus:
                    found = True
                    break

            if found:
                validated_count += 1
            else:
                within_gtdb_genus = False
                for rid in rids:
                    gtdb_genus = final_taxonomy[rid][Taxonomy.GENUS_INDEX]
                    root_genus = genus_map[gtdb_genus]
                    if root_genus.root_distance < median_red_genera:
                        # RED is to the left (smaller) of the median 
                        # RED for genera so the genus could be split
                        leaf_node = leaf_map[rid]
                        cur_node = leaf_node.parent_node
                        while cur_node != root_genus:
                            for leaf in cur_node.leaf_iter():
                                if leaf.taxon.label.startswith('D-'):
                                    continue
                                    
                                leaf_ncbi_genus = cur_genomes[leaf.taxon.label].ncbi_taxa.genus
                                if leaf_ncbi_genus == gtdb_genus:
                                    within_gtdb_genus = True

                            cur_node = cur_node.parent_node
                    else:
                        # flag genus as not suitable for splitting
                        within_gtdb_genus = True
                    
                if not within_gtdb_genus:
                    gtdb_genera = ['{}: {} (RED={:.2f})'.format(rid, final_taxonomy[rid][Taxonomy.GENUS_INDEX], red_genus[final_taxonomy[rid][Taxonomy.GENUS_INDEX]]) for rid in rids]
                    invalid_missing_ncbi_genera[ncbi_genus] = ', '.join(gtdb_genera)
                else:
                    validated_count += 1
                
        self.report_validation(
            invalid_missing_ncbi_genera, 
            validated_count, 
            'Identified {:,} NCBI genera without justification for being absent in GTDB (NCBI genus, GTDB genera):'.format(len(invalid_missing_ncbi_genera)))

    def validate_no_latin_genus_name(self, final_taxonomy, cur_genomes):
        """Validate placeholder genus names lacking potential Latin names."""
        
        # get genomes comprising GTDB and NCBI genera
        gtdb_placeholder_genus_gids = defaultdict(list)
        ncbi_genus_gids = defaultdict(set)
        gtdb_genera = set()
        for rid in final_taxonomy:
            gtdb_genus = final_taxonomy[rid][Taxonomy.GENUS_INDEX]
            gtdb_genera.add(gtdb_genus)
            
            if is_placeholder_taxon(gtdb_genus):
                gtdb_placeholder_genus_gids[gtdb_genus].append(rid)
                
            ncbi_genus = cur_genomes[rid].ncbi_taxa.genus
            if ncbi_genus != 'g__':
                ncbi_genus_gids[ncbi_genus].add(rid)
            
        # check if genus could be assigned a Latin name
        fout = open(os.path.join(self.output_dir, 'latin_names_for_gtdb_genera.tsv'), 'w')
        fout.write('Proposed GTDB genus\tNCBI genus\tPropostion of NCBI genus genomes\tNote\n')
        invalid_placeholder = {}
        validated_count = 0
        for gtdb_genus, rids in gtdb_placeholder_genus_gids.items():
            ncbi_genera = []
            for rid in rids:
                ncbi_genus = cur_genomes[rid].ncbi_taxa.genus
                if ncbi_genus != 'g__':
                    ncbi_genera.append(ncbi_genus)
                    
            if not ncbi_genera:
                validated_count += 1
            else:
                top_ncbi_genera, count = Counter(ncbi_genera).most_common(1)[0]
                perc_in_gtdb_genus = count*100.0/len(ncbi_genus_gids[top_ncbi_genera])
                
                # check if genomes from NCBI genus not in this GTDB genus, contain
                # the type species of the genus
                missing_type_species = False
                for rid in ncbi_genus_gids[top_ncbi_genera] - set(rids):
                    if cur_genomes[rid].is_gtdb_type_species():
                        missing_type_species = True
                        
                # check if there is an NCBI species with the NCBI genus name
                top_ncbi_sp_count = 0
                for rid in rids:
                    ncbi_sp = cur_genomes[rid].ncbi_taxa.species
                    if 'g__' + generic_name(ncbi_sp) == top_ncbi_genera:
                        top_ncbi_sp_count += 1

                if (not missing_type_species 
                    and perc_in_gtdb_genus >= 50
                    and canonical_taxon(gtdb_genus) != top_ncbi_genera
                    and (top_ncbi_genera not in gtdb_genera or top_ncbi_sp_count > 1)):
                    note = []
                    if top_ncbi_genera in gtdb_genera:
                        note.append('NCBI genus used elsewhere in GTDB')
                    if top_ncbi_sp_count:
                        note.append('Contains {:,} genomes with Latin NCBI species name.'.format(top_ncbi_sp_count))
                        
                    perc_str = '{} of {} ({:.1f}%)'.format(count, 
                                                            len(ncbi_genus_gids[top_ncbi_genera]), 
                                                            perc_in_gtdb_genus)
                    invalid_placeholder[gtdb_genus] = (top_ncbi_genera, 
                                                        perc_str,
                                                        '; '.join(note))
                                                        
                    fout.write('{}\t{}\t{}\t{}\n'.format(gtdb_genus,
                                                            top_ncbi_genera, 
                                                            perc_str,
                                                            '; '.join(note)))

                else:
                    validated_count += 1

        fout.close()

        self.report_validation(
            invalid_placeholder, 
            validated_count, 
            'Identified {:,} GTDB placeholder genera with candidate Latin names (Proposed GTDB genus, NCBI genus, Proportion of NCBI genus genomes):'.format(len(invalid_placeholder)))
            
    def validate_genus_placeholder_change(self, final_taxonomy, cur_genomes, prev_genomes):
        """Validate that genus names are not a replacement of one placeholder for another."""

        # get genomes that have changed from one placeholder to another placeholder
        changed_placeholder = defaultdict(lambda: defaultdict(list))
        for rid in final_taxonomy:
            if rid not in prev_genomes:
                continue
                
            prev_genus = prev_genomes[rid].gtdb_taxa.genus
            cur_genus = final_taxonomy[rid][Taxonomy.GENUS_INDEX]

            if (is_placeholder_taxon(prev_genus) 
                and is_placeholder_taxon(cur_genus)
                and prev_genus != cur_genus):
                changed_placeholder[prev_genus][cur_genus].append(rid)
                
        # get all genera in current and previous GTDB release
        cur_gtdb_genera = set([final_taxonomy[rid][Taxonomy.GENUS_INDEX] for rid in final_taxonomy])
        prev_gtdb_genera = set([prev_genomes[rid].gtdb_taxa.genus for rid in prev_genomes.sp_clusters])
        self.logger.info(' - identified {:,} current and {:,} previous GTDB genera.'.format(
                        len(cur_gtdb_genera),
                        len(prev_gtdb_genera)))

        # identified changed placeholders that look like they may be invalid
        invalid_placeholder = {}
        validated_count = 0
        for prev_genus in changed_placeholder:
            for cur_genus, rids in changed_placeholder[prev_genus].items():
                likely_valid_change = False
                if cur_genus in prev_gtdb_genera:
                    # no guarantee, but likely indicates species have simply moved
                    # into previous genus or two previous genera have been merged
                    likely_valid_change = True
                elif prev_genus in cur_gtdb_genera:
                    # no guarantee, but likely indicates species have simply moved
                    # out of previous genus and thus requires a new name
                    likely_valid_change = True
                elif is_suffixed_taxon(prev_genus) or is_suffixed_taxon(cur_genus):
                    # this sort of change is allowed since changing g__Riesia_B to g__NZ1215
                    # indicates the genus is no longer associated with Riesia; conversely,
                    # changing g__NZ1215 to g__Riesia_B indicates an assertion on our part 
                    # that this lineage is associated with Riesia at NCBI
                    likely_valid_change = True
                else:
                    for rid in rids:
                        ncbi_genus = 'g__' + generic_name(cur_genomes[rid].ncbi_taxa.species)
                        if canonical_taxon(cur_genus) == ncbi_genus:
                            # necessarily modified to reflect NCBI species assignment
                            # (e.g., FW-11 to Sphingomonas_H as representative is now S. oleivorans at NCBI)
                            likely_valid_change = True

                if likely_valid_change:
                    validated_count += 1
                else:
                    invalid_placeholder[prev_genus] = (cur_genus, len(rids))
                
        self.report_validation(
            invalid_placeholder, 
            validated_count, 
            'Identified {:,} genomes with invalid change to placeholder genus/generic name (Previous GTDB genus, Proposed GTDB genus, No. species changed):'.format(len(invalid_placeholder)))

    def validate_misclassified_genomes(self, final_taxonomy, cur_genomes, ncbi_synonyms):
        """Validating species names of genomes misclassified at NCBI."""
        
        # get NCBI species anchored by a type strain genome, or
        # being placed in a GTDB genus matching the NCBI generic name
        ncbi_anchored_species = {}
        for rid in final_taxonomy:
            ncbi_species = cur_genomes[rid].ncbi_taxa.species
            ncbi_species = ncbi_synonyms.get(ncbi_species, ncbi_species)
            
            if ncbi_species != 's__':
                gtdb_genus = final_taxonomy[rid][Taxonomy.GENUS_INDEX]
                ncbi_generic = generic_name(ncbi_species)
                anchored_sp = (cur_genomes[rid].is_effective_type_strain()
                                or gtdb_genus == 'g__' + ncbi_generic)
                                
                if anchored_sp:
                    gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
                    
                    if ncbi_species not in ncbi_anchored_species:
                        ncbi_anchored_species[ncbi_species] = (gtdb_species, rid)
                    elif cur_genomes[rid].is_effective_type_strain():
                        # favor having the GTDB species which is a type strain genome
                        ncbi_anchored_species[ncbi_species] = (gtdb_species, rid)

        # identify misclassified genomes and ensure they have a placeholder specific name
        invalid_misclassified = {}
        validated_count = 0
        misclassified_gids = set()
        for rid in final_taxonomy:
            if rid in self.mc_gids:
                continue
                
            ncbi_sp = cur_genomes[rid].ncbi_taxa.species
            ncbi_sp = ncbi_synonyms.get(ncbi_sp, ncbi_sp)
            
            gtdb_genus = final_taxonomy[rid][Taxonomy.GENUS_INDEX]
            gtdb_sp = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]

            if ncbi_sp in ncbi_anchored_species:
                gtdb_type_sp, type_gid = ncbi_anchored_species[ncbi_sp]
                
                if generic_name(gtdb_sp) != generic_name(gtdb_type_sp):
                    # genome is misclassified so should have a placeholder specific name
                    misclassified_gids.add(rid)
                    
                    gtdb_specific = specific_epithet(gtdb_sp)
                    if not is_placeholder_sp_epithet(gtdb_specific):
                        invalid_misclassified[rid] = (gtdb_sp, ncbi_sp)
                    else:
                        validated_count += 1

        self.report_validation(
            invalid_misclassified, 
            validated_count, 
            'Identified {:,} misclassified genomes with an invalid specific name (Proposed GTDB species, NCBI classification):'.format(len(invalid_misclassified)))
            
        return misclassified_gids
            
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
                final_scaled_tree,
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
                if cur_genomes[gid].gtdb_taxa.domain == prev_genomes[gid].gtdb_taxa.domain:
                    # there are known cases of genomes being transferred between domain which 
                    # necessitates a different GTDB assignment
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
        
        # validate manually-curated species names
        self.logger.info('Validating manually-curated species names.')
        self.mc_gids = self.validate_manually_curated_species(final_taxonomy, manual_sp_names, pmc_custom_species)
        
        # establish state of NCBI species
        ncbi_species_mngr = NCBI_SpeciesManager(cur_genomes, cur_clusters, self.output_dir)
        
        # get list of synonyms in order to restrict usage of species names
        self.logger.info('Reading GTDB synonyms.')
        ncbi_synonyms = ncbi_species_mngr.parse_synonyms_table(synonym_file)
        self.logger.info(' - identified {:,} synonyms from {:,} distinct species.'.format(
                            len(ncbi_synonyms),
                            len(set(ncbi_synonyms.values()))))

        # identified genomes with misclassified species assignments at NCBI
        self.logger.info('Identify genomes with misclassified NCBI species assignments.')
        misclassified_gids = ncbi_species_mngr.identify_misclassified_genomes(ncbi_synonyms, self.mc_gids)
        
        
        ncbi_sp_gids = defaultdict(list)
        for gid in cur_genomes:
            ncbi_species = cur_genomes[gid].ncbi_taxa.species
            if ncbi_species != 's__':
                ncbi_sp_gids[ncbi_species].append(gid)
        
        num_unclassified_ncbi_rep = 0
        for rid, cids in cur_clusters.items():
            rep_ncbi_species = cur_genomes[rid].ncbi_taxa.species
            
            if rep_ncbi_species == 's__':
                cluster_ncbi_species = []
                for cid in cids:
                    if cid == rid:
                        continue
                        
                    if cid in misclassified_gids:
                        continue
                        
                    ncbi_species = cur_genomes[cid].ncbi_taxa.species
                    if ncbi_species != 's__':
                        cluster_ncbi_species.append(ncbi_species)
                        
                if len(set(cluster_ncbi_species)) != 1:
                    continue
                    
                if len(cluster_ncbi_species) != len(ncbi_sp_gids[cluster_ncbi_species[0]]):
                    continue

                print(rid, len(cids), Counter(cluster_ncbi_species))
                num_unclassified_ncbi_rep += 1
                    
        print('num_unclassified_ncbi_rep', num_unclassified_ncbi_rep)
        
        

        # classify each Latin NCBI species as either unambiguously or ambiguously assigned
        self.logger.info('Classifying NCBI species as unambiguous, ambiguous, or synonym.')
        ncbi_species, unambiguous_ncbi_sp, ambiguous_ncbi_sp = ncbi_species_mngr.classify_ncbi_species(
                                                                    ncbi_synonyms, 
                                                                    misclassified_gids)
        
        unambiguous_type_strain_count = 0
        unambiguous_unanimous_consensus_count = 0
        for _, _, classification, _ in unambiguous_ncbi_sp.values():
            if classification == 'TYPE_STRAIN_GENOME':
                unambiguous_type_strain_count += 1
            elif classification == 'UNANIMOUS_CONSENSUS':
                unambiguous_unanimous_consensus_count += 1
            else:
                self.logger.error('Unknown classification type: {}'.format(classification))
                sys.exit(-1)
                
        self.logger.info(' - identified {:,} NCBI species.'.format(len(ncbi_species)))
        self.logger.info(' - identified {:,} ({:.2f}%) as unambiguous assignments.'.format(
                            len(unambiguous_ncbi_sp),
                            len(unambiguous_ncbi_sp) * 100.0 / len(ncbi_species)))
        self.logger.info('   - assigned {:,} ({:.2f}%) by type strain genome.'.format(
                            unambiguous_type_strain_count,
                            unambiguous_type_strain_count * 100.0 / len(ncbi_species)))
        self.logger.info('   - assigned {:,} ({:.2f}%) by unanimous consensus.'.format(
                            unambiguous_unanimous_consensus_count,
                            unambiguous_unanimous_consensus_count * 100.0 / len(ncbi_species)))
        self.logger.info(' - identified {:,} ({:.2f}%) as ambiguous assignments.'.format(
                            len(ambiguous_ncbi_sp),
                            len(ambiguous_ncbi_sp) * 100.0 / len(ncbi_species)))
        self.logger.info(' - identified {:,} ({:.2f}%) as GTDB synonyms.'.format(
                            len(ncbi_synonyms),
                            len(ncbi_synonyms) * 100.0 / len(ncbi_species)))
                            
        sys.exit(-1)

        # get species clusters with reassigned representative
        self.logger.info('Identifying species clusters with reassigned representatives.')
        reassigned_rid_mapping, prev_to_new_rids = prev_genomes.sp_clusters.reassigned_rids(cur_clusters)
        self.logger.info(' - identified {:,} species clusters with a reassigned representative.'.format(
                            len(reassigned_rid_mapping)))
                            
        # get polyphyletic suffixes for GTDB taxa
        self.taxon_suffix_manager = TaxonSuffixManager()
        
        # read scaled curation tree
        self.logger.info('Reading scaled and decorated tree.')
        curation_tree = dendropy.Tree.get_from_path(final_scaled_tree, 
                                            schema='newick', 
                                            rooting='force-rooted', 
                                            preserve_underscores=True)
        curation_tree.calc_node_root_distances()
        
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

        # validate ground truth test cases
        self.logger.info('Validating ground truth test cases.')
        self.validate_ground_truth(final_taxonomy, ground_truth_test_cases)
        
        # validate manually-curated species names
        self.logger.info('Validating manually-curated species names.')
        self.mc_gids = self.validate_manually_curated_species(final_taxonomy, manual_sp_names, pmc_custom_species)
        
        # validate manually established suffixes for specific names changed due to genus transfer
        self.logger.info('Validating manually-curated specific name suffix changes resulting from genus transfers.')
        self.validate_mannually_curated_specifix_suffix(final_taxonomy, cur_genomes, specific_epithet_ledger)
        
        # validate species names of genomes misclassified at NCBI
        self.logger.info('Validating species names of genomes misclassified at NCBI.')
        misclassified_gids = self.validate_misclassified_genomes(final_taxonomy, cur_genomes, ncbi_synonyms)

        # validate synonyms
        self.logger.info('Validating that GTDB synonyms are respected.')
        self.validate_synonyms(final_taxonomy, cur_genomes, ncbi_synonyms)
        
        # validate that each GTDB genus has a single type species
        self.logger.info('Validating that each GTDB genus has a single type species cluster.')
        self.validate_single_type_species(final_taxonomy, cur_genomes)

        # validate GTDB specific name assignments without a polyphyletic suffixes
        self.logger.info('Validating specific name assignments without a placeholder suffix.')
        self.validate_specific_no_suffix(final_taxonomy, cur_genomes, reassigned_rid_mapping)
        
        # validate that species clusters use same specific name placeholder suffixes as previously assigned
        self.logger.info('Validating specific names use same placeholder suffixes as assigned in previous release.')
        self.validate_prev_specific_placeholder_suffixes(final_taxonomy, 
                                                    cur_genomes, 
                                                    misclassified_gids,
                                                    reassigned_rid_mapping)
                                                    
        # validate suffix assignment rules for specific names
        self.logger.info('Validating specific names follow suffix assignment rules.')
        self.validate_specific_placeholder_suffixes(final_taxonomy, 
                                                    cur_genomes, 
                                                    misclassified_gids,
                                                    ncbi_synonyms)

        # validate placeholder suffixes for genera
        self.logger.info('Validating placeholder suffixes of genera.')
        self.validate_genera_placeholder_suffixes(final_taxonomy, cur_genomes, reassigned_rid_mapping)
        
        # validate spelling of genera
        self.logger.info('Validating spelling of genera.')
        self.validate_genera_spelling(final_taxonomy, cur_genomes, reassigned_rid_mapping)
        
        # validate that no forbidden specific names are used
        self.logger.info('Validating that no forbidden specific names were used.')
        self.validate_forbidden_specific_names(final_taxonomy)
       
        # validate that Latinized specific names are never used with non-suffixed, placeholder generic names
        self.logger.info('Validating that Latinized specific names are not used with non-suffixed, placeholder generic names.')
        self.validate_latin_sp_placeholder_generic(final_taxonomy, cur_genomes)
        
        # validate missing NCBI genus names
        self.logger.info('Validating missing NCBI-defined genus names.')
        self.validate_missing_ncbi_genera(final_taxonomy, cur_genomes, curation_tree)
        
        # validate placeholder genus names lack potential Latin names
        self.logger.info('Validating placeholder genus names lack potential Latin names.')
        self.validate_no_latin_genus_name(final_taxonomy, cur_genomes)
        
        # validate that genus names are not a replacement of one placeholder for another
        self.logger.info('Validating that genus names are not a change from one placeholder to another.')
        self.validate_genus_placeholder_change(final_taxonomy, cur_genomes, prev_genomes)

        # create table with modified GTDB species names for manual inspection
        self.logger.info('Creating table with modified GTDB species names.')
        self.report_modified_species_names(final_taxonomy, prev_genomes, prev_to_new_rids)
