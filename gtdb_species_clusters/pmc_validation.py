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
from itertools import product
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
from gtdb_species_clusters.genome_utils import canonical_gid
from gtdb_species_clusters.type_genome_utils import (read_clusters,
                                                        parse_updated_species_reps,
                                                        infer_prev_gtdb_reps)
from gtdb_species_clusters.taxon_utils import (generic_name,
                                                specific_epithet,
                                                canonical_taxon,
                                                canonical_species,
                                                taxon_suffix,
                                                is_placeholder_taxon,
                                                is_placeholder_sp_epithet,
                                                is_alphanumeric_taxon,
                                                is_alphanumeric_sp_epithet,
                                                is_suffixed_taxon,
                                                is_suffixed_sp_epithet,
                                                taxon_type,
                                                specific_epithet_type,
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
            if rid in self.mc_gids:
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
        
    def validate_type_strain_ledger(self, final_taxonomy, gtdb_type_strains_ledger):
        """Validating genomes in type strain ledger."""
        
        # read genomes in type strain ledger
        gtdb_type_strains = {}
        with open(gtdb_type_strains_ledger) as f:
            f.readline()
            
            for line in f:
                tokens = line.strip().split('\t')
                gid = canonical_gid(tokens[0])
                sp = tokens[1].strip()
                if not sp.startswith('s__'):
                    sp = 's__' + sp
                gtdb_type_strains[gid] = sp
                
        # validate names of genomes in type strain ledger
        invalid_type_strain = {}
        validated_count = 0
        for gid in gtdb_type_strains:
            if gid in final_taxonomy:
                gtdb_species = final_taxonomy[gid][Taxonomy.SPECIES_INDEX]
                gtdb_specific = specific_epithet(gtdb_species)
                
                ledger_specific = specific_epithet(gtdb_type_strains[gid])

                if not test_same_epithet(gtdb_specific, ledger_specific):
                    invalid_type_strain[gid] = (gtdb_species, gtdb_type_strains[gid])
                else:
                    validated_count += 1
                    
        self.report_validation(
            invalid_type_strain, 
            validated_count, 
            'Identified {:,} genomes that do not match name proposed in GTDB type strain ledger (Proposed GTDB species, Ledger species name):'.format(
                len(invalid_type_strain)))
                
        return gtdb_type_strains
                
    def validate_species_classification_ledger(self, final_taxonomy, cur_clusters, species_classification_ledger):
        """Validate genomes in GTDB type strain ledger."""
        
        # get map from genomes to their representative
        gtdb_gid_to_rid = {}
        for rid, cids in cur_clusters.items():
            for cid in cids:
                gtdb_gid_to_rid[cid] = rid
                
        # read genomes in species classification ledger
        gtdb_sp_ledger = {}
        with open(species_classification_ledger, encoding='utf-8') as f:
            f.readline()
            
            for line in f:
                tokens = line.strip().split('\t')
                gid = canonical_gid(tokens[0])
                sp = tokens[1].strip()
                if not sp.startswith('s__'):
                    sp = 's__' + sp
                    
                rid = gtdb_gid_to_rid[gid]
                gtdb_sp_ledger[rid] = sp

        # validate names of genomes in species classification ledger
        invalid_sp = {}
        validated_count = 0
        for gid in gtdb_sp_ledger:
            if gid in final_taxonomy:
                gtdb_species = final_taxonomy[gid][Taxonomy.SPECIES_INDEX]
                gtdb_specific = specific_epithet(gtdb_species)
                
                ledger_specific = specific_epithet(gtdb_sp_ledger[gid])
                
                if not test_same_epithet(gtdb_specific, ledger_specific):
                    invalid_sp[gid] = (gtdb_species, gtdb_sp_ledger[gid])
                else:
                    validated_count += 1
                    
        self.report_validation(
            invalid_sp, 
            validated_count, 
            'Identified {:,} genomes that do not match name proposed in GTDB species classification ledger (Proposed GTDB species, Ledger species name):'.format(
                len(invalid_sp)))
                
        return gtdb_sp_ledger
                
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
            
    def validate_prev_specific_placeholder_suffixes(self, final_taxonomy, cur_genomes, ncbi_misclassified_gids, new_to_prev_rid, unambiguous_ncbi_sp):
        """Validate that species clusters have same specific name placeholder suffixes as previously assigned."""
        
        unambiguous_rids = set()
        for ncbi_sp, (rid, assignment_type) in unambiguous_ncbi_sp.items():
            unambiguous_rids.add(rid)
            
        ambiguous_rids = set(final_taxonomy) - unambiguous_rids

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
            
            prev_rid = new_to_prev_rid.get(rid, rid)
            if prev_rid in cur_genomes:
                prev_species = cur_genomes[prev_rid].gtdb_taxa.species
                prev_specific = specific_epithet(prev_species)
                prev_suffix = taxon_suffix(prev_specific)
                
                if rid in ncbi_misclassified_gids:
                    # specific name should be a placeholder
                    if is_placeholder_sp_epithet(gtdb_specific):
                        validated_count += 1
                    else:
                        invalid_suffix[rid] = (gtdb_species, prev_species, 'MISCLASSIFIED')
                elif prev_species == 's__':
                    # since genome has no previous GTDB assignment it should have a
                    # suffix higher than observed in the previous taxonomy
                    highest_suffix = self.taxon_suffix_manager.highest_suffix(canonical_taxon(gtdb_species))
                    
                    if highest_suffix is None or self.taxon_suffix_manager.is_higher_suffix(suffix, highest_suffix):
                        validated_count += 1
                    else:
                        invalid_suffix[rid] = (gtdb_species, prev_species, 'NO_PRIOR_SPECIES')
                elif prev_suffix is None:
                    # this should only occur if a NCBI species name is now considered ambiguous,
                    # so all instances must be suffixed
                    if rid not in ambiguous_rids:
                        invalid_suffix[rid] = (gtdb_species, prev_species, 'NO_PRIOR_SUFFIX')
                    else:
                        validated_count += 1
                else:
                    # suffix should be the same between releases
                    if suffix != prev_suffix:
                        invalid_suffix[rid] = (gtdb_species, prev_species, 'SUFFIX DOES NOT MATCH')
                    else:
                        validated_count += 1

        print('Identified {:,} genomes with a suffixed specific name.'.format(suffix_count))
        self.report_validation(
            invalid_suffix, 
            validated_count, 
            'Identified {:,} genomes with an suffixed specific name that disagrees with previous GTDB assignment (Proposed GTDB species, Previous GTDB species, Case):'.format(len(invalid_suffix)))

    def validate_specific_placeholder_suffixes(self,
                                                final_taxonomy, 
                                                cur_genomes, 
                                                ncbi_misclassified_gids,
                                                unambiguous_ncbi_sp):
        """Validate suffix assignment rules for specific names."""
        
        unambiguous_rids = set()
        for ncbi_sp, (rid, assignment_type) in unambiguous_ncbi_sp.items():
            unambiguous_rids.add(rid)
            
        ambiguous_rids = set(final_taxonomy) - unambiguous_rids
            
        # validate that unambiguous NCBI names are not suffixed and that all
        # ambiguous NCBI names are placeholder names
        invalid_suffix = {}
        validated_count = 0
        for rid in final_taxonomy:
            if rid in self.mc_gids:
                continue
                
            gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            gtdb_specific = specific_epithet(gtdb_species)
            ncbi_species = cur_genomes[rid].ncbi_taxa.species
                
            if rid in unambiguous_rids:
                if taxon_suffix(gtdb_specific) is not None:
                    invalid_suffix[rid] = (gtdb_species, ncbi_species, 'UNAMBIGUOUS_SPECIES_WITH_SUFFIX')
                else:
                    validated_count += 1
            else:
                if not is_placeholder_sp_epithet(gtdb_specific):
                    invalid_suffix[rid] = (gtdb_species, ncbi_species, 'AMBIGUOUS_SPECIES_WITH_LATIN_SPECIFIC_NAME')
                else:
                    validated_count += 1

        self.report_validation(
            invalid_suffix, 
            validated_count, 
            'Identified {:,} genomes which violate the suffixing rules for specific names (Proposed GTDB species, NCBI species, Case):'.format(len(invalid_suffix)))
            
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
        
    def validate_specific_no_suffix(self, final_taxonomy, cur_genomes):
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
            if rid in self.mc_gids:
                continue
                
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

    def validate_ncbi_genera_monophyly(self, final_taxonomy, cur_genomes, cur_clusters):
        """Validate monophyly of NCBI genera."""
        
        # get NCBI genus assignments of genomes
        ncbi_genera_gids = defaultdict(set)
        gids_with_ncbi_genus = set()
        for gid in cur_genomes:
            ncbi_genus = cur_genomes[gid].ncbi_taxa.genus
            
            if ncbi_genus != 'g__':
                ncbi_genera_gids[ncbi_genus].add(gid)
                gids_with_ncbi_genus.add(gid)
                
        self.logger.info(' - identified {:,} NCBI genera.'.format(len(ncbi_genera_gids)))
        
        # get GTDB genus assignments of genomes
        gtdb_genera_gids = defaultdict(set)
        for rid, cids in cur_clusters.items():
            if rid not in final_taxonomy:
                continue
                
            gtdb_genus = final_taxonomy[rid][Taxonomy.GENUS_INDEX]
            for cid in cids:
                gtdb_genera_gids[gtdb_genus].add(cid)
                
        # identify NCBI genera that appear to have been merged with a GTDB genus
        gtdb_mergers = defaultdict(list)
        merged_ncbi_genera = 0
        for gtdb_genus, gtdb_gids in gtdb_genera_gids.items():
            canonical_gtdb_genus = canonical_taxon(gtdb_genus)
            
            for ncbi_genus, ncbi_gids in ncbi_genera_gids.items():
                if ncbi_genus == canonical_gtdb_genus:
                    continue
                    
                ncbi_genus_in_gtdb_genus = gtdb_gids.intersection(ncbi_gids)
                ncbi_genus_in_gtdb_genus_perc = len(ncbi_genus_in_gtdb_genus)*100.0/max(1, len(ncbi_gids))
                if ncbi_genus_in_gtdb_genus_perc > 50:
                    gtdb_mergers[gtdb_genus].append(ncbi_genus)
                    merged_ncbi_genera += 1
        self.logger.info(' - identified {:,} NCBI genera that were merged with GTDB genera.'.format(merged_ncbi_genera))
        
        # establish percentage of genomes descendant from each GTDB genus
        # with the same NCBI name
        fout = open(os.path.join(self.output_dir, 'gtdb_genus_monophyly_test.tsv'), 'w')
        fout.write('GTDB genera\tCanonical GTDB genera\tType species of genus\tMerged genera\t% cluster with NCBI genus\t% NCBI genus in cluster\tNCBI genus assignments\n')
        gtdb_genus_minority_support = 0
        for gtdb_genus, gtdb_gids in gtdb_genera_gids.items():
            canonical_gtdb_genus = canonical_taxon(gtdb_genus)
            
            ncbi_genus_gids = ncbi_genera_gids[canonical_gtdb_genus]

            ncbi_genus_in_gtdb_genus = gtdb_gids.intersection(ncbi_genus_gids)
            ncbi_genus_in_gtdb_genus_perc = len(ncbi_genus_in_gtdb_genus)*100.0/max(1, len(ncbi_genus_gids))
            
            gtdb_gids_with_ncbi_genus = gtdb_gids.intersection(gids_with_ncbi_genus)
            gtdb_gids_from_cur_ncbi_genus_perc = len(ncbi_genus_in_gtdb_genus)*100.0/max(1, len(gtdb_gids_with_ncbi_genus))
            
            if gtdb_gids_from_cur_ncbi_genus_perc < 50:
                gtdb_genus_minority_support += 1
            
            ncbi_genus_str = []
            for sp, count in Counter([cur_genomes[gid].ncbi_taxa.genus for gid in gtdb_gids]).most_common():
                ncbi_genus_str.append('{} ({:,})'.format(sp, count))
                
            type_species_of_genus = False
            for gtdb_gid in gtdb_gids:
                if cur_genomes[gtdb_gid].is_gtdb_type_species():
                    ncbi_genus = cur_genomes[gtdb_gid].ncbi_taxa.genus
                    if ncbi_genus == canonical_gtdb_genus:
                        type_species_of_genus = True
                    else:
                        if ncbi_genus not in gtdb_mergers[gtdb_genus]:
                            self.logger.warning('GTDB genus {} contains type species for {}.'.format(
                                                    gtdb_genus,
                                                    cur_genomes[gtdb_gid].ncbi_taxa.genus))

            fout.write('{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\t{}\n'.format(
                        gtdb_genus,
                        canonical_gtdb_genus,
                        type_species_of_genus,
                        ', '.join(sorted(gtdb_mergers[gtdb_genus])),
                        gtdb_gids_from_cur_ncbi_genus_perc,
                        ncbi_genus_in_gtdb_genus_perc,
                        ', '.join(ncbi_genus_str)))
            
        fout.close()
        
        self.logger.info(' - identified {:,} GTDB genera with <50% congruence with underlying NCBI genera.'.format(gtdb_genus_minority_support))

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

    def validate_unambiguous_ncbi_sp(self, final_taxonomy,
                                            unambiguous_ncbi_sp):
        """'Validating NCBI species names classified as having an unambiguous assignment within the GTDB."""
        
        invalid_unambiguous = {}
        validated_count = 0
        for ncbi_sp, (rid, assignment_type) in unambiguous_ncbi_sp.items():
            if rid in self.mc_gids:
                continue
                
            if rid not in final_taxonomy:
                # must be from other domain
                continue
                
            gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            gtdb_specific = specific_epithet(gtdb_species)
            
            ncbi_specific = specific_epithet(ncbi_sp)
            
            if not test_same_epithet(gtdb_specific, ncbi_specific):
                invalid_unambiguous[rid] = (gtdb_species, ncbi_sp)
            else:
                validated_count += 1
            
        self.report_validation(
            invalid_unambiguous, 
            validated_count, 
            'Identified {:,} unambiguous NCBI species with an invalid GTDB species name (Proposed GTDB species, NCBI classification):'.format(
                len(invalid_unambiguous)))
                
    def validate_ambiguous_ncbi_sp(self, final_taxonomy,
                                            cur_genomes,
                                            unambiguous_ncbi_sp):
        """'Validating NCBI species names classified as having an ambiguous assignment within the GTDB."""
        
        unambiguous_rids = set()
        for ncbi_sp, (rid, assignment_type) in unambiguous_ncbi_sp.items():
            unambiguous_rids.add(rid)

        invalid_ambiguous = {}
        validated_count = 0
        for gid in final_taxonomy:
            if gid in self.mc_gids:
                continue
                
            if gid in unambiguous_rids:
                continue
                
            gtdb_species = final_taxonomy[gid][Taxonomy.SPECIES_INDEX]
            gtdb_specific = specific_epithet(gtdb_species)

            if not is_placeholder_sp_epithet(gtdb_specific):
                invalid_ambiguous[gid] = (gtdb_species, ncbi_sp)
            else:
                validated_count += 1
            
        self.report_validation(
            invalid_ambiguous, 
            validated_count, 
            'Identified {:,} ambiguous NCBI species with an invalid GTDB species name (Proposed GTDB species, NCBI classification):'.format(
                len(invalid_ambiguous)))

    def validate_misclassified_genomes(self, final_taxonomy, cur_genomes, ncbi_synonyms, ncbi_misclassified_gids, unambiguous_ncbi_sp):
        """Validating species names of genomes misclassified at NCBI."""
        
        unambiguous_rids = set()
        for ncbi_sp, (rid, assignment_type) in unambiguous_ncbi_sp.items():
            unambiguous_rids.add(rid)

        # validate that any misclassified genomes that are (unfortunately)
        # GTDB species representatives have a placeholder name
        invalid_misclassified = {}
        validated_count = 0
        for gid in ncbi_misclassified_gids:
            if gid in self.mc_gids:
                continue
                
            if gid in final_taxonomy:
                gtdb_species = final_taxonomy[gid][Taxonomy.SPECIES_INDEX]
                gtdb_specific = specific_epithet(gtdb_species)
                    
                if gid in unambiguous_rids:
                    # this occurs when a misclassified genome is the representative of a
                    # unambiguous NCBI species (e.g., G003263425 is representative of 
                    # species cluster with the GTDB name s__Pseudomonas inefficax since this
                    # cluster has a unanimous consensus regarding this name once misclassification
                    # are taken into account). As such, in this case the cluster has a Latin
                    # name and this is considered correct.
                    if is_placeholder_sp_epithet(gtdb_specific): 
                        ncbi_species = cur_genomes[gid].ncbi_taxa.species
                        ncbi_species = ncbi_synonyms.get(ncbi_species, ncbi_species)
                        invalid_misclassified[gid] = (gtdb_species, ncbi_species)
                    else:
                        validated_count += 1
                elif not is_placeholder_sp_epithet(gtdb_specific):
                    ncbi_species = cur_genomes[gid].ncbi_taxa.species
                    ncbi_species = ncbi_synonyms.get(ncbi_species, ncbi_species)
                    invalid_misclassified[gid] = (gtdb_species, ncbi_species)
                else:
                    validated_count += 1

        self.report_validation(
            invalid_misclassified, 
            validated_count, 
            'Identified {:,} misclassified genomes with an invalid specific name (Proposed GTDB species, NCBI classification):'.format(
                len(invalid_misclassified)))
            
    def report_modified_species_names(self, 
                                        final_taxonomy, 
                                        cur_genomes, 
                                        prev_genomes, 
                                        cur_clusters,
                                        prev_to_new_rids,
                                        unambiguous_ncbi_sp):
        """Report GTDB species clusters with modified names."""
        
        unambiguous_rids = set()
        for ncbi_sp, (rid, assignment_type) in unambiguous_ncbi_sp.items():
            unambiguous_rids.add(rid)
        
        domains_of_interest = set([taxa[Taxonomy.DOMAIN_INDEX] for taxa in final_taxonomy.values()])
        
        fout = open(os.path.join(self.output_dir, 'modified_sp_names.tsv'), 'w')
        fout.write('Previous ID\tNew ID\tPrevious name\tNew name')
        fout.write('\tModified generic name\tModified specific name\tUpdated representative')
        fout.write('\tGenus name transition\tSpecific name transition\tSpecific name change category')
        fout.write('\tRepresentative is NCBI representative\tCluster contain NCBI representative')
        fout.write('\tUnambiguous NCBI species\tNCBI species assignments\n')
        
        test_reassigned_ncbi_rep = 0
        test_reassigned_ncbi_rep_in_cluster = 0
        
        modified_species_names = 0
        modified_generic_names = 0
        modified_specific_names = 0
        
        modified_generic_transition = defaultdict(lambda: defaultdict(int))
        modified_specific_transition = defaultdict(lambda: defaultdict(int))
        modified_specific_type = defaultdict(int)
        
        same_species_names = 0
        lost_sp_cluster = 0
        total_names = 0
        for prev_rid in prev_genomes.sp_clusters:
            prev_domain = prev_genomes[prev_rid].gtdb_taxa.domain
            if prev_domain not in domains_of_interest:
                continue
                
            total_names += 1

            prev_sp = prev_genomes[prev_rid].gtdb_taxa.species
            prev_generic = generic_name(prev_sp)
            prev_specific = specific_epithet(prev_sp)
            
            new_rid = prev_to_new_rids[prev_rid]
            if new_rid is None:
                lost_sp_cluster += 1
                fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            prev_rid,
                            'n/a',
                            prev_sp,
                            'n/a', 'n/a', 'n/a', 'n/a', 'n/a',
                            'n/a', 'n/a', 'n/a', 'n/a'))
            elif new_rid in final_taxonomy:
                new_sp = final_taxonomy[new_rid][Taxonomy.SPECIES_INDEX]
                new_generic = generic_name(new_sp)
                new_specific = specific_epithet(new_sp)
                
                if prev_sp == new_sp:
                    same_species_names += 1
                else:
                    modified_species_names += 1
                    
                    # check if GTDB representative is in cluster
                    ncbi_rep_in_cluster = False
                    for cid in cur_clusters[new_rid]:
                        ncbi_specific = specific_epithet(cur_genomes[new_rid].ncbi_taxa.species)
                        if cur_genomes[cid].is_ncbi_representative():
                            ncbi_rep_in_cluster = True
                
                    genus_transition = 'n/a'
                    if prev_generic != new_generic:
                        modified_generic_names += 1
                        modified_generic_transition[taxon_type('g__' + prev_generic)][taxon_type('g__' + new_generic)] += 1
                        genus_transition = '{} -> {}'.format(
                                                taxon_type('g__' + prev_generic), 
                                                taxon_type('g__' + new_generic))
                        
                    specific_transition = 'n/a'
                    specific_change_type = 'n/a'
                    if prev_specific != new_specific:
                        modified_specific_names += 1
                        modified_specific_transition[specific_epithet_type(prev_specific)][specific_epithet_type(new_specific)] += 1
                        specific_transition = '{} -> {}'.format(
                                                specific_epithet_type(prev_specific), 
                                                specific_epithet_type(new_specific))
                        
                        #***
                        if specific_epithet_type(prev_specific) == 'LATIN' and specific_epithet_type(new_specific) != 'LATIN':
                            if cur_genomes[new_rid].is_ncbi_representative():
                                test_reassigned_ncbi_rep += 1
                            if ncbi_rep_in_cluster:
                                test_reassigned_ncbi_rep_in_cluster += 1

                        if test_same_epithet(prev_specific, new_specific):
                            specific_change_type = 'GENUS_TRANSFER_SUFFIX_CHANGE'
                        elif cur_genomes[new_rid].is_effective_type_strain():
                            specific_change_type = 'ESTABLISH_BY_TYPE_STRAIN_GENOME'
                        else:
                            specific_change_type = 'OTHER'
                            
                        modified_specific_type[specific_change_type] += 1
                            
                    # get NCBI species assignments for genomes in cluster
                    ncbi_sp_str = []
                    for sp, count in Counter([cur_genomes[cid].ncbi_taxa.species for cid in cur_clusters[new_rid]]).most_common():
                        ncbi_sp_str.append('{} ({:,})'.format(sp, count))


                    fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                prev_rid,
                                new_rid,
                                prev_sp,
                                new_sp,
                                prev_generic != new_generic,
                                prev_specific != new_specific,
                                prev_rid != new_rid,
                                genus_transition,
                                specific_transition,
                                specific_change_type,
                                cur_genomes[new_rid].is_ncbi_representative(),
                                ncbi_rep_in_cluster,
                                new_rid in unambiguous_rids,
                                ', '.join(ncbi_sp_str)))
            else:
                # genome must be from other domain
                pass
            
        print('***test_reassigned_ncbi_rep', test_reassigned_ncbi_rep)
        print('***test_reassigned_ncbi_rep_in_cluster', test_reassigned_ncbi_rep_in_cluster)
        
        self.logger.info(' - total previous species names: {:,}'.format(total_names))
        self.logger.info(' - unmodified names: {:,} ({:.2f}%)'.format(
                        same_species_names,
                        same_species_names*100.0/total_names))
        self.logger.info(' - modified species names: {:,} ({:.2f}%)'.format(
                        modified_species_names,
                        modified_species_names*100.0/total_names))
                        
        self.logger.info(' - modified generic names: {:,} ({:.2f}%)'.format(
                        modified_generic_names,
                        modified_generic_names*100.0/total_names))
        for prev_type, new_type in product(['LATIN', 'SUFFIXED_LATIN', 'ALPHANUMERIC'], repeat=2):
            self.logger.info('   - modified generic names from {} to {}: {:,} ({:.2f}%)'.format(
                        prev_type,
                        new_type,
                        modified_generic_transition[prev_type][new_type],
                        modified_generic_transition[prev_type][new_type]*100.0/modified_generic_names))
                        
        self.logger.info(' - modified specific names: {:,} ({:.2f}%)'.format(
                        modified_specific_names,
                        modified_specific_names*100.0/total_names))
        for prev_type, new_type in product(['LATIN', 'SUFFIXED_LATIN', 'ALPHANUMERIC'], repeat=2):
            self.logger.info('   - modified specific names from {} to {}: {:,} ({:.2f}%)'.format(
                        prev_type,
                        new_type,
                        modified_specific_transition[prev_type][new_type],
                        modified_specific_transition[prev_type][new_type]*100.0/modified_specific_names))
        for specific_change_type, count in modified_specific_type.items():
            self.logger.info(' - {}: {:,} ({:.2f}%)'.format(
                                specific_change_type,
                                count,
                                count*100.0/modified_specific_names))

        self.logger.info(' - lost species clusters: {:,} ({:.2f}%)'.format(
                        lost_sp_cluster,
                        lost_sp_cluster*100.0/total_names))

        fout.close()
        
        # do a comparison between GTDB and NCBI in terms of number of modified names
        gtdbd_modified_sp_name = 0
        ncbi_modified_sp_name = 0
        gtdb_modified_genus_names = 0
        ncbi_modified_genus_names = 0
        gtdb_modified_specific_names = 0
        ncbi_modified_specific_names = 0
        total_names = 0
        for prev_rid in prev_genomes.sp_clusters:
            prev_domain = prev_genomes[prev_rid].gtdb_taxa.domain
            if prev_domain not in domains_of_interest:
                continue
                
            new_rid = prev_to_new_rids.get(prev_rid, prev_rid)
            if new_rid not in final_taxonomy:
                continue 

            prev_ncbi_sp = prev_genomes[prev_rid].ncbi_taxa.species
            if prev_ncbi_sp == 's__':
                continue
            prev_ncbi_genus = prev_genomes[prev_rid].ncbi_taxa.genus
            prev_ncbi_specific = specific_epithet(prev_ncbi_sp)

            new_ncbi_sp = cur_genomes[new_rid].ncbi_taxa.species
            new_ncbi_genus = cur_genomes[new_rid].ncbi_taxa.genus
            new_ncbi_specific = specific_epithet(new_ncbi_sp)
            
            prev_gtdb_sp = prev_genomes[prev_rid].gtdb_taxa.species
            prev_gtdb_genus = prev_genomes[prev_rid].gtdb_taxa.genus
            prev_gtdb_specific = specific_epithet(prev_gtdb_sp)
            
            new_gtdb_sp = final_taxonomy[new_rid][Taxonomy.SPECIES_INDEX]
            new_gtdb_genus = final_taxonomy[new_rid][Taxonomy.GENUS_INDEX]
            new_gtdb_specific = specific_epithet(new_gtdb_sp)
            
            total_names += 1
            
            if prev_ncbi_sp != new_ncbi_sp:
                ncbi_modified_sp_name += 1
                
            if prev_gtdb_sp != new_gtdb_sp:
                gtdbd_modified_sp_name += 1

            if prev_ncbi_genus != new_ncbi_genus:
                ncbi_modified_genus_names += 1
            
            if prev_gtdb_genus != new_gtdb_genus:
                gtdb_modified_genus_names += 1
                
            if prev_ncbi_specific != new_ncbi_specific:
                ncbi_modified_specific_names += 1
                
            if prev_gtdb_specific != new_gtdb_specific:
                gtdb_modified_specific_names += 1
                
        self.logger.info('Comparison of number of name changes in NCBI and GTDB.')
        self.logger.info(' - total previous GTDB representatives with NCBI species assignment: {:,}'.format(total_names))
        self.logger.info(' - modified NCBI species names: {:,} ({:.2f}%)'.format(
                        ncbi_modified_sp_name,
                        ncbi_modified_sp_name*100.0/total_names))
        self.logger.info(' - modified GTDB species names: {:,} ({:.2f}%)'.format(
                        gtdbd_modified_sp_name,
                        gtdbd_modified_sp_name*100.0/total_names))
        self.logger.info(' - modified NCBI genus names: {:,} ({:.2f}%)'.format(
                        ncbi_modified_genus_names,
                        ncbi_modified_genus_names*100.0/total_names))
        self.logger.info(' - modified GTDB genus names: {:,} ({:.2f}%)'.format(
                        gtdb_modified_genus_names,
                        gtdb_modified_genus_names*100.0/total_names))
        self.logger.info(' - modified NCBI specific names: {:,} ({:.2f}%)'.format(
                        ncbi_modified_specific_names,
                        ncbi_modified_specific_names*100.0/total_names))
        self.logger.info(' - modified GTDB specific names: {:,} ({:.2f}%)'.format(
                        gtdb_modified_specific_names,
                        gtdb_modified_specific_names*100.0/total_names))

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
                ncbi_misclassified_file,
                ncbi_genbank_assembly_file,
                untrustworthy_type_file,
                synonym_file,
                updated_species_reps,
                gtdb_type_strains_ledger,
                species_classification_ledger,
                sp_priority_ledger,
                genus_priority_ledger,
                specific_epithet_ledger,
                dsmz_bacnames_file,
                ground_truth_test_cases,
                skip_genus_checks):
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
        
        # get mappings between new and previous GTDB representatives
        self.logger.info('Mapping current GTDB representatives to previous representatives.')
        updated_gtdb_rids = parse_updated_species_reps(updated_species_reps)
        new_to_prev_rid = infer_prev_gtdb_reps(prev_genomes, cur_clusters, updated_gtdb_rids)
        self.logger.info(' - mapped {:,} current representatives to previous representatives.'.format(len(new_to_prev_rid)))
        
        # get species clusters with reassigned representative
        self.logger.info('Identifying species clusters with reassigned representatives.')
        prev_to_new_rids = prev_genomes.sp_clusters.updated_representatives(cur_clusters)
        self.logger.info(' - mapped {:,} previous representatives to current representatives.'.format(
                            len(prev_to_new_rids)))
                            
        # get polyphyletic suffixes for GTDB taxa
        self.taxon_suffix_manager = TaxonSuffixManager()
        
        # read scaled curation tree
        self.logger.info('Reading scaled and decorated tree.')
        curation_tree = dendropy.Tree.get_from_path(final_scaled_tree, 
                                            schema='newick', 
                                            rooting='force-rooted', 
                                            preserve_underscores=True)
        curation_tree.calc_node_root_distances()

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
        ncbi_misclassified_gids = ncbi_species_mngr.parse_ncbi_misclassified_table(ncbi_misclassified_file)
        self.logger.info(' - identified {:,} genomes with erroneous NCBI species assignments'.format(
                            len(ncbi_misclassified_gids)))

        # classify each Latin NCBI species as either unambiguously or ambiguously assigned
        self.logger.info('Classifying NCBI species as unambiguous, ambiguous, or synonym.')
        ncbi_classification_table = os.path.join(self.output_dir, 'ncbi_sp_classification.tsv')
        if os.path.exists(ncbi_classification_table):
            self.logger.warning('Reading classification of NCBI species from existing table: {}'.format(
                                    ncbi_classification_table))
            ncbi_species, unambiguous_ncbi_sp, ambiguous_ncbi_sp = ncbi_species_mngr.parse_ncbi_classification_table(
                                                                    ncbi_classification_table)
        else:
            ncbi_species, unambiguous_ncbi_sp, ambiguous_ncbi_sp = ncbi_species_mngr.classify_ncbi_species(
                                                                        ncbi_synonyms, 
                                                                        ncbi_misclassified_gids)
        
        unambiguous_classifications = [classification for rid, classification in unambiguous_ncbi_sp.values()]
        unambiguous_type_strain_count = unambiguous_classifications.count('TYPE_STRAIN_GENOME')
        unambiguous_unanimous_consensus_count = unambiguous_classifications.count('UNANIMOUS_CONSENSUS')
        assert unambiguous_type_strain_count + unambiguous_unanimous_consensus_count == len(unambiguous_ncbi_sp)
                
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
                            
        # validate manually-curated species names
        self.logger.info('Validating manually-curated species names.')
        self.mc_gids = self.validate_manually_curated_species(final_taxonomy, 
                                                                manual_sp_names, 
                                                                pmc_custom_species)
        
        # validate genomes in GTDB type strain ledger
        self.logger.info('Validating genomes in type strain ledger.')
        mc_type_strain = self.validate_type_strain_ledger(final_taxonomy, gtdb_type_strains_ledger)
        self.mc_gids.update(mc_type_strain)
        
        # validate genomes in GTDB type strain ledger
        self.logger.info('Validating genomes in species classification ledger.')
        mc_sp_classification = self.validate_species_classification_ledger(final_taxonomy, 
                                                                            cur_clusters, 
                                                                            species_classification_ledger)
        self.mc_gids.update(mc_sp_classification)

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

        # validate manually established suffixes for specific names changed due to genus transfer
        self.logger.info('Validating manually-curated specific name suffix changes resulting from genus transfers.')
        self.validate_mannually_curated_specifix_suffix(final_taxonomy, cur_genomes, specific_epithet_ledger)
        
        # validate NCBI species considered to have an unambiguous assignment within the GTDB
        self.logger.info('Validating NCBI species names classified as having an unambiguous assignment within the GTDB.')
        self.validate_unambiguous_ncbi_sp(final_taxonomy,
                                            unambiguous_ncbi_sp)
                                            
        # validate NCBI species not considered to have an unambiguous assignment within the GTDB
        self.logger.info('Validating NCBI species names classified as having an ambiguous assignment within the GTDB.')
        self.validate_ambiguous_ncbi_sp(final_taxonomy,
                                            cur_genomes,
                                            unambiguous_ncbi_sp)

        # validate species names of genomes misclassified at NCBI
        self.logger.info('Validating species names of genomes misclassified at NCBI.')
        self.validate_misclassified_genomes(final_taxonomy, 
                                            cur_genomes, 
                                            ncbi_synonyms,
                                            ncbi_misclassified_gids,
                                            unambiguous_ncbi_sp)

        # validate synonyms
        self.logger.info('Validating that GTDB synonyms are respected.')
        self.validate_synonyms(final_taxonomy, cur_genomes, ncbi_synonyms)
        
        # validate that each GTDB genus has a single type species
        self.logger.info('Validating that each GTDB genus has a single type species cluster.')
        self.validate_single_type_species(final_taxonomy, cur_genomes)

        # validate GTDB specific name assignments without a polyphyletic suffixes
        self.logger.info('Validating specific name assignments without a placeholder suffix.')
        self.validate_specific_no_suffix(final_taxonomy, cur_genomes)

        # validate that species clusters use same specific name placeholder suffixes as previously assigned
        self.logger.info('Validating specific names use same placeholder suffixes as assigned in previous release.')
        self.validate_prev_specific_placeholder_suffixes(final_taxonomy, 
                                                    cur_genomes, 
                                                    ncbi_misclassified_gids,
                                                    new_to_prev_rid,
                                                    unambiguous_ncbi_sp)

        # validate suffix assignment rules for specific names
        self.logger.info('Validating specific names follow suffix assignment rules.')
        self.validate_specific_placeholder_suffixes(final_taxonomy, 
                                                    cur_genomes, 
                                                    ncbi_misclassified_gids,
                                                    unambiguous_ncbi_sp)

        # validate spelling of genera
        self.logger.info('Validating spelling of genera.')
        self.validate_genera_spelling(final_taxonomy, cur_genomes, new_to_prev_rid)
        
        # validate that no forbidden specific names are used
        self.logger.info('Validating that no forbidden specific names were used.')
        self.validate_forbidden_specific_names(final_taxonomy)

        # validate that Latinized specific names are never used with non-suffixed, placeholder generic names
        self.logger.info('Validating that Latinized specific names are not used with non-suffixed, placeholder generic names.')
        self.validate_latin_sp_placeholder_generic(final_taxonomy, cur_genomes)

        if not skip_genus_checks:
            # validate placeholder suffixes for genera
            self.logger.info('Validating placeholder suffixes of genera.')
            self.validate_genera_placeholder_suffixes(final_taxonomy, cur_genomes, new_to_prev_rid)

            # validate missing NCBI genus names
            self.logger.info('Validating missing NCBI-defined genus names.')
            self.validate_missing_ncbi_genera(final_taxonomy, cur_genomes, curation_tree)

            # validate placeholder genus names lack potential Latin names
            self.logger.info('Validating placeholder genus names lack potential Latin names.')
            self.validate_no_latin_genus_name(final_taxonomy, cur_genomes)
            
            # validate that genus names are not a replacement of one placeholder for another
            self.logger.info('Validating that genus names are not a change from one placeholder to another.')
            self.validate_genus_placeholder_change(final_taxonomy, cur_genomes, prev_genomes)
            
            # validate monophyly of NCBI genera
            self.logger.info('Validating monophyly of NCBI genera.')
            self.validate_ncbi_genera_monophyly(final_taxonomy, 
                                                cur_genomes, 
                                                cur_clusters)

        # create table with modified GTDB species names for manual inspection
        self.logger.info('Creating table with modified GTDB species names.')
        self.report_modified_species_names(final_taxonomy, 
                                            cur_genomes, 
                                            prev_genomes, 
                                            cur_clusters,
                                            prev_to_new_rids,
                                            unambiguous_ncbi_sp)
