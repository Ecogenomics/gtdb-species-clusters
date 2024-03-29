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
import logging
from itertools import product
from collections import defaultdict, Counter

from numpy import median as np_median

import dendropy

from Levenshtein import (distance as lvn_distance)

from gtdblib.taxonomy.taxonomy import Taxonomy
from gtdblib.taxonomy.taxonomy import read_taxonomy
from gtdblib.taxonomy.validation import validate_taxonomy
from gtdblib.util.bio.newick import parse_label
from gtdblib.util.bio.accession import canonical_gid

from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.taxon_suffix_manager import TaxonSuffixManager
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters.ncbi_species_manager import NCBI_SpeciesManager
from gtdb_species_clusters.type_genome_utils import (read_clusters,
                                                     parse_updated_species_reps,
                                                     infer_prev_gtdb_reps,
                                                     parse_manual_sp_curation_files)
from gtdb_species_clusters.taxon_utils import (generic_name,
                                               specific_epithet,
                                               canonical_taxon,
                                               canonical_species,
                                               taxon_suffix,
                                               is_placeholder_taxon,
                                               is_placeholder_sp_epithet,
                                               is_latin_sp_epithet,
                                               is_suffixed_taxon,
                                               is_alphanumeric_taxon,
                                               is_suffixed_sp_epithet,
                                               taxon_type,
                                               specific_epithet_type,
                                               test_same_epithet,
                                               ncbi_to_gtdb_synonyms)

from gtdb_species_clusters.lpsn import LPSN


class PMC_Validation(object):
    """Validate final species names."""

    def __init__(self, output_dir):
        """Initialization."""

        self.output_dir = output_dir
        self.log = logging.getLogger('rich')

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
                print('{}'.format('\t'.join([str(f)
                                             for f in invalid_cases[gid]])))
            print('')

    def validate_genus_generic_match(self, final_taxonomy):
        """Validate that genus and generic names match."""

        invalid_names = {}
        validated_count = 0
        for gid in final_taxonomy:
            gtdb_genus = final_taxonomy[gid][Taxonomy.GENUS_INDEX]
            if gtdb_genus == 'g__':
                # GTDB genus assignment has not yet been performed,
                # which happens when evaluating preliminary taxonomies
                continue

            gtdb_species = final_taxonomy[gid][Taxonomy.SPECIES_INDEX]
            gtdb_generic = generic_name(gtdb_species)

            if gtdb_genus != 'g__' + gtdb_generic:
                invalid_names[gid] = (gid, gtdb_genus, gtdb_species)
            else:
                validated_count += 1

        self.report_validation(
            invalid_names,
            validated_count,
            ' - identified {:,} genomes with unmatched genus and generic names (Proposed GTDB genus, GTDB species):'.format(len(invalid_names)))

    def validate_type_species(self, final_taxonomy, cur_genomes, sp_priority_mngr):
        """Validate generic name of type species of genus."""

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
                if gtdb_genus == 'g__':
                    # GTDB genus assignment has not yet been performed,
                    # which happens when evaluating preliminary taxonomies
                    continue

                gtdb_species = final_taxonomy[gid][Taxonomy.SPECIES_INDEX]
                ncbi_genus = cur_genomes[gid].ncbi_taxa.genus

                if (gtdb_genus != ncbi_genus
                        and sp_priority_mngr.genus_priority(gtdb_genus, ncbi_genus) != gtdb_genus):
                    invalid_type_sp[gid] = (gid, gtdb_species, ncbi_genus)
                else:
                    validated_count += 1

        self.report_validation(
            invalid_type_sp,
            validated_count,
            ' - identified {:,} invalid type species genomes (Proposed GTDB species, NCBI genus):'.format(len(invalid_type_sp)))

    def validate_type_strain(self, final_taxonomy, cur_genomes):
        """Validate specific name of type strains."""

        # get all GTDB species represented by a type strain:
        gtdb_type_species = set()
        for rid in final_taxonomy:
            if cur_genomes[rid].is_effective_type_strain():
                gtdb_type_species.add(
                    final_taxonomy[rid][Taxonomy.SPECIES_INDEX])

        # check that species clusters represented by type strain genomes
        # have the specific name of the type strain
        invalid_type_strain = {}
        validated_count = 0
        for rid in final_taxonomy:
            if rid in self.mc_species:
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
                    invalid_type_strain[rid] = (
                        rid, gtdb_species, ncbi_species)
                else:
                    validated_count += 1

        self.report_validation(
            invalid_type_strain,
            validated_count,
            ' - identified {:,} invalid type strain genomes (Proposed GTDB species, NCBI species):'.format(len(invalid_type_strain)))

    def validate_lpsn_correct_names(self, final_taxonomy, lpsn):
        """Validate that GTDB is using the 'correct' names as specified at LPSN."""

        fout = open(os.path.join(self.output_dir,
                                 'validate_lpsn_correct_sp.tsv'), 'w')
        fout.write(
            'Genome ID\tGTDB species\tLPSN correct species\tLPSN synonyms\tGeneric change\tSpecific change\n')

        incorrect_names = {}
        validated_count = 0
        incorrect_synonym_count = 0
        for rid in final_taxonomy:
            gtdb_sp = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            if is_placeholder_taxon(gtdb_sp):
                continue

            gtdb_generic = generic_name(gtdb_sp)
            gtdb_specific = specific_epithet(gtdb_sp)

            if gtdb_sp in lpsn.sp_correct_names:
                # check if the GTDB species name is considered the correct name at LPSN
                lpsn_corr_taxon = lpsn.sp_correct_names[gtdb_sp]
                if 'subsp.' in lpsn_corr_taxon:
                    # remove subspecies portion
                    lpsn_corr_sp = ' '.join(lpsn_corr_taxon.split()[0:2])
                else:
                    lpsn_corr_sp = lpsn_corr_taxon

                if gtdb_sp != lpsn_corr_taxon:
                    synonym = lpsn_corr_taxon in lpsn.sp_synonyms[gtdb_sp]
                    row = (rid,
                           gtdb_sp,
                           lpsn_corr_taxon,
                           synonym,
                           gtdb_generic != generic_name(lpsn_corr_sp),
                           gtdb_specific != specific_epithet(lpsn_corr_sp))

                    fout.write('{}\n'.format('\t'.join([str(f) for f in row])))
                    if not synonym:
                        incorrect_names[rid] = row
                        incorrect_synonym_count += 1

                else:
                    validated_count += 1

        self.log.info(
            f' - identified {incorrect_synonym_count:,} cases where GTDB uses a name considered a synonym at LPSN (see validate_lpsn_correct_sp.tsv)')
        self.report_validation(
            incorrect_names,
            validated_count,
            ' - identified {:,} GTDB species that disagrees with correct species names at LPSN (GTDB species, LPSN correct species, LPSN synonyms, Generic change, Specific change):'.format(len(incorrect_names)))

        fout.close()

    def validate_genus_transfer_specific_epithet(self,
                                                 final_taxonomy,
                                                 cur_genomes,
                                                 lpsn):
        """Validate transferred species do not use a specific epithet already validated for the new species."""

        fout = open(os.path.join(self.output_dir,
                                 'validate_epithet_genus_transfer.tsv'), 'w')
        fout.write('Genome ID\tGTDB species\tNCBI species\tLPSN species\n')

        manual_check = {}
        validated_count = 0
        for rid in final_taxonomy:
            if rid in self.mc_species:
                continue

            gtdb_sp = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            if is_placeholder_taxon(gtdb_sp):
                # this is always acceptable as the GTDB species is either
                # from a placeholder genus or has a placeholder specific epithet
                continue

            gtdb_genus = final_taxonomy[rid][Taxonomy.GENUS_INDEX]
            gtdb_generic = generic_name(gtdb_sp)
            gtdb_specific = specific_epithet(gtdb_sp)

            ncbi_genus = cur_genomes[rid].ncbi_taxa.genus
            if ncbi_genus == 'g__':
                continue
            ncbi_sp = cur_genomes[rid].ncbi_taxa.species.replace(
                '[', '').replace(']', '')

            # check if genome has been transfered to a new genus
            potential_issue = False
            if gtdb_genus != ncbi_genus:
                if gtdb_sp == ncbi_sp:
                    # this *should* never happen, but NCBI doesn't always
                    # have consistent genus and generic names for genomes!
                    pass
                elif gtdb_sp == lpsn.sp_correct_names.get(ncbi_sp, None):
                    # GTDB assignment reflects the 'correct name' as
                    # indicated at LPSN so this is a valid (and expected)
                    # genus transfer
                    pass
                elif ncbi_sp == lpsn.sp_correct_names.get(gtdb_sp, None):
                    # It appears the GTDB is using an older name that isn't
                    # considered correct by LPSN, but that this is otherwise
                    # an acceptable genus transfer.
                    pass
                elif gtdb_sp in lpsn.sp_synonyms.get(ncbi_sp, []):
                    # GTDB assignment is using a recognized synonym at LPSN,
                    # though NCBI appears to be using the correct name according
                    # to LPSN
                    pass
                elif ncbi_sp in lpsn.sp_synonyms.get(gtdb_sp, []):
                    # GTDB assignment is the correct name at LPSN and
                    # NCBI is using a later synonym
                    pass
                else:
                    # check if GTDB species name has already been proposed
                    # and thus retaining the specific name after the genus
                    # transfer may be incorrect
                    # (e.g. G003007735 is Methylarcula marina at NCBI, but transferred into
                    # Paracoccus under the GTDB. However, P. marinus is already a validated
                    # name so setting the GTDB species to P. marinus would be incorrect.)

                    matched_lpsn_sp = None
                    if gtdb_sp in lpsn.sp_correct_names:
                        potential_issue = True
                        matched_lpsn_sp = gtdb_sp
                    else:
                        # need to check if a LPSN species name exists that only differs in terms
                        # of changes to the suffix of the specific epithet due to the gender of
                        # the genus
                        for lpsn_sp in lpsn.sp_correct_names:
                            lpsn_generic = generic_name(lpsn_sp)
                            lpsn_specific = specific_epithet(lpsn_sp)
                            if lpsn_generic == gtdb_generic and test_same_epithet(gtdb_specific, lpsn_specific):
                                potential_issue = True
                                matched_lpsn_sp = lpsn_sp
                                break

            if potential_issue:
                manual_check[rid] = (rid, gtdb_sp, ncbi_sp, matched_lpsn_sp)
                fout.write(f'{rid}\t{gtdb_sp}\t{ncbi_sp}\t{matched_lpsn_sp}\n')
            else:
                validated_count += 1

        self.report_validation(
            manual_check,
            validated_count,
            ' - identified {:,} genomes with GTDB assignments that may be invalid after genus transfer (Proposed GTDB species, NCBI species, LPSN species):'.format(len(manual_check)))

        fout.close()

    def validate_unique_species_names(self, final_taxonomy, cur_genomes):
        """Validate that species names are unique."""

        invalid_dup_name = {}
        validated_count = 0

        identified_gtdb_species = {}
        for rid in final_taxonomy:
            gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            if gtdb_species in identified_gtdb_species:
                prev_rid = identified_gtdb_species[gtdb_species]
                invalid_dup_name[rid] = (
                    rid, gtdb_species, cur_genomes[rid].ncbi_taxa.species)
                invalid_dup_name[prev_rid] = (
                    prev_rid, gtdb_species, cur_genomes[prev_rid].ncbi_taxa.species)
            else:
                identified_gtdb_species[gtdb_species] = rid
                validated_count += 1

        self.report_validation(
            invalid_dup_name,
            validated_count,
            ' - identified {:,} genomes with identical species names (Proposed GTDB species, NCBI species):'.format(len(invalid_dup_name)))

    def validate_manually_curated_species(self, final_taxonomy, mc_species):
        """Validate final assignment of species that were set by manual curation."""

        # read species names explicitly set via manual curation
        invalid_mc = {}
        validated_count = 0
        for gid, mc_sp in mc_species.items():
            if gid not in final_taxonomy:
                continue  # taxon from other domain

            if final_taxonomy[gid][Taxonomy.SPECIES_INDEX] != mc_sp:
                invalid_mc[gid] = (gid,
                                   final_taxonomy[gid]
                                   [Taxonomy.SPECIES_INDEX], mc_sp)
            else:
                validated_count += 1

        self.report_validation(
            invalid_mc,
            validated_count,
            ' - identified {:,} genomes that did not match manual curation (Proposed GTDB species, manually curated species):'.format(len(invalid_mc)))

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
                    invalid_type_strain[gid] = (
                        gid,
                        gtdb_species,
                        gtdb_type_strains[gid])
                else:
                    validated_count += 1

        self.report_validation(
            invalid_type_strain,
            validated_count,
            ' - identified {:,} genomes that do not match name proposed in GTDB type strain ledger (Proposed GTDB species, Ledger species name):'.format(
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

                if gid in gtdb_gid_to_rid:
                    # if gid is not in gtdb_gid_to_rid it indicates
                    # the genome has been removed from GTDB which does occur
                    rid = gtdb_gid_to_rid[gid]
                    gtdb_sp_ledger[rid] = sp

        # validate names of genomes in species classification ledger
        invalid_sp = {}
        validated_count = 0
        for gid in gtdb_sp_ledger:
            if gid in self.mc_species:
                continue

            if gid in final_taxonomy:
                gtdb_species = final_taxonomy[gid][Taxonomy.SPECIES_INDEX]
                gtdb_specific = specific_epithet(gtdb_species)

                ledger_specific = specific_epithet(gtdb_sp_ledger[gid])

                if not test_same_epithet(gtdb_specific, ledger_specific):
                    invalid_sp[gid] = (gid, gtdb_species, gtdb_sp_ledger[gid])
                else:
                    validated_count += 1

        self.report_validation(
            invalid_sp,
            validated_count,
            ' - identified {:,} genomes that do not match name proposed in GTDB species classification ledger (Proposed GTDB species, Ledger species name):'.format(
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
                    if (test_same_epithet(gtdb_canonical_specific, curated_epithet_changes[gtdb_generic][ncbi_specific])
                            and gtdb_canonical_specific != curated_epithet_changes[gtdb_generic][ncbi_specific]):
                        # epithets are highly similar, but not identical indicating a likely curation issue
                        invalid_suffix[rid] = (rid,
                                               gtdb_species, curated_epithet_changes[gtdb_generic][ncbi_specific])
                    else:
                        validated_count += 1

        self.report_validation(
            invalid_suffix,
            validated_count,
            ' - identified {:,} genomes that did not match manually curated changes to the suffix of specific name (Proposed GTDB species, manually curated specific name):'.format(len(invalid_suffix)))

    def validate_ground_truth(self, final_taxonomy, ground_truth_test_cases):
        """Validate final assignments against a set of ground truth assignments."""

        if ground_truth_test_cases.lower() == 'none':
            self.log.info(
                ' - no ground truth cases provided for validation.')
            return

        invalid_gt = {}
        validated_count = 0
        with open(ground_truth_test_cases) as f:
            for line in f:
                tokens = line.strip().split('\t')
                gid = tokens[0]
                gt_taxa_str = tokens[1]
                gt_taxa = [t.strip() for t in gt_taxa_str.split(';')]

                if final_taxonomy[gid] != gt_taxa:
                    invalid_gt[gid] = (
                        gid, final_taxonomy[gid][Taxonomy.SPECIES_INDEX], gt_taxa[Taxonomy.SPECIES_INDEX])
                else:
                    validated_count += 1

        self.report_validation(
            invalid_gt,
            validated_count,
            ' - identified {:,} genomes that did not match ground truth assignments (Proposed GTDB species, ground truth species):'.format(len(invalid_gt)))

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
                    invalid_synonym[rid] = (rid,
                                            gtdb_species,
                                            ncbi_species,
                                            '{} is a synonym of {}'.format(ncbi_species, ncbi_synonyms[ncbi_species]))
                else:
                    validated_count += 1

        self.report_validation(
            invalid_synonym,
            validated_count,
            ' - identified {:,} genomes that do not respect GTDB synonyms (Proposed GTDB species, NCBI species, NCBI synonym):'.format(len(invalid_synonym)))

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
                    genus_gtdb_type_species[gtdb_genus].add(
                        (gtdb_species, ncbi_species))

        # report genus with multiple type species
        invalid_num_type_species = {}
        validated_count = 0
        for genus, type_species in genus_gtdb_type_species.items():
            if len(type_species) > 1:
                invalid_num_type_species[genus] = (genus, [
                    ', '.join([str(gtdb_ncbi_sp) for gtdb_ncbi_sp in type_species])])
            else:
                validated_count += 1

        self.report_validation(
            invalid_num_type_species,
            validated_count,
            ' - identified {:,} genera with multiple type species (GTDB genera, (GTDB species, NCBI species)):'.format(len(invalid_num_type_species)))

    def validate_prev_specific_placeholder_suffixes(self, final_taxonomy, cur_genomes, new_to_prev_rid, unambiguous_ncbi_sp):
        """Validate that species clusters have same specific name placeholder suffixes as previously assigned."""

        unambiguous_rids = set()
        for rid, _assignment_type in unambiguous_ncbi_sp.values():
            unambiguous_rids.add(rid)

        ambiguous_rids = set(final_taxonomy) - unambiguous_rids

        invalid_suffix = {}
        validated_count = 0
        suffix_count = 0
        for rid in final_taxonomy:
            if rid in self.mc_species:
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

                if not test_same_epithet(canonical_taxon(gtdb_specific), canonical_taxon(prev_specific)):
                    # no need to match the placeholder suffix as the specific name
                    # has changed between releases and it isn't just a small change
                    # indicative of a genus transfer
                    validated_count += 1
                elif prev_species == 's__':
                    # since genome has no previous GTDB assignment it should have a
                    # suffix higher than observed in the previous taxonomy
                    highest_suffix = self.taxon_suffix_manager.highest_suffix(
                        canonical_taxon(gtdb_species))

                    if highest_suffix is None or self.taxon_suffix_manager.is_higher_suffix(suffix, highest_suffix):
                        validated_count += 1
                    else:
                        invalid_suffix[rid] = (
                            rid, gtdb_species, prev_species, 'NO_PRIOR_SPECIES')
                elif prev_suffix is None:
                    # this should only occur if a NCBI species name is now considered ambiguous,
                    # so all instances must be suffixed
                    if rid not in ambiguous_rids:
                        invalid_suffix[rid] = (
                            rid, gtdb_species, prev_species, 'NO_PRIOR_SUFFIX')
                    else:
                        validated_count += 1
                else:
                    # suffix should be the same between releases
                    if suffix != prev_suffix:
                        invalid_suffix[rid] = (
                            rid, gtdb_species, prev_species, 'SUFFIX DOES NOT MATCH')
                    else:
                        validated_count += 1

        print('Identified {:,} genomes with a suffixed specific name.'.format(
            suffix_count))
        self.report_validation(
            invalid_suffix,
            validated_count,
            ' - identified {:,} genomes with an suffixed specific name that disagrees with previous GTDB assignment (Proposed GTDB species, Previous GTDB species, Case):'.format(len(invalid_suffix)))

    def validate_specific_placeholder_suffixes(self,
                                               final_taxonomy,
                                               cur_genomes,
                                               unambiguous_ncbi_sp,
                                               unambiguous_ncbi_subsp):
        """Validate suffix assignment rules for specific names."""

        unambiguous_sp_rids = set()
        for rid, _assignment_type in unambiguous_ncbi_sp.values():
            unambiguous_sp_rids.add(rid)

        unambiguous_subsp_rids = set()
        for rid, _assignment_type in unambiguous_ncbi_subsp.values():
            unambiguous_subsp_rids.add(rid)

        # validate that unambiguous NCBI names are not suffixed and that all
        # ambiguous NCBI names are placeholder names
        invalid_suffix = {}
        validated_count = 0
        for rid in final_taxonomy:
            if rid in self.mc_species:
                continue

            gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            gtdb_specific = specific_epithet(gtdb_species)
            ncbi_species = cur_genomes[rid].ncbi_taxa.species

            if rid in unambiguous_sp_rids:
                if not is_latin_sp_epithet(gtdb_specific):
                    invalid_suffix[rid] = (
                        rid, gtdb_species, ncbi_species, 'UNAMBIGUOUS_SPECIES_WITH_SUFFIX')
                else:
                    validated_count += 1
            elif rid in unambiguous_subsp_rids:
                if not is_latin_sp_epithet(gtdb_specific):
                    invalid_suffix[rid] = (
                        rid, gtdb_species, ncbi_species, 'UNAMBIGUOUS_SUBSPECIES_WITH_SUFFIX')
                else:
                    validated_count += 1
            else:
                if not is_placeholder_sp_epithet(gtdb_specific):
                    invalid_suffix[rid] = (
                        rid, gtdb_species, ncbi_species, 'AMBIGUOUS_SPECIES_WITH_LATIN_SPECIFIC_NAME')
                else:
                    validated_count += 1

        self.report_validation(
            invalid_suffix,
            validated_count,
            ' - identified {:,} genomes which violate the suffixing rules for specific names (Proposed GTDB species, NCBI species, Case):'.format(len(invalid_suffix)))

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
                        invalid_suffix[rid] = (rid, gtdb_species, prev_species)
                else:
                    # suffix should be the same between releases
                    if suffix != prev_suffix:
                        invalid_suffix[rid] = (rid, gtdb_species, prev_species)
                    else:
                        validated_count += 1

        print('Identified {:,} genomes with a suffixed specific name.'.format(
            suffix_count))
        self.report_validation(
            invalid_suffix,
            validated_count,
            ' - identified {:,} genomes with a modified placeholder suffix which should be manually verified (Proposed GTDB species, Previous GTDB species):'.format(len(invalid_suffix)))

    def validate_genera_spelling(self, final_taxonomy, cur_genomes):
        """Validate spelling of genus names."""

        invalid_spelling = {}
        validated_count = 0
        for rid in final_taxonomy:
            gtdb_genus = final_taxonomy[rid][Taxonomy.GENUS_INDEX]
            canonical_genus = canonical_taxon(gtdb_genus)

            prev_gtdb_genus = cur_genomes[rid].gtdb_taxa.genus
            canonical_prev_genus = canonical_taxon(prev_gtdb_genus)

            if (canonical_genus != canonical_prev_genus
                    and lvn_distance(canonical_genus, canonical_prev_genus) <= 2
                    and not (canonical_genus.startswith('g__UBA') and canonical_prev_genus.startswith('g__UBA'))):
                # Names are highly similar, but not identical so may be a typo.
                # UBA genera are ignored since these appear similar to each other so are often flagged as false positives,
                # given that valid genus transfers do occur.
                invalid_spelling[gtdb_genus] = [gtdb_genus, prev_gtdb_genus]
            else:
                validated_count += 1

        self.report_validation(
            invalid_spelling,
            validated_count,
            ' - identified {:,} genera with potential spelling error (Proposed GTDB genus, Previous GTDB genus):'.format(len(invalid_spelling)))

    def validate_specific_no_suffix(self, final_taxonomy, cur_genomes, unambiguous_ncbi_sp, unambiguous_ncbi_subsp):
        """Validate GTDB specific name assignments without a polyphyletic suffixes."""

        unambiguous_sp_rids = set()
        for rid, _assignment_type in unambiguous_ncbi_sp.values():
            unambiguous_sp_rids.add(rid)

        unambiguous_subsp_rids = set()
        for rid, _assignment_type in unambiguous_ncbi_subsp.values():
            unambiguous_subsp_rids.add(rid)

        # determine number of GTDB species with same canonical name
        canonical_sp_count = defaultdict(set)
        for rid in final_taxonomy:
            gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            canonical_sp = canonical_taxon(gtdb_species)

            canonical_sp_count[canonical_sp].add(gtdb_species)

        invalid_names = {}
        validated_count = 0
        case_count = 0
        for rid in final_taxonomy:
            if rid in self.mc_species:
                continue

            gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            canonical_sp = canonical_taxon(gtdb_species)
            gtdb_specific = specific_epithet(gtdb_species)

            suffix = taxon_suffix(gtdb_specific)
            if suffix is not None:
                continue

            case_count += 1

            if (cur_genomes[rid].is_effective_type_strain()
                    or rid in unambiguous_sp_rids
                    or rid in unambiguous_subsp_rids
                    or len(canonical_sp_count[canonical_sp]) == 1):
                validated_count += 1
            else:
                invalid_names[rid] = (rid,
                                      gtdb_species,
                                      cur_genomes[rid].is_effective_type_strain(),
                                      ', '.join(canonical_sp_count[canonical_sp]))

        print(
            'Identified {:,} genomes without a suffixed specific name.'.format(case_count))
        self.report_validation(
            invalid_names,
            validated_count,
            ' - identified {:,} genomes which should have a polyphyletic suffix on specific name (Proposed GTDB species, effective type material, GTDB canonical species):'.format(len(invalid_names)))

    def validate_forbidden_specific_names(self, final_taxonomy):
        """Validate that no species names contain a forbidden specific name."""

        invalid_specific = {}
        validated_count = 0
        for gid, taxa in final_taxonomy.items():
            gtdb_species = taxa[Taxonomy.SPECIES_INDEX]
            specific = specific_epithet(gtdb_species)

            if specific in self.forbidden_specific_names:
                invalid_specific[gid] = gtdb_species
            else:
                validated_count += 1

        self.report_validation(
            invalid_specific,
            validated_count,
            ' - identified {:,} genomes with a forbidden specific name:'.format(len(invalid_specific)))

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
                invalid_sp_name[gid] = (gid, gtdb_species, ncbi_species)
            else:
                validated_count += 1

        self.report_validation(
            invalid_sp_name,
            validated_count,
            ' - identified {:,} genomes with a non-suffixed, placeholder generic name, but Latin specific name (Proposed GTDB species, NCBI species):'.format(len(invalid_sp_name)))

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

        self.log.info(
            ' - identified {:,} NCBI genera.'.format(len(ncbi_genera_gids)))

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
                ncbi_genus_in_gtdb_genus_perc = len(
                    ncbi_genus_in_gtdb_genus)*100.0/max(1, len(ncbi_gids))
                if ncbi_genus_in_gtdb_genus_perc > 50:
                    gtdb_mergers[gtdb_genus].append(ncbi_genus)
                    merged_ncbi_genera += 1
        self.log.info(
            ' - identified {:,} NCBI genera that were merged with GTDB genera.'.format(merged_ncbi_genera))

        # establish percentage of genomes descendant from each GTDB genus
        # with the same NCBI name
        fout = open(os.path.join(self.output_dir,
                                 'gtdb_genus_monophyly_test.tsv'), 'w')
        fout.write('GTDB genera\tCanonical GTDB genera\tType species of genus\tMerged genera\t% cluster with NCBI genus\t% NCBI genus in cluster\tNCBI genus assignments\n')
        gtdb_genus_minority_support = 0
        for gtdb_genus, gtdb_gids in gtdb_genera_gids.items():
            canonical_gtdb_genus = canonical_taxon(gtdb_genus)

            ncbi_genus_gids = ncbi_genera_gids[canonical_gtdb_genus]

            ncbi_genus_in_gtdb_genus = gtdb_gids.intersection(ncbi_genus_gids)
            ncbi_genus_in_gtdb_genus_perc = len(
                ncbi_genus_in_gtdb_genus)*100.0/max(1, len(ncbi_genus_gids))

            gtdb_gids_with_ncbi_genus = gtdb_gids.intersection(
                gids_with_ncbi_genus)
            gtdb_gids_from_cur_ncbi_genus_perc = len(
                ncbi_genus_in_gtdb_genus)*100.0/max(1, len(gtdb_gids_with_ncbi_genus))

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
                            self.log.warning('GTDB genus {} contains type species for {}.'.format(
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

        self.log.info(
            ' - identified {:,} GTDB genera with <50% congruence with underlying NCBI genera.'.format(gtdb_genus_minority_support))

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
        self.log.info(
            ' - median RED of genera: {:.3f}'.format(median_red_genera))

        # get all genera present in the GTDB
        gtdb_genera = set([final_taxonomy[rid][Taxonomy.GENUS_INDEX]
                           for rid in final_taxonomy])

        # get NCBI genera not present in the GTDB
        missing_ncbi_genera = defaultdict(list)
        for rid in final_taxonomy:
            ncbi_genus = cur_genomes[rid].ncbi_taxa.genus
            if ncbi_genus != 'g__' and ncbi_genus not in gtdb_genera:
                missing_ncbi_genera[ncbi_genus].append(rid)
        self.log.info(
            ' - found {:,} NCBI genera not represented in the GTDB.'.format(len(missing_ncbi_genera)))

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
                    gtdb_genera = ['{}: {} (RED={:.2f})'.format(
                        rid, final_taxonomy[rid][Taxonomy.GENUS_INDEX], red_genus[final_taxonomy[rid][Taxonomy.GENUS_INDEX]]) for rid in rids]
                    invalid_missing_ncbi_genera[ncbi_genus] = (
                        ncbi_genus, ', '.join(gtdb_genera))
                else:
                    validated_count += 1

        self.report_validation(
            invalid_missing_ncbi_genera,
            validated_count,
            ' - identified {:,} NCBI genera without justification for being absent in GTDB (NCBI genus, GTDB genera):'.format(len(invalid_missing_ncbi_genera)))

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
        fout = open(os.path.join(self.output_dir,
                                 'latin_names_for_gtdb_genera.tsv'), 'w')
        fout.write(
            'Proposed GTDB genus\tNCBI genus\tPropostion of NCBI genus genomes\tNote\n')
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
                perc_in_gtdb_genus = count*100.0 / \
                    len(ncbi_genus_gids[top_ncbi_genera])

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
                        note.append('Contains {:,} genomes with Latin NCBI species name.'.format(
                            top_ncbi_sp_count))

                    perc_str = '{} of {} ({:.1f}%)'.format(count,
                                                           len(
                                                               ncbi_genus_gids[top_ncbi_genera]),
                                                           perc_in_gtdb_genus)
                    invalid_placeholder[gtdb_genus] = (gtdb_genus,
                                                       top_ncbi_genera,
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
            ' - identified {:,} GTDB placeholder genera with candidate Latin names (Proposed GTDB genus, NCBI genus, Proportion of NCBI genus genomes):'.format(len(invalid_placeholder)))

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
        cur_gtdb_genera = set(
            [final_taxonomy[rid][Taxonomy.GENUS_INDEX] for rid in final_taxonomy])
        prev_gtdb_genera = set(
            [prev_genomes[rid].gtdb_taxa.genus for rid in prev_genomes.sp_clusters])
        self.log.info(' - identified {:,} current and {:,} previous GTDB genera.'.format(
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
                        ncbi_genus = 'g__' + \
                            generic_name(cur_genomes[rid].ncbi_taxa.species)
                        if canonical_taxon(cur_genus) == ncbi_genus:
                            # necessarily modified to reflect NCBI species assignment
                            # (e.g., FW-11 to Sphingomonas_H as representative is now S. oleivorans at NCBI)
                            likely_valid_change = True

                if likely_valid_change:
                    validated_count += 1
                else:
                    invalid_placeholder[prev_genus] = (
                        prev_genus, cur_genus, len(rids))

        self.report_validation(
            invalid_placeholder,
            validated_count,
            ' - identified {:,} genomes with invalid change to placeholder genus/generic name (Previous GTDB genus, Proposed GTDB genus, No. species changed):'.format(len(invalid_placeholder)))

    def validate_unambiguous_ncbi_sp(self, final_taxonomy,
                                     unambiguous_ncbi_sp,
                                     lpsn):
        """'Validating NCBI species names classified as having an unambiguous assignment within the GTDB."""

        invalid_unambiguous = {}
        validated_count = 0
        for ncbi_sp, (rid, _assignment_type) in unambiguous_ncbi_sp.items():
            if rid in self.mc_species:
                continue

            if rid not in final_taxonomy:
                # must be from other domain
                continue

            gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            gtdb_specific = specific_epithet(gtdb_species)

            ncbi_specific = specific_epithet(ncbi_sp)

            if not test_same_epithet(gtdb_specific, ncbi_specific):
                invalid_unambiguous[rid] = (
                    rid, gtdb_species, ncbi_sp, ','.join(lpsn.sp_synonyms[ncbi_sp]))
            else:
                validated_count += 1

        self.report_validation(
            invalid_unambiguous,
            validated_count,
            ' - identified {:,} unambiguous NCBI species with an invalid GTDB species name (Proposed GTDB species, NCBI classification, LPSN synonym):'.format(
                len(invalid_unambiguous)))

    def validate_ambiguous_ncbi_sp(self, final_taxonomy,
                                   cur_genomes,
                                   unambiguous_ncbi_sp,
                                   unambiguous_ncbi_subsp):
        """'Validating NCBI species names classified as having an ambiguous assignment within the GTDB."""

        unambiguous_sp_rids = set()
        for rid, _assignment_type in unambiguous_ncbi_sp.values():
            unambiguous_sp_rids.add(rid)

        unambiguous_subsp_rids = set()
        for rid, _assignment_type in unambiguous_ncbi_subsp.values():
            unambiguous_subsp_rids.add(rid)

        invalid_ambiguous = {}
        validated_count = 0
        for gid in final_taxonomy:
            if gid in self.mc_species:
                continue

            if (gid in unambiguous_sp_rids
                    or gid in unambiguous_subsp_rids):
                continue

            gtdb_species = final_taxonomy[gid][Taxonomy.SPECIES_INDEX]
            gtdb_specific = specific_epithet(gtdb_species)

            if not is_placeholder_sp_epithet(gtdb_specific):
                ncbi_sp = cur_genomes[gid].ncbi_taxa.species
                invalid_ambiguous[gid] = (gid, gtdb_species, ncbi_sp)
            else:
                validated_count += 1

        self.report_validation(
            invalid_ambiguous,
            validated_count,
            ' - identified {:,} ambiguous NCBI species with an invalid GTDB species name (Proposed GTDB species, NCBI classification):'.format(
                len(invalid_ambiguous)))

    def validate_misclassified_genomes(self,
                                       final_taxonomy,
                                       cur_genomes,
                                       ncbi_synonyms,
                                       ncbi_misclassified_gids,
                                       unambiguous_ncbi_sp,
                                       unambiguous_ncbi_subsp):
        """Validating species names of genomes misclassified at NCBI."""

        unambiguous_sp_rids = set()
        for rid, _assignment_type in unambiguous_ncbi_sp.values():
            unambiguous_sp_rids.add(rid)

        unambiguous_subsp_rids = set()
        for rid, _assignment_type in unambiguous_ncbi_subsp.values():
            unambiguous_subsp_rids.add(rid)

        # validate that any misclassified genomes that are (unfortunately)
        # GTDB species representatives have a placeholder name
        invalid_misclassified = {}
        validated_count = 0
        for gid in ncbi_misclassified_gids:
            if gid in self.mc_species:
                continue

            if gid in final_taxonomy:
                gtdb_species = final_taxonomy[gid][Taxonomy.SPECIES_INDEX]
                gtdb_specific = specific_epithet(gtdb_species)

                if gid in unambiguous_sp_rids or gid in unambiguous_subsp_rids:
                    # this occurs when a misclassified genome is the representative of a
                    # unambiguous NCBI species (e.g., G003263425 is representative of
                    # species cluster with the GTDB name s__Pseudomonas inefficax since this
                    # cluster has a unanimous consensus regarding this name once misclassification
                    # are taken into account). As such, in this case the cluster has a Latin
                    # name and this is considered correct.
                    if is_placeholder_sp_epithet(gtdb_specific):
                        ncbi_species = cur_genomes[gid].ncbi_taxa.species
                        ncbi_species = ncbi_synonyms.get(
                            ncbi_species, ncbi_species)
                        invalid_misclassified[gid] = (
                            gid, gtdb_species, ncbi_species)
                    else:
                        validated_count += 1
                elif not is_placeholder_sp_epithet(gtdb_specific):
                    ncbi_species = cur_genomes[gid].ncbi_taxa.species
                    ncbi_species = ncbi_synonyms.get(
                        ncbi_species, ncbi_species)
                    invalid_misclassified[gid] = (
                        gid, gtdb_species, ncbi_species)
                else:
                    validated_count += 1

        self.report_validation(
            invalid_misclassified,
            validated_count,
            ' - identified {:,} misclassified genomes with an invalid specific name (Proposed GTDB species, NCBI classification):'.format(
                len(invalid_misclassified)))

    def _get_genus_stem(self, genus):
        """Get stem of genus."""

        # Determining the stem for a genus is  non-trivial so we use a simple
        # heuristic which is to drop the last 2 letters, except for placeholder
        # names where the full name is the stem.

        if is_placeholder_taxon(genus):
            return genus[3:]
        else:
            return genus[3:-2]

    def validate_merged_taxa(self, final_taxonomy, cur_genomes, lpsn):
        """Validate name of merged taxa."""

        # get list of genomes assembled from the type species of the genus
        # for each GTDB taxon
        gtdb_taxon_gids = defaultdict(set)
        for rid in final_taxonomy:
            if cur_genomes[rid].is_gtdb_type_species():
                for taxon in final_taxonomy[rid]:
                    gtdb_taxon_gids[taxon].add(rid)

        # check if GTDB taxa is a merger of multiple NCBI taxa
        fout = open(os.path.join(self.output_dir,
                                 'validate_merged_taxa.tsv'), 'w')
        fout.write(
            'Case\tGTDB taxon\tPriority taxon\tYear of priority\tTaxa priority\tGTDB representatives\n')
        num_cases = 0
        for gtdb_taxon, type_gids in gtdb_taxon_gids.items():
            # Validate families and orders
            if gtdb_taxon[0:3] not in set(['f__', 'o__']):
                continue

            ncbi_taxa_with_type = {}
            for gid in type_gids:
                if gtdb_taxon[0:3] == 'f__':
                    ncbi_cur_taxon = cur_genomes[gid].ncbi_taxa.family
                elif gtdb_taxon[0:3] == 'o__':
                    ncbi_cur_taxon = cur_genomes[gid].ncbi_taxa.order

                ncbi_type_genus = cur_genomes[gid].ncbi_taxa.genus
                lpsn_type_genus = lpsn.type_material(ncbi_cur_taxon)

                if lpsn_type_genus == ncbi_type_genus:
                    ncbi_taxa_with_type[ncbi_cur_taxon] = gid

            if len(ncbi_taxa_with_type) >= 2:
                # check if GTDB name follows priority of merged taxa
                if gtdb_taxon not in ncbi_taxa_with_type:
                    print('GTDB taxon {} does not reflect descendant type material: {}'.format(
                        gtdb_taxon,
                        ', '.join(ncbi_taxa_with_type)))
                else:
                    priority_taxa = []
                    priority_year = 1e6
                    ncbi_taxa_priority = {}
                    for ncbi_taxon in ncbi_taxa_with_type:
                        if ncbi_taxon in lpsn.taxa:
                            priority = lpsn.taxa[ncbi_taxon].priority_year
                            ncbi_taxa_priority[ncbi_taxon] = priority
                            if priority < priority_year:
                                priority_taxa = [ncbi_taxon]
                                priority_year = priority
                            elif priority == priority_year:
                                priority_taxa.append(ncbi_taxon)
                                priority_year = priority

                    if len(priority_taxa) == 1:
                        if priority_taxa[0] != gtdb_taxon:
                            num_cases += 1
                            print('Incorrect GTDB name?',
                                  gtdb_taxon,
                                  priority_taxa[0],
                                  priority_year,
                                  ncbi_taxa_priority,
                                  ncbi_taxa_with_type)

                            fout.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                'Priority violation?',
                                gtdb_taxon,
                                priority_taxa[0],
                                priority_year,
                                ncbi_taxa_priority,
                                ncbi_taxa_with_type))

                    elif len(priority_taxa) > 1:
                        num_cases += 1
                        print('Type material has equal priority (manual curation required)',
                              gtdb_taxon,
                              priority_taxa,
                              priority_year,
                              ncbi_taxa_priority)

                        fout.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            'Equal priority (manual validation required)',
                            gtdb_taxon,
                            priority_taxa,
                            priority_year,
                            ncbi_taxa_priority,
                            ncbi_taxa_with_type))

        fout.close()

        self.log.info(
            f' - identified {num_cases:,} cases (see validate_merged_taxa.tsv)')

    def validate_parent_child_placeholder_names(self, final_taxonomy):
        """Validate parent-child placeholder names for simple, unnecessary changes.

        This test primarily aims to catch simply typos where the parent and child
        should in fact have identical names.

        Example: f__JDFR-13 conflicts with g__JdFR-13
        """

        invalid_name = {}
        validated_count = 0
        for gid, taxa in final_taxonomy.items():
            for rank_idx, parent in enumerate(taxa[0:Taxonomy.GENUS_INDEX]):
                child = taxa[rank_idx+1]

                if parent[3:] != child[3:] and parent[3:].lower() == child[3:].lower():
                    # identical except for a change in captialization
                    invalid_name[gid] = (gid, parent, child)

        self.report_validation(
            invalid_name,
            validated_count,
            ' - identified {:,} GTDB placeholder names that may conflict with child name (Genome ID, Parent name, Child name):'.format(
                len(invalid_name)))

    def validate_gtdb_spelling(self, final_taxonomy, cur_genomes):
        """Validate spelling of GTDB names by comparison with NCBI names."""

        invalid_spelling = {}
        validated_count = 0
        for gid, gtdb_taxa in final_taxonomy.items():
            if gid not in cur_genomes:
                continue

            for idx, gtdb_taxon in enumerate(gtdb_taxa):
                ncbi_taxon = cur_genomes[gid].ncbi_taxa.get_taxa(idx)
  
                if len(ncbi_taxon) == 3 or gtdb_taxon.startswith('s__'):
                    continue

                if (not gtdb_taxon.startswith(ncbi_taxon) 
                    and 0 < lvn_distance(gtdb_taxon, ncbi_taxon) <= 2):
                    # names are highly similar, but not identical so may be a typo
                    invalid_spelling[gid] = [gid, gtdb_taxon, ncbi_taxon]
            else:
                validated_count += 1

        derep_invalid_spelling = {}
        processed = set()
        for (gid, gtdb_taxon, ncbi_taxon) in invalid_spelling.values():
            pair = (gtdb_taxon, ncbi_taxon)
            if pair in processed:
                continue

            processed.add(pair)

            derep_invalid_spelling[gid] = [gid, gtdb_taxon, ncbi_taxon]

        self.report_validation(
            derep_invalid_spelling,
            validated_count,
            ' - identified {:,} potential spelling error (Genome ID, GTDB taxon, NCBI taxon):'.format(len(derep_invalid_spelling)))

    def validate_placeholder_names_not_latin(self, final_taxonomy, lpsn):
        """Validate placeholder names that might be mistaken as Latin names.

        Placeholder names should always have a number or additional capitalized letters
        to clear indicate they are not a Latin name. Catching cases violating these rules
        is challenging even though they are often obvious to humans (e.g. g__Notlatin).
        Right now, this test simply looks for names under 5 characters that appear to be
        Latin under the assumption this is too short for valid Latin names.

        Example: g__Gsub
        """

        invalid_name = {}
        validated_count = 0
        for gid, taxa in final_taxonomy.items():
            for taxon in taxa:
                if taxon in lpsn.taxa:
                    continue

                if not is_placeholder_taxon(taxon) and len(taxon[3:]) <= 5:
                    # Latin-looking name that is under 6 characters
                    invalid_name[taxon] = (gid, taxon)
                else:
                    validated_count += 1

        self.report_validation(
            invalid_name,
            validated_count,
            " - identified {:,} GTDB Latin-looking names that are short and may in fact be placehold names (GTDB Taxon):".format(
                len(invalid_name)))

    def validate_sp_with_alphanumeric_generic_names(self, final_taxonomy):
        """Validate species with alphanumeric generic names do not have suffixed Latin specific names.

        Examples:
          s__DTPE01 lauensis_A
          s__DRZC01 fontis_A
        """

        invalid_name = {}
        validated_count = 0
        for gid, taxa in final_taxonomy.items():
            genus = taxa[Taxonomy.GENUS_INDEX]
            species = taxa[Taxonomy.SPECIES_INDEX]
            specific = specific_epithet(species)

            if not species.startswith(genus.replace('g__', 's__')):
                print(f'Generic name does not match genus name: {gid}, {genus}, {species}')
                continue

            if is_alphanumeric_taxon(genus) and is_suffixed_sp_epithet(specific):
                # Latin-looking name that is under 6 characters
                invalid_name[gid] = (gid, species)
            else:
                validated_count += 1

        self.report_validation(
            invalid_name,
            validated_count,
            " - identified {:,} GTDB species names with an alphanumeric generic name and suffixed Latin specific name:".format(
                len(invalid_name)))

    def validate_suffix_of_specific_names(self, final_taxonomy, cur_genomes, lpsn):
        """Validate suffix of specific names by comparison with LPSN names.

        Example: Lacrimispora algidixylanolyticum vs. L. algidixylanolytica
        """

        # get LPSN names in genus
        lpsn_species_in_genus = defaultdict(set)
        for lpsn_taxon in lpsn.taxa:
            if lpsn_taxon.startswith('s__'):
                lpsn_genus = lpsn_taxon.split()[0].replace('s__', 'g__')
                lpsn_species_in_genus[lpsn_genus].add(lpsn_taxon)

        invalid_suffix = {}
        validated_count = 0
        for gid, taxa in final_taxonomy.items():
            gtdb_genus = taxa[Taxonomy.GENUS_INDEX]
            gtdb_sp = taxa[Taxonomy.SPECIES_INDEX]
            gtdb_specific = gtdb_sp.split()[1]

            # check if name is recognized at LPSN
            if gtdb_sp in lpsn_species_in_genus[gtdb_genus]:
                validated_count += 1
                continue

            # check if a name that only differs in the suffix of
            # the specific name is recognized at LPSN
            for lpsn_sp in lpsn_species_in_genus[gtdb_genus]:
                lpsn_specific = lpsn_sp.split()[1]

                if test_same_epithet(lpsn_specific, gtdb_specific):
                    invalid_suffix[gid] = (
                        gid, gtdb_sp, cur_genomes[gid].ncbi_taxa.species, lpsn_sp)

        self.report_validation(
            invalid_suffix,
            validated_count,
            ' - identified {:,} GTDB species with specific name with suffix that differs slightly from LPSN name (Genome ID, GTDB species, NCBI species, matched LPSN species):'.format(
                len(invalid_suffix)))

    def validate_taxa_by_lpsn_type_material(self, final_taxonomy, lpsn):
        """Validating placement of taxon based on LPSN type material."""

        # get all GTDB taxa and GTDB taxa descendant from each taxon
        gtdb_all_taxa = set()
        children_taxa = defaultdict(set)
        gtdb_lineage = {}
        for taxa in final_taxonomy.values():
            for idx, taxon in enumerate(taxa):
                gtdb_all_taxa.add(taxon)
                children_taxa[taxon].update(taxa[idx+1:])
                gtdb_lineage[taxon] = taxa[0:idx]

        # validate placement of GTDB taxa
        invalid_taxon = {}
        validated_count = 0

        for taxon in gtdb_all_taxa:
            if taxon.startswith('s__'):
                # can't validate species using this approach
                continue

            if taxon not in lpsn.taxa:
                # can't validate taxa without LPSN type material
                continue

            type_taxon = lpsn.taxa[taxon].type_material

            if type_taxon not in gtdb_all_taxa:
                # can't validate taxa where GTDB does not
                # define the type material
                continue

            if type_taxon not in children_taxa[taxon]:
                # type material is in a different lineage which
                # is not valid
                row = (taxon, type_taxon, '; '.join(gtdb_lineage[type_taxon]))
                invalid_taxon[taxon] = row
            else:
                # taxon contain type material as expected
                validated_count += 1

        self.report_validation(
            invalid_taxon,
            validated_count,
            ' - identified {:,} GTDB taxa with type material in a different lineage (Taxon, Type material, Type lineage):'.format(
                len(invalid_taxon)))

    def validate_taxa_by_genus_stem(self, final_taxonomy):
        """Validating placement of higher taxon names derived from stem of GTDB genera."""

        # Common stems that result in false positives
        # (e.g., Nitrospi from g__Nitrospina and g__Nitrospira)
        common_genus_stems = set(['Nitrospi'])

        # Find stem of all genus names and ensure higher order taxa derived from this
        # stem contain the genus.Example: stem of Abditibacterium is Abditibacteri
        # which matches with Abditibacteriales, so Abditibacteriales should contain
        # Abditibacterium. Note that determining the stem for a genus is actually
        # non-trivial so we just use a simple heuristic which is to drop the
        # last 2 letters and then to match this with the start of higher taxa.
        # We also deal with placeholder names, where the only exception is that
        # 2 letters are not dropped.
        genus_stems = {}
        genus_lineages = {}
        all_higher_taxa = defaultdict(set)
        for taxa in final_taxonomy.values():
            genus = taxa[Taxonomy.GENUS_INDEX]
            if genus == 'g__' or genus in genus_stems:
                continue

            if is_suffixed_taxon(genus):
                # suffixed genera such as Bacillus_A or BOG42_A,
                # are never the stems for higher taxon names
                continue

            for taxon in taxa[Taxonomy.PHYLUM_INDEX:Taxonomy.GENUS_INDEX]:
                if is_suffixed_taxon(taxon):
                    # suffixed taxon such as o__Thiohalomonadales_A or f__UBA11359_B
                    # do not speak to the correct placement of the stem genus so are
                    # skipped
                    continue
                all_higher_taxa[taxon].add(genus)

            stem = self._get_genus_stem(genus)
            if stem in common_genus_stems:
                continue

            genus_stems[genus] = stem
            genus_lineages[genus] = taxa[0:Taxonomy.SPECIES_INDEX]

        # test the placement of higher order taxa based on genus name
        taxon_suffixes = {'o__': 'ales', 'f__': 'aceae'}
        invalid_taxon = {}
        validated_count = 0
        for genus, genus_stem in genus_stems.items():
            passed_test = True
            for taxon, taxon_genera in all_higher_taxa.items():
                if genus in taxon_genera:
                    continue

                # check if taxon could be derived from one of the genera it contains
                valid = True
                for g in taxon_genera:
                    if g in genus_stems and taxon.startswith(genus_stems[g]):
                        valid = True
                        break

                if valid:
                    continue

                if is_placeholder_taxon(genus) and taxon[3:] == genus_stem:
                    # higher taxon appears to be derived for a placeholder GTDB genus
                    # name, but does not contain the genus
                    invalid_taxon[taxon] = (
                        taxon, genus, genus_stem, '; '.join(genus_lineages[genus]))
                    passed_test = False
                elif not is_placeholder_taxon(genus) and taxon[3:].startswith(genus_stem):
                    # higher taxon appears to be derived for a Latin GTDB genus
                    # name, but does not contain the genus. This results in a
                    # far number of false positives due to common stems (e.g.
                    # g__Thermus has the stem Therm and matches f__Thermopetrobacteraceae).
                    # To avoid these false positives a length check is done to ensure the
                    # higher taxon name appears to be derived from the genus.
                    if taxon[0:3] in ['p__', 'c__']:
                        # no set phylum or class suffix so chomp last 2 characters
                        # as a simple heuristic
                        taxon_stem = taxon[:-2]
                    else:
                        taxon_stem = taxon.replace(
                            taxon_suffixes[taxon[0:3]], '')

                    if 0 <= len(taxon_stem[3:]) - len(genus_stem) <= 3:
                        # taxon stem is only slightly longer than the
                        # genus stem we can't discount it being derived
                        # from the genus
                        row = (taxon, genus, genus_stem,
                               '; '.join(genus_lineages[genus]))
                        invalid_taxon[taxon] = row
                        passed_test = False

            if passed_test:
                validated_count += 1

        self.report_validation(
            invalid_taxon,
            validated_count,
            ' - identified {:,} taxon names that appear to be derived from a GTDB genus name in a different lineage (Taxon, Matching genus, Matching stem, Genus lineage):'.format(
                len(invalid_taxon)))

    def justify_specific_name_change(self,
                                     prev_rid,
                                     new_rid,
                                     prev_specific,
                                     new_specific,
                                     final_taxonomy,
                                     cur_genomes,
                                     prev_genomes,
                                     cur_clusters,
                                     unambiguous_sp_rids,
                                     unambiguous_ncbi_sp,
                                     unambiguous_subsp_rids,
                                     ncbi_misclassified_gids,
                                     gtdb_synonyms):
        """Provide justification for why specific name of species was changed."""

        specific_change_types = []
        reasons = []

        prev_rid_prev_ncbi_species = prev_genomes[prev_rid].ncbi_taxa.species
        prev_rid_prev_ncbi_specific = specific_epithet(
            prev_rid_prev_ncbi_species)

        prev_rid_new_ncbi_species = None
        if prev_rid in cur_genomes:
            prev_rid_new_ncbi_species = cur_genomes[prev_rid].ncbi_taxa.species

        prev_gtdb_species = prev_genomes[prev_rid].gtdb_taxa.species
        prev_canonical_gtdb_species = canonical_species(prev_gtdb_species)

        prev_rids_in_cluster = [
            cid for cid in cur_clusters[new_rid] if cid in prev_genomes.sp_clusters]

        new_ncbi_species = cur_genomes[new_rid].ncbi_taxa.species
        new_gtdb_species = final_taxonomy[new_rid][Taxonomy.SPECIES_INDEX]

        if new_rid in self.mc_species:
            specific_change_types.append('MANUAL_CURATION')
            reasons.append('Changed by curators')

        if (prev_rid_new_ncbi_species in gtdb_synonyms
                and gtdb_synonyms[prev_rid_new_ncbi_species] == new_gtdb_species):
            specific_change_types.append('GTDB_SYNONYM')
            reasons.append('Previous species is a synonym of new species')

        if test_same_epithet(prev_specific, new_specific):
            specific_change_types.append('GENUS_TRANSFER_SUFFIX_CHANGE')
            reasons.append('Changed suffix to reflect generic name of species')

        if cur_genomes[new_rid].is_effective_type_strain():
            specific_change_types.append('EFFECTIVE_TYPE_STRAIN')
            reasons.append(
                'Changed to reflect name of effective type strain of species')

        if (new_rid not in unambiguous_sp_rids
                and specific_epithet_type(prev_specific) == 'LATIN'
                and specific_epithet_type(new_specific) == 'SUFFIXED_LATIN'):
            specific_change_types.append('AMBIGUOUS_SPECIES')
            reasons.append(
                'Added alphabetic suffix as correct placement of species is ambiguous')

        if (new_rid not in unambiguous_sp_rids
                and specific_epithet_type(prev_specific) == 'LATIN'
                and specific_epithet_type(new_specific) == 'ALPHANUMERIC'
                and new_ncbi_species == 's__'):
            specific_change_types.append('AMBIGUOUS_SPECIES_UNCLASSIFIED_REP')
            reasons.append(
                'Changed to alphanumeric specific name as placement of species is ambiguous and GTDB representative has no NCBI species assignment')

        if prev_rid in ncbi_misclassified_gids:
            specific_change_types.append('PREVIOUS_REP_MISCLASSIFIED')
            reasons.append(
                'NCBI species assignment of previous GTDB representative considered misclassification under the GTDB')

        if new_rid in ncbi_misclassified_gids:
            specific_change_types.append('CURRENT_REP_MISCLASSIFIED')
            reasons.append(
                'NCBI species assignment of current GTDB representative considered misclassification under the GTDB')

        if (new_rid in unambiguous_sp_rids
                and unambiguous_sp_rids[new_rid] == NCBI_SpeciesManager.MAJORITY_VOTE
                and specific_epithet_type(new_specific) == 'LATIN'):
            specific_change_types.append('MAJORITY_VOTE_ASSIGNMENT')
            reasons.append(
                'Changed to reflect consensus species assignment of genomes in cluster')

        if (new_rid in unambiguous_subsp_rids
                and unambiguous_subsp_rids[new_rid] == new_specific):
            specific_change_types.append('PROMOTED_SUBSPECIES')
            reasons.append(
                'Changed to reflect consensus subspecies assignment of genomes in cluster')

        if (prev_canonical_gtdb_species in unambiguous_ncbi_sp
                and new_rid != unambiguous_ncbi_sp[prev_canonical_gtdb_species][0]
                and unambiguous_ncbi_sp[prev_canonical_gtdb_species][1] == NCBI_SpeciesManager.TYPE_STRAIN_GENOME
                and is_latin_sp_epithet(prev_specific)):
            specific_change_types.append(
                'CONFLICTS_WITH_EFFECTIVE_TYPE_STRAIN')
            reasons.append(
                'Effective type strain genome of previous species assignment identified in different GTDB species cluster')

        if (prev_rid_prev_ncbi_species != prev_rid_new_ncbi_species
                and prev_rid_new_ncbi_species is not None):
            specific_change_types.append('NCBI_REASSIGNED_SPECIES')
            reasons.append(
                'NCBI species assignment changed for previous GTDB representative')

        if prev_rid == new_rid:
            if (prev_rid_prev_ncbi_specific == prev_specific
                    and new_ncbi_species == 's__'
                    and specific_epithet_type(new_specific) == 'ALPHANUMERIC'):
                specific_change_types.append('NCBI_RETRACTED_SPECIES')
                reasons.append(
                    'Changed to reflect genome no longer having a species assignment at NCBI')

        if (len(prev_rids_in_cluster) > 1
                and specific_epithet_type(prev_specific) == 'ALPHANUMERIC'
                and specific_epithet_type(new_specific) == 'ALPHANUMERIC'):
            specific_change_types.append('MULTIPLE_PREVIOUS_GTDB_REPS')
            reasons.append(
                'Necessary name change as GTDB species cluster contains multiple previous representatives')

        if prev_specific in self.forbidden_specific_names:
            specific_change_types.append('INVALID_SPECIFIC_NAME')
            reasons.append(
                'Changed as previous assignment is not an effectively or validly published specific epithet')

        if len(specific_change_types) == 0:
            specific_change_types.append('OTHER')
            reasons.append(
                'Reason for change could not be automatically established')

        return '; '.join(specific_change_types), '; '.join(reasons)

    def report_modified_species_names(self,
                                      final_taxonomy,
                                      cur_genomes,
                                      prev_genomes,
                                      cur_clusters,
                                      prev_to_new_rids,
                                      unambiguous_ncbi_sp,
                                      unambiguous_ncbi_subsp,
                                      ncbi_misclassified_gids,
                                      gtdb_synonyms):
        """Report GTDB species clusters with modified names."""

        unambiguous_sp_rids = {}
        for rid, assignment_type in unambiguous_ncbi_sp.values():
            unambiguous_sp_rids[rid] = assignment_type

        unambiguous_subsp_rids = {}
        for ncbi_subsp, (rid, _assignment_type) in unambiguous_ncbi_subsp.items():
            unambiguous_subsp_rids[rid] = ncbi_subsp.split()[-1]

        domains_of_interest = set([taxa[Taxonomy.DOMAIN_INDEX]
                                   for taxa in final_taxonomy.values()])

        fout_sp = open(os.path.join(self.output_dir,
                                    'modified_specific_names.tsv'), 'w')
        fout_sp.write('Previous ID\tNew ID\tPrevious name\tNew name')
        fout_sp.write('\tUpdated generic name\tUpdated representative')
        fout_sp.write('\tEffective type strain')
        fout_sp.write(
            '\tRepresentative is NCBI representative\tCluster contain NCBI representative\tUnambiguous NCBI species')
        fout_sp.write('\tTransition\tCategory\tReason for change')
        fout_sp.write('\tNCBI species assignments\n')

        fout_generic = open(os.path.join(
            self.output_dir, 'modified_generic_names.tsv'), 'w')
        fout_generic.write('Previous ID\tNew ID\tPrevious name\tNew name')
        fout_generic.write('\tUpdated specific name\tUpdated representative')
        fout_generic.write('\tType species of genus')
        fout_generic.write('\tTransition\tCategory\tReason for change')
        fout_generic.write('\tNCBI species assignments\n')

        modified_species = set()
        modified_generic = set()
        modified_specific = set()

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
            elif new_rid in final_taxonomy:
                new_sp = final_taxonomy[new_rid][Taxonomy.SPECIES_INDEX]
                new_generic = generic_name(new_sp)
                new_specific = specific_epithet(new_sp)

                if prev_sp == new_sp:
                    same_species_names += 1
                else:
                    modified_species.add(new_rid)

                    # get NCBI species assignments for genomes in cluster
                    ncbi_sp_str = []
                    for sp, count in Counter([cur_genomes[cid].ncbi_taxa.species for cid in cur_clusters[new_rid]]).most_common():
                        ncbi_sp_str.append('{} ({:,})'.format(sp, count))

                    if prev_generic != new_generic:
                        modified_generic.add(new_rid)
                        modified_generic_transition[taxon_type(
                            'g__' + prev_generic)][taxon_type('g__' + new_generic)] += 1
                        generic_transition = '{} -> {}'.format(
                            taxon_type('g__' + prev_generic),
                            taxon_type('g__' + new_generic))

                        generic_change_category = 'TBD'
                        reason = 'TBD'

                        fout_generic.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            prev_rid, new_rid, prev_sp, new_sp,
                            prev_specific != new_specific, prev_rid != new_rid,
                            cur_genomes[new_rid].is_gtdb_type_species(),
                            generic_transition, generic_change_category, reason,
                            ', '.join(ncbi_sp_str)))

                    if prev_specific != new_specific:
                        modified_specific.add(new_rid)
                        modified_specific_transition[specific_epithet_type(
                            prev_specific)][specific_epithet_type(new_specific)] += 1
                        specific_transition = '{} -> {}'.format(
                            specific_epithet_type(prev_specific),
                            specific_epithet_type(new_specific))

                        specific_change_types, reason = self.justify_specific_name_change(
                            prev_rid,
                            new_rid,
                            prev_specific,
                            new_specific,
                            final_taxonomy,
                            cur_genomes,
                            prev_genomes,
                            cur_clusters,
                            unambiguous_sp_rids,
                            unambiguous_ncbi_sp,
                            unambiguous_subsp_rids,
                            ncbi_misclassified_gids,
                            gtdb_synonyms)

                        for change_type in [t.strip() for t in specific_change_types.split(';')]:
                            modified_specific_type[change_type] += 1

                        # check if GTDB representative is in cluster
                        ncbi_rep_in_cluster = False
                        for cid in cur_clusters[new_rid]:
                            if cur_genomes[cid].is_ncbi_representative():
                                ncbi_rep_in_cluster = True

                        fout_sp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            prev_rid, new_rid, prev_sp, new_sp,
                            prev_generic != new_generic, prev_rid != new_rid,
                            cur_genomes[new_rid].is_effective_type_strain(),
                            cur_genomes[new_rid].is_ncbi_representative(
                            ), ncbi_rep_in_cluster,
                            new_rid in unambiguous_sp_rids,
                            specific_transition, specific_change_types, reason,
                            ', '.join(ncbi_sp_str)))
            else:
                # genome must be from other domain
                pass

        self.log.info(
            ' - total previous species names: {:,}'.format(total_names))
        self.log.info(' - unmodified names: {:,} ({:.2f}%)'.format(
            same_species_names,
            same_species_names*100.0/total_names))
        self.log.info(' - modified species names: {:,} ({:.2f}%)'.format(
            len(modified_species),
            len(modified_species)*100.0/total_names))

        self.log.info(' - modified generic names: {:,} ({:.2f}%)'.format(
            len(modified_generic),
            len(modified_generic)*100.0/total_names))
        for prev_type, new_type in product(['LATIN', 'SUFFIXED_LATIN', 'ALPHANUMERIC'], repeat=2):
            self.log.info('   - modified generic names from {} to {}: {:,} ({:.2f}%)'.format(
                prev_type,
                new_type,
                modified_generic_transition[prev_type][new_type],
                modified_generic_transition[prev_type][new_type]*100.0/len(modified_generic)))

        self.log.info(' - modified specific names: {:,} ({:.2f}%)'.format(
            len(modified_specific),
            len(modified_specific)*100.0/total_names))
        for prev_type, new_type in product(['LATIN', 'SUFFIXED_LATIN', 'ALPHANUMERIC'], repeat=2):
            self.log.info('   - modified specific names from {} to {}: {:,} ({:.2f}%)'.format(
                prev_type,
                new_type,
                modified_specific_transition[prev_type][new_type],
                modified_specific_transition[prev_type][new_type]*100.0/len(modified_specific)))

        for specific_change_type, count in modified_specific_type.items():
            self.log.info(' - {}: {:,} ({:.2f}%)'.format(
                specific_change_type,
                count,
                count*100.0/len(modified_specific)))

        self.log.info(' - lost species clusters: {:,} ({:.2f}%)'.format(
            lost_sp_cluster,
            lost_sp_cluster*100.0/total_names))

        fout_sp.close()
        fout_generic.close()

        # generate curate trees to add manual inspection of modified species
        fout = open(os.path.join(self.output_dir,
                                 'modified_species.tree'), 'w')
        fout.write('({});\n'.format(','.join(modified_species)))
        fout.close()

        fout = open(os.path.join(self.output_dir,
                                 'modified_generic.tree'), 'w')
        fout.write('({});\n'.format(','.join(modified_generic)))
        fout.close()

        fout = open(os.path.join(self.output_dir,
                                 'modified_specific.tree'), 'w')
        fout.write('({});\n'.format(','.join(modified_specific)))
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

        self.log.info(
            'Comparison of number of name changes in NCBI and GTDB.')
        self.log.info(
            ' - total previous GTDB representatives with NCBI species assignment: {:,}'.format(total_names))
        self.log.info(' - modified NCBI species names: {:,} ({:.2f}%)'.format(
            ncbi_modified_sp_name,
            ncbi_modified_sp_name*100.0/total_names))
        self.log.info(' - modified GTDB species names: {:,} ({:.2f}%)'.format(
            gtdbd_modified_sp_name,
            gtdbd_modified_sp_name*100.0/total_names))
        self.log.info(' - modified NCBI genus names: {:,} ({:.2f}%)'.format(
            ncbi_modified_genus_names,
            ncbi_modified_genus_names*100.0/total_names))
        self.log.info(' - modified GTDB genus names: {:,} ({:.2f}%)'.format(
            gtdb_modified_genus_names,
            gtdb_modified_genus_names*100.0/total_names))
        self.log.info(' - modified NCBI specific names: {:,} ({:.2f}%)'.format(
            ncbi_modified_specific_names,
            ncbi_modified_specific_names*100.0/total_names))
        self.log.info(' - modified GTDB specific names: {:,} ({:.2f}%)'.format(
            gtdb_modified_specific_names,
            gtdb_modified_specific_names*100.0/total_names))

    def classify_ncbi_species(self,
                              ncbi_species_mngr,
                              ncbi_synonyms,
                              ncbi_misclassified_gids):
        """Classify each Latin NCBI species as either unambiguously or ambiguously assigned."""

        ncbi_classification_table = os.path.join(
            self.output_dir, 'ncbi_sp_classification.tsv')
        if os.path.exists(ncbi_classification_table):
            self.log.warning('Reading classification of NCBI species from existing table: {}'.format(
                ncbi_classification_table))
            ncbi_species, unambiguous_ncbi_sp, ambiguous_ncbi_sp = ncbi_species_mngr.parse_ncbi_classification_table(
                ncbi_classification_table)
        else:
            ncbi_species, unambiguous_ncbi_sp, ambiguous_ncbi_sp = ncbi_species_mngr.classify_ncbi_species(
                ncbi_synonyms,
                ncbi_misclassified_gids)

        unambiguous_classifications = [
            classification for rid, classification in unambiguous_ncbi_sp.values()]
        unambiguous_type_strain_count = unambiguous_classifications.count(
            NCBI_SpeciesManager.TYPE_STRAIN_GENOME)
        unambiguous_unanimous_consensus_count = unambiguous_classifications.count(
            NCBI_SpeciesManager.MAJORITY_VOTE)
        assert unambiguous_type_strain_count + \
            unambiguous_unanimous_consensus_count == len(unambiguous_ncbi_sp)

        self.log.info(
            ' - identified {:,} NCBI species.'.format(len(ncbi_species)))
        self.log.info(' - identified {:,} ({:.2f}%) as unambiguous assignments.'.format(
            len(unambiguous_ncbi_sp),
            len(unambiguous_ncbi_sp) * 100.0 / len(ncbi_species)))
        self.log.info('   - assigned {:,} ({:.2f}%) by type strain genome.'.format(
            unambiguous_type_strain_count,
            unambiguous_type_strain_count * 100.0 / len(ncbi_species)))
        self.log.info('   - assigned {:,} ({:.2f}%) by unanimous consensus.'.format(
            unambiguous_unanimous_consensus_count,
            unambiguous_unanimous_consensus_count * 100.0 / len(ncbi_species)))
        self.log.info(' - identified {:,} ({:.2f}%) as ambiguous assignments.'.format(
            len(ambiguous_ncbi_sp),
            len(ambiguous_ncbi_sp) * 100.0 / len(ncbi_species)))
        self.log.info(' - identified {:,} ({:.2f}%) as GTDB synonyms.'.format(
            len(ncbi_synonyms),
            len(ncbi_synonyms) * 100.0 / len(ncbi_species)))

        return ncbi_species, unambiguous_ncbi_sp, ambiguous_ncbi_sp

    def classify_ncbi_subspecies(self, ncbi_species_mngr):
        """Classify each NCBI subspecies as either unambiguously or ambiguously assigned."""

        (_ncbi_all_subspecies,
         unambiguous_ncbi_subsp,
         ambiguous_ncbi_subsp) = ncbi_species_mngr.classify_ncbi_subspecies()

        ncbi_subsp_classification_table = os.path.join(
            self.output_dir, 'ncbi_subsp_classification.tsv')
        if os.path.exists(ncbi_subsp_classification_table):
            self.log.warning('Reading classification of NCBI subspecies from existing table: {}'.format(
                ncbi_subsp_classification_table))
            ncbi_subspecies, unambiguous_ncbi_subsp, ambiguous_ncbi_subsp = ncbi_species_mngr.parse_ncbi_classification_table(
                ncbi_subsp_classification_table)
        else:
            ncbi_subspecies, unambiguous_ncbi_subsp, ambiguous_ncbi_subsp = ncbi_species_mngr.classify_ncbi_subspecies()

        unambiguous_classifications = [
            classification for rid, classification in unambiguous_ncbi_subsp.values()]
        unambiguous_type_strain_count = unambiguous_classifications.count(
            NCBI_SpeciesManager.TYPE_STRAIN_GENOME)
        unambiguous_unanimous_consensus_count = unambiguous_classifications.count(
            NCBI_SpeciesManager.MAJORITY_VOTE)
        assert unambiguous_type_strain_count + \
            unambiguous_unanimous_consensus_count == len(
                unambiguous_ncbi_subsp)

        self.log.info(
            ' - identified {:,} NCBI subspecies.'.format(len(ncbi_subspecies)))
        self.log.info(' - identified {:,} ({:.2f}%) as unambiguous assignments.'.format(
            len(unambiguous_ncbi_subsp),
            len(unambiguous_ncbi_subsp) * 100.0 / len(ncbi_subspecies)))
        self.log.info('   - assigned {:,} ({:.2f}%) by type strain genome.'.format(
            unambiguous_type_strain_count,
            unambiguous_type_strain_count * 100.0 / len(ncbi_subspecies)))
        self.log.info('   - assigned {:,} ({:.2f}%) by unanimous consensus.'.format(
            unambiguous_unanimous_consensus_count,
            unambiguous_unanimous_consensus_count * 100.0 / len(ncbi_subspecies)))
        self.log.info(' - identified {:,} ({:.2f}%) as ambiguous assignments.'.format(
            len(ambiguous_ncbi_subsp),
            len(ambiguous_ncbi_subsp) * 100.0 / len(ncbi_subspecies)))

        return ncbi_subspecies, unambiguous_ncbi_subsp, ambiguous_ncbi_subsp

    def run(self,
            final_taxonomy,
            final_scaled_tree,
            manual_sp_names,
            pmc_custom_species,
            gtdb_clusters_file,
            prev_gtdb_metadata_file,
            cur_gtdb_metadata_file,
            qc_passed_file,
            ncbi_misclassified_file,
            ncbi_genbank_assembly_file,
            untrustworthy_type_file,
            ncbi_synonym_file,
            updated_species_reps,
            gtdb_type_strains_ledger,
            species_classification_ledger,
            sp_priority_ledger,
            genus_priority_ledger,
            specific_epithet_ledger,
            ncbi_env_bioproject_ledger,
            lpsn_gss_metadata_file,
            lpsn_data,
            ground_truth_test_cases,
            skip_full_taxonomy_checks,
            skip_genus_checks):
        """Validate final species names."""

        # read LPSN metadata
        self.log.info('Reading LPSN metadata:')
        lpsn = LPSN(lpsn_data, lpsn_gss_metadata_file)
        self.log.info(' - identified {:,} correct LPSN species and {:,} synonyms.'.format(
            len(lpsn.sp_correct_names),
            len(lpsn.sp_synonyms)))
        self.log.info(' - identified {:,} taxa with defined type material at LPSN.'.format(
            len(lpsn.taxa)))

        # read manually-curated taxonomy
        self.log.info('Parsing final taxonomy:')
        final_taxonomy = read_taxonomy(final_taxonomy, use_canonical_gid=True)
        self.log.info(' - identified taxonomy strings for {:,} genomes.'.format(
            len(final_taxonomy)))

        # read species names explicitly set via manual curation
        self.mc_species = parse_manual_sp_curation_files(
            manual_sp_names, pmc_custom_species)

        # get priority of species and generic (genus) names
        sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger,
                                                  genus_priority_ledger,
                                                  lpsn_gss_metadata_file,
                                                  self.output_dir)

        # create previous and current GTDB genome sets
        self.log.info('Creating previous GTDB genome set.')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                             gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                             ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                             untrustworthy_type_ledger=untrustworthy_type_file,
                                             ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)
        self.log.info(' - previous genome set has {:,} species clusters spanning {:,} genomes.'.format(
            len(prev_genomes.sp_clusters),
            prev_genomes.sp_clusters.total_num_genomes()))

        self.log.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                            gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                            create_sp_clusters=False,
                                            qc_passed_file=qc_passed_file,
                                            ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                            untrustworthy_type_ledger=untrustworthy_type_file,
                                            ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)
        self.log.info(
            f' - current genome set contains {len(cur_genomes):,} genomes.')

        cur_genomes.set_prev_gtdb_classifications(prev_genomes)

        # sanity check
        self.log.info(
            'Checking that previous and current genomes have same GTDB assignments.')
        dup_genomes = 0
        for gid in set(cur_genomes).intersection(prev_genomes):
            if cur_genomes[gid].gtdb_taxa != prev_genomes[gid].gtdb_taxa:
                if cur_genomes[gid].gtdb_taxa.domain == prev_genomes[gid].gtdb_taxa.domain:
                    # there are known cases of genomes being transferred between domain which
                    # necessitates a different GTDB assignment
                    self.log.warning('Current GTDB assignments differ from previous assignment: {}, {}, {}'.format(
                        gid, cur_genomes[gid].gtdb_taxa, prev_genomes[gid].gtdb_taxa))
            dup_genomes += 1
        self.log.info(
            ' - validate GTDB assignments for {:,} genomes in common.'.format(dup_genomes))

        # read named GTDB species clusters
        self.log.info('Reading GTDB species clusters.')
        cur_clusters, _rep_radius = read_clusters(gtdb_clusters_file)
        self.log.info(' - identified {:,} clusters spanning {:,} genomes.'.format(
            len(cur_clusters),
            sum([len(gids) + 1 for gids in cur_clusters.values()])))
        assert len(set(final_taxonomy) - set(cur_clusters)) == 0

        # get mappings between new and previous GTDB representatives
        self.log.info(
            'Mapping current GTDB representatives to previous representatives.')
        updated_gtdb_rids = parse_updated_species_reps(updated_species_reps)
        new_to_prev_rid = infer_prev_gtdb_reps(
            prev_genomes, cur_clusters, updated_gtdb_rids)
        self.log.info(
            ' - mapped {:,} current representatives to previous representatives.'.format(len(new_to_prev_rid)))

        # get species clusters with reassigned representative
        self.log.info(
            'Identifying species clusters with reassigned representatives.')
        prev_to_new_rids = prev_genomes.sp_clusters.updated_representatives(
            cur_clusters)
        self.log.info(' - mapped {:,} previous representatives to current representatives.'.format(
            len(prev_to_new_rids)))

        # get polyphyletic suffixes for GTDB taxa
        self.taxon_suffix_manager = TaxonSuffixManager()

        # establish state of NCBI species
        ncbi_species_mngr = NCBI_SpeciesManager(cur_genomes,
                                                cur_clusters,
                                                self.mc_species,
                                                self.output_dir)
        ncbi_species_mngr.ncbi_sp_gtdb_cluster_table(final_taxonomy)
        self.forbidden_specific_names = ncbi_species_mngr.forbidden_specific_names

        # get list of synonyms in order to restrict usage of species names
        self.log.info('Reading NCBI synonyms.')
        ncbi_synonyms = ncbi_species_mngr.parse_synonyms_table(
            ncbi_synonym_file)
        self.log.info(' - identified {:,} synonyms from {:,} distinct NCBI species.'.format(
            len(ncbi_synonyms),
            len(set(ncbi_synonyms.values()))))

        # get synonyms as GTDB species
        self.log.info('Converting NCBI synonyms to GTDB species.')
        self.log.info('Reading NCBI synonyms.')
        gtdb_synonyms = ncbi_to_gtdb_synonyms(
            ncbi_synonym_file, final_taxonomy)
        self.log.info(' - identified {:,} synonyms from {:,} distinct GTDB species.'.format(
            len(gtdb_synonyms),
            len(set(gtdb_synonyms.values()))))

        # identified genomes with misclassified species assignments at NCBI
        self.log.info(
            'Identify genomes with misclassified NCBI species assignments.')
        ncbi_misclassified_gids = ncbi_species_mngr.parse_ncbi_misclassified_table(
            ncbi_misclassified_file)
        self.log.info(' - identified {:,} genomes with erroneous NCBI species assignments'.format(
            len(ncbi_misclassified_gids)))

        # classify each Latin NCBI species as either unambiguously or ambiguously assigned
        self.log.info(
            'Classifying NCBI species as unambiguous, ambiguous, or synonym.')
        _ncbi_species, unambiguous_ncbi_sp, _ambiguous_ncbi_sp = self.classify_ncbi_species(
            ncbi_species_mngr,
            ncbi_synonyms,
            ncbi_misclassified_gids)

        # classify each NCBI subspecies as either unambiguously or ambiguously assigned
        self.log.info(
            'Classifying NCBI subspecies as unambiguous or ambiguous.')
        _ncbi_all_subspecies, unambiguous_ncbi_subsp, _ambiguous_ncbi_subsp = self.classify_ncbi_subspecies(
            ncbi_species_mngr)

        # validate properly formatted species names
        if not skip_full_taxonomy_checks:
            self.log.info('Performing basic validation of taxonomy.')
            validate_taxonomy(final_taxonomy,
                              check_prefixes=True,
                              check_ranks=True,
                              check_hierarchy=True,
                              check_species=True,
                              check_group_names=True,
                              check_duplicate_names=True,
                              check_capitalization=True,
                              report_errors=True)

        # validate parent-child placeholder names (e.g. f__JDFR-13 conflicts with g__JdFR-13)
        # Note: It was decided that this capitalization inconsistency is acceptable in order
        # to avoid replacing long-standing placeholder names with slight modifications.
        #self.log.info('Validating parent-child placeholder names for simple typos.')
        # self.validate_parent_child_placeholder_names(final_taxonomy)

        # validate spelling of GTDB names by comparing to NCBI names
        self.log.info('Validating spelling of GTDB names by comparison to NCBI.')
        self.validate_gtdb_spelling(final_taxonomy, cur_genomes)

        # validate placeholder names that might be mistaken as Latin names (e.g., g__Gsub)
        self.log.info(
            'Validating placeholder names that might be mistaken as Latin names.')
        self.validate_placeholder_names_not_latin(final_taxonomy, lpsn)

        # validate alphanumeric generic names do not have suffixed Latin specific names
        # (e.g., s__DTPE01 lauensis_A is considered invalid)
        self.log.info(
            'Validating species with alphanumeric generic names do not have suffixed Latin specific names.')
        self.validate_sp_with_alphanumeric_generic_names(final_taxonomy)

        # validate suffix of specific names by looking for small deviations relative to LPSN names
        self.log.info(
            "Validating suffix of specific names by identifying small deviations from LPSN names.")
        self.validate_suffix_of_specific_names(
            final_taxonomy, cur_genomes, lpsn)

        # validate application of names for merged taxa
        self.log.info("Validating names of merged taxa.")
        self.validate_merged_taxa(final_taxonomy, cur_genomes, lpsn)

        # validate placement of higher taxon names derived from stem of GTDB genera
        self.log.info(
            "Validating placement of higher taxon names based on LPSN type material")
        self.validate_taxa_by_lpsn_type_material(final_taxonomy, lpsn)

        # validate that GTDB is using the 'correct' names as specified at LPSN
        self.log.info(
            "Validating that GTDB species names are using the 'correct' name according to LPSN.")
        self.validate_lpsn_correct_names(final_taxonomy, lpsn)

        # validate that transferred species do not use a specific epithet which
        # has already been validated for the species
        self.log.info(
            'Validating that specific epithet of GTDB species names are valid after genus transfers.')
        self.validate_genus_transfer_specific_epithet(final_taxonomy,
                                                      cur_genomes,
                                                      lpsn)

        # validate pplacement of higher taxon names derived from stem of GTDB genera
        # (test is very slow)
        #***self.log.info(
        #***    "Validating placement of higher taxon names derived from stem of GTDB genera.")
        #***self.validate_taxa_by_genus_stem(final_taxonomy)

        # validate that all species names are unique
        self.log.info('Validating that species names are unique.')
        self.validate_unique_species_names(final_taxonomy, cur_genomes)

        # validate manually-curated species names
        self.log.info('Validating manually-curated species names.')
        self.validate_manually_curated_species(final_taxonomy, self.mc_species)

        # validate genomes in GTDB type strain ledger
        self.log.info('Validating genomes in type strain ledger.')
        mc_type_strain = self.validate_type_strain_ledger(
            final_taxonomy, gtdb_type_strains_ledger)
        self.mc_species.update(mc_type_strain)

        # validate genomes in GTDB type strain ledger
        self.log.info(
            'Validating genomes in species classification ledger.')
        mc_sp_classification = self.validate_species_classification_ledger(final_taxonomy,
                                                                           cur_clusters,
                                                                           species_classification_ledger)
        self.mc_species.update(mc_sp_classification)

        # validate genus and generic names are the same
        self.log.info(
            'Validating that genus and generic names are the same.')
        self.validate_genus_generic_match(final_taxonomy)

        # validate type species
        self.log.info('Validating type species of genus.')
        self.validate_type_species(
            final_taxonomy, cur_genomes, sp_priority_mngr)

        # validate type strain
        self.log.info('Validating type strains of species.')
        self.validate_type_strain(final_taxonomy, cur_genomes)

        # validate ground truth test cases
        self.log.info('Validating ground truth test cases.')
        self.validate_ground_truth(final_taxonomy, ground_truth_test_cases)

        # validate manually established suffixes for specific names changed due to genus transfer
        self.log.info(
            'Validating manually-curated specific name suffix changes resulting from genus transfers.')
        self.validate_mannually_curated_specifix_suffix(
            final_taxonomy, cur_genomes, specific_epithet_ledger)

        # validate NCBI species considered to have an unambiguous assignment within the GTDB
        self.log.info(
            'Validating NCBI species names classified as having an unambiguous assignment within the GTDB.')
        self.validate_unambiguous_ncbi_sp(final_taxonomy,
                                          unambiguous_ncbi_sp,
                                          lpsn)

        # validate NCBI species not considered to have an unambiguous assignment within the GTDB
        self.log.info(
            'Validating NCBI species names classified as having an ambiguous assignment within the GTDB.')
        self.validate_ambiguous_ncbi_sp(final_taxonomy,
                                        cur_genomes,
                                        unambiguous_ncbi_sp,
                                        unambiguous_ncbi_subsp)

        # validate species names of genomes misclassified at NCBI
        self.log.info(
            'Validating species names of genomes misclassified at NCBI.')
        self.validate_misclassified_genomes(final_taxonomy,
                                            cur_genomes,
                                            ncbi_synonyms,
                                            ncbi_misclassified_gids,
                                            unambiguous_ncbi_sp,
                                            unambiguous_ncbi_subsp)

        # validate synonyms
        self.log.info('Validating that GTDB synonyms are respected.')
        self.validate_synonyms(final_taxonomy, cur_genomes, ncbi_synonyms)

        # validate that each GTDB genus has a single type species
        self.log.info(
            'Validating that each GTDB genus has a single type species cluster.')
        self.validate_single_type_species(final_taxonomy, cur_genomes)

        # validate GTDB specific name assignments without a polyphyletic suffixes
        self.log.info(
            'Validating specific name assignments without a placeholder suffix.')
        self.validate_specific_no_suffix(final_taxonomy,
                                         cur_genomes,
                                         unambiguous_ncbi_sp,
                                         unambiguous_ncbi_subsp)

        # validate that species clusters use same specific name placeholder suffixes as previously assigned
        self.log.info(
            'Validating specific names use same placeholder suffixes as assigned in previous release.')
        self.validate_prev_specific_placeholder_suffixes(final_taxonomy,
                                                         cur_genomes,
                                                         new_to_prev_rid,
                                                         unambiguous_ncbi_sp)

        # validate suffix assignment rules for specific names
        self.log.info(
            'Validating specific names follow suffix assignment rules.')
        self.validate_specific_placeholder_suffixes(final_taxonomy,
                                                    cur_genomes,
                                                    unambiguous_ncbi_sp,
                                                    unambiguous_ncbi_subsp)

        # validate spelling of genera
        self.log.info('Validating spelling of genera.')
        self.validate_genera_spelling(final_taxonomy, cur_genomes)

        # validate that no forbidden specific names are used
        self.log.info(
            'Validating that no forbidden specific names were used.')
        self.validate_forbidden_specific_names(final_taxonomy)

        # validate that Latinized specific names are never used with non-suffixed, placeholder generic names
        self.log.info(
            'Validating that Latinized specific names are not used with non-suffixed, placeholder generic names.')
        self.validate_latin_sp_placeholder_generic(final_taxonomy, cur_genomes)

        if not skip_genus_checks:
            # validate placeholder suffixes for genera
            self.log.info('Validating placeholder suffixes of genera.')
            self.validate_genera_placeholder_suffixes(
                final_taxonomy, cur_genomes, new_to_prev_rid)

            # read scaled curation tree
            if final_scaled_tree.lower() != 'none':
                self.log.info('Reading scaled and decorated tree.')
                curation_tree = dendropy.Tree.get_from_path(final_scaled_tree,
                                                            schema='newick',
                                                            rooting='force-rooted',
                                                            preserve_underscores=True)
                curation_tree.calc_node_root_distances()

                self.log.info(
                    'Validating missing NCBI-defined genus names.')
                self.validate_missing_ncbi_genera(
                    final_taxonomy, cur_genomes, curation_tree)

            # validate placeholder genus names lack potential Latin names
            self.log.info(
                'Validating placeholder genus names lack potential Latin names.')
            self.validate_no_latin_genus_name(final_taxonomy, cur_genomes)

            # validate that genus names are not a replacement of one placeholder for another
            self.log.info(
                'Validating that genus names are not a change from one placeholder to another.')
            self.validate_genus_placeholder_change(
                final_taxonomy, cur_genomes, prev_genomes)

            # validate monophyly of NCBI genera
            self.log.info('Validating monophyly of NCBI genera.')
            self.validate_ncbi_genera_monophyly(final_taxonomy,
                                                cur_genomes,
                                                cur_clusters)

        # create table with modified GTDB species names for manual inspection
        self.log.info('Creating table with modified GTDB species names.')
        self.report_modified_species_names(final_taxonomy,
                                           cur_genomes,
                                           prev_genomes,
                                           cur_clusters,
                                           prev_to_new_rids,
                                           unambiguous_ncbi_sp,
                                           unambiguous_ncbi_subsp,
                                           ncbi_misclassified_gids,
                                           gtdb_synonyms)
