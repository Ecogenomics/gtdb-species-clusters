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
import logging
from collections import defaultdict

import dendropy

from biolib.taxonomy import Taxonomy
from biolib.newick import parse_label, create_label

from gtdb_species_clusters.genomes import Genomes

from gtdb_species_clusters.species_name_manager import SpeciesNameManager
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters.specific_epithet_manager import SpecificEpithetManager
from gtdb_species_clusters.ncbi_species_manager import NCBI_SpeciesManager
from gtdb_species_clusters.genome_utils import canonical_gid
from gtdb_species_clusters.type_genome_utils import (read_clusters,
                                                     parse_updated_species_reps,
                                                     infer_prev_gtdb_reps,
                                                     parse_manual_sp_curation_files)
from gtdb_species_clusters.taxon_utils import (generic_name,
                                               specific_epithet,
                                               canonical_taxon,
                                               taxon_suffix,
                                               sort_by_naming_priority,
                                               is_latin_sp_epithet,
                                               is_placeholder_sp_epithet,
                                               test_same_epithet,
                                               longest_common_suffix)


class PMC_SpeciesNames(object):
    """Establish final species names based on manual curation."""

    def __init__(self, output_dir):
        """Initialization."""

        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')

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

        if is_placeholder_sp_epithet(specific):
            specific = self.sp_name_mngr.numeric_placeholder_sp_epithet(gid)
        else:
            specific = self.sp_name_mngr.suffixed_placeholder_sp_epithet(
                generic, specific)

        return 's__{} {}'.format(generic, specific)

    def create_numeric_sp_placeholder(self, gid, generic):
        """Explicitly create a numeric sp<accn> placeholder name."""

        specific = self.sp_name_mngr.numeric_placeholder_sp_epithet(gid)

        return 's__{} {}'.format(generic, specific)

    def set_type_species(self,
                         rid,
                         final_taxonomy,
                         cur_genomes,
                         unambiguous_ncbi_sp,
                         case_count):
        """Establish final species name for type species of genus genome."""

        (gtdb_genus,
         gtdb_species,
         _gtdb_generic,
         gtdb_specific) = self.key_taxon(rid, final_taxonomy)

        ncbi_genus = cur_genomes[rid].ncbi_taxa.genus
        ncbi_sp = cur_genomes[rid].ncbi_taxa.species

        if ncbi_sp not in unambiguous_ncbi_sp or rid != unambiguous_ncbi_sp[ncbi_sp][0]:
            self.logger.error('NCBI species {} represented by type species genome {} not designated as an unambiguous NCBI species.'.format(
                ncbi_sp,
                rid))
            sys.exit(-1)

        note = ''
        if gtdb_species != ncbi_sp:
            gtdb_priority = self.sp_priority_mngr.genus_priority(
                gtdb_genus, ncbi_genus)

            if gtdb_genus != ncbi_genus and (gtdb_priority is None or gtdb_priority == ncbi_genus):
                case_count['INCONGRUENT_TYPE_SPECIES_GENERIC'] += 1
                note = 'incongruent type species'
                print('INCONGRUENT_TYPE_SPECIES_GENERIC',
                      rid, gtdb_genus, ncbi_genus)

            if not test_same_epithet(gtdb_specific, specific_epithet(ncbi_sp)):
                case_count['INCONGRUENT_TYPE_SPECIES_SPECIFIC'] += 1
                note = 'incongruent type species'
                print('INCONGRUENT_TYPE_SPECIES_SPECIFIC',
                      rid, gtdb_species, ncbi_sp)

        return gtdb_species, note

    def set_type_strain(self, rid, final_taxonomy, cur_genomes, unambiguous_ncbi_sp, case_count):
        """Establish final species name for type strain of species genome."""

        (_gtdb_genus,
         gtdb_species,
         gtdb_generic,
         gtdb_specific) = self.key_taxon(rid, final_taxonomy)

        ncbi_sp = cur_genomes[rid].ncbi_taxa.species
        ncbi_specific = specific_epithet(ncbi_sp)

        if ncbi_sp not in unambiguous_ncbi_sp or rid != unambiguous_ncbi_sp[ncbi_sp][0]:
            self.logger.warning('NCBI species {} represented by type strain genome {} not designated as an unambiguous NCBI species.'.format(
                                ncbi_sp,
                                rid))

        note = ''
        final_sp = gtdb_species
        if ncbi_sp == 's__':
            # Genome lacks species assignment at NCBI so it being a type strain must
            # have been manually asserted
            pass
        elif not test_same_epithet(canonical_taxon(gtdb_specific), ncbi_specific):
            # specific epithets disagree which should never occur for type strains,
            # so changing this to the NCBI proposed specific epithet
            case_count['INCONGRUENT_TYPE_STRAIN'] += 1
            note = 'type strain genome had different GTDB and NCBI specific epithets'
            final_sp = 's__{} {}'.format(gtdb_generic, ncbi_specific)
            print('INCONGRUENT_TYPE_STRAIN', rid,
                  gtdb_species, ncbi_sp, final_sp)
        elif is_placeholder_sp_epithet(gtdb_specific):
            # Latin name with a suffix or sp<accn> specific name is by definition an
            # error since this is a type strain genome
            # [DHP (Mar. 22, 2021): actually this isn't always true since we have the
            # case of a type strain genome being transferred to a new genus where
            # the specific name needs a suffix, e.g.: G002240355 is the type strain
            # of Prauserella marina, but was transferred to the genus Saccharomonospora
            # and S. marina is already a recognized name so this genome should have
            # the name S. marina_A].
            case_count['ERRONEOUS_SUFFIX_TYPE_STRAIN'] += 1
            note = 'type strain genome had erroneous placeholder suffix'
            final_sp = 's__{} {}'.format(gtdb_generic, ncbi_specific)
            print('ERRONEOUS_SUFFIX_TYPE_STRAIN', rid,
                  gtdb_species, ncbi_sp, final_sp)

        return final_sp, note

    def generate_suffixed_final_name(self, gtdb_generic, prev_gtdb_specific, ncbi_specific):
        """Generate Latin name with suffix, using previously defined suffix where possible."""

        if (prev_gtdb_specific
                and canonical_taxon(prev_gtdb_specific) == ncbi_specific
                and taxon_suffix(prev_gtdb_specific)):
            # must recycle previously used suffix for this GTDB species cluster
            note = 'Using previously assigned suffix'
            prev_suffix = taxon_suffix(prev_gtdb_specific)
            final_species = 's__{} {}_{}'.format(
                gtdb_generic, ncbi_specific, prev_suffix)
        else:
            # determine next unused suffix associated with this GTDB species
            note = 'Generated new suffix'
            suffixed_specific_name = self.sp_name_mngr.suffixed_placeholder_sp_epithet(
                gtdb_generic, ncbi_specific)
            final_species = 's__{} {}'.format(
                gtdb_generic, suffixed_specific_name)

        return final_species, note

    def set_nontype_cluster_conservative(self,
                                         rid,
                                         final_taxonomy,
                                         cur_genomes,
                                         prev_genomes,
                                         new_to_prev_rid,
                                         ncbi_misclassified_gids,
                                         unambiguous_sp_rids,
                                         unambiguous_subsp_rids,
                                         ncbi_species_type_strain_rid,
                                         ncbi_synonyms,
                                         cur_gtdb_sp):
        """Establish final species name for non-type material genome.

        GTDB species clusters are assigned names that modify previously assigned
        Latin names as required for correctness (e.g., suffix changes due to genus transfers,
        fixing of previous errors, assigning Latin-suffix names to clarify valid genus transfers
        from erroneous NCBI species assignments, etc.). This method emphasizes names being static
        between releases and a placeholder name (Latin-suffixed or alphanumeric) is never replaced
        with another placeholder name. Consequently, it is entirely possible the placeholder name may
        no longer reflect the NCBI species assignments of genomes in the cluster. This is deemed an
        acceptable trade-off for the increased stability of names.
        """

        # get taxon names
        (_gtdb_genus,
         gtdb_species,
         gtdb_generic,
         _gtdb_specific) = self.key_taxon(rid, final_taxonomy)

        prev_gtdb_species = None
        prev_gtdb_specific = None
        if rid in new_to_prev_rid:
            prev_rid = new_to_prev_rid[rid]
            prev_gtdb_species = prev_genomes[prev_rid].gtdb_taxa.species
            prev_gtdb_specific = specific_epithet(prev_gtdb_species)

        final_species = gtdb_species
        note = ''
        if rid in unambiguous_sp_rids:
            # all genomes with this NCBI species assignment are contained in this
            # GTDB species cluster and all genomes in the cluster have this species
            # assignment (after discounting identified NCBI misclassifications),
            # so the cluster should be assigned the NCBI specific epithet without a suffix
            note = 'GTDB cluster has consensus for NCBI species'

            ncbi_species = unambiguous_sp_rids[rid]
            ncbi_species = ncbi_synonyms.get(ncbi_species, ncbi_species)
            ncbi_specific = specific_epithet(ncbi_species)

            final_species = 's__{} {}'.format(gtdb_generic, ncbi_specific)
            if final_species in cur_gtdb_sp:
                # name already exists so this collision needs to be resolved
                prev_sp_rid = cur_gtdb_sp[final_species][0]
                ncbi_generic = generic_name(ncbi_species)
                if (cur_genomes[prev_sp_rid].is_effective_type_strain()
                        and ncbi_generic != gtdb_generic):
                    # previous assignment should unambiguously have this species name,
                    # and this genome must have a conflict since it was transfer into
                    # this genus and happens to have the same specific name
                    # (e.g. Natronohydrobacter thiooxidans transferred into Roseinatronobacter,
                    #  while already contains the species R. thiooxidans)
                    final_species, note = self.generate_suffixed_final_name(
                        gtdb_generic, prev_gtdb_specific, ncbi_specific)
                    note = 'GTDB cluster has consensus NCBI species; Species transferred to genus with a species containing the same specific name; {}'.format(
                        note)
                else:
                    self.logger.error('Manual curation require to resolve name of genome {} as proposed name {} was assigned previously.'.format(
                        rid, final_species))

        elif rid in unambiguous_subsp_rids:
            # subspecies name can be promoted to the specific epithet of the assigned GTDB species name
            note = 'Subspecies promoted to specific epithet'
            ncbi_subsp = unambiguous_subsp_rids[rid].split()[-1]
            final_species = 's__{} {}'.format(gtdb_generic, ncbi_subsp)
        else:
            # species cluster must be assigned a placeholder name, which will be derived from
            # the selected GTDB representative
            if prev_gtdb_specific and is_placeholder_sp_epithet(prev_gtdb_specific):
                # re-use existing placeholder name whenever it exists
                note = 'Re-using previous placeholder name.'
                final_species = 's__{} {}'.format(
                    gtdb_generic, prev_gtdb_specific)
                self.sp_name_mngr.add_suffixed_species(final_species)
            else:
                # a new placeholder name is required
                ncbi_species = cur_genomes[rid].ncbi_taxa.species
                ncbi_specific = specific_epithet(ncbi_species)

                # if genome is considered to have an erroneous NCBI
                # species name, determine if it is at least in the
                # same GTDB family
                misclassified_diff_family = False
                if rid in ncbi_misclassified_gids:
                    gtdb_family = final_taxonomy[rid][Taxonomy.FAMILY_INDEX]
                    type_rid = ncbi_species_type_strain_rid[ncbi_species]
                    type_strain_gtdb_family = final_taxonomy[type_rid][Taxonomy.FAMILY_INDEX]
                    if gtdb_family != type_strain_gtdb_family:
                        misclassified_diff_family = True

                if (ncbi_specific
                        and ncbi_specific not in self.forbidden_specific_names
                        and not misclassified_diff_family):
                    # base placeholder off NCBI species assignment of representative
                    note = 'Generated Latin-suffix name from NCBI species assignment'
                    suffixed_specific_name = self.sp_name_mngr.suffixed_placeholder_sp_epithet(
                        gtdb_generic, ncbi_specific)
                    final_species = 's__{} {}'.format(
                        gtdb_generic, suffixed_specific_name)
                else:
                    # representative lacks a NCBI species assignment so create a numeric name
                    note = 'Generated alphanumeric name as representative lacks an NCBI species assignment or NCBI family considered misclassified'
                    final_species = self.create_numeric_sp_placeholder(
                        rid, gtdb_generic)

        return final_species, note

    def resolve_specific_epithet_suffixes(self, final_taxonomy):
        """Resolve cases where specific epithet needs to be modified to account for genus transfer."""

        for rid in final_taxonomy:
            (_gtdb_genus,
             gtdb_species,
             _gtdb_generic,
             gtdb_specific) = self.key_taxon(rid, final_taxonomy)

            canonical_gtdb_species = canonical_taxon(gtdb_species)
            gtdb_species_corr = self.sp_epithet_mngr.translate_species(
                canonical_gtdb_species)
            gtdb_specific_corr = specific_epithet(gtdb_species_corr)

            if canonical_taxon(gtdb_specific) != canonical_taxon(gtdb_specific_corr):
                suffix = taxon_suffix(gtdb_specific)
                if suffix is None:
                    final_taxonomy[rid][Taxonomy.SPECIES_INDEX] = gtdb_species_corr
                else:
                    final_taxonomy[rid][Taxonomy.SPECIES_INDEX] = '{}_{}'.format(
                        gtdb_species_corr, suffix)

    def resolve_species_classification_ledger(self, final_taxonomy, cur_genomes, cur_clusters, species_classification_ledger, manual_curation_rids):
        """Validate genomes in GTDB type strain ledger."""

        for gid in final_taxonomy:
            target_domain = final_taxonomy[gid][Taxonomy.DOMAIN_INDEX]
            break

        # read genomes in species classification ledger
        gtdb_sp_ledger = {}
        with open(species_classification_ledger, encoding='utf-8') as f:
            f.readline()

            for line in f:
                tokens = line.strip().split('\t')
                gid = canonical_gid(tokens[0])

                domain = tokens[2].strip()
                if not domain.startswith('d__'):
                    domain = 'd__' + domain

                if domain != target_domain:
                    continue

                sp = tokens[1].strip()
                if not sp.startswith('s__'):
                    sp = 's__' + sp

                gtdb_sp_ledger[gid] = sp

        if len(gtdb_sp_ledger) == 0:
            self.logger.error(
                'Failed to identify any entries in specifies classification ledger.')
            sys.exit(-1)

        # get map from genomes to their representative
        gtdb_gid_to_rid = {}
        for rid, cids in cur_clusters.items():
            for cid in cids:
                gtdb_gid_to_rid[cid] = rid

        # get all species defined in final taxonomy
        final_species = {}
        for rid, taxa in final_taxonomy.items():
            final_species[rid] = taxa[Taxonomy.SPECIES_INDEX]

        # defer to assignments in species classification ledger unless they
        # actively conflict with other assignments, in which case manual
        # curation is required
        fout = open(os.path.join(self.output_dir,
                                 'species_classification_ledger_updates.tsv'), 'w')
        fout.write('Genome ID in ledger\tRepresentive ID of cluster')
        fout.write(
            '\tLedger species\tGTDB species name after applying ledger\tSame as ledger\tGTDB species name before applying ledger')
        fout.write(
            '\tNCBI species of ledger genome\tNCBI species of representative\n')
        for gid, mc_species in gtdb_sp_ledger.items():
            if gid not in gtdb_gid_to_rid:
                self.logger.warning(
                    'Species classification genome {} no longer present in GTDB.'.format(gid))
                continue

            rid = gtdb_gid_to_rid[gid]

            gtdb_species = final_taxonomy[rid][Taxonomy.SPECIES_INDEX]
            gtdb_generic = generic_name(gtdb_species)
            gtdb_specific = specific_epithet(gtdb_species)

            mc_specific = specific_epithet(mc_species)

            final_sp = 's__{} {}'.format(gtdb_generic, mc_specific)
            if final_sp in final_species and rid != final_species[final_sp]:
                self.logger.error('Species {} proposed for {} in species classification ledger already exists in GTDB taxonomy: {}'.format(
                    final_sp, gid, final_species[final_sp]))

            if (rid in manual_curation_rids
                    and final_taxonomy[rid][Taxonomy.SPECIES_INDEX] != final_sp):
                self.logger.warning('Deferring to manual curation over update suggested by species classification ledger {}: {} {}'.format(
                    rid, final_taxonomy[rid][Taxonomy.SPECIES_INDEX], final_sp))
            else:
                final_taxonomy[rid][Taxonomy.SPECIES_INDEX] = final_sp

            if is_latin_sp_epithet(gtdb_specific) and not test_same_epithet(gtdb_specific, mc_specific):
                if gid in final_taxonomy:
                    self.logger.warning('Species classification ledger resulted in active change to Latin name of representative {}: {} -> {}'.format(
                        rid, gtdb_species, mc_species))
                else:
                    self.logger.warning('Species classification ledger entry for {} resulted in active change to Latin name of cluster now represented by {}: {} -> {}'.format(
                        gid, rid, gtdb_species, mc_species))

            fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                gid,
                rid,
                mc_species,
                final_taxonomy[rid][Taxonomy.SPECIES_INDEX],
                mc_species == final_taxonomy[rid][Taxonomy.SPECIES_INDEX],
                gtdb_species,
                cur_genomes[gid].ncbi_taxa.species,
                cur_genomes[rid].ncbi_taxa.species))

        fout.close()

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
                            new_sp = self.create_placeholder_species_name(
                                rid, generic_name(sp), specific_epithet(sp))

                            prev_gtdb_sp = cur_genomes[rid].gtdb_taxa.species
                            if (taxon_suffix(specific_epithet(prev_gtdb_sp))
                                    and canonical_taxon(prev_gtdb_sp) == canonical_taxon(new_sp)):
                                # should retain previously assign alphabetic suffixed
                                # (e.g., previously assigned s__Saccharomonospora marina_A, so should not be reassign as s__Saccharomonospora marina_B)
                                # (e.g., previously assigned as s__Dorea phocaeensis, so should be reassigned as s__Dorea phocaeensis_A)
                                new_sp = cur_genomes[rid].gtdb_taxa.species
                            mc_taxonomy[rid][Taxonomy.SPECIES_INDEX] = new_sp

                            print('Modified {} from {} to {} as multiple type strain species had the proposed name and this genome was transferred from a different NCBI genera, {}.'.format(
                                rid,
                                sp,
                                mc_taxonomy[rid][Taxonomy.SPECIES_INDEX],
                                cur_genomes[rid].ncbi_taxa.species))

        return resolved_gids

    def finalize_species_name(self,
                              rid,
                              final_sp,
                              case,
                              note,
                              final_taxonomy,
                              cur_gtdb_sp):
        "Set final species name for GTDB representative and report potential issues."

        if final_sp in cur_gtdb_sp:
            prev_rid, prev_case = cur_gtdb_sp[final_sp]
            self.logger.warning('Finalized GTDB species name {} already exists: rid={} prev_rid={} prev_case={}'.format(
                                final_sp, rid, prev_rid, prev_case))

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
                               new_to_prev_rid,
                               prev_genomes,
                               cur_genomes,
                               ncbi_misclassified_gids,
                               unambiguous_ncbi_sp,
                               unambiguous_ncbi_subsp,
                               ncbi_synonyms,
                               species_classification_ledger):
        """Establish final species names."""

        self.final_name_log = open(os.path.join(
            self.output_dir, 'final_sp_name.log'), 'w')
        self.final_name_log.write('Genome ID\tGTDB species\tCase\tNote\n')

        final_taxonomy = copy.deepcopy(mc_taxonomy)

        # order representatives by naming priority
        rtn = sort_by_naming_priority(final_taxonomy,
                                      cur_clusters,
                                      prev_genomes,
                                      cur_genomes,
                                      mc_species)
        manual_curation_rids, type_species_rids, type_strains_rids, binomial_rids, placeholder_rids = rtn

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
                                                   unambiguous_ncbi_sp,
                                                   case_count)

            self.finalize_species_name(rid,
                                       final_sp,
                                       case,
                                       note,
                                       final_taxonomy,
                                       cur_gtdb_sp)

        for rid in type_strains_rids:
            case = 'TYPE_STRAIN_OF_SPECIES'
            final_sp, note = self.set_type_strain(rid,
                                                  final_taxonomy,
                                                  cur_genomes,
                                                  unambiguous_ncbi_sp,
                                                  case_count)

            self.finalize_species_name(rid,
                                       final_sp,
                                       case,
                                       note,
                                       final_taxonomy,
                                       cur_gtdb_sp)

        # resolve type strain genomes with the same name
        self.resolve_duplicate_type_strain_species(final_taxonomy, cur_genomes)

        # establish final names for non-type genomes
        unambiguous_sp_rids = {}
        for ncbi_species, (rid, _assignment_type) in unambiguous_ncbi_sp.items():
            unambiguous_sp_rids[rid] = ncbi_species

        unambiguous_subsp_rids = {}
        for ncbi_subsp, (rid, _assignment_type) in unambiguous_ncbi_subsp.items():
            unambiguous_subsp_rids[rid] = ncbi_subsp

        ncbi_species_type_strain_rid = {}
        for rid, cids in cur_clusters.items():
            for cid in cids:
                if cur_genomes[cid].is_effective_type_strain():
                    ncbi_sp = cur_genomes[cid].ncbi_taxa.species
                    ncbi_species_type_strain_rid[ncbi_sp] = rid

        for rid in binomial_rids + placeholder_rids:
            case = 'NONTYPE_SPECIES_CLUSTER'
            final_sp, note = self.set_nontype_cluster_conservative(rid,
                                                                   final_taxonomy,
                                                                   cur_genomes,
                                                                   prev_genomes,
                                                                   new_to_prev_rid,
                                                                   ncbi_misclassified_gids,
                                                                   unambiguous_sp_rids,
                                                                   unambiguous_subsp_rids,
                                                                   ncbi_species_type_strain_rid,
                                                                   ncbi_synonyms,
                                                                   cur_gtdb_sp)

            self.finalize_species_name(rid,
                                       final_sp,
                                       case,
                                       note,
                                       final_taxonomy,
                                       cur_gtdb_sp)

        self.final_name_log.close()

        for case, count in case_count.items():
            print('{}\t{}'.format(case, count))

        # resolve specific epithets requiring changes due to genus transfers
        self.logger.info(
            'Resolving changes to suffix of specific epithets due to genus transfers.')
        self.resolve_specific_epithet_suffixes(final_taxonomy)

        # Incorporate species assignments in species classification ledger. This
        # is done last since conflicts with the ledger need to be identified and
        # manually validated.
        self.logger.info(
            'Resolving species assignments indicated in species classification ledger.')
        self.resolve_species_classification_ledger(final_taxonomy,
                                                   cur_genomes,
                                                   cur_clusters,
                                                   species_classification_ledger,
                                                   manual_curation_rids)

        return final_taxonomy

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
                        gid,
                        gtdb_sp_name,
                        ncbi_sp))

        return gtdb_type_strain_ledger

    def classify_ncbi_species(self,
                              ncbi_species_mngr,
                              ncbi_synonyms,
                              ncbi_misclassified_gids):
        """Classify each Latin NCBI species as either unambiguously or ambiguously assigned."""

        ncbi_classification_table = os.path.join(
            self.output_dir, 'ncbi_sp_classification.tsv')
        if os.path.exists(ncbi_classification_table):
            self.logger.warning('Reading classification of NCBI species from existing table: {}'.format(
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
        unambiguous_majority_vote_count = unambiguous_classifications.count(
            NCBI_SpeciesManager.MAJORITY_VOTE)
        assert unambiguous_type_strain_count + \
            unambiguous_majority_vote_count == len(unambiguous_ncbi_sp)

        self.logger.info(
            ' - identified {:,} NCBI species.'.format(len(ncbi_species)))
        self.logger.info(' - resolved {:,} ({:.2f}%) of species names.'.format(
            len(unambiguous_ncbi_sp),
            len(unambiguous_ncbi_sp) * 100.0 / len(ncbi_species)))
        self.logger.info('   - assigned {:,} ({:.2f}%) by type strain genome.'.format(
            unambiguous_type_strain_count,
            unambiguous_type_strain_count * 100.0 / len(ncbi_species)))
        self.logger.info('   - assigned {:,} ({:.2f}%) by majority vote.'.format(
            unambiguous_majority_vote_count,
            unambiguous_majority_vote_count * 100.0 / len(ncbi_species)))
        self.logger.info(' - identified {:,} ({:.2f}%) as ambiguous assignments.'.format(
            len(ambiguous_ncbi_sp),
            len(ambiguous_ncbi_sp) * 100.0 / len(ncbi_species)))
        self.logger.info(' - identified {:,} ({:.2f}%) as GTDB synonyms.'.format(
            len(ncbi_synonyms),
            len(ncbi_synonyms) * 100.0 / len(ncbi_species)))

        return ncbi_species, unambiguous_ncbi_sp, ambiguous_ncbi_sp

    def classify_ncbi_subspecies(self, ncbi_species_mngr):
        """Classify each NCBI subspecies as either unambiguously or ambiguously assigned."""

        ncbi_subsp_classification_table = os.path.join(
            self.output_dir, 'ncbi_subsp_classification.tsv')
        if os.path.exists(ncbi_subsp_classification_table):
            self.logger.warning('Reading classification of NCBI subspecies from existing table: {}'.format(
                ncbi_subsp_classification_table))
            ncbi_subspecies, unambiguous_ncbi_subsp, ambiguous_ncbi_subsp = ncbi_species_mngr.parse_ncbi_classification_table(
                ncbi_subsp_classification_table)
        else:
            ncbi_subspecies, unambiguous_ncbi_subsp, ambiguous_ncbi_subsp = ncbi_species_mngr.classify_ncbi_subspecies()

        unambiguous_classifications = [
            classification for rid, classification in unambiguous_ncbi_subsp.values()]
        unambiguous_type_strain_count = unambiguous_classifications.count(
            NCBI_SpeciesManager.TYPE_STRAIN_GENOME)
        unambiguous_majority_vote_count = unambiguous_classifications.count(
            NCBI_SpeciesManager.MAJORITY_VOTE)
        assert unambiguous_type_strain_count + \
            unambiguous_majority_vote_count == len(unambiguous_ncbi_subsp)

        self.logger.info(
            ' - identified {:,} NCBI subspecies.'.format(len(ncbi_subspecies)))
        self.logger.info(' - resolved {:,} ({:.2f}%) of species names.'.format(
            len(unambiguous_ncbi_subsp),
            len(unambiguous_ncbi_subsp) * 100.0 / len(ncbi_subspecies)))
        self.logger.info('   - assigned {:,} ({:.2f}%) by type strain genome.'.format(
            unambiguous_type_strain_count,
            unambiguous_type_strain_count * 100.0 / len(ncbi_subspecies)))
        self.logger.info('   - assigned {:,} ({:.2f}%) by majority vote.'.format(
            unambiguous_majority_vote_count,
            unambiguous_majority_vote_count * 100.0 / len(ncbi_subspecies)))
        self.logger.info(' - identified {:,} ({:.2f}%) as ambiguous assignments.'.format(
            len(ambiguous_ncbi_subsp),
            len(ambiguous_ncbi_subsp) * 100.0 / len(ncbi_subspecies)))

        return ncbi_subspecies, unambiguous_ncbi_subsp, ambiguous_ncbi_subsp

    def write_gtdb_synonym_table(self, cur_genomes, final_taxonomy, ncbi_synonym_file, out_file):
        """Write synonyms table with GTDB species names."""

        fout = open(out_file, 'w')
        fout.write(
            'Synonym type\tGTDB species\tNCBI species\tGTDB representative\tStrain IDs\tType sources\tPriority year')
        fout.write('\tGTDB type species\tGTDB type strain\tNCBI assembly type')
        fout.write(
            '\tSynonym\tHighest-quality synonym genome\tSynonym strain IDs\tSynonym type sources\tSynonym priority year')
        fout.write(
            '\tSynonym GTDB type species\tSynonym GTDB type strain\tSynonym NCBI assembly type')
        fout.write('\tANI\tAF\tNotes\n')

        with open(ncbi_synonym_file) as f:
            header = f.readline().strip().split('\t')

            type_idx = header.index('Synonym type')
            assert type_idx == 0

            rid_idx = header.index('GTDB representative')
            synonym_genome_idx = header.index('Highest-quality synonym genome')
            warning_idx = header.index('Warnings')

            for line in f:
                tokens = line.strip().split('\t')

                synonym_type = tokens[type_idx]
                if synonym_type == 'CONSENSUS_SYNONYM':
                    # HACK: renaming this as we eventually settled on using a
                    # strict majority vote rule
                    synonym_type = 'MAJORITY_VOTE_SYNONYM'

                rid = tokens[rid_idx]
                if rid in final_taxonomy:  # other representative must be from other domain
                    tokens[rid_idx] = cur_genomes.full_gid[tokens[rid_idx]]
                    tokens[synonym_genome_idx] = cur_genomes.full_gid[tokens[synonym_genome_idx]]

                    if len(tokens) > warning_idx and tokens[warning_idx]:
                        self.logger.warning('Assuming synonym priority for {} ({}) has been fixed manually.'.format(
                            final_taxonomy[rid][Taxonomy.SPECIES_INDEX],
                            rid))
                        tokens[warning_idx] = 'Incorrect priority establish manually and fixed accordingly'

                    fout.write('{}\t{}\t{}\n'.format(
                        synonym_type,
                        final_taxonomy[rid][Taxonomy.SPECIES_INDEX],
                        '\t'.join(tokens[1:])))

        fout.close()

    def update_curation_tree(self, curation_tree, final_taxonomy):
        """Update curation tree wtih finalized species names."""

        self.logger.info('Reading curation tree.')
        curation_tree = dendropy.Tree.get_from_path(curation_tree,
                                                    schema='newick',
                                                    rooting='force-rooted',
                                                    preserve_underscores=True)

        # strip species labels
        for node in curation_tree.postorder_node_iter():
            support, taxon_label, aux_info = parse_label(node.label)
            if taxon_label:
                none_sp_taxa = [t.strip() for t in taxon_label.split(
                    ';') if not t.strip().startswith('s__')]
                node.label = create_label(
                    support, '; '.join(none_sp_taxa), aux_info)

        # put new species labels in tree
        processed_rids = set()
        for leaf in curation_tree.leaf_node_iter():
            rid = leaf.taxon.label
            if not rid.startswith('D-'):
                sp_curation_node = leaf.parent_node
                support, taxon_label, aux_info = parse_label(
                    sp_curation_node.label)

                taxa = []
                if taxon_label:
                    if 's__' in taxon_label:
                        print(rid, taxon_label)
                    taxa = [t.strip() for t in taxon_label.split(';')]

                node_taxa = taxa + [final_taxonomy[rid]
                                    [Taxonomy.SPECIES_INDEX]]
                sp_curation_node.label = create_label(
                    support, '; '.join(node_taxa), aux_info)
                processed_rids.add(rid)

        assert processed_rids == set(final_taxonomy)

        # write new tree to file
        output_tree = os.path.join(
            self.output_dir, 'curation_tree.final_species.tree')
        curation_tree.write_to_path(output_tree,
                                    schema='newick',
                                    suppress_rooting=True,
                                    suppress_leaf_node_labels=False,
                                    unquoted_underscores=True)
        self.logger.info(
            'Updated curation tree written to: {}'.format(output_tree))

    def write_taxonomy(self, 
                        final_taxonomy, 
                        cur_genomes, 
                        cur_clusters,
                        unambiguous_ncbi_sp,
                        ambiguous_ncbi_sp,
                        ncbi_synonyms, 
                        unambiguous_ncbi_subsp, 
                        ambiguous_ncbi_subsp):
        """Write taxonomy information to file."""

        # write out standard taxonomy file
        final_taxonomy_file = os.path.join(
            self.output_dir, 'final_taxonomy.tsv')
        fout = open(final_taxonomy_file, 'w')
        for rid, taxa in final_taxonomy.items():
            fout.write('{}\t{}\n'.format(
                rid,
                ';'.join(taxa)))
        fout.close()
        self.logger.info(
            'Final taxonomy written to: {}'.format(final_taxonomy_file))

        # write out taxonomy file with NCBI species classification for
        # all genomes in GTDB species clusters
        sp_classification_file = os.path.join(self.output_dir, 'gtdb_sp_clusters.ncbi_sp.tsv')
        fout = open(sp_classification_file, 'w')
        fout.write('Genome ID\tGTDB taxonomy\tGTDB species\tNo. clustered\tNCBI species\n')
        for rid, taxa in final_taxonomy.items():
            ncbi_sp = []
            for cid in cur_clusters[rid]:
                ncbi_sp.append(f'{cid}: {cur_genomes[cid].ncbi_taxa.species}')

            fout.write('{}\t{}\t{}\t{}\t{}\n'.format(
                rid,
                ';'.join(taxa),
                taxa[Taxonomy.SPECIES_INDEX],
                len(ncbi_sp),
                ','.join(ncbi_sp)
            ))
        fout.close()

        # write out final taxonomy with additional metadata
        fout = open(os.path.join(self.output_dir,
                                 'final_taxonomy_metadata.tsv'), 'w')
        fout.write(
            'Genome ID\tGTDB taxonomy\tType strain\tEffective type strain\tType species')
        fout.write(
            '\tNCBI species\tNCBI representative\tNCBI species classification\tNCBI subspecies classification\n')
        for rid, taxa in final_taxonomy.items():
            ncbi_species = cur_genomes[rid].ncbi_taxa.species
            ncbi_sp_classification = 'n/a'
            if ncbi_species in unambiguous_ncbi_sp:
                ncbi_sp_classification = 'UNAMBIGUOUS'
            elif ncbi_species in ambiguous_ncbi_sp:
                ncbi_sp_classification = 'AMBIGUOUS'
            elif ncbi_species in ncbi_synonyms:
                ncbi_sp_classification = 'SYNONYM: {}'.format(
                    ncbi_synonyms[ncbi_species])

            ncbi_subspecies = cur_genomes[rid].ncbi_unfiltered_taxa.subspecies
            ncbi_subsp_classification = 'n/a'
            if ncbi_subspecies in unambiguous_ncbi_subsp:
                ncbi_subsp_classification = 'UNAMBIGUOUS'
            elif ncbi_subspecies in ambiguous_ncbi_subsp:
                ncbi_subsp_classification = 'AMBIGUOUS'

            fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                rid,
                ';'.join(taxa),
                cur_genomes[rid].is_gtdb_type_strain(),
                cur_genomes[rid].is_effective_type_strain(),
                cur_genomes[rid].is_gtdb_type_species(),
                ncbi_species,
                cur_genomes[rid].is_ncbi_representative(),
                ncbi_sp_classification,
                ncbi_subsp_classification))
        fout.close()

    def run(self,
            curation_tree,
            manual_taxonomy,
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
            ncbi_untrustworthy_sp_ledger,
            ncbi_env_bioproject_ledger,
            lpsn_gss_file):
        """Finalize species names based on results of manual curation."""

        # read manually-curated taxonomy
        self.logger.info('Parsing manually-curated taxonomy.')
        mc_taxonomy = Taxonomy().read(manual_taxonomy, use_canonical_gid=True)
        self.logger.info(' - identified taxonomy strings for {:,} genomes.'.format(
            len(mc_taxonomy)))

        # read species names explicitly set via manual curation
        mc_species = parse_manual_sp_curation_files(manual_sp_names,
                                                    pmc_custom_species)

        # sanity check that manually curated species names don't conflict with ledger
        with open(species_classification_ledger) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')

                gid = canonical_gid(tokens[0])
                if gid in mc_species and mc_species[gid].replace('s__', '') != tokens[1].replace('s__', ''):
                    self.logger.error(
                        f'Manually-curated name for {gid} conflicts with {species_classification_ledger}: {mc_species[gid]} {tokens[1]}')
                    sys.exit(-1)

        # create previous and current GTDB genome sets
        self.logger.info('Creating previous GTDB genome set.')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                             gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                             ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                             untrustworthy_type_ledger=untrustworthy_type_file,
                                             ncbi_untrustworthy_sp_ledger=ncbi_untrustworthy_sp_ledger,
                                             ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)
        self.logger.info(' - previous genome set has {:,} species clusters spanning {:,} genomes.'.format(
            len(prev_genomes.sp_clusters),
            prev_genomes.sp_clusters.total_num_genomes()))

        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                            gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                            create_sp_clusters=False,
                                            qc_passed_file=qc_passed_file,
                                            ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                            untrustworthy_type_ledger=untrustworthy_type_file,
                                            ncbi_untrustworthy_sp_ledger=ncbi_untrustworthy_sp_ledger,
                                            ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)
        self.logger.info(
            f' - current genome set contains {len(cur_genomes):,} genomes.')

        cur_genomes.set_prev_gtdb_classifications(prev_genomes)

        # read ledger with explicit names for type strain genomes
        self.logger.info(
            'Validating species names given in GTDB type strain ledger.')
        gtdb_type_strain_ledger = self.parse_gtdb_type_strain_ledger(
            gtdb_type_strains_ledger, cur_genomes)
        self.logger.info(
            ' - identified {:,} genomes in ledger.'.format(len(gtdb_type_strain_ledger)))
        for gid in gtdb_type_strain_ledger:
            if gid not in mc_taxonomy:
                continue

            cur_gtdb_sp = mc_taxonomy[gid][Taxonomy.SPECIES_INDEX]
            ledger_sp = gtdb_type_strain_ledger[gid]
            if cur_gtdb_sp != ledger_sp:
                self.logger.warning('Assigned species name for {} conflicts with GTDB type strain ledger: {} {}'.format(
                    gid,
                    cur_gtdb_sp,
                    ledger_sp))

        # read named GTDB species clusters
        self.logger.info('Reading GTDB species clusters.')
        cur_clusters, _rep_radius = read_clusters(gtdb_clusters_file)
        self.logger.info(' - identified {:,} clusters spanning {:,} genomes.'.format(
            len(cur_clusters),
            sum([len(gids) + 1 for gids in cur_clusters.values()])))
        assert len(set(mc_taxonomy) - set(cur_clusters)) == 0

        # create specific epithet manager
        self.sp_epithet_mngr = SpecificEpithetManager()
        self.sp_epithet_mngr.parse_specific_epithet_ledger(
            specific_epithet_ledger)
        self.sp_epithet_mngr.infer_epithet_map(mc_taxonomy,
                                               mc_species,
                                               cur_genomes)
        self.sp_epithet_mngr.write_diff_epithet_map(
            os.path.join(self.output_dir, 'specific_epithet_diff_map.tsv'))
        self.sp_epithet_mngr.write_epithet_map(
            os.path.join(self.output_dir, 'specific_epithet_map.tsv'))
        self.sp_epithet_mngr.write_epithet_map(
            os.path.join(self.output_dir,
                         'specific_epithet_map.new_cases.tsv'),
            filtered_previously_checked=True)

        # create species name manager
        self.logger.info('Initializing species name manager.')
        self.sp_name_mngr = SpeciesNameManager(prev_genomes, cur_genomes)

        # initialize species priority manager
        self.sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger,
                                                       genus_priority_ledger,
                                                       lpsn_gss_file,
                                                       self.output_dir)

        # establish state of NCBI species
        ncbi_species_mngr = NCBI_SpeciesManager(cur_genomes,
                                                cur_clusters,
                                                mc_species,
                                                self.output_dir)
        self.forbidden_specific_names = ncbi_species_mngr.forbidden_specific_names

        # get list of synonyms in order to restrict usage of species names
        self.logger.info('Reading GTDB synonyms.')
        ncbi_synonyms = ncbi_species_mngr.parse_synonyms_table(
            ncbi_synonym_file)
        self.logger.info(' - identified {:,} synonyms from {:,} distinct species.'.format(
            len(ncbi_synonyms),
            len(set(ncbi_synonyms.values()))))

        # identified genomes with misclassified species assignments at NCBI
        self.logger.info(
            'Identify genomes with misclassified NCBI species assignments.')
        ncbi_misclassified_gids = ncbi_species_mngr.parse_ncbi_misclassified_table(
            ncbi_misclassified_file)
        self.logger.info(' - identified {:,} genomes with erroneous NCBI species assignments'.format(
            len(ncbi_misclassified_gids)))

        # classify each Latin NCBI species as either unambiguously or ambiguously assigned
        self.logger.info(
            'Classifying NCBI species as unambiguous, ambiguous, or synonym.')
        ncbi_species, unambiguous_ncbi_sp, ambiguous_ncbi_sp = self.classify_ncbi_species(
            ncbi_species_mngr,
            ncbi_synonyms,
            ncbi_misclassified_gids)

        # classify each NCBI subspecies as either unambiguously or ambiguously assigned
        self.logger.info(
            'Classifying NCBI subspecies as unambiguous or ambiguous.')
        (_ncbi_all_subspecies,
         unambiguous_ncbi_subsp,
         ambiguous_ncbi_subsp) = self.classify_ncbi_subspecies(ncbi_species_mngr)

        # get mapping between new and old representatives. This can't just be taken
        # from the `updated_species_reps` file since the final results of de novo
        # cluster can, in some rare cases, cause small movements in the genomes
        # associated with each species clusters and the generation or loss of
        # a representative
        self.logger.info(
            'Mapping current GTDB representatives to previous representatives.')
        updated_gtdb_rids = parse_updated_species_reps(updated_species_reps)
        new_to_prev_rid = infer_prev_gtdb_reps(
            prev_genomes, cur_clusters, updated_gtdb_rids)
        self.logger.info(
            ' - mapped {:,} current representatives to previous representatives.'.format(len(new_to_prev_rid)))

        # establish appropriate species names for GTDB clusters with new representatives
        final_taxonomy = self.finalize_species_names(mc_taxonomy,
                                                     mc_species,
                                                     cur_clusters,
                                                     new_to_prev_rid,
                                                     prev_genomes,
                                                     cur_genomes,
                                                     ncbi_misclassified_gids,
                                                     unambiguous_ncbi_sp,
                                                     unambiguous_ncbi_subsp,
                                                     ncbi_synonyms,
                                                     species_classification_ledger)

        # sanity check that all species clusters have a unique name
        sp_names = {}
        for rid, taxa in final_taxonomy.items():
            sp = taxa[Taxonomy.SPECIES_INDEX]
            if sp in sp_names:
                self.logger.error('Species name {} assigned to at least 2 GTDB representatives: {} {}'.format(
                    sp, sp_names[sp], rid))
            sp_names[sp] = rid

        # write out final taxonomy
        self.write_taxonomy(final_taxonomy, 
                            cur_genomes, 
                            cur_clusters,
                            unambiguous_ncbi_sp, 
                            ambiguous_ncbi_sp,
                            ncbi_synonyms, 
                            unambiguous_ncbi_subsp, 
                            ambiguous_ncbi_subsp)

        # re-decorate curation tree
        if curation_tree is not None:
            self.update_curation_tree(curation_tree, final_taxonomy)

        # write synonym table with GTDB species names
        out_file = os.path.join(self.output_dir, 'gtdb_synonyms_final.tsv')
        self.write_gtdb_synonym_table(
            cur_genomes, final_taxonomy, ncbi_synonym_file, out_file)
