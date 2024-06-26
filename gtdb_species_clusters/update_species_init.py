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
from collections import defaultdict

from numpy import mean as np_mean

from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_name_manager import SpeciesNameManager
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters.ncbi_species_manager import NCBI_SpeciesManager
from gtdb_species_clusters.type_genome_utils import (read_clusters,
                                                     parse_gtdb_type_strain_ledger)
from gtdb_species_clusters.taxon_utils import (generic_name,
                                               specific_epithet,
                                               canonical_taxon,
                                               longest_common_prefix,
                                               is_placeholder_taxon)


class UpdateSpeciesInit(object):
    """Produce initial best guess at GTDB species clusters."""

    def __init__(self, output_dir):
        """Initialization."""

        self.output_dir = output_dir
        self.log = logging.getLogger('rich')

        self.sp_name_log = open(os.path.join(
            self.output_dir, 'sp_name_log.tsv'), 'w')
        self.sp_name_log.write(
            'GTDB domain\tGenome ID\tPrevious GTDB species\tNew GTDB species\tAction\n')

        self.sp_curation_log = open(os.path.join(
            self.output_dir, 'sp_curation_log.tsv'), 'w')
        self.sp_curation_log.write(
            'GTDB domain\tGenome ID\tPrevious NCBI species\tCurrent NCBI species')
        self.sp_curation_log.write(
            '\tPrevious GTDB genus\tPrevious GTDB species\tProposed GTDB species')
        self.sp_curation_log.write(
            '\tModified generic name\tMandatory fix\tGTDB type strain ledger\tCase\tIssue\n')

    def update_log(self, gid, cur_genomes, prev_genomes, sp, action):
        """Add entry to update log."""

        prev_gtdb_sp = 'n/a'
        if gid in prev_genomes:
            prev_gtdb_sp = prev_genomes[gid].gtdb_taxa.species

        self.sp_name_log.write('{}\t{}\t{}\t{}\t{}\n'.format(
            cur_genomes[gid].gtdb_taxa.domain,
            gid,
            prev_gtdb_sp,
            sp,
            action))

    def curation_log(self,
                     gid,
                     cur_genomes,
                     prev_genomes,
                     proposed_gtdb_sp,
                     mandatory,
                     case,
                     issue,
                     in_gtdb_type_strain_ledger):
        """Add entry to curation log."""

        prev_ncbi_sp = 'n/a'
        if gid in prev_genomes:
            prev_ncbi_sp = prev_genomes[gid].ncbi_taxa.species

        prev_gtdb_sp = 'n/a'
        prev_gtdb_genus = 'n/a'
        if gid in cur_genomes:
            prev_gtdb_sp = cur_genomes[gid].gtdb_taxa.species
            prev_gtdb_genus = cur_genomes[gid].gtdb_taxa.genus

        cur_ncbi_sp = cur_genomes[gid].ncbi_taxa.species

        self.sp_curation_log.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            cur_genomes[gid].gtdb_taxa.domain,
            gid,
            prev_ncbi_sp,
            cur_ncbi_sp,
            prev_gtdb_genus,
            prev_gtdb_sp,
            proposed_gtdb_sp,
            generic_name(proposed_gtdb_sp) != generic_name(proposed_gtdb_sp),
            mandatory,
            in_gtdb_type_strain_ledger,
            case,
            issue))

    def _test_same_epithet(self, epithet1, epithet2):
        """Test if species epithet are the same, except for changes due to difference in the gender of the genus."""

        if epithet1 == epithet2:
            return True

        lcp = longest_common_prefix(epithet1, epithet2)
        if len(lcp) >= min(len(epithet1), len(epithet2)) - 3:
            # a small change to the suffix presumably reflecting
            # a change in the gender of the genus
            return True

        return False

    def _determine_sp_epithet(self, rid, cur_genomes, synonyms):
        """Determine most appropriate epithet for species name.

        This method assumes the epithet of the NCBI species name
        should be favoured; presumable because the genome is
        assembled from the type strain of the species.
        """

        ncbi_sp = cur_genomes[rid].ncbi_taxa.species
        if ncbi_sp in synonyms:
            ncbi_sp = synonyms[ncbi_sp]
        ncbi_genus = 'g__' + generic_name(ncbi_sp)

        gtdb_sp = cur_genomes[rid].gtdb_taxa.species
        if gtdb_sp in synonyms:
            gtdb_sp = synonyms[gtdb_sp]
        gtdb_genus = 'g__' + generic_name(gtdb_sp)

        gtdb_sp_epithet = specific_epithet(gtdb_sp)
        ncbi_sp_epithet = specific_epithet(ncbi_sp)
        if gtdb_sp_epithet != ncbi_sp_epithet:
            if gtdb_genus == ncbi_genus or is_placeholder_taxon(gtdb_genus):
                # take NCBI specific name since there is no concern about
                # a suffix change due to a change in the gender of the genus
                sp_epithet = specific_epithet(ncbi_sp)
            else:
                # genera do not agree so the suffix of the
                # specific name may have been changed to
                # reflect the gender of the GTDB genus
                if self._test_same_epithet(gtdb_sp_epithet, ncbi_sp_epithet):
                    # should keep GTDB specific name since it is only
                    # a small change to the suffix presumably reflecting
                    # a change in the gender of the genus
                    sp_epithet = gtdb_sp_epithet
                else:
                    # looks like a real conflict so favour the NCBI
                    # specific name since this is a type strain genome;
                    # there may be an issue with the gender of the GTDB
                    # genus, but this will be resolved later since the
                    # GTDB genus names are ultimately set by manual curation
                    sp_epithet = ncbi_sp_epithet
        else:
            sp_epithet = gtdb_sp_epithet

        return sp_epithet

    def name_type_species_clusters(self,
                                   rids_by_naming_priority,
                                   cluster_sp_names,
                                   used_sp_names,
                                   clusters,
                                   prev_genomes,
                                   cur_genomes,
                                   gtdb_type_strain_ledger,
                                   synonyms):
        """Assign names to clusters represented by a type species of genus genome."""

        num_updates = 0
        conflicting_type_sp = defaultdict(set)
        for rid in rids_by_naming_priority:
            if not cur_genomes[rid].is_gtdb_type_species():
                continue

            ncbi_sp = cur_genomes[rid].ncbi_taxa.species
            if ncbi_sp == 's__':
                self.log.error(
                    'Type species of genus genome has no NCBI species assignment: {}'.format(rid))
                sys.exit(-1)

            if ncbi_sp in synonyms:
                ncbi_sp = synonyms[ncbi_sp]

            if rid in gtdb_type_strain_ledger:
                proposed_gtdb_sp = gtdb_type_strain_ledger[rid]
            else:
                # update GTDB assignment only if it is unassigned or a placeholder name
                gtdb_genus = cur_genomes[rid].gtdb_taxa.genus
                gtdb_sp = cur_genomes[rid].gtdb_taxa.species
                if gtdb_sp in synonyms:
                    gtdb_sp = synonyms[gtdb_sp]

                genus = gtdb_genus
                if is_placeholder_taxon(genus):
                    genus = 'g__' + generic_name(ncbi_sp)

                sp_epithet = self._determine_sp_epithet(
                    rid, cur_genomes, synonyms)

                proposed_gtdb_sp = '{} {}'.format(
                    genus.replace('g__', 's__'), sp_epithet)

                assert genus != 'g__' and sp_epithet != ''

            ncbi_genus = 'g__' + generic_name(ncbi_sp)
            proposed_gtdb_genus = 'g__' + generic_name(proposed_gtdb_sp)
            if proposed_gtdb_genus != ncbi_genus:
                conflicting_type_sp[ncbi_genus].add(proposed_gtdb_genus)

                ncbi_year = self.sp_priority_mngr.genus_priority_year(
                    ncbi_genus)
                gtdb_year = self.sp_priority_mngr.genus_priority_year(
                    proposed_gtdb_genus)

                if gtdb_year < ncbi_year:
                    self.curation_log(rid,
                                      cur_genomes,
                                      prev_genomes,
                                      proposed_gtdb_sp,
                                      False,
                                      "TYPE_SPECIES_OF_GENUS",
                                      "GTDB genus assignment does not reflect type species of genus",
                                      rid in gtdb_type_strain_ledger)
                else:
                    self.curation_log(rid,
                                      cur_genomes,
                                      prev_genomes,
                                      proposed_gtdb_sp,
                                      True,
                                      "TYPE_SPECIES_OF_GENUS",
                                      "GTDB genus assignment does not reflect type species of genus and may violate naming priority ({}: {}, {}: {})".format(
                                          ncbi_genus, ncbi_year, proposed_gtdb_genus, gtdb_year),
                                      rid in gtdb_type_strain_ledger)

            if proposed_gtdb_sp != cur_genomes[rid].gtdb_taxa.species:
                if rid in prev_genomes:
                    num_updates += 1
                    self.update_log(rid, cur_genomes, prev_genomes, proposed_gtdb_sp,
                                    "Updated to reflect genome assembled from type species of genus.")
                for cid in clusters[rid]:
                    # ***cur_genomes[cid].gtdb_taxa.genus = 'g__' + generic_name(proposed_gtdb_sp)
                    cur_genomes[cid].gtdb_taxa.species = proposed_gtdb_sp

            if proposed_gtdb_sp in used_sp_names:
                self.log.warning('NCBI species name of type species representative already in use: {}: {}, {}'.format(
                    ncbi_sp, rid, used_sp_names[proposed_gtdb_sp]))

            if proposed_gtdb_sp in synonyms:
                self.log.error('Assigning synonym to type species representative cluster: {}: {}'.format(
                    proposed_gtdb_sp, rid))

            cluster_sp_names[rid] = proposed_gtdb_sp
            used_sp_names[proposed_gtdb_sp] = rid

        self.log.info(f' - updated name of {num_updates:,} GTDB clusters')
        num_conflicts = sum([len(v) for v in conflicting_type_sp.values()])
        self.log.warning(
            f' - proposed GTDB genus did not reflect type species of genus for {num_conflicts:,} clusters (occurs due to merging of genera)')

        fout = open(os.path.join(self.output_dir,
                                 'conflicting_type_species_of_genus.tsv'), 'w')
        fout.write('NCBI genus\tProposed GTDB genera\n')
        for ncbi_genus, gtdb_genera in conflicting_type_sp.items():
            fout.write('{}\t{}\n'.format(ncbi_genus, ', '.join(gtdb_genera)))
        fout.close()

    def name_type_strain_clusters(self,
                                  rids_by_naming_priority,
                                  cluster_sp_names,
                                  used_sp_names,
                                  clusters,
                                  prev_genomes,
                                  cur_genomes,
                                  gtdb_type_strain_ledger,
                                  synonyms):
        """Assign names to clusters represented by a validally or effectively published type strain genome."""

        num_updates = 0
        num_conflicting_specific = 0
        num_conflicting_generic = 0
        dup_name_count = 0
        for rid in rids_by_naming_priority:
            if rid in cluster_sp_names:
                continue

            if not cur_genomes[rid].is_effective_type_strain():
                continue

            ncbi_sp = cur_genomes[rid].ncbi_taxa.species
            if ncbi_sp == 's__':
                # NCBI does define type material that lacks a binomial species name
                continue

            if ncbi_sp in synonyms:
                ncbi_sp = synonyms[ncbi_sp]

            if rid in gtdb_type_strain_ledger:
                proposed_gtdb_sp = gtdb_type_strain_ledger[rid]
            else:
                # update GTDB generic name only if it is unassigned, and
                # the GTDB specific name unless it looks like only a
                # small change to the suffix which likely reflect a
                # change in the gender of the genus
                gtdb_genus = cur_genomes[rid].gtdb_taxa.genus
                gtdb_sp = cur_genomes[rid].gtdb_taxa.species
                if gtdb_sp in synonyms:
                    gtdb_sp = synonyms[gtdb_sp]

                genus = gtdb_genus
                if genus == 'g__':
                    genus = 'g__' + generic_name(ncbi_sp)

                sp_epithet = self._determine_sp_epithet(
                    rid, cur_genomes, synonyms)

                proposed_gtdb_sp = '{} {}'.format(
                    genus.replace('g__', 's__'), sp_epithet)

                assert genus != 'g__' and sp_epithet != ''

            ncbi_specific = specific_epithet(ncbi_sp)
            proposed_specific = specific_epithet(proposed_gtdb_sp)
            if not self._test_same_epithet(ncbi_specific, proposed_specific):
                num_conflicting_specific += 1
                self.curation_log(rid,
                                  cur_genomes,
                                  prev_genomes,
                                  proposed_gtdb_sp,
                                  True,
                                  "TYPE_STRAIN_OF_SPECIES",
                                  "GTDB specific species name does not reflect type strain of species",
                                  rid in gtdb_type_strain_ledger)

            if generic_name(proposed_gtdb_sp) != generic_name(ncbi_sp):
                mandatory = False
                issue = "GTDB generic name does not reflect type strain of species"
                if rid in prev_genomes and generic_name(prev_genomes[rid].ncbi_taxa.species) != generic_name(ncbi_sp):
                    issue += '; NCBI generic name has been updated'
                    mandatory = True
                    num_conflicting_generic += 1
                    self.curation_log(rid,
                                      cur_genomes,
                                      prev_genomes,
                                      proposed_gtdb_sp,
                                      mandatory,
                                      "TYPE_STRAIN_OF_SPECIES",
                                      issue,
                                      rid in gtdb_type_strain_ledger)

            if proposed_gtdb_sp != cur_genomes[rid].gtdb_taxa.species:
                if rid in prev_genomes:
                    num_updates += 1
                    self.update_log(rid, cur_genomes, prev_genomes, proposed_gtdb_sp,
                                    "Updated to reflect genome assembled from type strain of species.")
                for cid in clusters[rid]:
                    # ***cur_genomes[cid].gtdb_taxa.genus = 'g__' + generic_name(proposed_gtdb_sp)
                    cur_genomes[cid].gtdb_taxa.species = proposed_gtdb_sp

            if proposed_gtdb_sp in used_sp_names:
                dup_name_count += 1
                self.curation_log(rid,
                                  cur_genomes,
                                  prev_genomes,
                                  proposed_gtdb_sp,
                                  True,
                                  "ALREADY_USED_SPECIES_NAME",
                                  "Assigned species name supported by multiple type strain genomes: {} {}".format(
                                      rid, used_sp_names[proposed_gtdb_sp]),
                                  rid in gtdb_type_strain_ledger)

                self.log.warning('NCBI species name of type strain representative already in use: {}: {}, {}'.format(
                    proposed_gtdb_sp, rid, used_sp_names[proposed_gtdb_sp]))

            if proposed_gtdb_sp in synonyms:
                self.log.error('Assigning synonym to type strain representative cluster: {}: {}'.format(
                    proposed_gtdb_sp, rid))
                #***sys.exit(-1)

            cluster_sp_names[rid] = proposed_gtdb_sp
            used_sp_names[proposed_gtdb_sp] = rid

        self.log.info(f' - updated name of {num_updates:,} GTDB clusters')

        if num_conflicting_generic > 0:
            self.log.warning(
                f' - GTDB generic name does not reflect type strain of species and NCBI genus assignment changed this release for {num_conflicting_generic:,} species clusters')

        if num_conflicting_specific > 0:
            self.log.warning(
                f' - GTDB specific name does not reflect type strain of species for {num_conflicting_specific:,} species clusters')

        if dup_name_count > 0:
            self.log.warning(
                f' - same name assigned to {dup_name_count:,} species clusters')

    def name_binomial_clusters(self,
                               rids_by_naming_priority,
                               cluster_sp_names,
                               used_sp_names,
                               clusters,
                               prev_genomes,
                               cur_genomes,
                               synonyms):
        """Assign names to clusters with suitable binomial Latin names which lack type strain genomes."""

        # get species assignment of previous GTDB species clusters
        gtdb_canonical_sp_epithets = defaultdict(set)
        for rid in clusters:
            gtdb_genus = cur_genomes[rid].gtdb_taxa.genus
            gtdb_sp_epithet = cur_genomes[rid].gtdb_taxa.specific_epithet
            gtdb_canonical_sp_epithets[gtdb_genus].add(
                canonical_taxon(gtdb_sp_epithet))

        num_updates = 0
        dup_name_count = 0
        for rid in rids_by_naming_priority:
            if rid in cluster_sp_names:
                continue

            ncbi_sp = cur_genomes[rid].ncbi_taxa.species
            if is_placeholder_taxon(ncbi_sp):
                continue

            if ncbi_sp in synonyms:
                ncbi_sp = synonyms[ncbi_sp]

            # update GTDB generic name only if unassigned, and
            # GTDB specific name only if no other species has
            # this assignment in the genus
            gtdb_genus = cur_genomes[rid].gtdb_taxa.genus
            gtdb_sp = cur_genomes[rid].gtdb_taxa.species
            if gtdb_sp in synonyms:
                gtdb_sp = synonyms[gtdb_sp]

            genus = gtdb_genus
            if genus == 'g__':
                genus = 'g__' + generic_name(ncbi_sp)

            gtdb_sp_epithet = specific_epithet(gtdb_sp)
            ncbi_sp_epithet = specific_epithet(ncbi_sp)

            # *** HACK: this isn't stable since a given epithet could
            # occur twice
            if ncbi_sp_epithet not in gtdb_canonical_sp_epithets[gtdb_genus]:
                # should update to NCBI epithet in order to reflect
                # recognized species name
                sp_epithet = ncbi_sp_epithet
            else:
                sp_epithet = gtdb_sp_epithet
                if not gtdb_sp_epithet:
                    sp = '{} {}'.format(genus.replace(
                        'g__', 's__'), ncbi_sp_epithet)
                    next_suffix = self.sp_name_mngr.taxon_suffix_manager.next_suffix(
                        sp)
                    sp_epithet = '{}_{}'.format(ncbi_sp_epithet, next_suffix)

            proposed_gtdb_sp = '{} {}'.format(
                genus.replace('g__', 's__'), sp_epithet)
            if proposed_gtdb_sp in used_sp_names:
                # better suffix name
                next_suffix = self.sp_name_mngr.taxon_suffix_manager.next_suffix(
                    proposed_gtdb_sp)
                sp_epithet = '{}_{}'.format(canonical_taxon(
                    specific_epithet(proposed_gtdb_sp)), next_suffix)
                proposed_gtdb_sp = '{} {}'.format(
                    genus.replace('g__', 's__'), sp_epithet)

            assert genus != 'g__' and sp_epithet != ''

            if proposed_gtdb_sp != cur_genomes[rid].gtdb_taxa.species:
                if rid in prev_genomes:
                    num_updates += 1
                    self.update_log(rid, cur_genomes, prev_genomes, proposed_gtdb_sp,
                                    "Updated to reflect bionomial NCBI species name.")
                for cid in clusters[rid]:
                    # ***cur_genomes[cid].gtdb_taxa.genus = 'g__' + generic_name(proposed_gtdb_sp)
                    cur_genomes[cid].gtdb_taxa.species = proposed_gtdb_sp

            if proposed_gtdb_sp in used_sp_names:
                dup_name_count += 1
                self.curation_log(rid,
                                  cur_genomes,
                                  prev_genomes,
                                  proposed_gtdb_sp,
                                  True,
                                  "ALREADY_USED_SPECIES_NAME",
                                  "Assigned species name supported by NCBI binomial species name: {} {}".format(
                                      rid, used_sp_names[proposed_gtdb_sp]),
                                  False)
                self.log.warning('NCBI binomial species name already in use: {}: {}, {}'.format(
                    proposed_gtdb_sp, rid, used_sp_names[proposed_gtdb_sp]))

            if proposed_gtdb_sp in synonyms:
                self.log.error(
                    'Assigning synonym to non-type strain species cluster: {}: {}'.format(ncbi_sp, rid))

            cluster_sp_names[rid] = proposed_gtdb_sp
            used_sp_names[proposed_gtdb_sp] = rid

        self.log.info(f' - updated name of {num_updates:,} GTDB clusters')

        if dup_name_count > 0:
            self.log.warning(
                f' - same name assigned to {dup_name_count:,} species clusters')

    def name_placeholder_clusters(self,
                                  rids_by_naming_priority,
                                  cluster_sp_names,
                                  used_sp_names,
                                  clusters,
                                  prev_genomes,
                                  cur_genomes,
                                  synonyms):
        """Assign placeholder names to clusters that could not be assigned a binomial Latin name."""

        num_updates = 0
        dup_name_count = 0
        for rid in rids_by_naming_priority:
            if rid in cluster_sp_names:
                continue

            gtdb_sp = cur_genomes[rid].gtdb_taxa.species
            if gtdb_sp != 's__' and gtdb_sp not in used_sp_names:
                proposed_gtdb_sp = gtdb_sp
            else:
                genus = 'g__{unresolved}'
                gtdb_genus = cur_genomes[rid].gtdb_taxa.genus
                if gtdb_genus != 'g__':
                    genus = gtdb_genus

                proposed_gtdb_sp = 's__{} sp{}'.format(
                    genus.replace('g__', ''), rid[1:])

            if proposed_gtdb_sp != cur_genomes[rid].gtdb_taxa.species:
                if rid in prev_genomes:
                    num_updates += 1
                    self.update_log(rid, cur_genomes, prev_genomes, proposed_gtdb_sp,
                                    "Updated to reflect binomial NCBI species name.")

                for cid in clusters[rid]:
                    # ***cur_genomes[cid].gtdb_taxa.genus = 'g__' + generic_name(proposed_gtdb_sp)
                    cur_genomes[cid].gtdb_taxa.species = proposed_gtdb_sp

            if proposed_gtdb_sp in used_sp_names:
                dup_name_count += 1
                self.curation_log(rid,
                                  cur_genomes,
                                  prev_genomes,
                                  proposed_gtdb_sp,
                                  True,
                                  "ALREADY_USED_SPECIES_NAME",
                                  "Placeholder name assigned to multiple genomes: {} {}".format(
                                      rid, used_sp_names[proposed_gtdb_sp]),
                                  False)
                self.log.warning('Placeholder name already in use: {}: {}, {}'.format(
                    proposed_gtdb_sp, rid, used_sp_names[proposed_gtdb_sp]))

            if proposed_gtdb_sp in synonyms:
                self.log.error(
                    'Assigning synonym to non-type strain species cluster: {}: {}'.format(gtdb_sp, rid))

            cluster_sp_names[rid] = proposed_gtdb_sp
            used_sp_names[proposed_gtdb_sp] = rid

        self.log.info(f' - updated name of {num_updates:,} GTDB clusters')
        self.log.warning(
            f' - same name assigned to {dup_name_count:,} species clusters')

    def _sort_by_naming_priority(self,
                                 clusters,
                                 prev_genomes,
                                 cur_genomes,
                                 synonyms):
        """Sort representatives by naming priority."""

        # group by naming priority
        type_species = []
        type_strains = []
        binomial = []
        placeholder = []
        for rid in clusters:
            ncbi_sp = cur_genomes[rid].ncbi_taxa.species
            if ncbi_sp in synonyms:
                ncbi_sp = synonyms[ncbi_sp]

            gtdb_sp = cur_genomes[rid].gtdb_taxa.species
            if gtdb_sp in synonyms:
                gtdb_sp = synonyms[gtdb_sp]

            if cur_genomes[rid].is_gtdb_type_species():
                type_species.append(rid)
            elif cur_genomes[rid].is_effective_type_strain():
                type_strains.append(rid)
            elif not is_placeholder_taxon(ncbi_sp):
                binomial.append(rid)
            else:
                placeholder.append(rid)

        assert len(clusters) == len(type_species) + \
            len(type_strains) + len(binomial) + len(placeholder)

        # sort groups so previous GTDB representatives
        # are processed first
        rids_by_naming_priority = []
        for d in [type_species, type_strains, binomial, placeholder]:
            prev_reps = []
            new_reps = []
            for rid in d:
                if rid in prev_genomes.sp_clusters:
                    prev_reps.append(rid)
                else:
                    new_reps.append(rid)

            rids_by_naming_priority.extend(prev_reps)
            rids_by_naming_priority.extend(new_reps)

        assert len(rids_by_naming_priority) == len(clusters)

        return rids_by_naming_priority

    def update_species_names(self,
                             clusters,
                             prev_genomes,
                             cur_genomes,
                             gtdb_type_strain_ledger,
                             synonyms):
        """Determine most suitable name for each GTDB species cluster."""

        # identify NCBI species with multiple genomes assembled from type strain of species
        self.log.info(
            'Determining type species genomes in each NCBI species.')
        ncbi_sp_type_strain_genomes = defaultdict(set)
        for gid in cur_genomes:
            if cur_genomes[gid].is_effective_type_strain():
                if gid in gtdb_type_strain_ledger:
                    ncbi_sp = gtdb_type_strain_ledger[gid]
                else:
                    ncbi_sp = cur_genomes[gid].ncbi_taxa.species

                if ncbi_sp != 's__':
                    # yes, NCBI has genomes marked as assembled from type material
                    # that do not actually have a binomial species name
                    ncbi_sp_type_strain_genomes[ncbi_sp].add(gid)

        self.log.info(' - identified effective type strain genomes for {:,} NCBI species'.format(
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
                self.log.warning('Type strain genomes from NCBI species {} were assigned to {:,} GTDB species clusters: {} {}.'.format(
                    ncbi_sp,
                    len(gtdb_rids),
                    [(gid, rid_map[gid]) for gid in type_gids],
                    [(rid, cur_genomes[rid].gtdb_taxa.species, cur_genomes[rid].ncbi_taxa.species) for rid in gtdb_rids]))

        # order representatives by naming priority
        rids_by_naming_priority = self._sort_by_naming_priority(clusters,
                                                                prev_genomes,
                                                                cur_genomes,
                                                                synonyms)

        # determine names for species clusters represented by a validly
        # published species name with a designated type species of genus
        self.log.info(
            'Assigning names to cluster represented by type species genomes:')
        cluster_sp_names = {}
        used_sp_names = {}
        self.name_type_species_clusters(rids_by_naming_priority,
                                        cluster_sp_names,
                                        used_sp_names,
                                        clusters,
                                        prev_genomes,
                                        cur_genomes,
                                        gtdb_type_strain_ledger,
                                        synonyms)

        num_type_species = len(cluster_sp_names)
        self.log.info(
            f' - assigned {num_type_species:,} name from type species')

        # determine names for species clusters represented by a validly
        # or effectively published species name with a designated type strain
        self.log.info(
            'Assigning names to cluster represented by type strain genomes:')
        self.name_type_strain_clusters(rids_by_naming_priority,
                                       cluster_sp_names,
                                       used_sp_names,
                                       clusters,
                                       prev_genomes,
                                       cur_genomes,
                                       gtdb_type_strain_ledger,
                                       synonyms)
        num_type_strains = len(cluster_sp_names) - num_type_species
        total_sp_names = len(cluster_sp_names)
        self.log.info(
            f' - assigned {num_type_strains:,} names from type strains and {total_sp_names:,} total names')

        # determine name for species clusters where a NCBI-define
        # binomial latin name exists, but there is no designated type strain
        # (e.g., Candidatus species)
        self.log.info(
            'Assigning names to clusters with binomial NCBI species names:')
        self.name_binomial_clusters(rids_by_naming_priority,
                                    cluster_sp_names,
                                    used_sp_names,
                                    clusters,
                                    prev_genomes,
                                    cur_genomes,
                                    synonyms)
        num_assigned = len(cluster_sp_names) - total_sp_names
        total_sp_names = len(cluster_sp_names)
        self.log.info(
            f' - assigned {num_assigned:,} additional names and {total_sp_names:,} total names')

        # determine names for remaining clusters which should all have a
        # placeholder name
        self.log.info(
            'Assigning names to clusters without NCBI species names:')
        self.name_placeholder_clusters(rids_by_naming_priority,
                                       cluster_sp_names,
                                       used_sp_names,
                                       clusters,
                                       prev_genomes,
                                       cur_genomes,
                                       synonyms)
        num_assigned = len(cluster_sp_names) - total_sp_names
        total_sp_names = len(cluster_sp_names)
        self.log.info(
            f' - assigned {num_assigned:,} additional names and {total_sp_names:,} total names')

        self.sp_curation_log.close()

        # sanity check
        for rid, sp in cluster_sp_names.items():
            assert cur_genomes[rid].gtdb_taxa.species == sp

    def run(self,
            gtdb_clusters_file,
            prev_gtdb_metadata_file,
            cur_gtdb_metadata_file,
            qc_passed_file,
            gtdbtk_classify_file,
            ncbi_genbank_assembly_file,
            untrustworthy_type_file,
            synonym_file,
            gtdb_type_strains_ledger,
            sp_priority_ledger,
            genus_priority_ledger,
            ncbi_untrustworthy_sp_ledger,
            ncbi_env_bioproject_ledger,
            lpsn_gss_file):
        """Produce initial best guess at GTDB species clusters."""

        # create previous and current GTDB genome sets
        self.log.info('Creating previous GTDB genome set.')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                             gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                             ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                             untrustworthy_type_ledger=untrustworthy_type_file,
                                             ncbi_untrustworthy_sp_ledger=ncbi_untrustworthy_sp_ledger,
                                             ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)
        self.log.info(' - previous genome set has {:,} species clusters spanning {:,} genomes'.format(
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
                                            ncbi_untrustworthy_sp_ledger=ncbi_untrustworthy_sp_ledger,
                                            ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)

        cur_genomes.set_prev_gtdb_classifications(prev_genomes)

        # update current genomes with GTDB-Tk classifications
        self.log.info(
            'Updating current genomes with GTDB-Tk classifications.')
        num_updated, num_ncbi_sp = cur_genomes.set_gtdbtk_classification(
            gtdbtk_classify_file, prev_genomes)
        self.log.info(
            f' - set GTDB taxa for {num_updated:,} genomes with {num_ncbi_sp:,} genomes using NCBI genus and species name')

        # read named GTDB species clusters
        self.log.info('Reading GTDB species clusters.')
        cur_clusters, rep_radius = read_clusters(gtdb_clusters_file)
        self.log.info(' - identified {:,} clusters spanning {:,} genomes'.format(
            len(cur_clusters),
            sum([len(gids) + 1 for gids in cur_clusters.values()])))

        # write out archaeal and bacterial genome files
        self.log.info(
            'Creating genomic files for archaeal and bacterial reference genomes.')
        for prefix, domain in [('ar', 'd__Archaea'), ('bac', 'd__Bacteria')]:
            fout = open(os.path.join(self.output_dir,
                                     '{}_genomes.lst'.format(prefix)), 'w')
            genome_count = 0
            for rid in cur_clusters:
                gtdb_domain = cur_genomes[rid].gtdb_taxa.domain
                if gtdb_domain == domain:
                    fout.write('{}\n'.format(cur_genomes[rid].ncbi_accn))
                    genome_count += 1

            self.log.info(
                f' - identified {genome_count:,} {domain} reference genomes')
            fout.close()

        # read ledger with explicit names for type strain genomes
        self.log.info('Parsing GTDB type strain ledger.')
        gtdb_type_strain_ledger = parse_gtdb_type_strain_ledger(
            gtdb_type_strains_ledger, cur_genomes)
        self.log.info(
            ' - identified {:,} genomes in ledger'.format(len(gtdb_type_strain_ledger)))

        # get list of synonyms in order to restrict usage of species names
        self.log.info('Reading GTDB synonyms:')
        ncbi_species_mngr = NCBI_SpeciesManager(cur_genomes,
                                                cur_clusters,
                                                {},
                                                self.output_dir)
        synonyms = ncbi_species_mngr.parse_synonyms_table(synonym_file)
        self.log.info(' - identified {:,} synonyms from {:,} distinct species'.format(
            len(synonyms),
            len(set(synonyms.values()))))

        # create species name manager
        self.log.info('Initializing species name manager.')
        self.sp_name_mngr = SpeciesNameManager(prev_genomes, cur_genomes)

        # initialize species priority manager
        self.sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger,
                                                       genus_priority_ledger,
                                                       lpsn_gss_file,
                                                       self.output_dir)

        # establish appropriate species names for GTDB clusters with new representatives
        self.update_species_names(cur_clusters,
                                  prev_genomes,
                                  cur_genomes,
                                  gtdb_type_strain_ledger,
                                  synonyms)

        # write out taxonomy files
        bac_taxonomy_out = open(os.path.join(
            self.output_dir, 'gtdb_bac_taxonomy.tsv'), 'w')
        ar_taxonomy_out = open(os.path.join(
            self.output_dir, 'gtdb_ar_taxonomy.tsv'), 'w')
        for rid in cur_clusters:
            gtdb_domain = cur_genomes[rid].gtdb_taxa.domain
            fout = bac_taxonomy_out
            if gtdb_domain == 'd__Archaea':
                fout = ar_taxonomy_out
            fout.write('{}\t{}\n'.format(rid, cur_genomes[rid].gtdb_taxa))

        bac_taxonomy_out.close()
        ar_taxonomy_out.close()
