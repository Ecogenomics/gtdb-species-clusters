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
import csv
import logging
import functools
from collections import defaultdict

from gtdb_species_clusters.genome import Genome
from gtdb_species_clusters.genome_utils import select_highest_quality
from gtdb_species_clusters.lpsn import LPSN


class SpeciesPriorityManager(object):
    """Resolve naming priority of species names."""

    def __init__(self,
                 species_priority_ledger,
                 genus_priority_ledger,
                 lpsn_gss_file,
                 output_dir):
        """Initialization."""

        self.log = logging.getLogger('timestamp')

        self.output_dir = output_dir

        self._fout_ambiguous_priority = open(os.path.join(
            output_dir, 'ambiguous_sp_priority.tsv'), 'w')
        self._fout_ambiguous_priority.write(
            'NCBI species A\tNCBI species B\tPriority year\n')

        self._parse_sp_priority_ledger(species_priority_ledger)
        self._parse_genus_priority_ledger(genus_priority_ledger)
        self._parse_lpsn_gss(lpsn_gss_file)

    def __del__(self):
        """Destructor."""

        self._fout_ambiguous_priority.close()

    def _parse_sp_priority_ledger(self, species_priority_ledger):
        """Parse manually resolved priority cases."""

        self.log.info('Parsing species priority ledger:')
        self.manual_sp_priority = defaultdict(lambda: {})
        self.manual_species_priority_year = {}

        num_cases = 0
        with open(species_priority_ledger, encoding='utf-8') as f:
            header = f.readline().strip().split('\t')

            spA_idx = header.index('NCBI species A')
            spB_idx = header.index('NCBI species B')
            sp_year_idx = header.index('Priority year')
            priority_sp_idx = header.index('Species with priority')

            for line in f:
                tokens = line.strip().split('\t')

                spA = tokens[spA_idx].strip()
                spB = tokens[spB_idx].strip()
                year = int(tokens[sp_year_idx].strip())
                priority_sp = tokens[priority_sp_idx].strip()

                assert spA.startswith('s__')
                assert spB.startswith('s__')
                assert priority_sp.startswith('s__')
                if priority_sp not in [spA, spB]:
                    self.log.error('Error in species priority ledger. Species {} cannot have priority for {} and {}.'.format(
                        priority_sp, spA, spB))
                    sys.exit(-1)

                self.manual_sp_priority[spA][spB] = priority_sp
                self.manual_sp_priority[spB][spA] = priority_sp

                self.manual_species_priority_year[spA] = year
                self.manual_species_priority_year[spB] = year
                num_cases += 1

        self.log.info(
            f' - identified {num_cases:,} manually resolved cases')

    def _parse_genus_priority_ledger(self, genus_priority_ledger):
        """Parse manually resolved priority cases."""

        self.log.info('Parsing genus priority ledger:')
        self.manual_genus_priority = defaultdict(lambda: {})

        num_cases = 0
        with open(genus_priority_ledger, encoding='utf-8') as f:
            header = [f.strip() for f in f.readline().strip().split('\t')]

            genusA_idx = header.index('Genus A')
            genusB_idx = header.index('Genus B')
            priority_idx = header.index('Priority')

            for line in f:
                tokens = [v.strip() for v in line.strip().split('\t')]

                genusA = tokens[genusA_idx]
                genusB = tokens[genusB_idx]
                priority_genus = tokens[priority_idx]

                assert genusA.startswith('g__')
                assert genusB.startswith('g__')
                assert priority_genus.startswith('g__')
                assert priority_genus in [genusA, genusB]

                self.manual_genus_priority[genusA][genusB] = priority_genus
                self.manual_genus_priority[genusB][genusA] = priority_genus
                num_cases += 1

        self.log.info(
            f' - identified {num_cases:,} manually resolved cases')

    def _parse_lpsn_gss(self, lpsn_gss_file):
        """Parse genus priority information and validity of names from LPSN GSS metadata file."""

        fout = open(os.path.join(self.output_dir,
                                 'lpsn_genus_priorities.tsv'), 'w')
        fout.write('Genus\tLPSN authors\tParse priority\n')

        genus_priority_count = 0
        self.lpsn_genus_priority = {}
        self.lpsn_valid_names = set()
        valid_genera_count = 0
        valid_species_count = 0
        valid_subsp_count = 0
        illegitimate_names = set()
        with open(lpsn_gss_file, encoding='utf-8', errors='ignore') as f:
            csv_reader = csv.reader(f)

            for line_num, tokens in enumerate(csv_reader):
                if line_num == 0:
                    genus_idx = tokens.index('genus_name')
                    species_idx = tokens.index('sp_epithet')
                    subsp_idx = tokens.index('subsp_epithet')
                    status_idx = tokens.index('status')
                    author_idx = tokens.index('authors')
                else:
                    # get taxon name
                    generic = tokens[genus_idx].strip().replace('"', '')
                    specific = tokens[species_idx].strip().replace('"', '')
                    subsp = tokens[subsp_idx].strip().replace('"', '')
                    if subsp:
                        taxon = 'sb__{} {} subsp. {}'.format(
                            generic, specific, subsp)
                        valid_subsp_count += 1
                    elif specific:
                        taxon = 's__{} {}'.format(generic, specific)
                        valid_species_count += 1
                    else:
                        taxon = 'g__{}'.format(generic)
                        valid_genera_count += 1

                    # parse status of taxon
                    status = tokens[status_idx].strip().replace('"', '')
                    status_tokens = [t.strip() for t in status.split(';')]
                    status_tokens = [tt.strip()
                                     for t in status_tokens for tt in t.split(',')]

                    if 'illegitimate name' in status_tokens:
                        illegitimate_names.add(taxon)
                        if taxon in self.lpsn_genus_priority:
                            continue
                    else:
                        self.lpsn_valid_names.add(taxon)

                    # get priority references, ignoring references if they are
                    # marked as being a revied name as indicated by a 'ex' or 'emend'
                    # (e.g. Holospora (ex Hafkine 1890) Gromov and Ossipov 1981)
                    ref_str = tokens[author_idx]
                    priority_year = LPSN.parse_lpsn_priority_year(ref_str)

                    # if 'gen. nov.' in status_tokens:
                    if taxon.startswith('g__'):
                        if (taxon in self.lpsn_genus_priority
                                and taxon not in illegitimate_names
                                and self.lpsn_genus_priority[taxon] != priority_year):
                            # conflict that can't be attributed to one of the entries being
                            # considered an illegitimate name
                            self.log.error('Conflicting genus priority for {}: {} {}'.format(
                                taxon, priority_year, self.lpsn_genus_priority[taxon]))

                        self.lpsn_genus_priority[taxon] = priority_year
                        genus_priority_count += 1
                        fout.write('{}\t{}\t{}\n'.format(taxon,
                                                         tokens[author_idx],
                                                         priority_year))

        self.log.info(
            f' - established genus priority for {genus_priority_count:,} genera using LPSN GSS metadata')
        self.log.info(' - identified {:,} genera, {:,} species and {:,} subpecies with validly published names in LPSN GSS metadata'.format(
            valid_genera_count,
            valid_species_count,
            valid_subsp_count))

    def genus_priority_year(self, genus):
        """Determine year of priority for genus."""

        if genus in self.lpsn_genus_priority:
            return self.lpsn_genus_priority[genus]

        return Genome.NO_PRIORITY_YEAR

    def genus_priority(self, genus1, genus2):
        """Resolve priority between two genera."""

        year1 = self.genus_priority_year(genus1)
        year2 = self.genus_priority_year(genus2)

        if genus1 in self.manual_genus_priority and genus2 in self.manual_genus_priority:
            return self.manual_genus_priority[genus1][genus2]
        elif year1 < year2:
            return genus1
        elif year2 < year1:
            return genus2

        return None

    def species_priority_year(self, cur_genomes, gid):
        """Get year of priority for genome."""

        ncbi_sp = cur_genomes[gid].ncbi_taxa.species
        if (cur_genomes[gid].is_effective_type_strain()
                and ncbi_sp in self.manual_species_priority_year):
            return self.manual_species_priority_year[ncbi_sp]

        return cur_genomes[gid].year_of_priority()

    def consensus_gtdb_genus(self, cur_genomes, gid1, gid2):
        """Get consensus GTDB genus taking into account a genome may not be classified."""

        g1 = cur_genomes[gid1].gtdb_taxa.genus
        g2 = cur_genomes[gid2].gtdb_taxa.genus
        if g1 == g2:
            return g1
        elif g1 == 'g__':
            return g2
        elif g2 == 'g__':
            return g1

        return 'g__'

    def species_priority(self, cur_genomes, gid1, gid2):
        """Resolve species priority of genomes."""

        ncbi_sp1 = cur_genomes[gid1].ncbi_taxa.species
        ncbi_sp2 = cur_genomes[gid2].ncbi_taxa.species

        ncbi_genus1 = cur_genomes[gid1].ncbi_taxa.genus
        ncbi_genus2 = cur_genomes[gid2].ncbi_taxa.genus

        # give priority to genomes assembled from the
        # type strain of the species
        if (cur_genomes[gid1].is_effective_type_strain()
                and not cur_genomes[gid2].is_effective_type_strain()):
            return gid1, 'effective type strain'
        if (not cur_genomes[gid1].is_effective_type_strain()
                and cur_genomes[gid2].is_effective_type_strain()):
            return gid2, 'effective type strain'

        # give priority to valid species names
        if (ncbi_sp1 in self.lpsn_valid_names
                and ncbi_sp2 not in self.lpsn_valid_names):
            return gid1, 'selected valid species name'
        if (ncbi_sp1 not in self.lpsn_valid_names
                and ncbi_sp2 in self.lpsn_valid_names):
            return gid2, 'selected valid species name'

        # give priority to type species of genus first
        if (cur_genomes[gid1].is_gtdb_type_species()
                and not cur_genomes[gid2].is_gtdb_type_species()):
            return gid1, 'type species of genus'
        if (not cur_genomes[gid1].is_gtdb_type_species()
                and cur_genomes[gid2].is_gtdb_type_species()):
            return gid2, 'type species of genus'

        # give priority to NCBI species with generic name matching
        # GTDB genus if only one NCBI species matches GTDB genus
        consensus_gtdb_genus = self.consensus_gtdb_genus(
            cur_genomes, gid1, gid2)
        if (ncbi_genus1 == consensus_gtdb_genus
                and ncbi_genus2 != consensus_gtdb_genus):
            return gid1, 'match to expected GTDB genus assignment'
        if (ncbi_genus1 != consensus_gtdb_genus
                and ncbi_genus2 == consensus_gtdb_genus):
            return gid2, 'match to expected GTDB genus assignment'

        # give priority based on date of valid publication
        if self.species_priority_year(cur_genomes, gid1) < self.species_priority_year(cur_genomes, gid2):
            return gid1, 'year of priority'
        if self.species_priority_year(cur_genomes, gid1) > self.species_priority_year(cur_genomes, gid2):
            return gid2, 'year of priority'

        # either neither genome is type material, or both genomes are type material with the same
        # priority year (type material could also be the type strain of a subspecies)

        # base priority on genome quality if genomes are not
        # assembled from type strain of species
        elif (not cur_genomes[gid1].is_effective_type_strain()
              and not cur_genomes[gid2].is_effective_type_strain()):
            # neither genome is effective type strain
            hq_gid = select_highest_quality([gid1, gid2], cur_genomes)
            return hq_gid, 'highest quality genome'
        elif (cur_genomes[gid1].is_exclusively_effective_type_strain()
              and cur_genomes[gid2].is_exclusively_effective_type_strain()):
            # both genomes are only effective, but not validated type material
            hq_gid = select_highest_quality([gid1, gid2], cur_genomes)
            return hq_gid, 'highest quality genome'

        # both genomes are at least the effective type strain of the species
        assert cur_genomes[gid1].is_effective_type_strain(
        ) and cur_genomes[gid2].is_effective_type_strain()
        if ncbi_sp1 == ncbi_sp2:
            # genomes are from the same NCBI species so have the same naming priority,
            # so resolve ordering using genome quality
            hq_gid = select_highest_quality([gid1, gid2], cur_genomes)
            return hq_gid, 'highest quality genome'

        # both genomes are type strain of species, but priority
        # is ambiguous using just publication date so species need
        # to be in the priority ledger
        if ncbi_sp1 not in self.manual_sp_priority or ncbi_sp2 not in self.manual_sp_priority[ncbi_sp1]:
            self.log.error('Ambiguous priority based on publication date.')
            self.log.error('Species need to be manually resolved in priority ledger: {}: {} / {} / {}, {}: {} / {} / {}'.format(
                gid1, ncbi_sp1, cur_genomes[gid1].is_gtdb_type_strain(
                ), cur_genomes[gid1].is_effective_type_strain(),
                gid2, ncbi_sp2, cur_genomes[gid2].is_gtdb_type_strain(), cur_genomes[gid2].is_effective_type_strain()))

            self._fout_ambiguous_priority.write('{}\t{}\t{}\n'.format(
                ncbi_sp1, ncbi_sp2,
                self.species_priority_year(cur_genomes, gid1)))

            # arbitrarily resolve naming priority based on genome quality
            hq_gid = select_highest_quality([gid1, gid2], cur_genomes)
            return hq_gid, 'error: ambiguous priority date'  # ***

        priority_sp = self.manual_sp_priority[ncbi_sp1][ncbi_sp2]
        if ncbi_sp1 == priority_sp:
            return gid1, 'manual curation'

        return gid2, 'manual curation'

    def test_species_priority(self, cur_genomes, gid1, gid2):
        """Check if species of genome 1 has priority over species of genome 2."""

        if (cur_genomes[gid1].is_gtdb_type_strain()
                and not cur_genomes[gid2].is_gtdb_type_strain()):
            return True
        elif (not cur_genomes[gid1].is_gtdb_type_strain()
              and cur_genomes[gid2].is_gtdb_type_strain()):
            return False

        priority_gid, _note = self.species_priority(cur_genomes, gid1, gid2)

        if priority_gid == gid1:
            return True

        return False

    def _sp_priority_sort(self, gid1, gid2):
        """Sort genomes by priority."""

        priority_gid, _note = self.species_priority(
            self.cur_genomes, gid1, gid2)
        if priority_gid == gid1:
            return -1

        return 1

    def sort_by_sp_priority(self, cur_genomes, gids, reverse=False):
        """Sort set of genomes by priority with highest priority genome first."""

        self.cur_genomes = cur_genomes
        sorted_gids = sorted(gids,
                             key=functools.cmp_to_key(self._sp_priority_sort),
                             reverse=reverse)

        return sorted_gids
