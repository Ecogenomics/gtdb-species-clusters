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
from gtdb_species_clusters.genome_utils import select_highest_quality, parse_seqcode_table
from gtdb_species_clusters.lpsn import LPSN
from gtdb_species_clusters.taxon_utils import generic_name, canonical_species, specific_epithet, test_same_epithet


class SpeciesPriorityManager(object):
    """Resolve naming priority of species names."""

    def __init__(self,
                 species_priority_ledger,
                 genus_priority_ledger,
                 lpsn_gss_file,
                 output_dir,
                 seqcode_file=None):
        """Initialization.

        Species names are the species genomes are type material of, which are their
        reference (SeqCode or NCBI) species unless they are in the GTDB type strains
        ledger (see priority_names). Names validly published under the SeqCode, read
        from the SeqCode table when given, are valid names alongside those
        validly published under the ICNP according to LPSN. The year of priority
        of a genus without an entry in LPSN is the SeqCode year of priority of
        its type species, when the SeqCode table gives one.
        """

        self.log = logging.getLogger('rich')

        self.output_dir = output_dir

        self._fout_ambiguous_priority = open(os.path.join(
            output_dir, 'ambiguous_sp_priority.tsv'), 'w')
        self._fout_ambiguous_priority.write(
            'NCBI species A\tNCBI species B\tPriority year\n')
        self._reported_ambiguous_pairs = set()

        self._parse_sp_priority_ledger(species_priority_ledger)
        self._parse_genus_priority_ledger(genus_priority_ledger)
        self._parse_lpsn_gss(lpsn_gss_file)

        self.valid_names = set(self.lpsn_valid_names)
        self.seqcode_genus_priority = {}
        if seqcode_file:
            seqcode_records = parse_seqcode_table(seqcode_file).values()
            seqcode_species = {'s__' + r.name for r in seqcode_records}
            self.valid_names.update(seqcode_species)
            self.log.info(f' - identified {len(seqcode_species):,} species with validly published names in SeqCode')

            self.seqcode_genus_priority = self._seqcode_genus_priority(seqcode_records)
            self.log.info(' - established genus priority for {:,} genera using SeqCode type species'.format(
                len(self.seqcode_genus_priority)))

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
                if priority_genus not in [genusA, genusB]:
                    self.log.error(f"Priority genus {priority_genus} no in set of genera: {genusA} {genusB}")
                    sys.exit(1)

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

    @staticmethod
    def _seqcode_genus_priority(seqcode_records):
        """Year of priority of each genus from the SeqCode year of priority of its type species."""

        genus_priority = {}
        for r in seqcode_records:
            if not r.type_species_of_genus or r.priority_year == Genome.NO_PRIORITY_YEAR:
                continue

            genus = 'g__' + r.name.split()[0]
            genus_priority[genus] = min(r.priority_year,
                                        genus_priority.get(genus, Genome.NO_PRIORITY_YEAR))

        return genus_priority

    def genus_priority_year(self, genus):
        """Determine year of priority for genus.

        The LPSN year is used when LPSN has an entry for the genus, else
        the SeqCode year of priority of its type species.
        """

        if genus in self.lpsn_genus_priority:
            return self.lpsn_genus_priority[genus]

        if genus in self.seqcode_genus_priority:
            return self.seqcode_genus_priority[genus]

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

    @staticmethod
    def priority_names(genome):
        """Get species and genus used to establish the naming priority of a genome.

        Genomes are identified by the species they are type material of (see Genome.type_species), which is
        their reference (SeqCode or NCBI) species unless they are in the GTDB type strains ledger. The genus of a
        ledger species keeps any GTDB suffix so it can match the GTDB genus (e.g., g__Ruminococcus_D).
        """

        species = genome.type_species()
        if species == genome.ref_taxa.species:
            return species, genome.ref_taxa.genus

        return species, 'g__' + generic_name(species)

    def has_valid_name(self, genome):
        """Check if the species a genome is type material of has a validly published name (ICNP or SeqCode).

        A GTDB type strains ledger species is looked up without suffixes (e.g., s__Ruminococcus_D bicirculans
        as s__Ruminococcus bicirculans). A ledger species that places the validly published, trusted reference
        species of the genome in another genus is a new combination of that name, so is also considered valid
        (e.g., s__Globicatella yangjingiae for s__Tuanshanicoccus yangjingiae); the specific epithet may differ
        in gender (see test_same_epithet).
        """

        species = genome.type_species()
        ref_species = genome.ref_taxa.species
        if species == ref_species:
            return species in self.valid_names

        species = canonical_species(species)
        if species in self.valid_names:
            return True

        return (genome.is_ref_species_trusted()
                and ref_species in self.valid_names
                and test_same_epithet(specific_epithet(species), specific_epithet(ref_species)))

    def species_priority_year(self, cur_genomes, gid):
        """Get year of priority for genome."""

        species, _genus = self.priority_names(cur_genomes[gid])
        if (cur_genomes[gid].is_effective_type_strain()
                and species in self.manual_species_priority_year):
            return self.manual_species_priority_year[species]

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

        # species and genera are the names of the species the genomes are type material of,
        # which are their reference (SeqCode or NCBI) names unless they are in the GTDB type
        # strains ledger (see priority_names)
        ncbi_sp1, ncbi_genus1 = self.priority_names(cur_genomes[gid1])
        ncbi_sp2, ncbi_genus2 = self.priority_names(cur_genomes[gid2])
        valid_name1 = self.has_valid_name(cur_genomes[gid1])
        valid_name2 = self.has_valid_name(cur_genomes[gid2])

        # give priority to genomes assembled from the
        # type strain of the species
        if (cur_genomes[gid1].is_effective_type_strain()
                and not cur_genomes[gid2].is_effective_type_strain()):
            return gid1, 'effective type strain'
        if (not cur_genomes[gid1].is_effective_type_strain()
                and cur_genomes[gid2].is_effective_type_strain()):
            return gid2, 'effective type strain'

        # give priority to valid species names (validly published under the ICNP or SeqCode)
        if valid_name1 and not valid_name2:
            return gid1, 'selected valid species name'
        if not valid_name1 and valid_name2:
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

            # report each species pair only once, regardless of the order
            # in which the species are compared
            sp_pair = tuple(sorted([ncbi_sp1, ncbi_sp2]))
            if sp_pair not in self._reported_ambiguous_pairs:
                self._reported_ambiguous_pairs.add(sp_pair)
                self._fout_ambiguous_priority.write('{}\t{}\t{}\n'.format(
                    sp_pair[0], sp_pair[1],
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
