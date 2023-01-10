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

import logging
from collections import defaultdict, Counter

from biolib.taxonomy import Taxonomy

from gtdb_species_clusters.taxon_utils import (canonical_taxon,
                                               test_same_epithet,
                                               generic_name,
                                               specific_epithet,
                                               is_placeholder_taxon,
                                               longest_common_suffix)


class SpecificEpithetManager():
    """Manage specific epithets, including changing suffixes due to transfers to genera with a different gender."""

    def __init__(self):
        """Initialize."""

        self.log = logging.getLogger('rich')

        self.curated_epithet_changes = {}
        self.previously_curated = {}
        self.sp_epithet_map = defaultdict(lambda: {})
        self.gtdb_ncbi_generic_map = defaultdict(lambda: defaultdict(list))

    def parse_specific_epithet_ledger(self, specific_epithet_ledger):
        """Parse manually curated changes to specific names resulting from genus transfers."""

        self.log.info('Reading manually-curate specific epithet changes.')
        mc_count = 0
        with open(specific_epithet_ledger) as f:
            header = f.readline().strip().split('\t')

            generic_idx = header.index('GTDB generic')
            sp_corr_idx = header.index('GTDB specific_corrected')
            sp_idx = header.index('GTDB specific')
            ncbi_generic_idx = header.index('NCBI generic')
            ncbi_sp_idx = header.index('NCBI specific')

            for line in f:
                tokens = [t.strip() for t in line.strip().split('\t')]

                gtdb_generic = tokens[generic_idx]
                specific_corr = tokens[sp_corr_idx]
                gtdb_specific = tokens[sp_idx]
                ncbi_generic = tokens[ncbi_generic_idx]
                ncbi_specific = tokens[ncbi_sp_idx]

                self.previously_curated[f's__{gtdb_generic} {gtdb_specific}'] = f's__{ncbi_generic} {ncbi_specific}'

                if specific_corr:
                    self.previously_curated[f's__{gtdb_generic} {specific_corr}'] = f's__{ncbi_generic} {ncbi_specific}'

                    species_orig = 's__{} {}'.format(
                        gtdb_generic, gtdb_specific)
                    species_corr = 's__{} {}'.format(
                        gtdb_generic, specific_corr)
                    self.curated_epithet_changes[species_orig] = species_corr
                    mc_count += 1
                elif gtdb_specific != ncbi_specific:
                    species_orig = 's__{} {}'.format(
                        gtdb_generic, ncbi_specific)
                    species_corr = 's__{} {}'.format(
                        gtdb_generic, gtdb_specific)
                    self.curated_epithet_changes[species_orig] = species_corr
                    mc_count += 1

        self.log.info(
            ' - identified {:,} specified cases.'.format(mc_count))

    def translate_species(self, gtdb_species):
        """Translate GTDB species name to name with corrected specific epithet."""

        return self.curated_epithet_changes.get(gtdb_species, gtdb_species)

    def infer_epithet_map(self,
                          cur_taxonomy,
                          mc_species,
                          cur_genomes):
        """Infer mapping of NCBI epithet to GTDB epithet which may be different due to gender of genus."""

        # **************************************
        # This should be updated so it only includes valid transfers, and not
        # results due to misclassifications at NCBI. For example, right now this
        # code reports Enterobacter cancerogenus being transferred to Pantoea, but
        # really this is just a misclassified NCBI genome.

        # get species in GTDB genus
        generic_rids = defaultdict(list)
        for rid in cur_taxonomy:
            gtdb_generic = cur_taxonomy[rid][Taxonomy.GENUS_INDEX].replace(
                'g__', '')
            if rid in mc_species:
                gtdb_generic = generic_name(mc_species[rid])

            generic_rids[gtdb_generic].append(rid)

        # establish epithets that are nearly identical
        # except for small change to suffix which is
        # assumed to be due to a gender change
        for gtdb_generic, rids in generic_rids.items():
            ncbi_sp_epithet_list = defaultdict(list)
            for rid in rids:
                ncbi_species = cur_genomes[rid].ncbi_taxa.species
                if ncbi_species == 's__':
                    continue

                ncbi_generic = generic_name(ncbi_species)
                ncbi_specific = specific_epithet(ncbi_species)

                if rid in mc_species:
                    gtdb_species = mc_species[rid]
                else:
                    gtdb_species = cur_genomes[rid].gtdb_taxa.species

                gtdb_specific = canonical_taxon(specific_epithet(gtdb_species))

                self.gtdb_ncbi_generic_map[gtdb_generic][gtdb_specific].append(
                    ncbi_generic)

                if test_same_epithet(ncbi_specific, gtdb_specific):
                    ncbi_sp_epithet_list[ncbi_specific].append(gtdb_specific)

            for ncbi_specific, gtdb_specific_list in ncbi_sp_epithet_list.items():
                gtdb_specific_counter = Counter(gtdb_specific_list)

                top_gtdb_specific, count = gtdb_specific_counter.most_common(1)[
                    0]

                map_perc = count*100.0/len(gtdb_specific_list)
                if map_perc >= 50:
                    self.sp_epithet_map[gtdb_generic][ncbi_specific] = top_gtdb_specific

                if map_perc != 100:
                    self.log.warning('Imperfect suffix mapping between {} {} to {} at {:.1f}%.'.format(
                        gtdb_generic,
                        top_gtdb_specific,
                        ncbi_specific,
                        count*100.0/len(gtdb_specific_list)))

    def write_diff_epithet_map(self, output_file):
        """Write out epithet map for specific names that differ between GTDB and NCBI."""

        fout = open(output_file, 'w')

        fout.write('GTDB generic\tGTDB specific\tNCBI generic\tNCBI specific\n')

        for gtdb_generic in self.sp_epithet_map:
            for ncbi_specific, gtdb_specific in self.sp_epithet_map[gtdb_generic].items():
                if gtdb_specific != ncbi_specific:
                    ncbi_generic_list = self.gtdb_ncbi_generic_map[gtdb_generic][gtdb_specific]

                    ncbi_generic_counter = Counter(ncbi_generic_list)

                    top_ncbi_generic, count = ncbi_generic_counter.most_common(1)[
                        0]
                    map_perc = count*100.0/len(ncbi_generic_list)

                    if map_perc == 100:
                        fout.write('{}\t{}\t{}\t{}\n'.format(
                            gtdb_generic,
                            gtdb_specific,
                            top_ncbi_generic,
                            ncbi_specific))
                    else:
                        self.log.warning('Imperfect GTDB to NCBI genus mapping for {} {} -> {}'.format(
                            gtdb_generic,
                            gtdb_specific,
                            ncbi_generic_counter))

        fout.close()

    def write_epithet_map(self, output_file, filtered_previously_checked=False):
        """Write out epithet map."""

        fout = open(output_file, 'w')

        fout.write('GTDB generic\tGTDB specific\tNCBI generic\tNCBI specific\n')

        for gtdb_generic in self.sp_epithet_map:
            if is_placeholder_taxon('g__' + gtdb_generic):
                continue

            for ncbi_specific, gtdb_specific in self.sp_epithet_map[gtdb_generic].items():
                ncbi_generic_list = self.gtdb_ncbi_generic_map[gtdb_generic][gtdb_specific]

                ncbi_generic_counter = Counter(ncbi_generic_list)
                top_ncbi_generic, _count = ncbi_generic_counter.most_common(1)[
                    0]

                if top_ncbi_generic != gtdb_generic:
                    if filtered_previously_checked:
                        gtdb_sp = f's__{gtdb_generic} {gtdb_specific}'
                        ncbi_sp = f's__{top_ncbi_generic} {ncbi_specific}'
                        if gtdb_sp in self.previously_curated and ncbi_sp in self.previously_curated[gtdb_sp]:
                            continue

                        # also skip the following as recommended by Masha:
                        # genera end in -ium, -um, -ella, or -iella
                        suffixes = ('ium', 'um', 'ella', 'iella')
                        if gtdb_generic.endswith(suffixes) and top_ncbi_generic.endswith(suffixes):
                            continue

                        # genera that end in the same suffix can also be skipped,
                        # but there isn't an easy way to establish the suffix so
                        # here we just skip cases where the last 56+ characters are
                        # the same to catch cases like 'bacter', 'vibrio', 'plasma',
                        # 'spora', 'monas', etc.
                        lcs = longest_common_suffix(
                            gtdb_generic, top_ncbi_generic)
                        if len(lcs) >= 5:
                            continue

                    fout.write('{}\t{}\t{}\t{}\n'.format(
                        gtdb_generic,
                        gtdb_specific,
                        top_ncbi_generic,
                        ncbi_specific))

        fout.close()
