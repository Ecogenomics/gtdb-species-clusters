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

import sys
import logging
from collections import defaultdict

from gtdb_species_clusters.taxon_utils import canonical_taxon, is_placeholder_sp_epithet
from gtdb_species_clusters.taxon_suffix_manager import TaxonSuffixManager


class SpeciesNameManager():
    """Manage assignment and updating of GTDB species names.

    This manager handle the following:
     - tracking valid suffixes
     - determining presence of species name in GTDB
     - determining presence of species epithet in a given genus
     - tracking if a species name has been assigned multiple times
    """

    def __init__(self, prev_genomes, cur_genomes):
        """Initialization."""

        self.log = logging.getLogger('timestamp')

        self.taxon_suffix_manager = TaxonSuffixManager()

        self.prev_genomes = prev_genomes
        self.cur_genomes = cur_genomes

        # get species assignment of previous GTDB species clusters
        self.gtdb_sp_epithets = defaultdict(set)
        self.gtdb_canonical_sp_epithets = defaultdict(set)
        self.sp_epithets_rid = defaultdict(lambda: {})
        for rid, _sp in prev_genomes.sp_clusters.species():
            gtdb_genus = prev_genomes[rid].gtdb_taxa.genus
            gtdb_sp_epithet = prev_genomes[rid].gtdb_taxa.specific_epithet
            self.gtdb_sp_epithets[gtdb_genus].add(gtdb_sp_epithet)
            self.gtdb_canonical_sp_epithets[gtdb_genus].add(
                canonical_taxon(gtdb_sp_epithet))
            self.sp_epithets_rid[gtdb_genus][gtdb_sp_epithet] = rid

    def numeric_placeholder_sp_epithet(self, rid):
        """Generate placeholder species epithet from representative accession."""

        if rid.startswith('G'):
            sp_id = rid.replace('G', '')
        else:
            print('Unrecognized genome ID: {}'.format(rid))
            sys.exit(-1)

        return 'sp{}'.format(sp_id)

    def suffixed_placeholder_sp_epithet(self, generic, specific):
        """Determine suffix for species."""

        if is_placeholder_sp_epithet(specific):
            self.log.error(
                f'Only Latin specific names should have a suffix: generic = {generic}; specific = {specific}')
            raise Exception('Invalid species epithet')

        sp = 's__{} {}'.format(generic, specific)
        next_suffix = self.taxon_suffix_manager.next_suffix(sp)

        return '{}_{}'.format(specific, next_suffix)

    def add_suffixed_species(self, species):
        """Account for species names when generating further placeholder names.

        This is required as species can be transferred between genera and will
        retain their existing suffix. As such, we must track that this genera
        now has a species name with a given suffix.

        (e.g., Lactobacillus_G kunkeei_A is transferred to Apilactobacillus
               as A. kunkeei_A, so the next suffixed A. kunkeei representative
               must have a 'B' suffix.)
        """

        self.taxon_suffix_manager.add_suffixed_species(species)
