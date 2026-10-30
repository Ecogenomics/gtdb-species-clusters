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

from functools import lru_cache

from gtdblib.taxonomy.taxonomy import Taxonomy


RANK_INDEX = {rank_prefix: rank_idx for rank_idx, rank_prefix in enumerate(Taxonomy.RANK_PREFIXES)}


@lru_cache(maxsize=None)
def _parse_taxa_str(taxa_str, filtered):
    """Convert taxonomy string to taxa and standard ranks.

    Many genomes have the same taxonomy string, so results are cached. Immutable
    tuples are returned as each Taxa object needs its own copy it can modify.
    """

    if taxa_str and taxa_str != 'none':
        if filtered:
            taxa_str = taxa_str.replace('Candidatus ', '')
            taxa_str = taxa_str.replace('[', '')
            taxa_str = taxa_str.replace(']', '')
        taxa = tuple(t.strip() for t in taxa_str.split(';'))
    else:
        taxa = tuple(Taxonomy.RANK_PREFIXES)

    standard_taxa = []
    for taxon in taxa:
        rank_idx = RANK_INDEX.get(taxon[0:3])
        if rank_idx is not None:
            standard_taxa.append((rank_idx, taxon))

    return taxa, tuple(standard_taxa)


class Taxa(object):
    """Taxa for organism."""

    def __init__(self, taxa_str, filtered=True):
        """Initialization."""

        self._parse_taxa(taxa_str, filtered)

    def _parse_taxa(self, taxa_str, filtered):
        """Convert taxonomy string to taxa list."""

        taxa, standard_taxa = _parse_taxa_str(taxa_str, filtered)
        self.taxa = list(taxa)
        self.standard_taxa = dict(standard_taxa)

    def __str__(self):
        """User-friendly string representation."""

        return '; '.join(self.taxa)

    def __iter__(self):
        """Iterate over taxa."""

        for taxon in self.taxa:
            yield taxon

    def __getitem__(self, idx):
        """Get genome."""

        return self.taxa[idx]

    def __contains__(self, taxon):
        """Check if genome is in genome set."""

        return taxon in self.taxa

    def __len__(self):
        """Number of taxa."""

        return len(self.taxa)

    def __eq__(self, other):
        """Override the default Equals behavior"""
        return self.taxa == other.taxa

    def __ne__(self, other):
        """Override the default Unequal behavior"""
        return self.taxa != other.taxa

    def get_taxa(self, rank_idx):
        """Get taxon at specified rank."""

        return self.standard_taxa[rank_idx]

    def set_taxa(self, rank_idx, new_taxon):
        """Set taxon for specified rank."""

        self.standard_taxa[rank_idx] = new_taxon

        for idx, taxon in enumerate(self.taxa):
            if taxon[0:3] == Taxonomy.RANK_PREFIXES[rank_idx]:
                self.taxa[idx] = new_taxon
                break

    def update_taxa(self, taxa):
        """Update taxa."""

        self.taxa = taxa.taxa.copy()
        self.standard_taxa = taxa.standard_taxa.copy()

    @property
    def subspecies(self):
        """Get subspecies classification."""

        for taxon in self.taxa:
            if taxon.startswith('sb__'):
                return taxon

        return None

    @property
    def species(self):
        return self.standard_taxa[6]

    @species.setter
    def species(self, taxon):
        self.set_taxa(6, taxon)

    @property
    def specific_epithet(self):
        """Get specific epithet."""

        sp = self.species
        if sp == 's__':
            return ''

        _generic, specific = sp.split()

        return specific

    @specific_epithet.setter
    def specific_epithet(self, specific):
        """Set specific epithet."""

        sp = self.species
        generic, _specific = sp.split()

        new_sp = f'{generic} {specific}'

        self.species = new_sp

    @property
    def canonical_specific_epithet(self):
        """Get canonical version of specific_ pithet."""

        if '_' in self.specific_epithet:
            canonicalized = self.specific_epithet[0:self.specific_epithet.rfind(
                '_')]

        return canonicalized

    @property
    def genus(self):
        return self.standard_taxa[5]

    @genus.setter
    def genus(self, taxon):
        self.set_taxa(5, taxon)

    @property
    def family(self):
        return self.standard_taxa[4]

    @family.setter
    def family(self, taxon):
        self.set_taxa(4, taxon)

    @property
    def order(self):
        return self.standard_taxa[3]

    @order.setter
    def order(self, taxon):
        self.set_taxa(3, taxon)

    @property
    def class_taxon(self):
        return self.standard_taxa[2]

    @class_taxon.setter
    def class_taxon(self, taxon):
        self.set_taxa(2, taxon)

    @property
    def phylum(self):
        return self.standard_taxa[1]

    @phylum.setter
    def phylum(self, taxon):
        self.set_taxa(1, taxon)

    @property
    def domain(self):
        return self.standard_taxa[0]

    @domain.setter
    def domain(self, taxon):
        self.set_taxa(0, taxon)
