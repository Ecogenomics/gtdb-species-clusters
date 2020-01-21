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
import re
from collections import defaultdict
from collections import Counter

from biolib.taxonomy import Taxonomy

from gtdb_species_clusters.taxon_utils import canonical_taxon


class TaxonSuffixManager():
    """Manage last polyphyletic suffix used for each taxon."""

    def __init__(self):
        """Initialize."""
        
        self.logger = logging.getLogger('timestamp')
        self.prev_taxonomy_dir = '/srv/projects/gtdb/data/taxonomy_gtdb'

        # get all previous taxonomy files
        self.logger.info('Reading previous GTDB taxonomy files in {}:'.format(self.prev_taxonomy_dir))
        taxonomies = defaultdict(lambda: {})
        for f in os.listdir(self.prev_taxonomy_dir):
            if f.endswith('.tsv') and 'gtdb' in f:
                self.logger.info('  %s' % f)
                taxonomy_file = os.path.join(self.prev_taxonomy_dir, f)
                
                taxonomy_id ='_'.join(f.split('_')[0:2])
                taxonomies[taxonomy_id].update(Taxonomy().read(taxonomy_file))
            
        self.logger.info('Considering taxonomy from {:,} previous releases.'.format(len(taxonomies)))

        # get highest alphabetic suffix for each taxon
        self.logger.info('Determining highest polyphyletic alphabetic suffix for each taxon.')

        self.taxon_suffix = {}
        for taxonomy in taxonomies.values():
            for taxa in taxonomy.values():
                for taxon in taxa:
                    rank_prefix = taxon[0:3]
                    taxon_name = taxon[3:]

                    if '_' in taxon_name:
                        if rank_prefix != 's__':
                            taxon_name, suffix = taxon_name.rsplit('_', 1)
                        else:
                            # check if the specific name has a suffix
                            generic_name, specific_name = taxon_name.split()
                            if '_' in specific_name:
                                canonical_specific_name, suffix = specific_name.rsplit('_', 1)
                                taxon_name = '{} {}'.format(generic_name, canonical_specific_name)
                            else:
                                continue

                        canonical_taxon = '{}{}'.format(rank_prefix, taxon_name)
                        cur_suffix = self.taxon_suffix.get(canonical_taxon, 'A')
                        if self._suffix_value(suffix) >= self._suffix_value(cur_suffix):
                            self.taxon_suffix[canonical_taxon] = suffix

    def _suffix_value(self, suffix):
        """Get value of taxa suffix.
        
        A = 65, B = 66, ..., Z = 90
        AA = 16705, AB = 16705, BA = 16961
        """
        
        v = 0
        for idx, ch in enumerate(suffix):
            v += 256**(len(suffix)-idx-1) * ord(ch)
        
        return v
    
    def _increment_suffix(self, suffix):
        """Increment suffix by one character."""
        
        if len(suffix) > 2:
            print('[Error] Unable to handle suffixes with >2 characters.')
            sys.exit(-1)
        
        if len(suffix) == 1:
            if suffix != 'Z':
                suffix = chr(ord(suffix) + 1)
            else:
                suffix = 'AA'
        else:
            last_ch = suffix[1]
            if last_ch != 'Z':
                suffix = suffix[0] + chr(ord(last_ch) + 1)
            else:
                first_ch = suffix[0]
                first_ch = chr(ord(first_ch) + 1)
                suffix = first_ch + 'A'

        return suffix
        
    def next_suffix(self, taxon):
        """Get next suffix for taxon."""
        
        ct = canonical_taxon(taxon)
        if ct in self.taxon_suffix:
            cur_suffix = self.taxon_suffix[ct]
            next_suffix = self._increment_suffix(cur_suffix)
        else:
            next_suffix = 'A'
            
        self.taxon_suffix[ct] = next_suffix
        
        return next_suffix
