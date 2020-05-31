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
from collections import defaultdict, Counter

from biolib.taxonomy import Taxonomy

from gtdb_species_clusters.taxon_utils import (canonical_taxon,
                                                test_same_epithet,
                                                generic_name,
                                                specific_epithet,
                                                is_placeholder_taxon)


class SpecificEpithetManager():
    """Manage specific epithets, including changing suffixes due to transfers to genera with a different gender."""

    def __init__(self):
        """Initialize."""
        
        self.logger = logging.getLogger('timestamp')
        
        self.sp_epithet_map = defaultdict(lambda: {})
        self.gtdb_ncbi_generic_map = defaultdict(lambda: defaultdict(list))
        
    def infer_epithet_map(self, gids_of_interest, mc_species, cur_genomes, cur_clusters):
        """Infer mapping of NCBI epithet to GTDB epithet which may be different due to gender of genus."""
        
        # get species in GTDB genus
        generic_rids = defaultdict(list)
        for rid in cur_clusters:
            if rid not in gids_of_interest:
                continue
                
            gtdb_generic = cur_genomes[rid].gtdb_taxa.genus.replace('g__', '')
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

                self.gtdb_ncbi_generic_map[gtdb_generic][gtdb_specific].append(ncbi_generic)

                if test_same_epithet(ncbi_specific, gtdb_specific):
                    ncbi_sp_epithet_list[ncbi_specific].append(gtdb_specific)
                    
            for ncbi_specific, gtdb_specific_list in ncbi_sp_epithet_list.items():
                gtdb_specific_counter = Counter(gtdb_specific_list)
                
                top_gtdb_specific, count = gtdb_specific_counter.most_common(1)[0]

                map_perc = count*100.0/len(gtdb_specific_list)
                if map_perc >= 50:
                    self.sp_epithet_map[gtdb_generic][ncbi_specific] = top_gtdb_specific
                    
                if map_perc != 100:
                    self.logger.warning('Imperfect suffix mapping between from {} {} to {} at {:.1f}%.'.format(
                                            gtdb_generic, 
                                            top_gtdb_specific, 
                                            ncbi_specific,
                                            count*100.0/len(gtdb_specific_list)))
                                            
    def translate_epithet(self, gtdb_generic, ncbi_specific):
        """Translate NCBI specific name to GTDB equivalence."""
        
        if gtdb_generic in self.sp_epithet_map:
            return self.sp_epithet_map[gtdb_generic].get(ncbi_specific, ncbi_specific)
            
        return ncbi_specific
    
    def write_diff_epithet_map(self, output_file):
        """Write out epithet map for specific names that differ between GTDB and NCBI."""
        
        fout = open(output_file, 'w')
        
        fout.write('GTDB generic\tGTDB specific\tNCBI generic\tNCBI specific\n')
        
        for gtdb_generic in self.sp_epithet_map:
            for ncbi_specific, gtdb_specific in self.sp_epithet_map[gtdb_generic].items():
                if gtdb_specific != ncbi_specific:
                    ncbi_generic_list = self.gtdb_ncbi_generic_map[gtdb_generic][gtdb_specific]
                    ncbi_generic_counter = Counter(ncbi_generic_list)
                
                    top_ncbi_generic, count = ncbi_generic_counter.most_common(1)[0]
                    map_perc = count*100.0/len(ncbi_generic_list)
                
                    if map_perc == 100:
                        fout.write('{}\t{}\t{}\t{}\n'.format(
                                    gtdb_generic, 
                                    gtdb_specific,
                                    top_ncbi_generic,
                                    ncbi_specific))
                    else:
                        self.logger.warning('Imperfect GTDB to NCBI genus mapping for {} {} -> {}'.format(
                                                gtdb_generic, 
                                                gtdb_specific, 
                                                ncbi_generic_counter))
        
        fout.close()
        
        
    def write_epithet_map(self, output_file):
        """Write out epithet map."""
        
        fout = open(output_file, 'w')
        
        fout.write('GTDB generic\tGTDB specific\tNCBI generic\tNCBI specific\n')
        
        for gtdb_generic in self.sp_epithet_map:
            if is_placeholder_taxon('g__' + gtdb_generic):
                continue
                
            for ncbi_specific, gtdb_specific in self.sp_epithet_map[gtdb_generic].items():
                ncbi_generic_list = self.gtdb_ncbi_generic_map[gtdb_generic][gtdb_specific]
                ncbi_generic_counter = Counter(ncbi_generic_list)
                top_ncbi_generic, count = ncbi_generic_counter.most_common(1)[0]
                
                if top_ncbi_generic != gtdb_generic:
                    fout.write('{}\t{}\t{}\t{}\n'.format(
                                gtdb_generic, 
                                gtdb_specific,
                                top_ncbi_generic,
                                ncbi_specific))
        
        fout.close()
