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
import csv
import re
import datetime
import logging
from collections import defaultdict
from dataclasses import dataclass

from gtdb_species_clusters.genome import Genome


@dataclass
class LPSN_Taxon(object):
    """Metadata for single LPSN taxon."""
    
    # Ideally this would be a full record with the
    # majority of data at LPSN. However, this information
    # is not currently provided to the public except
    # on the web pages. We have requested that a data file
    # similar to the GSS file be provided.
    
    name: str
    priority_year: int
    type_material: str
    type_material_priority_year : int


class LPSN(object):
    """LPSN metadata for taxa."""
    
    def parse_lpsn_priority_year(reference_str):
        """Parse year of priority from LPSN reference string."""

        references = reference_str.replace('(', '').replace(')', '')
        years = re.sub(r'emend\.[^\d]*\d{4}', '', references)
        years = re.sub(r'ex [^\d]*\d{4}', ' ', years)
        years = re.findall('[1-3][0-9]{3}', years, re.DOTALL)
        years = [int(y) for y in years if int(y) <= datetime.datetime.now().year]
        
        if len(years) == 0:
            # assume this name is validated through ICN and just take the first 
            # date given as the year of priority
            years = re.findall('[1-3][0-9]{3}', references, re.DOTALL)
            years = [int(y) for y in years if int(y) <= datetime.datetime.now().year]
        
        return years[0]

    def __init__(self, lpsn_type_material_file, lpsn_gss_metadata_file):
        """Initialization."""
        
        self.logger = logging.getLogger('timestamp')
        
        if lpsn_type_material_file:
            self.taxa = self.parse_lpsn_type_material_file(lpsn_type_material_file)
            
        if lpsn_gss_metadata_file:
            self.sp_correct_names, self.sp_synonyms = self.parse_lpsn_gss_file(lpsn_gss_metadata_file)

    def type_material(self, taxon):
        """Get type material of taxon."""
        
        if taxon in self.taxa:
            return self.taxa[taxon].type_material
            
        return None
        
    def parse_lpsn_type_material_file(self, lpsn_type_material_file):
        """Parse type material for taxa as defined at LPSN."""
        
        taxa = {}
        with open(lpsn_type_material_file) as f:
            f.readline()
            for line in f:
                tokens = [t.strip() for t in line.split('\t')]
                
                rank = tokens[0]
                taxon = tokens[1].replace('Candidatus ', '')
                type_taxon = tokens[2].replace('Candidatus ', '')
                
                if type_taxon:
                    taxon_priority_year = LPSN.parse_lpsn_priority_year(tokens[3])
                    
                    if rank in ['Species', 'Subspecies']:
                        # the type material for species and subspecies is a specific strain,
                        # and this strain has the same priority as the species or subspecies
                        type_taxon_priority_year = taxon_priority_year
                    else:
                        if tokens[4]:
                            type_taxon_priority_year = LPSN.parse_lpsn_priority_year(tokens[4])
                        else:
                            self.logger.warning('Type material has no priority information: {}'.format(type_taxon))
                            type_taxon_priority_year = Genome.NO_PRIORITY_YEAR

                    taxa[taxon] = LPSN_Taxon(taxon, taxon_priority_year, type_taxon, type_taxon_priority_year)
                    
        return taxa

    def parse_lpsn_gss_file(self, lpsn_gss_metadata_file):
        """Parse LPSN GSS (genus-species-subspecies) metadata."""
        
        # get mapping between record numbers and LPSN taxa
        lpsn_record_taxon_map = {}
        with open(lpsn_gss_metadata_file, encoding='utf-8') as f:
            csv_reader = csv.reader(f)

            for line_num, tokens in enumerate(csv_reader):
                if line_num == 0:
                    genus_idx = tokens.index('genus_name')
                    sp_idx = tokens.index('sp_epithet')
                    subsp_idx = tokens.index('subsp_epithet')
                    record_no_idx = tokens.index('record_no')
                else:
                    if tokens[genus_idx] != '' and tokens[sp_idx] != '' and tokens[subsp_idx] != '':
                        taxon = 's__{} {} subsp. {}'.format(tokens[genus_idx].strip(), tokens[sp_idx].strip(), tokens[subsp_idx].strip())
                    elif tokens[genus_idx] != '' and  tokens[sp_idx] != '':
                        taxon = 's__{} {}'.format(tokens[genus_idx].strip(), tokens[sp_idx].strip())
                    else:
                        taxon = 'g__{}'.format(tokens[genus_idx].strip())

                    lpsn_record_taxon_map[tokens[record_no_idx]] = taxon
            
        # get synonyms and names considered correct at LPSN 
        lpsn_sp_correct_names = {}
        lpsn_sp_synonyms = defaultdict(set)
        with open(lpsn_gss_metadata_file, encoding='utf-8') as f:
            csv_reader = csv.reader(f)

            for line_num, tokens in enumerate(csv_reader):
                if line_num == 0:
                    genus_idx = tokens.index('genus_name')
                    sp_idx = tokens.index('sp_epithet')
                    subsp_idx = tokens.index('subsp_epithet')
                    status_idx = tokens.index('status')
                    type_idx = tokens.index('nomenclatural_type')
                    record_no_idx = tokens.index('record_no')
                    record_lnk_idx = tokens.index('record_lnk')
                else:
                    if tokens[genus_idx] != '' and  tokens[sp_idx] != '' and tokens[subsp_idx] != '':
                        # process subspecies
                        pass
                    elif tokens[genus_idx] != '' and  tokens[sp_idx] != '':
                        # process species
                        taxon = 's__{} {}'.format(tokens[genus_idx].strip(), tokens[sp_idx].strip())
                        
                        if 'correct name' in tokens[status_idx]:
                            assert tokens[record_lnk_idx] == ''
                            lpsn_sp_correct_names[taxon] = taxon
                        elif tokens[record_lnk_idx] != '':
                            # for species, this should link to the correct name for
                            # a species according to LPSN
                            lpsn_sp_correct_names[taxon] = lpsn_record_taxon_map[tokens[record_lnk_idx]]
                            
                        if 'synonym' in tokens[status_idx] and tokens[record_lnk_idx] != '':
                            synonym_sp = lpsn_record_taxon_map[tokens[record_lnk_idx]]
                            lpsn_sp_synonyms[taxon].add(synonym_sp)
                    else:
                        # process genus
                        pass
                        
        return lpsn_sp_correct_names, lpsn_sp_synonyms