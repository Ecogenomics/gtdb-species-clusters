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
    type_material_priority_year: int


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

    def __init__(self, lpsn_data, lpsn_gss_metadata_file):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        if lpsn_data:
            self.taxa = self.parse_lpsn_data_file(lpsn_data)

        if lpsn_gss_metadata_file:
            self.sp_correct_names, self.sp_synonyms = self.parse_lpsn_gss_file(lpsn_gss_metadata_file)

    def type_material(self, taxon):
        """Get type material of taxon."""

        if taxon in self.taxa:
            return self.taxa[taxon].type_material

        return None

    def parse_lpsn_data_file(self, lpsn_data):
        """Parse type material for taxa as defined at LPSN."""

        rank_prefixes = {'domain': 'd__',
                         'phylum': 'p__',
                         'class': 'c__',
                         'order': 'o__',
                         'family': 'f__',
                         'genus': 'g__',
                         'species': 's__'}

        # family is never type material
        type_prefixes = {0: 'c__', 1: 'o__', 2: 'g__', 3: 's__', 4: ''}

        taxa = {}
        with open(lpsn_data) as f:
            header = f.readline().strip().split('\t')

            rank_idx = header.index('Rank')
            name_idx = header.index('Name')
            priority_idx = header.index('Priority')
            type_class_idx = header.index('Type class')
            type_strain_idx = header.index('Type strain')

            for line in f:
                tokens = [t.strip() for t in line.split('\t')]

                rank = tokens[rank_idx]
                if rank == 'subspecies':
                    # not currently considering subspecies
                    continue

                prefix = rank_prefixes[rank]
                taxon = f"{prefix}{tokens[name_idx].replace('Candidatus ', '')}"

                priority_citations = tokens[priority_idx]
                if priority_citations.lower() == 'n/a':
                    continue
                taxon_priority_year = LPSN.parse_lpsn_priority_year(priority_citations)

                type_taxa = []
                for idx, type_idx in enumerate(range(type_class_idx, type_strain_idx+1)):
                    type_taxon = tokens[type_idx]
                    if type_taxon.lower() != 'n/a':
                        type_taxa.append(f"{type_prefixes[idx]}{type_taxon.replace('Candidatus ', '')}")

                if len(type_taxa) == 0:
                    # No type material for this LPSN taxon which occurs
                    # for placeholder taxa and candidatus taxa
                    type_taxa = ['n/a']
                elif len(type_taxa) > 1:
                    self.logger.error(f'Identified multiple type taxa for {taxon}: {type_taxa}')
                    sys.exit(-1)

                type_taxon = type_taxa[0]
                taxa[taxon] = LPSN_Taxon(taxon, taxon_priority_year, type_taxon, None)

        # fill in priority year information for type taxa
        missing_type_count = 0
        for taxon in taxa:
            if taxon.startswith('s__'):
                # the type material for species is a specific strain,
                # and this strain has the same priority as the species
                type_priority = taxa[taxon].priority_year
            else:
                type_taxon = taxa[taxon].type_material
                if type_taxon not in taxa:
                    missing_type_count += 1
                    continue
                type_priority = taxa[type_taxon].priority_year

            taxa[taxon].type_material_priority_year = type_priority

        if missing_type_count:
            self.logger.error(
                f'Identified {missing_type_count:,} taxa without entries for type material. This is currently a known issue that needs to be addressed.')

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
                        taxon = 's__{} {} subsp. {}'.format(
                            tokens[genus_idx].strip(), tokens[sp_idx].strip(), tokens[subsp_idx].strip())
                    elif tokens[genus_idx] != '' and tokens[sp_idx] != '':
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
                    if tokens[genus_idx] != '' and tokens[sp_idx] != '' and tokens[subsp_idx] != '':
                        # process subspecies
                        pass
                    elif tokens[genus_idx] != '' and tokens[sp_idx] != '':
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
