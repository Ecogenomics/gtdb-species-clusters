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
import functools
from collections import defaultdict

from gtdb_species_clusters.genome import Genome
from gtdb_species_clusters.genome_utils import select_highest_quality


class SpeciesPriorityManager(object):
    """Resolve naming priority of species names."""

    def __init__(self, species_priority_ledger, dsmz_bacnames_file):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')
        
        self._parse_sp_priority_ledger(species_priority_ledger)
        self._parse_dsmz_bacnames(dsmz_bacnames_file)

    def _parse_sp_priority_ledger(self, species_priority_ledger):
        """Parse manually resolved priority cases."""
        
        self.logger.info('Parsing species priority ledger.')
        self.manual_priority = defaultdict(lambda: {})
        
        num_cases = 0
        with open(species_priority_ledger, encoding='utf-8') as f:
            header = f.readline().strip().split('\t')
            
            spA_index = header.index('NCBI species A')
            spB_index = header.index('NCBI species B')
            priority_sp_index = header.index('Priority')
            
            for line in f:
                tokens = line.strip().split('\t')
                
                spA = tokens[spA_index].strip()
                spB = tokens[spB_index].strip()
                priority_sp = tokens[priority_sp_index].strip()
                
                assert spA.startswith('s__')
                assert spB.startswith('s__')
                assert priority_sp.startswith('s__')
                assert priority_sp in [spA, spB]
                
                self.manual_priority[spA][spB] = priority_sp
                self.manual_priority[spB][spA] = priority_sp
                num_cases += 1
                
        self.logger.info(f' ... identified {num_cases:,} manually resolved cases.')
        
    def _parse_dsmz_bacnames(self, dsmz_bacnames_file):
        """Parse priority information from LPSN at DSMZ."""
        
        self._genus_priority = {}
        with open(dsmz_bacnames_file, encoding='utf-8') as f:
            header = f.readline().strip().split('\t')
            
            genus_index = header.index('GENUS')
            sp_index = header.index('SPECIES')
            status_index = header.index('STATUS')
            author_index = header.index('AUTHORS')
            
            for line in f:
                tokens = line.strip().split('\t')
                
                genus = tokens[genus_index].strip().replace('"', '')
                status = tokens[status_index].strip().replace('"', '')
                authors = tokens[author_index].strip().replace('"', '')
                if authors.startswith('('):
                    authors = authors[authors.find(')')+1:].strip()

                for idx, ch in enumerate(authors): 
                    if ch.isdigit():
                        year = int(authors[idx:idx+4])
                        break

                if status.startswith('gen. nov.') or status.startswith('genus'):
                    self._genus_priority['g__' + genus] = year

    def genus_priority(self, genus):
        """Determine year of priority for genus."""
        
        if genus in self._genus_priority:
            return self._genus_priority[genus]
                    
        return Genome.NO_PRIORITY_YEAR

    def priority(self, cur_genomes, gid1, gid2):
        """Resolve priority between genomes."""

        if (cur_genomes[gid1].is_gtdb_type_species()
            and not cur_genomes[gid2].is_gtdb_type_species()):
            return gid1, 'type species of genus'
        elif (not cur_genomes[gid1].is_gtdb_type_species()
            and cur_genomes[gid2].is_gtdb_type_species()):
            return gid2, 'type species of genus'
        elif cur_genomes[gid1].year_of_priority() < cur_genomes[gid2].year_of_priority():
            return gid1, 'year of priority'
        elif cur_genomes[gid1].year_of_priority() > cur_genomes[gid2].year_of_priority():
            return gid2, 'year of priority'
        elif (cur_genomes[gid1].is_effective_type_strain() 
                and not cur_genomes[gid2].is_effective_type_strain()):
            return gid1, 'effective type strain'
        elif (not cur_genomes[gid1].is_effective_type_strain() 
                and cur_genomes[gid2].is_effective_type_strain()):
            return gid2, 'effective type strain'
        elif (not cur_genomes[gid1].is_effective_type_strain() 
                and not cur_genomes[gid2].is_effective_type_strain()):
            # neither genome it effective type strain
            hq_gid = select_highest_quality([gid1, gid2], cur_genomes)
            return hq_gid, 'highest quality genome'
        elif (cur_genomes[gid1].is_exclusively_effective_type_strain() 
                and cur_genomes[gid2].is_exclusively_effective_type_strain()):
            # both genomes are only effective, but not validated type material
            hq_gid = select_highest_quality([gid1, gid2], cur_genomes)
            return hq_gid, 'highest quality genome'

        # both genomes are type strain of species, but priority
        # is ambiguous using just publication date so species need
        # to be in the priority ledger
        sp1 = cur_genomes[gid1].ncbi_taxa.species
        sp2 = cur_genomes[gid2].ncbi_taxa.species
        if sp1 not in self.manual_priority or sp2 not in self.manual_priority[sp1]:
            self.logger.error('Ambiguous priority based on publication date.')
            self.logger.error('Species need to be manually resolved in priority ledger: {}: {} / {} / {}, {}: {} / {} / {}'.format(
                                gid1, sp1, cur_genomes[gid1].is_gtdb_type_strain(), cur_genomes[gid1].is_effective_type_strain(),
                                gid2, sp2, cur_genomes[gid2].is_gtdb_type_strain(), cur_genomes[gid2].is_effective_type_strain()))
            #***sys.exit(-1)
            return gid1, 'error: ambiguous priority date' #***
            
        priority_sp = self.manual_priority[sp1][sp2]
        if sp1 == priority_sp:
            return gid1, 'manual curation'
            
        return gid2, 'manual curation'
        
    def has_priority(self, cur_genomes, gid1, gid2):
        """Check if genome 1 has priority over genome 2."""
        
        if (cur_genomes[gid1].is_gtdb_type_strain() 
            and not cur_genomes[gid2].is_gtdb_type_strain()):
            return True
        elif (not cur_genomes[gid1].is_gtdb_type_strain() 
            and cur_genomes[gid2].is_gtdb_type_strain()):
            return False
            
        priority_gid, note = self.priority(cur_genomes, gid1, gid2)
        
        if priority_gid == gid1:
            return True
            
        return False
        
    def _priority_sort(self, gid1, gid2):
        """Sort genomes by priority."""
        
        priority_gid, note = self.priority(self.cur_genomes, gid1, gid2)
        if priority_gid == gid1:
            return -1
            
        return 1
        
    def sort_by_priority(self, cur_genomes, gids, reverse=False):
        """Sort set of genomes by priority with highest priority genome first."""
        
        self.cur_genomes = cur_genomes
        sorted_gids = sorted(gids, 
                                key=functools.cmp_to_key(self._priority_sort), 
                                reverse=reverse)
        
        return sorted_gids
