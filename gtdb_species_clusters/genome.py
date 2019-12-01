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
from collections import defaultdict


class Genome(object):
    """Single genome."""
    
    NCBI_TYPE_SPECIES = set(['assembly from type material', 
                        'assembly from neotype material',
                        'assembly designated as neotype'])
    NCBI_PROXYTYPE = set(['assembly from proxytype material'])
    NCBI_TYPE_SUBSP = set(['assembly from synonym type material'])

    GTDB_TYPE_SPECIES = set(['type strain of species', 'type strain of neotype'])
    GTDB_TYPE_SUBSPECIES = set(['type strain of subspecies', 'type strain of heterotypic synonym'])
    GTDB_NOT_TYPE_MATERIAL = set(['not type material'])

    def __init__(self,
                    gid,
                    ncbi_accn,
                    gtdb_rid,
                    gtdb_taxonomy,
                    ncbi_taxonomy,
                    ncbi_taxonomy_unfiltered,
                    gtdb_type_designation,
                    ncbi_type_material,
                    ncbi_assembly_level,
                    ncbi_genome_representation,
                    ncbi_refseq_category,
                    ncbi_genome_category,
                    comp,
                    cont,
                    gs,
                    contig_count,
                    n50,
                    scaffold_count,
                    ambiguous_bases,
                    total_gap_len,
                    ssu_count,
                    ssu_length,
                    ncbi_molecule_count,
                    ncbi_unspanned_gaps,
                    ncbi_spanned_gaps):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')
        
        self.gid = gid
        self.ncbi_accn = ncbi_accn
        self.gtdb_taxonomy = gtdb_taxonomy
        self.ncbi_taxonomy = ncbi_taxonomy
        self.ncbi_taxonomy_unfiltered = ncbi_taxonomy_unfiltered
        self.gtdb_type_designation = gtdb_type_designation
        self.ncbi_type_material = ncbi_type_material
        self.ncbi_assembly_level = ncbi_assembly_level
        self.ncbi_genome_representation = ncbi_genome_representation
        self.ncbi_refseq_category = ncbi_refseq_category
        self.ncbi_genome_category = ncbi_genome_category
        if self.gid.startswith('UBA'):
            self.ncbi_genome_category = 'metagenome'
            
        self.comp = comp
        self.cont = cont
        self.gs = gs
        self.contig_count = contig_count
        self.n50 = n50
        self.scaffold_count = scaffold_count
        self.ambiguous_bases = ambiguous_bases
        self.total_gap_len = total_gap_len
        self.ssu_count = ssu_count
        self.ssu_length = ssu_length
        self.ncbi_molecule_count = ncbi_molecule_count
        self.ncbi_unspanned_gaps = ncbi_unspanned_gaps
        self.ncbi_spanned_gaps = ncbi_spanned_gaps
        
        self.is_gtdb_sp_rep = (gtdb_rid == gid)
        self.is_ncbi_subspecies = self._is_ncbi_subspecies()
        self.is_isolate = self._is_isolate()
        self.is_mag = self._is_mag()
        self.is_sag = self._is_sag()
        assert(self.is_isolate ^ self.is_mag ^ self.is_sag)
        self.is_gtdb_type_strain = self._is_gtdb_type_strain()
        self.is_ncbi_type_strain = self._is_gtdb_type_strain()
        self.is_complete_genome = self._is_complete_genome()
        
        self.genomic_file = None
        
    def __str__(self):
        """User-friendly string representation."""

        return (
            f'{{gid:{self.gid}, '
            f'gtdb_sp:{self.gtdb_taxonomy[6]}, '
            f'ncbi_sp:{self.ncbi_taxonomy[6]}}}'
        )
        
    def _is_ncbi_subspecies(self):
        """Check if genome is a subspecies at NCBI."""

        ncbi_subspecies = None
        for taxon in self.ncbi_taxonomy_unfiltered:
            if taxon.startswith('sb__'):
                ncbi_subspecies = taxon[4:]
                
                # fix odd designation impacting less than a dozen genomes 
                ncbi_subspecies = ncbi_subspecies.replace(' pv. ', ' subsp. ') 

        if not ncbi_subspecies:
            return False

        if 'subsp.' not in ncbi_subspecies:
            self.logger.warning(f'Genome {self.gid} has invalid subspecies: {ncbi_subspecies}')
            return False

        tokens = ncbi_subspecies.split()
        subsp_index = tokens.index('subsp.')
        if tokens[subsp_index - 1] == tokens[subsp_index + 1]:
            return False

        return True
        
    def _is_isolate(self):
        """Check if genome is an isolate."""
        
        if self.ncbi_genome_category:
            if ('metagenome' in self.ncbi_genome_category.lower()
                or 'environmental' in self.ncbi_genome_category.lower()
                or 'single cell' in self.ncbi_genome_category.lower()):
                return False
                
        return True
        
    def _is_mag(self):
        """Check if genome is a MAG."""
        
        if not self.is_isolate:
            if ('metagenome' in self.ncbi_genome_category.lower()
                or 'environmental' in self.ncbi_genome_category.lower()):
                return True
                
        return False
        
    def _is_sag(self):
        """Check if genome is a SAG."""
        
        if not self.is_isolate:
            if 'single cell' in self.ncbi_genome_category.lower():
                return True
                
        return False
    
    def _is_gtdb_type_strain(self):
        """Check if genome is a type strain genome in GTDB."""
        
        return self.gtdb_type_designation in Genome.GTDB_TYPE_SPECIES
        
    def _is_ncbi_type_strain(self):
        """Check if genome is a type strain genome at NCBI."""
        
        if self.ncbi_type_material:
            if self.ncbi_type_material.lower() in Genome.NCBI_TYPE_SPECIES:
                if not self.is_ncbi_subspecies:
                    return True
                else:
                    # genome is marked as 'assembled from type material', but
                    # is a subspecies according to the NCBI taxonomy so is
                    # the type strain of a subspecies
                    # (e.g. Alpha beta subsp. gamma)
                    return False
                    
        return False
        
    def _is_complete_genome(self):
        """Check if genome is a complete assembly."""
        
        return (self.ncbi_assembly_level 
                and self.ncbi_assembly_level.lower() in ['complete genome', 'chromosome']
                and self.ncbi_genome_representation
                and self.ncbi_genome_representation.lower() == 'full'
                and self.scaffold_count == self.ncbi_molecule_count
                and self.ncbi_unspanned_gaps == 0
                and self.ncbi_spanned_gaps <= 10
                and self.ambiguous_bases <= 1e4
                and self.total_gap_len <= 1e4
                and self.ssu_count >= 1)
                
    def ncbi_subspecies(self):
        """Get NCBI subspecies classification of genome."""
        
        for taxon in self.ncbi_taxonomy_unfiltered:
            if taxon.startswith('sb__'):
                return taxon
                    
        return None
        
    def ncbi_species(self):
        """Get NCBI species classification."""
        
        return self.ncbi_taxonomy[6]
        
    def gtdb_species(self):
        """Get GTDB species classification."""
        
        return self.gtdb_taxonomy[6]
        
    def ncbi_specific_epithet(self):
        """Get NCBI specific epithet."""
        
        ncbi_sp = self.ncbi_species()
        if ncbi_sp == 's__':
            return ''
        
        s = ncbi_sp.replace('Candidatus ', '')
        generic, specific = s.split()
        
        return specific
        
    def gtdb_specific_epithet(self):
        """Get GTDB specific epithet."""
        
        gtdb_sp = self.gtdb_species()
        if gtdb_sp == 's__':
            return ''

        generic, specific = gtdb_sp.split()
        
        return specific
        
    def ncbi_genus(self):
        """Get NCBI genus classification."""
        
        return self.ncbi_taxonomy[5]

    def gtdb_genus(self):
        """Get GTDB genus classification."""
        
        return self.gtdb_taxonomy[5]

    def score_update(self):
        """"Calculate score of genomes for updating representatives."""

        # set base quality so genomes have the following priority order:
        #  type strain genome
        #  NCBI assembled from type material
        #  NCBI reference or representative genome
        #  type subspecies genome
        q = 0
        if self.is_gtdb_type_strain:
            q = 1e4
        
        if self.is_ncbi_type_strain:
            q = 1e3
            
        q += self.score_assembly()
            
        return q

    def score_assembly(self):
        """Calculate score indicating quality of genome assembly."""
        
        q = 0
        
        # check if genome appears to complete consist of only an unspanned
        # chromosome and unspanned plasmids and thus should be considered
        # very high quality
        
        if self.is_complete_genome:
            q += 100
            
        q += self.comp - 5*self.cont
        q -= 5*float(self.contig_count)/100
        q -= 5*float(self.ambiguous_bases)/1e5
        
        if self.is_sag:
            q -= 100 # genome is a SAG
        elif self.is_mag: 
            q -= 200 # genome is a MAG
        
        # check for near-complete 16S rRNA gene
        gtdb_domain = self.gtdb_taxonomy[0]
        min_ssu_len = 1200
        if gtdb_domain == 'd__Archaea':
            min_ssu_len = 900
            
        if self.ssu_length and self.ssu_length >= min_ssu_len:
            q += 10
            
        return q
