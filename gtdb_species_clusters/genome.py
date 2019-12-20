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
from dataclasses import dataclass


@dataclass
class Genome(object):
    """Single genome."""
    
    gid: str
    ncbi_accn: str
    gtdb_rid: str
    gtdb_taxonomy: list
    ncbi_taxonomy: list
    ncbi_taxonomy_unfiltered: list
    gtdb_type_designation: str
    gtdb_type_designation_sources: str
    ncbi_type_material: str
    ncbi_strain_identifiers: str
    ncbi_assembly_level: str
    ncbi_genome_representation: str
    ncbi_refseq_category: str
    ncbi_genome_category: str
    comp: float
    cont: float
    strain_heterogeneity_100: float
    length: int
    contig_count: int
    contig_n50: int
    scaffold_count: int
    ambiguous_bases: int
    total_gap_len: int
    ssu_count: int
    ssu_length: int
    ncbi_molecule_count: int
    ncbi_unspanned_gaps: int
    ncbi_spanned_gaps: int
    
    NCBI_TYPE_SPECIES = set(['assembly from type material', 
                                'assembly from neotype material',
                                'assembly designated as neotype',
                                'assembly designated as reftype'])
    NCBI_PROXYTYPE = set(['assembly from proxytype material'])
    NCBI_TYPE_SUBSP = set(['assembly from synonym type material'])

    GTDB_TYPE_SPECIES = set(['type strain of species', 'type strain of neotype'])
    GTDB_TYPE_SUBSPECIES = set(['type strain of subspecies', 'type strain of heterotypic synonym'])
    GTDB_NOT_TYPE_MATERIAL = set(['not type material'])

    def __post_init__(self):
        """Post data initialization."""

        self.logger = logging.getLogger('timestamp')
        
        if self.gid.startswith('UBA'):
            self.ncbi_genome_category = 'metagenome'
            
        self.genomic_file = None
        
    def __str__(self):
        """User-friendly string representation."""

        return (
            f'{{gid:{self.gid}, '
            f'gtdb_sp:{self.gtdb_taxonomy[6]}, '
            f'ncbi_sp:{self.ncbi_taxonomy[6]}}}'
        )

    def is_gtdb_sp_rep(self):
        """Check if genome is representative of species."""

        return self.gtdb_rid == self.gid
        
    def is_ncbi_subspecies(self):
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
        
    def is_isolate(self):
        """Check if genome is an isolate."""
        
        if self.ncbi_genome_category:
            if ('metagenome' in self.ncbi_genome_category.lower()
                or 'environmental' in self.ncbi_genome_category.lower()
                or 'single cell' in self.ncbi_genome_category.lower()):
                return False
                
        return True
        
    def is_mag(self):
        """Check if genome is a MAG."""
        
        if not self.is_isolate():
            if ('metagenome' in self.ncbi_genome_category.lower()
                or 'environmental' in self.ncbi_genome_category.lower()):
                return True
                
        return False
        
    def is_sag(self):
        """Check if genome is a SAG."""
        
        if not self.is_isolate():
            if 'single cell' in self.ncbi_genome_category.lower():
                return True
                
        return False
    
    def is_gtdb_type_strain(self):
        """Check if genome is a type strain genome in GTDB."""
        
        return self.gtdb_type_designation in Genome.GTDB_TYPE_SPECIES
        
    def is_ncbi_type_strain(self):
        """Check if genome is a type strain genome at NCBI."""
        
        if self.ncbi_type_material:
            if self.ncbi_type_material.lower() in Genome.NCBI_TYPE_SPECIES:
                if not self.is_ncbi_subspecies():
                    return True
                else:
                    # genome is marked as 'assembled from type material', but
                    # is a subspecies according to the NCBI taxonomy so is
                    # the type strain of a subspecies
                    # (e.g. Alpha beta subsp. gamma)
                    return False
                    
        return False
        
    def is_complete_genome(self):
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

    def ncbi_domain(self):
        """Get NCBI domain classification."""
        
        return self.ncbi_taxonomy[0]

    def gtdb_domain(self):
        """Get GTDB domain classification."""
        
        return self.gtdb_taxonomy[0]
        
    def score_ani(self, ani):
        """Calculate balanced score that accounts for ANI to previous representative."""
        
        ani_score = 100 - 20*(100-ani)
        return 0.5*ani_score + 0.5*self.score_type_strain()

    def score_type_strain(self):
        """"Calculate score of genomes with preference to type strain genomes."""

        # set base quality so genomes have the following priority order:
        #  type strain genome
        #  NCBI assembled from type material
        q = 0
        if self.is_gtdb_type_strain():
            q = 1e5
        
        if self.is_ncbi_type_strain():
            q = 1e4

        q += self.score_assembly()
            
        return q

    def score_assembly(self):
        """Calculate score indicating quality of genome assembly."""
        
        q = 0
        
        # check if genome appears to complete consist of only an unspanned
        # chromosome and unspanned plasmids and thus should be considered
        # very high quality
        
        if self.is_complete_genome():
            q += 100
            
        q += self.comp - 5*self.cont
        q -= 5*float(self.contig_count)/100
        q -= 5*float(self.ambiguous_bases)/1e5
        
        if self.is_sag():
            q -= 100 # genome is a SAG
        elif self.is_mag(): 
            q -= 100 # genome is a MAG
        
        # check for near-complete 16S rRNA gene
        gtdb_domain = self.gtdb_taxonomy[0]
        min_ssu_len = 1200
        if gtdb_domain == 'd__Archaea':
            min_ssu_len = 900
            
        if self.ssu_length and self.ssu_length >= min_ssu_len:
            q += 10
            
        return q
        
    def pass_qc(self,
                marker_perc,
                min_comp,
                max_cont,
                min_quality,
                sh_exception,
                min_perc_markers,
                max_contigs,
                min_N50,
                max_ambiguous,
                failed_tests):
        """Check if genome passes QC."""
        
        failed = False
        if self.comp < min_comp:
            failed_tests['comp'] += 1
            failed = True
        
        if self.strain_heterogeneity_100 >= sh_exception:
            if self.cont > 20:
                failed_tests['cont'] += 1
                failed = True
            q = self.comp - 5*self.cont*(1.0 - self.strain_heterogeneity_100/100.0)
            if q < min_quality:
                failed_tests['qual'] += 1
                failed = True
        else:
            if self.cont > max_cont:
                failed_tests['cont'] += 1
                failed = True
            q = self.comp - 5*self.cont
            if q < min_quality:
                failed_tests['qual'] += 1
                failed = True
                
        if marker_perc < min_perc_markers:
            failed_tests['marker_perc'] += 1
            failed = True
                
        if self.contig_count > max_contigs:
            failed_tests['contig_count'] += 1
            failed = True
        if self.contig_n50 < min_N50:
            failed_tests['N50'] += 1
            failed = True
        if self.ambiguous_bases > max_ambiguous:
            failed_tests['ambig'] += 1
            failed = True
        
        return not failed
