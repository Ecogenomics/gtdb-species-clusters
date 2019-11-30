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

from biolib.taxonomy import Taxonomy

from gtdb_species_clusters.common import canonical_gid


class Genomes(object):
    """A collection of genomes."""
    
    def __init__(self):
        """Initialization."""

        self.genomes = {}
        self.ncbi_accn_map = {}
        
        self.logger = logging.getLogger('timestamp')
        
    def __str__(self):
        """User-friendly string representation."""
        
        return '{{num_genomes:{}}}'.format(len(self.genomes))
        
    def __getitem__(self, gid):
        """Get genome."""
        
        return self.genomes[gid]
        
    def __contains__(self, gid):
        """Check if genome is in genome set."""
        
        return gid in self.genomes

    def __len__(self):
        """Size of genome set."""
        
        return len(self.genomes)

    def _get_taxa(self, taxonomy_str):
        """Convert taxonomy string to taxa list."""
        
        if taxonomy_str and taxonomy_str != 'none':
            return [t.strip() for t in taxonomy_str.split(';')]

        return list(Taxonomy.rank_prefixes)
        
    def _convert_int(self, value):
        """Convert database value to integer."""
        
        return int(value) if value and value != 'none' else 0
        
    def get_gid(self, idx):
        """Get ID of genome at specific index."""
        
        return list(self.genomes)[idx]
        
    def gtdb_type_strain_genomes(self):
        """Get genomes considered type strain of species by GTDB."""
        
        type_strain_gids = set()
        for gid, genome in self.genomes.items():
            if genome.is_gtdb_type_strain:
                type_strain_gids.add(gid)
                
        return type_strain_gids

    def load_from_metadata_file(self, metadata_file):
        """Load genomes from metadata file."""

        with open(metadata_file, encoding='utf-8') as f:
            headers = f.readline().strip().split('\t')

            genome_index = headers.index('accession')

            gtdb_taxonomy_index = headers.index('gtdb_taxonomy')
            ncbi_taxonomy_index = headers.index('ncbi_taxonomy')
            ncbi_taxonomy_unfiltered_index = headers.index('ncbi_taxonomy_unfiltered')
            
            gtdb_type_index = headers.index('gtdb_type_designation')
            ncbi_type_index = headers.index('ncbi_type_material_designation')
            ncbi_asm_level_index = headers.index('ncbi_assembly_level')
            ncbi_genome_representation_index = headers.index('ncbi_genome_representation')
            ncbi_refseq_cat_index = headers.index('ncbi_refseq_category')
            ncbi_genome_cat_index = headers.index('ncbi_genome_category')
            
            comp_index = headers.index('checkm_completeness')
            cont_index = headers.index('checkm_contamination')
            gs_index = headers.index('genome_size')
            contig_count_index = headers.index('contig_count')
            n50_index = headers.index('n50_contigs')
            scaffold_count_index = headers.index('scaffold_count')
            ambiguous_bases_index = headers.index('ambiguous_bases')
            total_gap_len_index = headers.index('total_gap_length')
            ssu_count_index = headers.index('ssu_count')
            ssu_length_index = headers.index('ssu_length')
            ncbi_molecule_count_index = headers.index('ncbi_molecule_count')
            ncbi_unspanned_gaps_index = headers.index('ncbi_unspanned_gaps')
            ncbi_spanned_gaps_index = headers.index('ncbi_spanned_gaps')

            for line in f:
                line_split = line.strip().split('\t')
                
                gid = canonical_gid(line_split[genome_index])
                
                if gid.startswith('U_'):
                    # check if genome has a UBA identifier
                    org_name_index = headers.index('organism_name')
                    org_name = line_split[org_name_index]
                    if '(UBA' in org_name:
                        gid = org_name[org_name.find('(')+1:-1]
                    else:
                        continue # skip non-UBA user genomes
                        
                self.ncbi_accn_map[gid] = line_split[genome_index]
                
                gtdb_taxonomy = self._get_taxa(line_split[gtdb_taxonomy_index])
                ncbi_taxonomy = self._get_taxa(line_split[ncbi_taxonomy_index])
                ncbi_taxonomy_unfiltered = self._get_taxa(line_split[ncbi_taxonomy_unfiltered_index])
                
                gtdb_type = line_split[gtdb_type_index]
                ncbi_type = line_split[ncbi_type_index]
                ncbi_asm_level = line_split[ncbi_asm_level_index]
                ncbi_genome_representation = line_split[ncbi_genome_representation_index]
                ncbi_refseq_cat = line_split[ncbi_refseq_cat_index]
                ncbi_genome_cat = line_split[ncbi_genome_cat_index]
                
                comp = float(line_split[comp_index])
                cont = float(line_split[cont_index])
                gs = int(line_split[gs_index])
                contig_count = int(line_split[contig_count_index])
                n50 = int(line_split[n50_index])
                scaffold_count = int(line_split[scaffold_count_index])
                ambiguous_bases = int(line_split[ambiguous_bases_index])
                total_gap_len = int(line_split[total_gap_len_index])
                ssu_count = int(line_split[ssu_count_index])
                ssu_length = self._convert_int(line_split[ssu_length_index])
                ncbi_molecule_count = self._convert_int(line_split[ncbi_molecule_count_index])
                ncbi_unspanned_gaps = self._convert_int(line_split[ncbi_unspanned_gaps_index])
                ncbi_spanned_gaps = self._convert_int(line_split[ncbi_spanned_gaps_index])

                self.genomes[gid] = Genome(gid,
                                            gtdb_taxonomy,
                                            ncbi_taxonomy,
                                            ncbi_taxonomy_unfiltered,
                                            gtdb_type,
                                            ncbi_type,
                                            ncbi_asm_level,
                                            ncbi_genome_representation,
                                            ncbi_refseq_cat,
                                            ncbi_genome_cat,
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
                                            ncbi_spanned_gaps)

class Genome(object):
    """Genome."""
    
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
        self.gtdb_taxonomy = gtdb_taxonomy
        self.ncbi_taxonomy = ncbi_taxonomy
        self.ncbi_taxonomy_unfiltered = ncbi_taxonomy_unfiltered
        self.gtdb_type_designation = gtdb_type_designation
        self.ncbi_type_material = ncbi_type_material
        self.ncbi_assembly_level = ncbi_assembly_level
        self.ncbi_genome_representation = ncbi_genome_representation
        self.ncbi_refseq_category = ncbi_refseq_category
        self.ncbi_genome_category = ncbi_genome_category
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
        
        self.is_ncbi_subspecies = self._is_ncbi_subspecies()
        self.is_isolate = self._is_isolate()
        self.is_gtdb_type_strain = self._is_gtdb_type_strain()
        self.is_ncbi_type_strain = self._is_gtdb_type_strain()
        self.is_complete_genome = self._is_complete_genome()
        
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
