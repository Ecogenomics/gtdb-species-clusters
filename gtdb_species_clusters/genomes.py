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

from gtdb_species_clusters.genome_utils import canonical_gid
from gtdb_species_clusters.genome import Genome
from gtdb_species_clusters.species_clusters import SpeciesClusters


class Genomes(object):
    """A collection of genomes."""

    def __init__(self):
        """Initialization."""

        self.genomes = {}
        self.sp_clusters = SpeciesClusters()
        
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
        
    def _apply_ncbi_taxonomy_ledgers(self,
                                        species_exception_file, 
                                        genus_exception_file):
        """Apply corrections to NCBI taxonomy."""

        species_updates = {}
        if species_exception_file:
            with open(species_exception_file, encoding='utf-8') as f:
                f.readline()
                for line in f:
                    line_split = [token.strip() for token in line.strip().split('\t')]
                    gid = canonical_gid(line_split[0])

                    sp = line_split[1].replace('Candidatus ', '')
                    if gid not in self.genomes:
                        self.logger.warning(f'Genome {gid} in species exception list not defined in genome set.')
                        continue
                        
                    if not sp.startswith('s__'):
                        sp = 's__' + sp
                        
                    self.genomes[gid].ncbi_taxonomy[6] = sp
                    species_updates[gid] = sp
        
        if genus_exception_file:
            with open(genus_exception_file, encoding='utf-8') as f:
                f.readline()
                for line in f:
                    line_split = [token.strip() for token in line.strip().split('\t')]
                    gid = canonical_gid(line_split[0])
                    genus = line_split[1]
                    if gid not in self.genomes:
                        self.logger.warning(f'Genome {gid} in genus exception list not defined in genome set.')
                        continue
                        
                    if genus.startswith('g__'):
                        genus = genus[3:]
                        
                    self.genomes[gid].ncbi_taxonomy[5] = f'g__{genus}'
                    
                    species = self.genomes[gid].ncbi_species()
                    if species != 's__':
                        specific = self.genomes[gid].ncbi_specific_epithet()
                        self.genomes[gid].ncbi_taxonomy[6] = f's__{genus} {specific}'

                    # sanity check ledgers
                    if gid in species_updates and genus not in species_updates[gid]:
                        self.logger.error(f'Species and genus ledgers have conflicting assignments for {gid}.')
                        sys.exit(-1)

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
        
    def load_genomic_file_paths(self, genome_path_file):
        """Determine path to genomic FASTA file for each genome."""

        for line in open(genome_path_file):
            line_split = line.strip().split('\t')
            
            gid = line_split[0]
            gid = canonical_gid(gid)
            if gid in self.genomes:
                genome_path = line_split[1]
                accession = os.path.basename(os.path.normpath(genome_path))
                genomic_file = os.path.join(genome_path, accession + '_genomic.fna')
                self.genomes[gid].genomic_file = genomic_file

    def load_from_metadata_file(self, 
                                metadata_file,
                                species_exception_file,
                                genus_exception_file,
                                create_sp_clusters=True):
        """Create genome set from file(s)."""

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
            
            gtdb_genome_rep_index = headers.index('gtdb_genome_representative')

            for line in f:
                line_split = line.strip().split('\t')
                
                ncbi_accn = line_split[genome_index]
                gid = canonical_gid(ncbi_accn)

                if gid.startswith('U_'):
                    # check if genome has a UBA identifier
                    org_name_index = headers.index('organism_name')
                    org_name = line_split[org_name_index]
                    if '(UBA' in org_name:
                        gid = org_name[org_name.find('(')+1:-1]
                    else:
                        continue # skip non-UBA user genomes
                
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
                
                gtdb_rid = canonical_gid(line_split[gtdb_genome_rep_index])
                if create_sp_clusters:
                    self.sp_clusters.update_sp_cluster(gtdb_rid, gid, gtdb_taxonomy[6])

                self.genomes[gid] = Genome(gid,
                                            ncbi_accn,
                                            gtdb_rid,
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
                                            
        self._apply_ncbi_taxonomy_ledgers(species_exception_file,
                                            genus_exception_file)
