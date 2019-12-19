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

from gtdb_species_clusters.genome import Genome
from gtdb_species_clusters.species_clusters import SpeciesClusters
from gtdb_species_clusters.genome_utils import canonical_gid
from gtdb_species_clusters.taxon_utils import is_placeholder_taxon


class Genomes(object):
    """A collection of genomes."""

    def __init__(self):
        """Initialization."""

        self.genomes = {}
        self.sp_clusters = SpeciesClusters()
        
        self.genomic_files = {} # this should be removed, but the FastANI interface
                                # currently requires data to be passed in as a dictionary
                                
        self.user_uba_id_map = {}   # can be removed once UBA genomes are no longer part
                                    # of the archaeal genome set
        
        self.logger = logging.getLogger('timestamp')

    def __str__(self):
        """User-friendly string representation."""
        
        return '{{num_genomes:{}}}'.format(len(self.genomes))
        
    def __iter__(self):
        """Iterate over genome IDs comprising genome set."""
        
        for gid in self.genomes:
            yield gid
        
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
            taxonomy_str = taxonomy_str.replace('Candidatus ', '')
            return [t.strip() for t in taxonomy_str.split(';')]

        return list(Taxonomy.rank_prefixes)
        
    def _convert_int(self, value):
        """Convert database value to integer."""
        
        return int(value) if value and value != 'none' else 0
        
    def _convert_float(self, value):
        """Convert database value to float."""
        
        return float(value) if value and value != 'none' else 0
        
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
                    gid = canonical_gid(line_split[0].strip())

                    sp = line_split[1].strip().replace('Candidatus ', '')
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
                    gid = canonical_gid(line_split[0].strip())
                    genus = line_split[1].strip()
                    if gid not in self.genomes:
                        self.logger.warning(f'Genome {gid} in genus exception list not defined in genome set.')
                        continue
                        
                    if genus.startswith('g__'):
                        genus = genus[3:]
                        
                    self.genomes[gid].ncbi_taxonomy[5] = f'g__{genus}'
                    
                    species = self.genomes[gid].ncbi_species
                    if species != 's__':
                        specific = self.genomes[gid].ncbi_specific_epithet
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
        
    def named_ncbi_species(self):
        """Get genomes in valid or effectively published, including Candidatus, species in NCBI taxonomy."""
        
        named_ncbi_sp = defaultdict(set)
        for gid in self.genomes:
            if not is_placeholder_taxon(self.genomes[gid].ncbi_species):
                named_ncbi_sp[self.genomes[gid].ncbi_species].add(gid)

        return named_ncbi_sp
        
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
                self.genomic_files[gid] = genomic_file

    def load_from_metadata_file(self, 
                                metadata_file,
                                species_exception_file=None,
                                genus_exception_file=None,
                                gtdb_type_strains_ledger=None,
                                create_sp_clusters=True,
                                uba_genome_file=None,
                                qc_passed_file=None):
        """Create genome set from file(s)."""
        
        valid_uba_ids = set()
        if uba_genome_file:
            with open(uba_genome_file) as f:
                for line in f:
                    line_split = line.strip().split('\t')
                    valid_uba_ids.add(line_split[0].strip())

        pass_qc_gids = set()
        if qc_passed_file:
            with open(qc_passed_file) as f:
                f.readline()
                for line in f:
                    line_split = line.strip().split('\t')
                    pass_qc_gids.add(line_split[0].strip())
                    
        gtdb_type_strains = set()
        if gtdb_type_strains_ledger:
            with open(gtdb_type_strains_ledger) as f:
                f.readline()
                for line in f:
                    tokens = line.strip().split('\t')
                    gid = canonical_gid(tokens[0].strip())
                    gtdb_type_strains.add(gid)

        with open(metadata_file, encoding='utf-8') as f:
            headers = f.readline().strip().split('\t')

            genome_index = headers.index('accession')

            gtdb_taxonomy_index = headers.index('gtdb_taxonomy')
            ncbi_taxonomy_index = headers.index('ncbi_taxonomy')
            ncbi_taxonomy_unfiltered_index = headers.index('ncbi_taxonomy_unfiltered')
            
            gtdb_type_index = headers.index('gtdb_type_designation')
            gtdb_type_sources_index = headers.index('gtdb_type_designation_sources')
            ncbi_strain_identifiers_index = headers.index('ncbi_strain_identifiers')
            ncbi_type_index = headers.index('ncbi_type_material_designation')
            ncbi_asm_level_index = headers.index('ncbi_assembly_level')
            ncbi_genome_representation_index = headers.index('ncbi_genome_representation')
            ncbi_refseq_cat_index = headers.index('ncbi_refseq_category')
            ncbi_genome_cat_index = headers.index('ncbi_genome_category')
            
            comp_index = headers.index('checkm_completeness')
            cont_index = headers.index('checkm_contamination')
            sh_100_index = None
            if 'checkm_strain_heterogeneity_100' in headers:
                sh_100_index = headers.index('checkm_strain_heterogeneity_100')
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
                        uba_id = org_name[org_name.find('(')+1:-1]
                        if uba_id in valid_uba_ids:
                            self.user_uba_id_map[gid] = uba_id
                            gid = uba_id
                    else:
                        continue # skip non-UBA user genomes
                        
                if pass_qc_gids and gid not in pass_qc_gids:
                    continue
                
                gtdb_taxonomy = self._get_taxa(line_split[gtdb_taxonomy_index])
                ncbi_taxonomy = self._get_taxa(line_split[ncbi_taxonomy_index])
                ncbi_taxonomy_unfiltered = self._get_taxa(line_split[ncbi_taxonomy_unfiltered_index])
                
                gtdb_type = line_split[gtdb_type_index]
                gtdb_type_sources = line_split[gtdb_type_sources_index]
                if gid in gtdb_type_strains:
                    gtdb_type = 'type strain of species'
                    gtdb_type_sources = 'GTDB curator'
                
                ncbi_type = line_split[ncbi_type_index]
                ncbi_strain_identifiers = line_split[ncbi_strain_identifiers_index]
                ncbi_asm_level = line_split[ncbi_asm_level_index]
                ncbi_genome_representation = line_split[ncbi_genome_representation_index]
                ncbi_refseq_cat = line_split[ncbi_refseq_cat_index]
                ncbi_genome_cat = line_split[ncbi_genome_cat_index]
                
                comp = float(line_split[comp_index])
                cont = float(line_split[cont_index])
                sh_100 = 0
                if sh_100_index:
                    sh_100 = self._convert_float(line_split[sh_100_index])
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
                                            gtdb_type_sources,
                                            ncbi_type,
                                            ncbi_strain_identifiers,
                                            ncbi_asm_level,
                                            ncbi_genome_representation,
                                            ncbi_refseq_cat,
                                            ncbi_genome_cat,
                                            comp,
                                            cont,
                                            sh_100,
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
