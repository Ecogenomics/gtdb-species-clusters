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

from gtdb_species_clusters.genome import Genome
from gtdb_species_clusters.taxa import Taxa
from gtdb_species_clusters.species_clusters import SpeciesClusters
from gtdb_species_clusters.genome_utils import canonical_gid, exclude_from_refseq, read_gtdbtk_classifications
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
                                    
        self.uba_user_id_map = {}
        
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
        
        gid = self.user_uba_id_map.get(gid, gid)
        return self.genomes[gid]
        
    def __contains__(self, gid):
        """Check if genome is in genome set."""
        
        gid = self.user_uba_id_map.get(gid, gid)
        return gid in self.genomes

    def __len__(self):
        """Size of genome set."""
        
        return len(self.genomes)

    def _convert_int(self, value, default_value=0):
        """Convert database value to integer."""
        
        return int(value) if value and value != 'none' else default_value
        
    def _convert_float(self, value, default_value=0.0):
        """Convert database value to float."""
        
        return float(value) if value and value != 'none' else default_value
        
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
                        
                    self.genomes[gid].ncbi_taxa.species = sp
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
                        
                    self.genomes[gid].ncbi_taxa.genus = f'g__{genus}'
                    
                    species = self.genomes[gid].ncbi_taxa.species
                    if species != 's__':
                        specific = self.genomes[gid].ncbi_taxa.specific_epithet
                        self.genomes[gid].ncbi_taxa.species = f's__{genus} {specific}'

                    # sanity check ledgers
                    if gid in species_updates and genus not in species_updates[gid]:
                        self.logger.error(f'Species and genus ledgers have conflicting assignments for {gid}.')
                        sys.exit(-1)
                        
    def gtdb_type_species_of_genus(self, gtdb_genus):
        """Get genome assembled from type species of genus."""

        for gid in self.genomes:
            if not self.genomes[gid].is_gtdb_type_species():
                continue
                
            if self.genomes[gid].gtdb_taxa.genus == gtdb_genus:
                return gid
                
        return None
    
    def gtdb_sp_rep(self, gtdb_sp):
        """Get representative genome for GTDB species cluster."""

        for gid in self.genomes:
            if not self.genomes[gid].is_gtdb_sp_rep():
                continue
                
            if self.genomes[gid].gtdb_taxa.species == gtdb_sp:
                return gid
                
        self.logger.error(f'Failed to find representative of GTDB species for {gtdb_sp}.')
        sys.exit(-1)

    def ncbi_sp_effective_type_genomes(self):
        """Get effect type genomes for each NCBI species."""
        
        ncbi_sp_type_strain_genomes = defaultdict(set)
        for gid in self.genomes:
            if self.genomes[gid].is_effective_type_strain():
                ncbi_sp = self.genomes[gid].ncbi_taxa.species
                if ncbi_sp != 's__':
                    # yes, NCBI has genomes marked as assembled from type material
                    # that do not actually have a binomial species name
                    ncbi_sp_type_strain_genomes[ncbi_sp].add(gid)
                    
        return ncbi_sp_type_strain_genomes

    def sort_by_assembly_score(self):
        """Return genomes sorted by their assembly score."""
        
        for gid in sorted(self.genomes.keys(), 
                            key=lambda gid: self.genomes[gid].score_assembly(), 
                            reverse=True):
            yield gid

    def get_gid(self, idx):
        """Get ID of genome at specific index."""
        
        return list(self.genomes)[idx]
        
    def gtdb_type_strain_genomes(self):
        """Get genomes considered type strain of species by GTDB."""
        
        type_strain_gids = set()
        for gid, genome in self.genomes.items():
            if genome.is_gtdb_type_strain():
                type_strain_gids.add(gid)
                
        return type_strain_gids
        
    def get_ncbi_type_strain_genomes(self):
        """Get type strain genomes for NCBI species."""
        
        type_strain_genomes = defaultdict(set)
        for gid, genome in self.genomes.items():
            if genome.is_effective_type_strain():
                type_strain_genomes[genome.ncbi_taxa.species].add(gid)
                
        return type_strain_genomes
        
    def named_ncbi_species(self):
        """Get genomes in valid or effectively published, including Candidatus, species in NCBI taxonomy."""
        
        named_ncbi_sp = defaultdict(set)
        for gid in self.genomes:
            if not is_placeholder_taxon(self.genomes[gid].ncbi_taxa.species):
                named_ncbi_sp[self.genomes[gid].ncbi_taxa.species].add(gid)

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
                
                #*** Need to handle UBA genomes that can reference via their U_ ID when run through MASH
                self.genomic_files[accession] = genomic_file 
                
    def parse_untrustworthy_type_ledger(self, untrustworth_type_ledger):
        """Determine genomes that should be considered untrustworthy as type material."""
        
        untrustworthy_as_type = set()
        with open(untrustworth_type_ledger) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                
                untrustworthy_as_type.add(canonical_gid(tokens[0]))
            
        return untrustworthy_as_type

    def load_from_metadata_file(self, 
                                metadata_file,
                                species_exception_file=None,
                                genus_exception_file=None,
                                gtdb_type_strains_ledger=None,
                                create_sp_clusters=True,
                                uba_genome_file=None,
                                qc_passed_file=None,
                                ncbi_genbank_assembly_file=None,
                                untrustworthy_type_ledger=None,
                                gtdbtk_classify_file=None):
        """Create genome set from file(s)."""
        
        pass_qc_gids = set()
        if qc_passed_file:
            with open(qc_passed_file) as f:
                f.readline()
                for line in f:
                    line_split = line.strip().split('\t')
                    pass_qc_gids.add(line_split[0].strip())
            self.logger.info(f' ... identified {len(pass_qc_gids):,} genomes passing QC.')
                    
        valid_uba_ids = set()
        if uba_genome_file:
            with open(uba_genome_file) as f:
                for line in f:
                    line_split = line.strip().split('\t')
                    valid_uba_ids.add(line_split[0].strip())
            self.logger.info(f' ... identified {len(valid_uba_ids):,} UBA genomes to retain.')

        gtdb_type_strains = set()
        if gtdb_type_strains_ledger:
            with open(gtdb_type_strains_ledger) as f:
                f.readline()
                for line in f:
                    tokens = line.strip().split('\t')
                    gid = canonical_gid(tokens[0].strip())
                    gtdb_type_strains.add(gid)
            self.logger.info(f' ... identified {len(gtdb_type_strains):,} manually annotated as type strain genomes.')
                    
        excluded_from_refseq_note = {}
        if ncbi_genbank_assembly_file:
            excluded_from_refseq_note = exclude_from_refseq(ncbi_genbank_assembly_file)
            
        untrustworthy_as_type = set()
        if untrustworthy_type_ledger:
            untrustworthy_as_type = self.parse_untrustworthy_type_ledger(untrustworthy_type_ledger)
            self.logger.info(f' ... identified {len(untrustworthy_as_type):,} genomes annotated as untrustworthy as type.')
            
        gtdbtk_classifications = {}
        if gtdbtk_classify_file:
            gtdbtk_classifications = read_gtdbtk_classifications(gtdbtk_classify_file)
            self.logger.info(f' ... identified {len(gtdbtk_classifications):,} GTDB-Tk classifications.')

        with open(metadata_file, encoding='utf-8') as f:
            headers = f.readline().strip().split('\t')

            genome_index = headers.index('accession')

            gtdb_taxonomy_index = headers.index('gtdb_taxonomy')
            ncbi_taxonomy_index = headers.index('ncbi_taxonomy')
            ncbi_taxonomy_unfiltered_index = headers.index('ncbi_taxonomy_unfiltered')
            
            gtdb_type_index = headers.index('gtdb_type_designation')
            gtdb_type_sources_index = headers.index('gtdb_type_designation_sources')
            gtdb_type_species_of_genus_index = headers.index('gtdb_type_species_of_genus')
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
            gtdb_rep_index = headers.index('gtdb_representative')
            
            if 'lpsn_priority_year' in headers:
                # this information will be missing from the previous
                # GTDB metadata file as we strip this out due to 
                # concerns over republishing this information
                lpsn_priority_index = headers.index('lpsn_priority_year')
                dsmz_priority_index = headers.index('dsmz_priority_year')
                straininfo_priority_index = headers.index('straininfo_priority_year')

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
                            self.uba_user_id_map[uba_id] = gid
                            gid = uba_id
                        else:
                            continue # retain only valid UBA genomes
                    else:
                        continue # skip non-UBA user genomes
                        
                if pass_qc_gids and gid not in pass_qc_gids:
                    continue
                
                gtdb_taxonomy = Taxa(line_split[gtdb_taxonomy_index])
                if gid in gtdbtk_classifications:
                    gtdb_taxonomy = Taxa(';'.join(gtdbtk_classifications[gid]))
                
                ncbi_taxonomy = Taxa(line_split[ncbi_taxonomy_index])
                ncbi_taxonomy_unfiltered = Taxa(line_split[ncbi_taxonomy_unfiltered_index])
                
                gtdb_type = line_split[gtdb_type_index]
                gtdb_type_sources = line_split[gtdb_type_sources_index]
                if gid in gtdb_type_strains:
                    gtdb_type = 'type strain of species'
                    gtdb_type_sources = 'GTDB curator'
                gtdb_type_species_of_genus = line_split[gtdb_type_species_of_genus_index] == 't'
                
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
                
                gtdb_is_rep = line_split[gtdb_rep_index] == 't'
                gtdb_rid = canonical_gid(line_split[gtdb_genome_rep_index])
                if create_sp_clusters:
                    self.sp_clusters.update_sp_cluster(gtdb_rid, gid, gtdb_taxonomy.species)
                
                if 'lpsn_priority_year' in headers:
                    lpsn_priority_year = self._convert_int(line_split[lpsn_priority_index], Genome.NO_PRIORITY_YEAR)
                    dsmz_priority_year = self._convert_int(line_split[dsmz_priority_index], Genome.NO_PRIORITY_YEAR)
                    straininfo_priority_year = self._convert_int(line_split[straininfo_priority_index], Genome.NO_PRIORITY_YEAR)
                else:
                    lpsn_priority_year = Genome.NO_PRIORITY_YEAR
                    dsmz_priority_year = Genome.NO_PRIORITY_YEAR
                    straininfo_priority_year = Genome.NO_PRIORITY_YEAR

                self.genomes[gid] = Genome(gid,
                                            ncbi_accn,
                                            gtdb_rid,
                                            gtdb_is_rep,
                                            gtdb_taxonomy,
                                            ncbi_taxonomy,
                                            ncbi_taxonomy_unfiltered,
                                            gtdb_type,
                                            gtdb_type_sources,
                                            gtdb_type_species_of_genus,
                                            gid in untrustworthy_as_type,
                                            ncbi_type,
                                            ncbi_strain_identifiers,
                                            ncbi_asm_level,
                                            ncbi_genome_representation,
                                            ncbi_refseq_cat,
                                            ncbi_genome_cat,
                                            excluded_from_refseq_note.get(gid, ''),
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
                                            ncbi_spanned_gaps,
                                            lpsn_priority_year,
                                            dsmz_priority_year,
                                            straininfo_priority_year)
                                            
        self._apply_ncbi_taxonomy_ledgers(species_exception_file,
                                            genus_exception_file)
