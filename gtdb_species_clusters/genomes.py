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

from gtdblib.util.bio.accession import canonical_gid

from gtdb_species_clusters.genome import Genome
from gtdb_species_clusters.taxa import Taxa
from gtdb_species_clusters.species_clusters import SpeciesClusters
from gtdb_species_clusters.genome_utils import (exclude_from_refseq,
                                                parse_ncbi_env_bioproject_ledger,
                                                parse_ncbi_bioproject,
                                                parse_seqcode_table,
                                                read_gtdbtk_classifications,
                                                open_file,
                                                MultiFileReader)
from gtdb_species_clusters.taxon_utils import is_placeholder_taxon


class Genomes(object):
    """A collection of genomes."""

    def __init__(self):
        """Initialization."""

        self.genomes = {}
        self.sp_clusters = SpeciesClusters()

        # this could be removed, but the skani interface
        # currently requires data to be passed in as a dictionary
        self.genomic_files = {}

        self.full_gid = {}  # translate from canonical GID to NCBI accession

        self.log = logging.getLogger('rich')

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
                    line_split = [token.strip()
                                  for token in line.strip().split('\t')]
                    gid = canonical_gid(line_split[0].strip())

                    sp = line_split[1].strip().replace('Candidatus ', '')
                    if gid not in self.genomes:
                        self.log.warning(
                            f'Genome {gid} in species exception list not defined in genome set.')
                        continue

                    if not sp.startswith('s__'):
                        sp = 's__' + sp

                    self.genomes[gid].ncbi_taxa.species = sp
                    species_updates[gid] = sp

        if genus_exception_file:
            with open(genus_exception_file, encoding='utf-8') as f:
                f.readline()
                for line in f:
                    line_split = [token.strip()
                                  for token in line.strip().split('\t')]
                    gid = canonical_gid(line_split[0].strip())
                    genus = line_split[1].strip()
                    if gid not in self.genomes:
                        self.log.warning(
                            f'Genome {gid} in genus exception list not defined in genome set.')
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
                        self.log.error(
                            f'Species and genus ledgers have conflicting assignments for {gid}.')
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

        self.log.error(
            f'Failed to find representative of GTDB species for {gtdb_sp}.')
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

    def get_ref_type_strain_genomes(self):
        """Get type strain genomes for each reference (SeqCode or NCBI) species."""

        type_strain_genomes = defaultdict(set)
        for gid, genome in self.genomes.items():
            if genome.is_effective_type_strain():
                type_strain_genomes[genome.ref_taxa.species].add(gid)

        return type_strain_genomes

    def named_ncbi_species(self):
        """Get genomes with valid or effectively published, including Candidatus, trusted species names in NCBI taxonomy."""

        named_ncbi_sp = defaultdict(set)
        for gid in self.genomes:
            if self.genomes[gid].ncbi_untrustworthy_sp:
                continue

            if not is_placeholder_taxon(self.genomes[gid].ncbi_taxa.species):
                named_ncbi_sp[self.genomes[gid].ncbi_taxa.species].add(gid)

        return named_ncbi_sp

    def named_ref_species(self):
        """Get genomes with valid or effectively published, including Candidatus, trusted species names in the reference (SeqCode or NCBI) taxonomy."""

        named_ref_sp = defaultdict(set)
        for gid in self.genomes:
            if self.genomes[gid].ncbi_untrustworthy_sp:
                continue

            if not is_placeholder_taxon(self.genomes[gid].ref_taxa.species):
                named_ref_sp[self.genomes[gid].ref_taxa.species].add(gid)

        return named_ref_sp

    def type_strain_genomes_by_species(self):
        """Get effective type strain genomes grouped by the species they are type material of (see Genome.type_species)."""

        type_strain_genomes = defaultdict(set)
        for gid, genome in self.genomes.items():
            if genome.is_effective_type_strain():
                type_sp = genome.type_species()
                if type_sp != 's__':
                    type_strain_genomes[type_sp].add(gid)

        return type_strain_genomes

    def load_genomic_file_paths(self, genome_path_file):
        """Determine path to genomic FASTA file for each genome."""

        with open_file(genome_path_file) as f:
            for line in f:
                line_split = line.strip().split('\t')

                gid = line_split[0]
                gid = canonical_gid(gid)
                if gid in self.genomes:
                    genome_path = line_split[1]
                    accession = os.path.basename(os.path.normpath(genome_path))
                    genomic_file = os.path.join(genome_path, accession + '_genomic.fna.gz')
                    self.genomes[gid].genomic_file = genomic_file
                    self.genomic_files[gid] = genomic_file

    def parse_untrustworthy_type_ledger(self, untrustworth_type_ledger):
        """Determine genomes that should be considered untrustworthy as type material."""

        untrustworthy_as_type = set()
        with open(untrustworth_type_ledger) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                untrustworthy_as_type.add(canonical_gid(tokens[0]))

        return untrustworthy_as_type

    def parse_ncbi_untrustworthy_sp_ledger(self, ncbi_untrustworthy_sp_ledger):
        """Determine genomes that should be considered as having untrustworthy NCBI species assignments."""

        untrustworthy_ncbi_sp = set()
        with open(ncbi_untrustworthy_sp_ledger) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                untrustworthy_ncbi_sp.add(canonical_gid(tokens[0]))

        return untrustworthy_ncbi_sp

    def set_gtdbtk_classification(self, gtdbtk_classify_file, prev_genomes):
        """Update classification of genomes based on GTDB-Tk results."""

        # get species names in previous GTDB release
        prev_ncbi_sp = set()
        prev_gtdb_sp = set()
        for gid in prev_genomes:
            prev_ncbi_sp.add(prev_genomes[gid].ncbi_taxa.species)
            prev_gtdb_sp.add(prev_genomes[gid].gtdb_taxa.species)

        # set new genomes to the predicted GTDB-Tk classification, but
        # change genus and species classifications of a genome if it has
        # a previously unseen NCBI species assignment. This is a problematic
        # case for curation as GTDB-Tk is unaware of these assignments.
        # A common example is a new genome being the most basal
        # member of a genus, this genome being classified to this genus by
        # GTDB-Tk, but this genome being from a newly proposed genera which
        # should be favored in order to have the GTDB reflect the opinion
        # of the community. New NCBI taxa above the rank of genus are not considered
        # as these are relatively uncommon and are picked up for manual curation
        # using curation trees that specifically highlight genomes with previously
        # unseen NCBI taxon names (see update_curation_trees)
        gtdbtk_classifications = read_gtdbtk_classifications(
            gtdbtk_classify_file)
        num_updated = 0
        num_ncbi_sp = 0
        for gid in self.genomes:
            if gid in gtdbtk_classifications:
                num_updated += 1

                gtdbtk_taxa = Taxa(';'.join(gtdbtk_classifications[gid]))
                self.genomes[gid].gtdb_taxa.update_taxa(gtdbtk_taxa)

                ncbi_sp = self.genomes[gid].ncbi_taxa.species
                if (ncbi_sp == 's__'
                    or ncbi_sp in prev_ncbi_sp
                        or ncbi_sp in prev_gtdb_sp):
                    continue

                self.genomes[gid].gtdb_taxa.set_taxa(
                    5, self.genomes[gid].ncbi_taxa.genus)
                self.genomes[gid].gtdb_taxa.set_taxa(6, ncbi_sp)
                num_ncbi_sp += 1

        return num_updated, num_ncbi_sp

    def set_prev_gtdb_classifications(self, prev_genomes):
        """Set genomes to GTDB assignments in previous release."""

        # set current genomes to have same GTDB assignments as in previous
        # GTDB release. This is necessary since genomes may have different
        # NCBI accession numbers between releases and thus the previous GTDB
        # taxonomy will not be reflected in the latest GTDB database. The
        # exception is if a genome has changed domains, in which case the
        # previous assignment is invalid.
        self.log.info(
            'Setting GTDB taxonomy of genomes in current genome set.')

        update_count = 0
        conflicting_domain_count = 0
        for prev_gid in prev_genomes:
            if prev_gid in self.genomes:
                if prev_genomes[prev_gid].gtdb_taxa != self.genomes[prev_gid].gtdb_taxa:
                    if prev_genomes[prev_gid].gtdb_taxa.domain == self.genomes[prev_gid].gtdb_taxa.domain:
                        # verify updating genomes without assigned GTDB taxa below domain
                        assert self.genomes[prev_gid].gtdb_taxa.phylum == 'p__'

                        update_count += 1
                        self.genomes[prev_gid].gtdb_taxa.update_taxa(
                            prev_genomes[prev_gid].gtdb_taxa)
                    else:
                        conflicting_domain_count += 1

        self.log.info(f' - updated {update_count:,} genomes')
        self.log.info(
            f' - identified {conflicting_domain_count:,} genomes with conflicting domain assignments')

    def load_from_metadata_file(self,
                                metadata_file,
                                species_exception_file=None,
                                genus_exception_file=None,
                                gtdb_type_strains_ledger=None,
                                create_sp_clusters=True,
                                qc_passed_file=None,
                                ncbi_genbank_assembly_file=None,
                                untrustworthy_type_ledger=None,
                                ncbi_untrustworthy_sp_ledger=None,
                                ncbi_env_bioproject_ledger=None,
                                seqcode_file=None):
        """Create genome set from file(s).

        The metadata_file and ncbi_genbank_assembly_file arguments may each be a
        single file or a list of files (e.g., separate archaeal and bacterial files),
        and files may be gzip compressed.

        The seqcode_file gives the type genome of each validly published SeqCode species.
        It sets the SeqCode taxonomy (Genome.seqcode_taxa) of these genomes, which takes
        priority over the NCBI taxonomy (see Genome.ref_taxa). The same file is used for
        the previous and current genome sets so names are comparable between releases;
        only the reference taxonomy of previous genomes is affected, never their GTDB taxonomy.
        """

        metadata_name = os.path.basename(metadata_file if isinstance(metadata_file, str) else metadata_file[0])

        seqcode = {}
        if seqcode_file:
            seqcode = parse_seqcode_table(seqcode_file)
            self.log.info(
                f' - identified {len(seqcode):,} type genomes of validly published SeqCode species')

        pass_qc_gids = set()
        if qc_passed_file:
            with open(qc_passed_file) as f:
                f.readline()
                for line in f:
                    line_split = line.strip().split('\t')
                    pass_qc_gids.add(line_split[0].strip())
            self.log.info(
                f' - identified {len(pass_qc_gids):,} genomes passing QC')

        gtdb_type_strains = {}
        if gtdb_type_strains_ledger:
            with open(gtdb_type_strains_ledger) as f:
                f.readline()
                for line in f:
                    tokens = line.strip().split('\t')
                    gid = canonical_gid(tokens[0].strip())
                    sp = tokens[1].strip()
                    if not sp.startswith('s__'):
                        sp = 's__' + sp
                    gtdb_type_strains[gid] = sp
            self.log.info(
                f' - identified {len(gtdb_type_strains):,} manually annotated as type strain genomes')

        untrustworthy_as_type = set()
        if untrustworthy_type_ledger:
            untrustworthy_as_type = self.parse_untrustworthy_type_ledger(
                untrustworthy_type_ledger)
            self.log.info(
                f' - identified {len(untrustworthy_as_type):,} genomes annotated as untrustworthy as type by GTDB')

            # the untrustworthy file may be a curated ledger or the output of
            # u_resolve_types so is reported explicitly (relative to the current
            # directory); genomes considered untrustworthy are never considered type material
            untrustworthy_type_file = os.path.relpath(untrustworthy_type_ledger)
            for gid in gtdb_type_strains:
                if gid in untrustworthy_as_type:
                    self.log.warning(f"Genome {gid} is in the GTDB type strains ledger, but is untrustworthy as type material "
                                     f"according to {untrustworthy_type_file}; it will not be considered type material.")

        excluded_from_refseq_note = {}
        ncbi_bioproject = {}
        if ncbi_genbank_assembly_file:
            ncbi_bioproject = parse_ncbi_bioproject(ncbi_genbank_assembly_file)
            excluded_from_refseq_note = exclude_from_refseq(
                ncbi_genbank_assembly_file)

        # only BioProjects determined to contain environmental genomes are
        # relevant here; BioProjects of isolates are only used by u_rep_actions
        ncbi_env_bioproject = set()
        if ncbi_env_bioproject_ledger:
            ncbi_env_bioproject, _isolate_bioprojects = parse_ncbi_env_bioproject_ledger(
                ncbi_env_bioproject_ledger)

        untrustworthy_ncbi_sp = set()
        if ncbi_untrustworthy_sp_ledger:
            untrustworthy_ncbi_sp = self.parse_ncbi_untrustworthy_sp_ledger(
                ncbi_untrustworthy_sp_ledger)
            self.log.info(
                f' - identified {len(untrustworthy_ncbi_sp):,} genomes annotated as having untrustworthy NCBI species assignments')

        with MultiFileReader(metadata_file) as f:
            headers = f.readline().strip().split('\t')

            genome_index = headers.index('accession')

            gtdb_taxonomy_index = headers.index('gtdb_taxonomy')
            ncbi_taxonomy_index = headers.index('ncbi_taxonomy')
            ncbi_taxonomy_unfiltered_index = headers.index(
                'ncbi_taxonomy_unfiltered')

            if 'gtdb_type_designation' in headers:
                gtdb_type_index = headers.index('gtdb_type_designation')
            elif 'gtdb_type_designation_ncbi_taxa' in headers:
                # Note: name of field changed in R214
                gtdb_type_index = headers.index('gtdb_type_designation_ncbi_taxa')
            else:
                self.log.error("Unable to find GTDB type designation field.")
                sys.exit(1)

            if 'gtdb_type_designation_sources' in headers:
                gtdb_type_sources_index = headers.index('gtdb_type_designation_sources')
            elif 'gtdb_type_designation_ncbi_taxa_sources' in headers:
                # Note: name of field changed in R214
                gtdb_type_sources_index = headers.index('gtdb_type_designation_ncbi_taxa_sources')
            else:
                self.log.error("Unable to find GTDB type designation sources field.")
                sys.exit(1)

            gtdb_type_species_of_genus_index = headers.index(
                'gtdb_type_species_of_genus')
            ncbi_strain_identifiers_index = headers.index(
                'ncbi_strain_identifiers')
            ncbi_type_index = headers.index('ncbi_type_material_designation')
            ncbi_asm_level_index = headers.index('ncbi_assembly_level')
            ncbi_genome_representation_index = headers.index(
                'ncbi_genome_representation')
            ncbi_refseq_cat_index = headers.index('ncbi_refseq_category')
            ncbi_genome_cat_index = headers.index('ncbi_genome_category')

            cm1_comp_index = headers.index('checkm_completeness')
            cm1_cont_index = headers.index('checkm_contamination')
            cm1_sh_100_index = None
            if 'checkm_strain_heterogeneity_100' in headers:
                cm1_sh_100_index = headers.index('checkm_strain_heterogeneity_100')

            cm2_comp_index = headers.index('checkm2_completeness')
            cm2_cont_index = headers.index('checkm2_contamination')

            gs_index = headers.index('genome_size')
            contig_count_index = headers.index('contig_count')
            n50_index = headers.index('n50_contigs')
            scaffold_count_index = headers.index('scaffold_count')
            ambiguous_bases_index = headers.index('ambiguous_bases')
            total_gap_len_index = headers.index('total_gap_length')
            coding_density_index = headers.index('coding_density')
            ssu_count_index = headers.index('ssu_count')
            ssu_length_index = headers.index('ssu_length')
            ncbi_molecule_count_index = headers.index('ncbi_molecule_count')
            ncbi_unspanned_gaps_index = headers.index('ncbi_unspanned_gaps')
            ncbi_spanned_gaps_index = headers.index('ncbi_spanned_gaps')

            gtdb_genome_rep_index = headers.index('gtdb_genome_representative')
            gtdb_rep_index = headers.index('gtdb_representative')

            lpsn_priority_index = None
            if 'lpsn_priority_year' in headers:
                # this information will be missing from the previous
                # GTDB metadata file as we strip this out due to
                # concerns over republishing this information
                lpsn_priority_index = headers.index('lpsn_priority_year')

            seqcode_name_index = None
            if seqcode and 'seqcode_name' in headers:
                # the SeqCode table is the source of SeqCode names, but should
                # agree with the SeqCode fields in the current GTDB metadata
                seqcode_name_index = headers.index('seqcode_name')

            num_seqcode = 0
            num_seqcode_ncbi_sp = 0
            num_seqcode_ncbi_sp_diff = 0
            for line in f:
                tokens = line.strip().split('\t')

                ncbi_accn = tokens[genome_index]
                gid = canonical_gid(ncbi_accn)
                self.full_gid[gid] = ncbi_accn

                if pass_qc_gids and gid not in pass_qc_gids:
                    continue

                gtdb_taxonomy = Taxa(tokens[gtdb_taxonomy_index])
                ncbi_taxonomy = Taxa(tokens[ncbi_taxonomy_index])
                ncbi_taxonomy_unfiltered = Taxa(
                    tokens[ncbi_taxonomy_unfiltered_index], filtered=False)

                seqcode_taxonomy = None
                seqcode_priority_year = Genome.NO_PRIORITY_YEAR
                if gid in seqcode:
                    seqcode_taxonomy = Taxa(seqcode[gid].classification)
                    seqcode_taxonomy.species = 's__' + seqcode[gid].name
                    seqcode_priority_year = seqcode[gid].priority_year

                    num_seqcode += 1
                    if ncbi_taxonomy.species != 's__':
                        num_seqcode_ncbi_sp += 1
                        if ncbi_taxonomy.species != seqcode_taxonomy.species:
                            num_seqcode_ncbi_sp_diff += 1

                if seqcode_name_index is not None:
                    metadata_seqcode_name = tokens[seqcode_name_index]
                    if metadata_seqcode_name not in ['', 'none']:
                        if gid not in seqcode:
                            self.log.warning(f'Genome {gid} has the SeqCode name {metadata_seqcode_name} in {metadata_name}, '
                                             'but is not in the SeqCode table.')
                        elif metadata_seqcode_name.replace('Candidatus ', '') != seqcode[gid].name:
                            self.log.warning(f'Genome {gid} has the SeqCode name {metadata_seqcode_name} in {metadata_name}, '
                                             f'but {seqcode[gid].name} in the SeqCode table.')
                    elif gid in seqcode:
                        self.log.warning(f'Genome {gid} has the SeqCode name {seqcode[gid].name} in the SeqCode table, '
                                         f'but no SeqCode name in {metadata_name}.')

                gtdb_type = tokens[gtdb_type_index]
                gtdb_type_sources = tokens[gtdb_type_sources_index]
                type_strains_ledger_sp = ''
                if gid in gtdb_type_strains:
                    type_strains_ledger_sp = gtdb_type_strains[gid]

                    # the ledger takes precedence over the SeqCode and NCBI names, so report
                    # ledger entries that conflict with SeqCode or that change nothing; an entry
                    # that sets a different GTDB species than the metadata is never redundant
                    same_gtdb_sp = gtdb_taxonomy.species == type_strains_ledger_sp
                    if seqcode_taxonomy is not None:
                        if seqcode_taxonomy.species != type_strains_ledger_sp:
                            self.log.warning(f'Genome {gid} is in the GTDB type strains ledger as {type_strains_ledger_sp}, '
                                             f'but is {seqcode_taxonomy.species} in SeqCode; the ledger entry takes precedence.')
                        elif same_gtdb_sp:
                            self.log.warning(f'Genome {gid} is in the GTDB type strains ledger as {type_strains_ledger_sp}, '
                                             'which is its SeqCode and GTDB species; the ledger entry is redundant.')
                    elif (gtdb_type in Genome.GTDB_TYPE_SPECIES
                            and ncbi_taxonomy.species == type_strains_ledger_sp
                            and same_gtdb_sp):
                        self.log.warning(f'Genome {gid} is in the GTDB type strains ledger as {type_strains_ledger_sp}, '
                                         f'which is its NCBI and GTDB species and it is already a type strain in {metadata_name}; '
                                         'the ledger entry is redundant.')

                    gtdb_type = 'type strain of species'
                    gtdb_type_sources = 'GTDB curator'

                    # force specified species name into GTDB taxonomy string
                    gtdb_taxonomy.species = gtdb_type_strains[gid]

                gtdb_type_species_of_genus = tokens[gtdb_type_species_of_genus_index] == 't'

                ncbi_type = tokens[ncbi_type_index]
                ncbi_strain_identifiers = tokens[ncbi_strain_identifiers_index]
                ncbi_asm_level = tokens[ncbi_asm_level_index]
                ncbi_genome_representation = tokens[ncbi_genome_representation_index]
                ncbi_refseq_cat = tokens[ncbi_refseq_cat_index]
                ncbi_genome_cat = tokens[ncbi_genome_cat_index]

                if ncbi_bioproject.get(gid, None) in ncbi_env_bioproject:  # ***
                    # HACK to force genomes from MAG mining projects
                    # to be indicated as MAGs which are currently
                    # not correctly annotated at NCBI
                    ncbi_genome_cat = 'derived from environmental source'

                # missing CheckM estimates are set to the worst possible values
                # (0% completeness, 100% contamination) so genomes fail CheckM QC
                cm1_comp = self._convert_float(tokens[cm1_comp_index], 0.0)
                cm1_cont = self._convert_float(tokens[cm1_cont_index], 100.0)
                cm1_sh_100 = 0
                if cm1_sh_100_index:
                    cm1_sh_100 = self._convert_float(tokens[cm1_sh_100_index], 0.0)

                cm2_comp = self._convert_float(tokens[cm2_comp_index], 0.0)
                cm2_cont = self._convert_float(tokens[cm2_cont_index], 100.0)

                gs = self._convert_int(tokens[gs_index], 0)
                contig_count = int(tokens[contig_count_index])
                n50 = int(tokens[n50_index])
                scaffold_count = int(tokens[scaffold_count_index])
                ambiguous_bases = int(tokens[ambiguous_bases_index])
                total_gap_len = int(tokens[total_gap_len_index])
                coding_density = self._convert_float(tokens[coding_density_index], 0.0)
                ssu_count = int(tokens[ssu_count_index])
                ssu_length = self._convert_int(tokens[ssu_length_index])
                ncbi_molecule_count = self._convert_int(
                    tokens[ncbi_molecule_count_index])
                ncbi_unspanned_gaps = self._convert_int(
                    tokens[ncbi_unspanned_gaps_index])
                ncbi_spanned_gaps = self._convert_int(
                    tokens[ncbi_spanned_gaps_index])

                gtdb_is_rep = tokens[gtdb_rep_index] == 't'
                gtdb_rid = canonical_gid(tokens[gtdb_genome_rep_index])
                if create_sp_clusters and gtdb_rid != 'none':
                    self.sp_clusters.update_sp_cluster(
                        gtdb_rid, gid, gtdb_taxonomy.species)

                lpsn_priority_year = Genome.NO_PRIORITY_YEAR
                if lpsn_priority_index is not None:
                    lpsn_priority_year = self._convert_int(
                        tokens[lpsn_priority_index], Genome.NO_PRIORITY_YEAR)

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
                                           gid in untrustworthy_ncbi_sp,
                                           ncbi_type,
                                           ncbi_strain_identifiers,
                                           ncbi_asm_level,
                                           ncbi_genome_representation,
                                           ncbi_refseq_cat,
                                           ncbi_genome_cat,
                                           excluded_from_refseq_note.get(gid, ''),
                                           cm1_comp,
                                           cm1_cont,
                                           cm1_sh_100,
                                           cm2_comp,
                                           cm2_cont,
                                           gs,
                                           contig_count,
                                           n50,
                                           scaffold_count,
                                           ambiguous_bases,
                                           total_gap_len,
                                           coding_density,
                                           ssu_count,
                                           ssu_length,
                                           ncbi_molecule_count,
                                           ncbi_unspanned_gaps,
                                           ncbi_spanned_gaps,
                                           lpsn_priority_year,
                                           seqcode_taxa=seqcode_taxonomy,
                                           seqcode_priority_year=seqcode_priority_year,
                                           type_strains_ledger_sp=type_strains_ledger_sp)

        if seqcode:
            self.log.info(f' - identified {num_seqcode:,} genomes with a SeqCode taxonomy, of which {num_seqcode_ncbi_sp:,} '
                          f'also have an NCBI species and {num_seqcode_ncbi_sp_diff:,} have an NCBI species that differs')

        self._apply_ncbi_taxonomy_ledgers(species_exception_file,
                                          genus_exception_file)
