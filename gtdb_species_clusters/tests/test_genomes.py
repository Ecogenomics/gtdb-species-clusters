# Test methods in Genomes class

import os
from collections import defaultdict

from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.update_qc_genomes import QcCriteria
from gtdb_species_clusters import defaults as Defaults

# GTDB R237 metadata for a genome with complete metadata (G003983815),
# a genome missing CheckM, coding density, and other values (G001578415),
# and a genome missing CheckM v2 and genome size values (G964261755)
METADATA_FILE = os.path.join(os.path.dirname(__file__), 'data', 'metadata_none_values.tsv')

qc_criteria = QcCriteria(min_comp=Defaults.QC_MIN_COMP,
                         max_cont=Defaults.QC_MAX_CONT,
                         min_quality=Defaults.QC_MIN_QUALITY,
                         sh_exception=Defaults.QC_SH_EXCEPTION,
                         min_perc_markers=Defaults.QC_MIN_PERC_MARKERS,
                         max_contigs=Defaults.QC_MAX_CONTIGS,
                         min_N50=Defaults.QC_MIN_N50,
                         max_ambiguous=Defaults.QC_MAX_AMBIGUOUS,
                         abnormal_genome_size_mb=Defaults.QC_ABNORMAL_GENOME_SIZE_MB,
                         abnormal_genome_size_coding_density=Defaults.QC_ABNORMAL_GENOME_SIZE_CODING_DENSITY,
                         frameshift_cm2_comp=Defaults.QC_FRAMESHIFT_CM2_COMP,
                         frameshift_coding_density=Defaults.QC_FRAMESHIFT_CODING_DENSITY,
                         frameshift_contigs=Defaults.QC_FRAMESHIFT_CONTIGS,
                         ncbi_frameshift_cm2_comp=Defaults.QC_NCBI_FRAMESHIFT_CM2_COMP,
                         ncbi_frameshift_coding_density=Defaults.QC_NCBI_FRAMESHIFT_CODING_DENSITY,
                         max_cm_xor_contigs=Defaults.MAX_CHECKM_XOR_CONTIGS)


class TestGenomes:
    """Test Genomes class."""

    def test_load_metadata_none_values(self):
        """Test default values assigned to missing ('none') metadata values."""

        genomes = Genomes()
        genomes.load_from_metadata_file(METADATA_FILE, create_sp_clusters=False)
        assert len(genomes) == 3

        # genome with complete metadata retains its values
        g = genomes['G003983815']
        assert g.cm1_comp > 0 and g.cm1_cont < 100
        assert g.cm2_comp > 0 and g.cm2_cont < 100
        assert g.length > 0 and g.coding_density > 0

        # missing CheckM v1 and v2 completeness and contamination are set to
        # the worst possible values, with missing coding density set to 0
        g = genomes['G001578415']
        assert (g.cm1_comp, g.cm1_cont, g.cm1_strain_heterogeneity_100) == (0.0, 100.0, 0.0)
        assert (g.cm2_comp, g.cm2_cont) == (0.0, 100.0)
        assert g.coding_density == 0.0

        # missing CheckM v2 estimates and genome size
        g = genomes['G964261755']
        assert g.cm1_comp == 100.0
        assert (g.cm2_comp, g.cm2_cont) == (0.0, 100.0)
        assert g.length == 0

    def test_qc_metadata_none_values(self):
        """Test genomes with missing CheckM estimates fail QC."""

        genomes = Genomes()
        genomes.load_from_metadata_file(METADATA_FILE, create_sp_clusters=False)

        for gid in ['G001578415', 'G964261755']:
            failed_tests = defaultdict(int)
            assert not genomes[gid].pass_qc(100, qc_criteria, failed_tests)
            assert 'cm2_comp' in failed_tests

    def test_type_strain_ledger_untrustworthy(self, tmp_path, caplog, monkeypatch):
        """Test genomes in the GTDB type strains ledger and untrustworthy file are not type material."""

        # untrustworthy file is reported relative to the current directory
        monkeypatch.chdir(tmp_path)
        (tmp_path / '5_u_resolve_types').mkdir()

        gtdb_type_strains_ledger = tmp_path / 'gtdb_type_strains.tsv'
        gtdb_type_strains_ledger.write_text('Genome ID\tProposed GTDB species name\n'
                                            'G003983815\tTestus ledgerensis\n')

        untrustworthy_file = tmp_path / '5_u_resolve_types' / 'untrustworthy_type_material.tsv'
        untrustworthy_file.write_text('Genome ID\tNCBI species\tGTDB species\tLTP species\tReason for declaring untrustworthy\n'
                                      'G003983815\ts__Testus ncbiensis\ts__Testus ledgerensis\t\t93.50% ANI to other type strain genomes\n')

        genomes = Genomes()
        with caplog.at_level('WARNING'):
            genomes.load_from_metadata_file(METADATA_FILE,
                                            gtdb_type_strains_ledger=str(gtdb_type_strains_ledger),
                                            create_sp_clusters=False,
                                            untrustworthy_type_ledger=str(untrustworthy_file))

        assert [r.getMessage() for r in caplog.records if r.levelname == 'WARNING'] == [
            f'Genome G003983815 is in the GTDB type strains ledger, but is untrustworthy as type material '
            'according to 5_u_resolve_types/untrustworthy_type_material.tsv; it will not be considered type material.']

        # species name from the ledger is retained, but genome is not considered type material
        g = genomes['G003983815']
        assert g.gtdb_taxa.species == 's__Testus ledgerensis'
        assert g.gtdb_type_designation == 'type strain of species'
        assert not g.is_gtdb_type_strain()
        assert not g.is_effective_type_strain()

    def test_seqcode_taxa(self, tmp_path, caplog):
        """Test SeqCode taxonomy is set from the SeqCode table and checked against the ledger and metadata."""

        # GCF accession in the table resolves to the GCA accession in the metadata;
        # G003983815 is in the ledger with its SeqCode species (redundant),
        # G001578415 is in the ledger with a different species (conflict)
        seqcode_file = tmp_path / 'seqcode_table.tsv'
        seqcode_file.write_text(
            'seqcode_type_material_accn\tseqcode_name\tseqcode_rank\tseqcode_species_status\tseqcode_priority_date'
            '\tseqcode_classification\tseqcode_type_species_of_genus\n'
            'GCF_001578415.1\tPseudomonas nova\tspecies\tValid (SeqCode)\t2025-03-01T00:00:00.000Z'
            '\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__\tTrue\n'
            'GCA_003983815.1\tCandidatus Escherichia nova\tspecies\tValid (SeqCode)\t2024-01-01T00:00:00.000Z'
            '\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__\tFalse\n'
            'GCA_000000001.1\tGenusia nova\tgenus\tValid (SeqCode)\t2024-01-01T00:00:00.000Z\td__Bacteria;p__;c__;o__;f__;g__Genusia;s__\tTrue\n')

        gtdb_type_strains_ledger = tmp_path / 'gtdb_type_strains.tsv'
        gtdb_type_strains_ledger.write_text('Genome ID\tProposed GTDB species name\n'
                                            'G003983815\tEscherichia nova\n'
                                            'G001578415\tPseudomonas_A nova\n')

        genomes = Genomes()
        with caplog.at_level('WARNING'):
            genomes.load_from_metadata_file(METADATA_FILE,
                                            gtdb_type_strains_ledger=str(gtdb_type_strains_ledger),
                                            create_sp_clusters=False,
                                            seqcode_file=str(seqcode_file))

        warnings = [r.getMessage() for r in caplog.records if r.levelname == 'WARNING']
        assert warnings == [
            'Skipped 1 rows in SeqCode table that are not validly published species.',
            'Genome G003983815 has the SeqCode name Escherichia nova in the SeqCode table, but no SeqCode name in metadata_none_values.tsv.',
            'Genome G003983815 is in the GTDB type strains ledger as s__Escherichia nova, which is its SeqCode species; the ledger entry is redundant.',
            'Genome G001578415 has the SeqCode name Pseudomonas nova in the SeqCode table, but no SeqCode name in metadata_none_values.tsv.',
            'Genome G001578415 is in the GTDB type strains ledger as s__Pseudomonas_A nova, but is s__Pseudomonas nova in SeqCode; the ledger entry takes precedence.']

        g = genomes['G001578415']
        assert g.seqcode_taxa.species == 's__Pseudomonas nova'
        assert g.seqcode_taxa.genus == 'g__Pseudomonas'
        assert g.seqcode_taxa.family == 'f__Pseudomonadaceae'
        assert g.seqcode_priority_year == 2025
        assert g.ref_taxa_source == 'SeqCode'
        assert g.ncbi_taxa.species == 's__Pseudomonas aeruginosa'
        assert g.type_strains_ledger_sp == 's__Pseudomonas_A nova'
        assert g.type_species() == 's__Pseudomonas_A nova'
        assert g.gtdb_taxa.species == 's__Pseudomonas_A nova'

        g = genomes['G003983815']
        assert g.seqcode_taxa.species == 's__Escherichia nova'
        assert g.seqcode_priority_year == 2024
        assert g.type_species() == 's__Escherichia nova'

        g = genomes['G964261755']
        assert g.seqcode_taxa is None
        assert g.ref_taxa_source == 'NCBI'
        assert g.type_strains_ledger_sp == ''
        assert g.type_species() == 's__'

        assert genomes.named_ref_species() == {'s__Pseudomonas nova': {'G001578415'},
                                               's__Escherichia nova': {'G003983815'}}
        assert genomes.type_strain_genomes_by_species() == {'s__Pseudomonas_A nova': {'G001578415'},
                                                            's__Escherichia nova': {'G003983815'}}

    def test_seqcode_taxa_not_loaded(self):
        """Test genomes have no SeqCode taxonomy when no SeqCode table is given."""

        genomes = Genomes()
        genomes.load_from_metadata_file(METADATA_FILE, create_sp_clusters=False)
        assert all(genomes[gid].seqcode_taxa is None for gid in genomes)
        assert all(genomes[gid].ref_taxa is genomes[gid].ncbi_taxa for gid in genomes)
