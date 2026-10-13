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
