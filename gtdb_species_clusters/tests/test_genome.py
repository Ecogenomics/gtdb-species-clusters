# Test methods in Genome class

from collections import defaultdict
import copy

from gtdb_species_clusters.genome import Genome
from gtdb_species_clusters.taxa import Taxa
from gtdb_species_clusters.update_qc_genomes import QcCriteria
from gtdb_species_clusters import defaults as Defaults

test_genome = Genome(gid="G012345678",
                     ncbi_accn="GCA_012345678.1",
                     gtdb_rid="G012345678",
                     gtdb_is_rep=True,
                     gtdb_taxa=Taxa(
                         "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri"),
                     ncbi_taxa=Taxa(
                         "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri"),
                     ncbi_unfiltered_taxa=Taxa(
                         "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri"),
                     gtdb_type_designation="type strain of species",
                     gtdb_type_designation_sources="LPSN",
                     gtdb_type_species_of_genus=True,
                     gtdb_untrustworthy_as_type=False,
                     ncbi_untrustworthy_sp=False,
                     ncbi_type_material='assembly from type material',
                     ncbi_strain_identifiers="K-12",
                     ncbi_assembly_level="Complete genome",
                     ncbi_genome_representation="Full",
                     ncbi_refseq_category="Reference genome",
                     ncbi_genome_category=None,
                     excluded_from_refseq_note="",
                     cm1_comp=100,
                     cm1_cont=0,
                     cm1_strain_heterogeneity_100=0,
                     cm2_comp=100,
                     cm2_cont=0,
                     length=3500000,
                     contig_count=1,
                     contig_n50=3500000,
                     scaffold_count=1,
                     ambiguous_bases=0,
                     total_gap_len=0,
                     coding_density=88.0,
                     ssu_count=1,
                     ssu_length=1600,
                     ncbi_molecule_count=1,
                     ncbi_unspanned_gaps=0,
                     ncbi_spanned_gaps=0,
                     lpsn_priority_year=1900)

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


class TestGenome:
    """Test Genome class."""

    def test_ncbi_subspecies(self):
        """Test identification of genomes classified as a subspecies at NCBI."""

        g = copy.copy(test_genome)

        # no subspecies defined
        g.ncbi_unfiltered_taxa = Taxa(
            "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri")
        assert not g.is_ncbi_subspecies()

        # subspecies defined
        g.ncbi_unfiltered_taxa = Taxa(
            "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri subsp. testus")
        assert g.is_ncbi_subspecies()

        # a variety or str is not a subspecies
        g.ncbi_unfiltered_taxa = Taxa(
            "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri var. testus")
        assert not g.is_ncbi_subspecies()

        g.ncbi_unfiltered_taxa = Taxa(
            "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri str testus")
        assert not g.is_ncbi_subspecies()

        # a variety is not a subspecies
        g.ncbi_unfiltered_taxa = Taxa(
            "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri var. testus")
        assert not g.is_ncbi_subspecies()

        # genome is not from a subspecies if the subspecies and specific name are the same since it should be considered from the species itself;
        # this is critical for identifying a genome as the type strain of the species and not the type strain of a subspecies
        g.ncbi_unfiltered_taxa = Taxa(
            "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri subsp. flexneri")
        assert not g.is_ncbi_subspecies()

    def test_is_isolate_mag_sag(self):
        """Test identification of a genome being a MAG, SAG, or isolate."""

        g = copy.copy(test_genome)

        g.ncbi_genome_category = "derived from metagenome"
        assert not g.is_isolate()
        assert g.is_mag()
        assert not g.is_sag()

        g.ncbi_genome_category = "derived from environmental sample"
        assert not g.is_isolate()
        assert g.is_mag()
        assert not g.is_sag()

        g.ncbi_genome_category = "derived from single cell"
        assert not g.is_isolate()
        assert g.is_sag()
        assert not g.is_mag()

        g.ncbi_genome_category = None
        assert g.is_isolate()
        assert not g.is_mag()
        assert not g.is_sag()

        g.ncbi_genome_category = "isolate"
        assert g.is_isolate()
        assert not g.is_mag()
        assert not g.is_sag()

        g.ncbi_genome_category = ""
        assert g.is_isolate()
        assert not g.is_mag()
        assert not g.is_sag()

    def test_gtdb_type_material(self):
        """Test identification of genome being type material according to GTDB."""

        g = copy.copy(test_genome)

        # untrusted type material should always return False
        g.gtdb_untrustworthy_as_type = True
        g.gtdb_type_species_of_genus = True
        assert not g.is_gtdb_type_species()
        g.gtdb_type_designation = 'type strain of species'
        assert not g.is_gtdb_type_strain()
        g.gtdb_type_designation = 'type strain of subspecies'
        assert not g.is_gtdb_type_subspecies()

        # test type material
        g.gtdb_untrustworthy_as_type = False
        assert g.is_gtdb_type_species()
        g.gtdb_type_designation = 'type strain of species'
        assert g.is_gtdb_type_strain()
        g.gtdb_type_designation = 'type strain of subspecies'
        assert g.is_gtdb_type_subspecies()

    def test_ncbi_type_material(self):
        """Test identification of genome being type material according to NCBI."""

        g = copy.copy(test_genome)

        # test type material flagged as untrustworthy within the GTDB
        g.excluded_from_refseq_note = ""
        g.gtdb_untrustworthy_as_type = True
        g.ncbi_type_material = 'assembly from type material'
        assert not g.is_ncbi_type_strain()
        g.ncbi_unfiltered_taxa = Taxa(
            "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri subsp. testus")
        assert not g.is_ncbi_type_subspecies()

        # test type material flagged as untrustworthy at NCBI
        g.excluded_from_refseq_note = "untrustworthy as type"
        g.gtdb_untrustworthy_as_type = False
        g.ncbi_type_material = 'assembly from type material'
        assert not g.is_ncbi_type_strain()
        g.ncbi_unfiltered_taxa = Taxa(
            "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri subsp. testus")
        assert not g.is_ncbi_type_subspecies()

        # test valid type material
        g.excluded_from_refseq_note = ""
        g.gtdb_untrustworthy_as_type = False
        g.ncbi_unfiltered_taxa = Taxa(
            "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri")
        g.ncbi_type_material = 'assembly from type material'
        assert g.is_ncbi_type_strain()
        assert not g.is_ncbi_type_subspecies()

        # test type material for NCBI subspecies
        g.ncbi_unfiltered_taxa = Taxa(
            "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri subsp. testus")
        g.ncbi_type_material = 'assembly from type material'
        assert not g.is_ncbi_type_strain()
        assert g.is_ncbi_type_subspecies()

    def test_is_complete_genome(self):
        """Test identifying a genome as a complete assembly."""

        g = copy.copy(test_genome)

        assert g.is_complete_genome()

        g.ssu_count = 0
        assert not g.is_complete_genome()

    def test_scoring(self):
        """Test assembly scoring."""

        g = copy.copy(test_genome)

        # highest possible score for isolate
        assembly_score = g.score_assembly()
        assert assembly_score == 210

        # highest possible score for isolate,
        # type strain genome
        type_score = g.score_type_strain()
        assert type_score == 1e6 + 210

        # highest possible score for MAG or SAG
        g.ncbi_genome_category = "derived from metagenome"
        assembly_score = g.score_assembly()
        assert assembly_score == 110

        # highest possible score for MAG/SAG
        # that is considered a type strain genome
        type_score = g.score_type_strain()
        assert type_score == 1e6 + 110

    def test_pass_qc(self):
        """Test reporting of genome passing QC."""

        def run_qc(marker_perc=100, **attrs):
            g = copy.copy(test_genome)
            for k, v in attrs.items():
                setattr(g, k, v)

            failed_tests = defaultdict(int)
            passed = g.pass_qc(marker_perc, qc_criteria, failed_tests)

            return passed, set(failed_tests)

        # high quality genome passes QC
        assert run_qc() == (True, set())

        # fails if both CheckM v1 and v2 estimates fail
        assert run_qc(cm1_comp=49, cm2_comp=49) == (False, {'cm1_comp', 'cm2_comp'})

        # passes if only one of CheckM v1 or v2 estimates fail,
        # provided the genome has few contigs
        assert run_qc(cm1_cont=11) == (True, set())
        assert run_qc(cm2_cont=11) == (True, set())
        assert run_qc(cm1_cont=11, contig_count=10) == (False, {'cm1_cont'})
        assert run_qc(cm2_cont=11, contig_count=10) == (False, {'cm2_cont'})

        # CheckM v1 contamination is discounted for genomes with high strain heterogeneity
        # (up to 20%), so long as the estimated quality is sufficient
        assert run_qc(cm1_cont=15, contig_count=10) == (False, {'cm1_cont'})
        assert run_qc(cm1_cont=15, cm1_strain_heterogeneity_100=90, contig_count=10) == (True, set())
        assert run_qc(cm1_cont=21, cm1_strain_heterogeneity_100=90, contig_count=10) == (False, {'cm1_cont'})
        assert run_qc(cm1_comp=60, cm1_cont=15, cm1_strain_heterogeneity_100=80,
                      contig_count=10) == (False, {'cm1_qual'})

        # assembly and marker gene criteria
        assert run_qc(marker_perc=39) == (False, {'marker_perc'})
        assert run_qc(contig_count=2001) == (False, {'contig_count'})
        assert run_qc(contig_n50=4999) == (False, {'N50'})
        assert run_qc(ambiguous_bases=100001) == (False, {'ambig'})

        # large genomes with low coding density
        assert run_qc(length=15e6, coding_density=69) == (False, {'abnormal_genome_size'})
        assert run_qc(length=14e6, coding_density=69) == (True, set())

        # genomes likely to have frameshift errors
        assert run_qc(cm2_comp=69, coding_density=69) == (False, {'frameshift'})
        assert run_qc(cm2_comp=69, coding_density=69, contig_count=20) == (True, set())

        # genomes marked by NCBI as having many frameshifted proteins
        assert run_qc(excluded_from_refseq_note='many frameshifted proteins',
                      cm2_comp=79, coding_density=79) == (False, {'ncbi_frameshift'})
        assert run_qc(excluded_from_refseq_note='many frameshifted proteins',
                      cm2_comp=80, coding_density=79) == (True, set())
