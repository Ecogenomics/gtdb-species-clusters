# Test methods in Genome class

from collections import defaultdict
import copy

import pytest

from gtdb_species_clusters.genome import Genome
from gtdb_species_clusters.taxa import Taxa

test_genome = Genome("G012345678",
                                "GCA_012345678.1",
                                "G012345678",
                                True,
                                Taxa("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri"),
                                Taxa("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri"),
                                Taxa("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri"),
                                "type strain of species",
                                "LPSN",
                                True,
                                False,
                                False,
                                'assembly from type material',
                                "K-12",
                                "Complete genome",
                                "Full",
                                "Reference genome",
                                None,
                                "",
                                100,
                                0,
                                0,
                                3500000,
                                1,
                                3500000,
                                1,
                                0,
                                0,
                                1,
                                1600,
                                1,
                                0,
                                0,
                                "1900")

class TestGenome:

    def test_ncbi_subspecies(self):
        g = copy.copy(test_genome)
        
        # no subspecies defined
        g.ncbi_unfiltered_taxa = Taxa("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri")
        assert not g.is_ncbi_subspecies()
        
        # subspecies defined
        g.ncbi_unfiltered_taxa = Taxa("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri subsp. testus")
        assert g.is_ncbi_subspecies()
        
        # a variety or str is not a subspecies
        g.ncbi_unfiltered_taxa = Taxa("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri var. testus")
        assert not g.is_ncbi_subspecies()
        
        g.ncbi_unfiltered_taxa = Taxa("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri str testus")
        assert not g.is_ncbi_subspecies()
        
        # a variety is not a subspecies
        g.ncbi_unfiltered_taxa = Taxa("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri var. testus")
        assert not g.is_ncbi_subspecies()
        
        # genome is not from a subspecies if the subspecies and specific name are the same since it should be considered from the species itself;
        # this is critical for identifying a genome as the type strain of the species and not the type strain of a subspecies
        g.ncbi_unfiltered_taxa = Taxa("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri subsp. flexneri")
        assert not g.is_ncbi_subspecies()

    def test_is_isolate_mag_sag(self):
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
        g = copy.copy(test_genome)
        
        # untrusted type material should always return False
        g.gtdb_untrustworthy_as_type = True
        g.gtdb_type_species_of_genus = True
        assert not g.is_gtdb_type_species()
        g.gtdb_type_designation == 'type strain of species'
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
        g = copy.copy(test_genome)

        # test type material flagged as untrustworthy within the GTDB
        g.excluded_from_refseq_note = ""
        g.gtdb_untrustworthy_as_type = True
        g.ncbi_type_material = 'assembly from type material'
        assert not g.is_ncbi_type_strain()
        g.ncbi_unfiltered_taxa = Taxa("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri subsp. testus")
        assert not g.is_ncbi_type_subspecies()
        
        # test type material flagged as untrustworthy at NCBI
        g.excluded_from_refseq_note = "untrustworthy as type"
        g.gtdb_untrustworthy_as_type = False
        g.ncbi_type_material = 'assembly from type material'
        assert not g.is_ncbi_type_strain()
        g.ncbi_unfiltered_taxa = Taxa("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri subsp. testus")
        assert not g.is_ncbi_type_subspecies()
        
        # test valid type material
        g.excluded_from_refseq_note = ""
        g.gtdb_untrustworthy_as_type = False
        g.ncbi_unfiltered_taxa = Taxa("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri")
        g.ncbi_type_material = 'assembly from type material'
        assert g.is_ncbi_type_strain()
        assert not g.is_ncbi_type_subspecies()
        
        # test type material for NCBI subspecies
        g.ncbi_unfiltered_taxa = Taxa("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri;sb__Escherichia flexneri subsp. testus")
        g.ncbi_type_material = 'assembly from type material'
        assert not g.is_ncbi_type_strain()
        assert g.is_ncbi_type_subspecies()
        
    def test_is_complete_genome(self):
        g = copy.copy(test_genome)
        
        assert g.is_complete_genome()
        
        g.ssu_count = 0
        assert not g.is_complete_genome()
        
    def test_scoring(self):
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
        g = copy.copy(test_genome)
        
        failed_tests = defaultdict(int)
        assert g.pass_qc(100,       #marker_perc
                            50,     #min_comp
                            10,     #max_cont
                            50,     #min_quality
                            80,     #sh_exception
                            40,     #min_perc_markers
                            1000,   #max_contigs
                            5000,   #min_N50
                            100000, #max_ambiguous
                            failed_tests)
        
        g.comp = 49
        assert not g.pass_qc(100,       #marker_perc
                                50,     #min_comp
                                10,     #max_cont
                                50,     #min_quality
                                80,     #sh_exception
                                40,     #min_perc_markers
                                1000,   #max_contigs
                                5000,   #min_N50
                                100000, #max_ambiguous
                                failed_tests)
                                
        g.comp = 51
        g.cont = 11
        assert not g.pass_qc(100,       #marker_perc
                                50,     #min_comp
                                10,     #max_cont
                                50,     #min_quality
                                80,     #sh_exception
                                40,     #min_perc_markers
                                1000,   #max_contigs
                                5000,   #min_N50
                                100000, #max_ambiguous
                                failed_tests)
        

