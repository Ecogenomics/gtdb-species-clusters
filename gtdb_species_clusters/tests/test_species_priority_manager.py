# Test naming priority of genomes in SpeciesPriorityManager class

from collections import defaultdict
from types import SimpleNamespace

import pytest

from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager


class StubGenome:
    def __init__(self, ref_sp, ledger_sp='', type_strain=True, year=2000, gtdb_genus='g__'):
        self.ref_taxa = SimpleNamespace(species=ref_sp, genus='g__' + ref_sp[3:].split()[0])
        self.gtdb_taxa = SimpleNamespace(genus=gtdb_genus)
        self.type_strains_ledger_sp = ledger_sp
        self.type_strain = type_strain
        self.year = year

    def type_species(self):
        # mirrors Genome.type_species: ledger species, else reference species
        return self.type_strains_ledger_sp or self.ref_taxa.species

    def is_effective_type_strain(self):
        return self.type_strain

    def is_exclusively_effective_type_strain(self):
        return False

    def is_gtdb_type_strain(self):
        return self.type_strain

    def is_gtdb_type_species(self):
        return False

    def year_of_priority(self):
        return self.year


@pytest.fixture
def sp_priority_mngr(tmp_path):
    mngr = SpeciesPriorityManager.__new__(SpeciesPriorityManager)
    mngr.valid_names = {'s__Catalinimonas alkaloidigena', 's__Mooreia alkaloidigena', 's__Ruminococcus bicirculans'}
    mngr.manual_sp_priority = defaultdict(dict)
    mngr.manual_species_priority_year = {}
    mngr._reported_ambiguous_pairs = set()
    mngr._fout_ambiguous_priority = open(tmp_path / 'ambiguous_sp_priority.tsv', 'w')
    yield mngr
    mngr._fout_ambiguous_priority.close()


def test_priority_names():
    """Verify genomes are identified by the species they are type material of, with suffixes ignored for valid names."""

    genome = StubGenome('s__Ruminococcus bicirculans', ledger_sp='s__Ruminococcus_D bicirculans')
    assert SpeciesPriorityManager.priority_names(genome) == ('s__Ruminococcus_D bicirculans',
                                                             'g__Ruminococcus_D',
                                                             's__Ruminococcus bicirculans')

    genome = StubGenome('s__Catalinimonas alkaloidigena', ledger_sp='s__Mooreia alkaloidigena')
    assert SpeciesPriorityManager.priority_names(genome) == ('s__Mooreia alkaloidigena',
                                                             'g__Mooreia',
                                                             's__Mooreia alkaloidigena')

    genome = StubGenome('s__Catalinimonas alkaloidigena')
    assert SpeciesPriorityManager.priority_names(genome) == ('s__Catalinimonas alkaloidigena',
                                                             'g__Catalinimonas',
                                                             's__Catalinimonas alkaloidigena')


def test_species_priority_ledger_species(sp_priority_mngr):
    """Verify a type strains ledger genome is not of the same species as a type strain of its NCBI species,
    so priority is resolved by the species priority ledger using its ledger species."""

    cur_genomes = {'L1': StubGenome('s__Catalinimonas alkaloidigena', ledger_sp='s__Mooreia alkaloidigena'),
                   'T1': StubGenome('s__Catalinimonas alkaloidigena')}
    sp_priority_mngr.manual_sp_priority['s__Mooreia alkaloidigena']['s__Catalinimonas alkaloidigena'] = 's__Mooreia alkaloidigena'
    sp_priority_mngr.manual_sp_priority['s__Catalinimonas alkaloidigena']['s__Mooreia alkaloidigena'] = 's__Mooreia alkaloidigena'

    assert sp_priority_mngr.species_priority(cur_genomes, 'T1', 'L1') == ('L1', 'manual curation')
    assert sp_priority_mngr.species_priority(cur_genomes, 'L1', 'T1') == ('L1', 'manual curation')


def test_species_priority_suffixed_ledger_species(sp_priority_mngr):
    """Verify a suffixed ledger species is a valid name and its suffixed genus matches the GTDB genus."""

    cur_genomes = {'L1': StubGenome('s__Ruminococcus bicirculans', ledger_sp='s__Ruminococcus_D bicirculans'),
                   'T1': StubGenome('s__Alpha beta')}
    assert sp_priority_mngr.species_priority(cur_genomes, 'T1', 'L1') == ('L1', 'selected valid species name')

    sp_priority_mngr.valid_names.add('s__Alpha beta')
    cur_genomes = {'L1': StubGenome('s__Ruminococcus bicirculans', ledger_sp='s__Ruminococcus_D bicirculans',
                                    gtdb_genus='g__Ruminococcus_D'),
                   'T1': StubGenome('s__Alpha beta', gtdb_genus='g__Ruminococcus_D')}
    assert sp_priority_mngr.species_priority(cur_genomes, 'T1', 'L1') == ('L1', 'match to expected GTDB genus assignment')


def test_species_priority_year_ledger_species(sp_priority_mngr):
    """Verify the manual year of priority of a type strains ledger genome is that of its ledger species."""

    cur_genomes = {'L1': StubGenome('s__Catalinimonas alkaloidigena', ledger_sp='s__Mooreia alkaloidigena', year=2013)}
    assert sp_priority_mngr.species_priority_year(cur_genomes, 'L1') == 2013

    sp_priority_mngr.manual_species_priority_year['s__Catalinimonas alkaloidigena'] = 2010
    assert sp_priority_mngr.species_priority_year(cur_genomes, 'L1') == 2013

    sp_priority_mngr.manual_species_priority_year['s__Mooreia alkaloidigena'] = 2020
    assert sp_priority_mngr.species_priority_year(cur_genomes, 'L1') == 2020
