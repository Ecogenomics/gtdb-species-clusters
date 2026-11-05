# Test identification of synonyms in NCBI_SpeciesManager class

from collections import defaultdict
from types import SimpleNamespace

import pytest

from gtdb_species_clusters.ncbi_species_manager import NCBI_SpeciesManager


class StubGenome:
    def __init__(self, ncbi_sp='s__', seqcode_sp=None, ledger_sp='', type_strain=False, score=0):
        self.ncbi_taxa = SimpleNamespace(species=ncbi_sp)
        self.seqcode_taxa = SimpleNamespace(species=seqcode_sp) if seqcode_sp else None
        self.gtdb_taxa = SimpleNamespace(species='s__')
        self.type_strains_ledger_sp = ledger_sp
        self.type_strain = type_strain
        self.score = score

    @property
    def ref_taxa(self):
        return self.seqcode_taxa if self.seqcode_taxa is not None else self.ncbi_taxa

    @property
    def ref_taxa_source(self):
        return 'SeqCode' if self.seqcode_taxa is not None else 'NCBI'

    def type_species(self):
        # mirrors Genome.type_species: ledger species, else SeqCode species, else NCBI species
        if self.type_strains_ledger_sp:
            return self.type_strains_ledger_sp
        return self.ref_taxa.species

    def is_gtdb_type_strain_ledger(self):
        return bool(self.type_strains_ledger_sp)

    def is_effective_type_strain(self):
        return self.type_strain or self.is_gtdb_type_strain_ledger()

    def score_type_strain(self):
        return self.score


class StubGenomes(dict):
    def type_strain_genomes_by_species(self):
        # mirrors Genomes.type_strain_genomes_by_species
        type_strain_genomes = defaultdict(set)
        for gid, genome in self.items():
            if genome.is_effective_type_strain() and genome.type_species() != 's__':
                type_strain_genomes[genome.type_species()].add(gid)
        return type_strain_genomes


@pytest.fixture
def cur_genomes():
    return StubGenomes({
        # SeqCode type genome with an NCBI species that differs from its SeqCode species
        'R1': StubGenome(ncbi_sp='s__Pelagibacter communis', seqcode_sp='s__Pelagibacter ubiqueversans', type_strain=True),
        # SeqCode type genome without an NCBI species is a type strain synonym under its SeqCode species
        'S1': StubGenome(seqcode_sp='s__Pelagibacter bermudensis', type_strain=True),
        # type strain genome without a species name is never a synonym
        'T0': StubGenome(type_strain=True),
        # genomes keeping the NCBI species of the SeqCode type genome are a majority vote synonym
        'G1': StubGenome(ncbi_sp='s__Pelagibacter communis', score=1),
        'G2': StubGenome(ncbi_sp='s__Pelagibacter communis', score=2),
        # ledger genome: its NCBI species is neither a type strain synonym nor a majority vote synonym
        'L1': StubGenome(ncbi_sp='s__Synechococcus lacustris', ledger_sp='s__Synechococcus_D lacustris'),
        'T1': StubGenome(ncbi_sp='s__Synechococcus lacustris', type_strain=True),
        'G3': StubGenome(ncbi_sp='s__Synechococcus lacustris'),
        # misclassified genomes are ignored
        'M1': StubGenome(ncbi_sp='s__Synechococcus rubescens', type_strain=True)})


@pytest.fixture
def ncbi_species_mngr(cur_genomes, tmp_path):
    cur_clusters = {'R1': {'R1', 'S1', 'T0', 'G1', 'G2'},
                    'L1': {'L1', 'T1', 'G3', 'M1'}}

    untrustworthy_type_file = tmp_path / 'untrustworthy_type_material.tsv'
    untrustworthy_type_file.write_text('Genome ID\n')

    return NCBI_SpeciesManager(cur_genomes, cur_clusters, {}, str(untrustworthy_type_file), str(tmp_path))


def test_identify_type_strain_synonyms(ncbi_species_mngr):
    """Verify type strain synonyms are identified by the species genomes are type material of."""

    synonyms = ncbi_species_mngr.identify_type_strain_synonyms({'M1'})
    assert synonyms == {'R1': ['S1']}


def test_identify_consensus_synonyms(ncbi_species_mngr):
    """Verify majority vote synonyms are identified by the SeqCode or NCBI species of genomes."""

    synonyms = ncbi_species_mngr.identify_consensus_synonyms({'M1'})
    assert synonyms == {'R1': ['G2']}


def test_synonym_species(ncbi_species_mngr):
    """Verify species names and sources written to the synonyms table."""

    assert ncbi_species_mngr.synonym_species('R1') == ('s__Pelagibacter ubiqueversans', 'SeqCode')
    assert ncbi_species_mngr.synonym_species('G1') == ('s__Pelagibacter communis', 'NCBI')
    assert ncbi_species_mngr.synonym_species('L1') == ('s__Synechococcus_D lacustris', 'GTDB curator')


def test_parse_synonyms_table(ncbi_species_mngr, tmp_path):
    """Verify synonyms are read from the reference species columns and rows without a species name are skipped."""

    synonym_file = tmp_path / 'synonyms.tsv'
    synonym_file.write_text('Synonym type\tNCBI species\tNCBI synonym\tWarnings\t'
                            'Reference species\tReference source used\tSynonym reference species\tSynonym reference source used\n'
                            'TYPE_STRAIN_SYNONYM\ts__Pelagibacter communis\ts__\t\t'
                            's__Pelagibacter ubiqueversans\tSeqCode\ts__Pelagibacter bermudensis\tSeqCode\n'
                            'MAJORITY_VOTE_SYNONYM\ts__\ts__Kariarchaeum pelagius\t\t'
                            's__\tNCBI\ts__Kariarchaeum pelagius\tNCBI\n')
    assert ncbi_species_mngr.parse_synonyms_table(str(synonym_file)) == {
        's__Pelagibacter bermudensis': 's__Pelagibacter ubiqueversans'}

    # tables written before the reference species columns were added
    synonym_file.write_text('Synonym type\tNCBI species\tNCBI synonym\tWarnings\n'
                            'TYPE_STRAIN_SYNONYM\ts__Escherichia coli\ts__Shigella flexneri\t\n')
    assert ncbi_species_mngr.parse_synonyms_table(str(synonym_file)) == {
        's__Shigella flexneri': 's__Escherichia coli'}
