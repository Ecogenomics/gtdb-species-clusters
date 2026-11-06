# Test identification of synonyms in NCBI_SpeciesManager class

import logging
from collections import defaultdict
from types import SimpleNamespace

import pytest

from gtdb_species_clusters.ncbi_species_manager import NCBI_SpeciesManager


class StubGenome:
    def __init__(self, ncbi_sp='s__', seqcode_sp=None, ledger_sp='', type_strain=False, score=0, untrustworthy_sp=False):
        self.ncbi_taxa = SimpleNamespace(species=ncbi_sp)
        self.ncbi_untrustworthy_sp = untrustworthy_sp
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
        # mirrors Genome.is_effective_type_strain: the type strains ledger takes precedence over an untrustworthy NCBI species
        if self.ncbi_untrustworthy_sp and not self.is_gtdb_type_strain_ledger():
            return False
        return self.type_strain or self.is_gtdb_type_strain_ledger()

    def is_ref_species_trusted(self):
        return self.ref_taxa_source != 'NCBI' or not self.ncbi_untrustworthy_sp

    def trusted_species_names(self):
        names = {self.ref_taxa.species} if self.is_ref_species_trusted() else set()
        if not self.ncbi_untrustworthy_sp:
            names.add(self.ncbi_taxa.species)
        return names

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


def test_classify_ncbi_species(ncbi_species_mngr, caplog):
    """Verify species are classified by the species type strain genomes are type material of and the reference species of other genomes.

    A type strain representative also anchors its reference and NCBI species when no type strain genome
    of these species is a representative (e.g., the NCBI species of a ledger genome whose other type strain
    genomes are in its cluster), and both names of the representative being assigned to its cluster is not an error.
    """

    with caplog.at_level(logging.ERROR, logger='rich'):
        _all_species, unambiguous, ambiguous = ncbi_species_mngr.classify_ncbi_species({}, {'M1'})

    assert unambiguous == {'s__Pelagibacter ubiqueversans': ('R1', NCBI_SpeciesManager.TYPE_STRAIN_GENOME),
                           's__Pelagibacter communis': ('R1', NCBI_SpeciesManager.TYPE_STRAIN_GENOME),
                           's__Synechococcus_D lacustris': ('L1', NCBI_SpeciesManager.TYPE_STRAIN_GENOME),
                           's__Synechococcus lacustris': ('L1', NCBI_SpeciesManager.TYPE_STRAIN_GENOME)}
    assert ambiguous == {'s__Pelagibacter bermudensis'}
    assert 'assigned to multiple unambiguous' not in caplog.text


def test_type_strain_fallback_species(tmp_path):
    """Verify a species is not anchored by type strain representatives of different clusters or when it has a type strain representative."""

    cur_genomes = StubGenomes({
        # two ledger genomes in different clusters with the same NCBI species
        'L1': StubGenome(ncbi_sp='s__Thioglobus perditus', ledger_sp='s__Thioglobus_B perditus'),
        'L2': StubGenome(ncbi_sp='s__Thioglobus perditus', ledger_sp='s__Thioglobus_C perditus'),
        'G1': StubGenome(ncbi_sp='s__Thioglobus perditus'),
        # ledger genome whose NCBI species has its own type strain representative
        'L3': StubGenome(ncbi_sp='s__Solincola tengchongensis', ledger_sp='s__Solincola_A tengchongensis'),
        'T3': StubGenome(ncbi_sp='s__Solincola tengchongensis', type_strain=True),
        # ledger genome whose NCBI species has no other genomes
        'L4': StubGenome(ncbi_sp='s__Aquicultor secundus', ledger_sp='s__Aquicultor_A secundus')})
    cur_clusters = {'L1': {'L1', 'G1'}, 'L2': {'L2'}, 'L3': {'L3'}, 'T3': {'T3'}, 'L4': {'L4'}}
    untrustworthy_type_file = tmp_path / 'untrustworthy_type_material.tsv'
    untrustworthy_type_file.write_text('Genome ID\n')
    mngr = NCBI_SpeciesManager(cur_genomes, cur_clusters, {}, str(untrustworthy_type_file), str(tmp_path))

    ncbi_species_gids = defaultdict(list)
    for gid in cur_genomes:
        ncbi_species_gids[mngr.classification_species(gid)].append(gid)

    assert mngr.type_strain_fallback_species(ncbi_species_gids) == {}


def untrustworthy_ncbi_sp_mngr(tmp_path):
    """Type strains ledger genomes whose NCBI species are in the NCBI untrustworthy species ledger, and so
    are misclassified genomes (see u_ncbi_erroneous), e.g. G900100765 is the type strain of s__Mooreia alkaloidigena
    and G002240565 is the type strain of s__Clavibacter insidiosus."""

    cur_genomes = StubGenomes({
        'L1': StubGenome(ncbi_sp='s__Catalinimonas alkaloidigena', ledger_sp='s__Mooreia alkaloidigena',
                         untrustworthy_sp=True, score=3),
        'T1': StubGenome(ncbi_sp='s__Catalinimonas alkaloidigena', type_strain=True, score=2),
        'G1': StubGenome(ncbi_sp='s__Catalinimonas alkaloidigena'),
        'R2': StubGenome(ncbi_sp='s__Clavibacter nebraskensis', type_strain=True, score=3),
        'L2': StubGenome(ncbi_sp='s__Clavibacter michiganensis', ledger_sp='s__Clavibacter insidiosus',
                         untrustworthy_sp=True, score=2)})
    cur_clusters = {'L1': {'L1'}, 'T1': {'T1', 'G1'}, 'R2': {'R2', 'L2'}}

    untrustworthy_type_file = tmp_path / 'untrustworthy_type_material.tsv'
    untrustworthy_type_file.write_text('Genome ID\n')

    return NCBI_SpeciesManager(cur_genomes, cur_clusters, {}, str(untrustworthy_type_file), str(tmp_path))


def test_classify_ledger_type_strain_with_untrustworthy_ncbi_species(tmp_path, caplog):
    """Verify a type strains ledger genome anchors its ledger species even though its NCBI species is untrustworthy,
    and does not anchor its NCBI species."""

    mngr = untrustworthy_ncbi_sp_mngr(tmp_path)
    with caplog.at_level(logging.ERROR, logger='rich'):
        _all_species, unambiguous, ambiguous = mngr.classify_ncbi_species({}, {'L1', 'L2'})

    assert unambiguous == {'s__Mooreia alkaloidigena': ('L1', NCBI_SpeciesManager.TYPE_STRAIN_GENOME),
                           's__Catalinimonas alkaloidigena': ('T1', NCBI_SpeciesManager.TYPE_STRAIN_GENOME),
                           's__Clavibacter nebraskensis': ('R2', NCBI_SpeciesManager.TYPE_STRAIN_GENOME)}
    assert ambiguous == {'s__Clavibacter insidiosus'}
    assert 'ERROR' not in caplog.text


def test_type_strain_synonyms_untrustworthy_ncbi_species(tmp_path):
    """Verify type strains ledger genomes with an untrustworthy NCBI species are type strain synonyms, and the
    untrustworthy NCBI species of a representative does not prevent type strain genomes of that species being synonyms."""

    mngr = untrustworthy_ncbi_sp_mngr(tmp_path)
    mngr.cur_clusters = {'L1': {'L1', 'T1', 'G1'}, 'R2': {'R2', 'L2'}}

    assert mngr.identify_type_strain_synonyms({'L1', 'L2'}) == {'L1': ['T1'], 'R2': ['L2']}


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


def test_write_synonym_table_seqcode_valid_names(ncbi_species_mngr, tmp_path, monkeypatch):
    """Verify names validly published under the SeqCode are valid names when checking synonym priority."""

    created = {}

    class StubSpeciesPriorityManager:
        def __init__(self, *args, seqcode_file=None):
            created['seqcode_file'] = seqcode_file

    monkeypatch.setattr('gtdb_species_clusters.ncbi_species_manager.SpeciesPriorityManager',
                        StubSpeciesPriorityManager)

    ncbi_species_mngr.write_synonym_table({}, {}, {},
                                          'sp_priority.tsv', 'genus_priority.tsv', 'lpsn_gss.csv',
                                          'seqcode_table.tsv')
    assert created['seqcode_file'] == 'seqcode_table.tsv'
    assert (tmp_path / 'synonyms.tsv').read_text().startswith('Synonym type\tNCBI species')
