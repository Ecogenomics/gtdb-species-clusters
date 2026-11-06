# Test naming of species clusters after manual curation in PMC_SpeciesNames class

import logging
from collections import defaultdict
from types import SimpleNamespace

import pytest

from gtdb_species_clusters.pmc_species_names import PMC_SpeciesNames
from gtdb_species_clusters.taxon_utils import sort_by_naming_priority

BURK = ['d__Bacteria', 'p__Pseudomonadota', 'c__Gammaproteobacteria', 'o__Burkholderiales', 'f__Burkholderiaceae']


class StubGenome:
    def __init__(self, gtdb_sp, ncbi_sp='s__', seqcode_sp=None, ledger_sp='',
                 type_species=False, type_strain=False):
        self.gtdb_taxa = SimpleNamespace(genus='g__' + gtdb_sp[3:].split()[0], species=gtdb_sp)
        self.ncbi_taxa = SimpleNamespace(species=ncbi_sp)
        self.seqcode_taxa = SimpleNamespace(species=seqcode_sp) if seqcode_sp else None
        self.type_strains_ledger_sp = ledger_sp
        self.gtdb_type_species = type_species
        self.type_strain = type_strain or type_species or bool(ledger_sp)

    @property
    def ref_taxa(self):
        return self.seqcode_taxa if self.seqcode_taxa is not None else self.ncbi_taxa

    def type_species(self):
        # mirrors Genome.type_species: ledger species, else SeqCode species, else NCBI species
        if self.type_strains_ledger_sp:
            return self.type_strains_ledger_sp
        return self.ref_taxa.species

    def is_gtdb_type_species(self):
        return self.gtdb_type_species

    def is_effective_type_strain(self):
        return self.type_strain

    def is_gtdb_type_strain(self):
        return self.type_strain


def taxonomy(genomes):
    return {gid: BURK + [g.gtdb_taxa.genus, g.gtdb_taxa.species] for gid, g in genomes.items()}


@pytest.fixture
def pmc(tmp_path):
    p = PMC_SpeciesNames(str(tmp_path))
    p.sp_priority_mngr = SimpleNamespace(genus_priority=lambda g1, g2: None)
    return p


def test_sort_by_naming_priority():
    """Verify SeqCode type genomes are named as type material rather than in a separate group."""

    cur_genomes = {
        # SeqCode type species in GTDB genus of its SeqCode name
        'S1': StubGenome('s__Acidiparvus lacustris', seqcode_sp='s__Acidiparvus lacustris', type_species=True),
        # SeqCode type species in another GTDB genus is named as a type strain
        'S2': StubGenome('s__DAHUXN01 sp000002', seqcode_sp='s__Allingihabitans vicinus', type_species=True),
        # SeqCode type strain
        'T1': StubGenome('s__Limnohabitans sp000003', seqcode_sp='s__Limnohabitans lacus', type_strain=True),
        # ledger genome in GTDB genus of its ledger name is a type species
        'L1': StubGenome('s__Ruminococcus_D bicirculans', ncbi_sp='s__Ruminococcus bicirculans',
                         ledger_sp='s__Ruminococcus_D bicirculans', type_species=True),
        # SeqCode genome that is not effective type material (e.g. untrustworthy as type)
        'B1': StubGenome('s__Rhodoferax sp000004', seqcode_sp='s__Rhodoferax antarcticus'),
        'P1': StubGenome('s__Rhodoferax sp000005')}
    prev_genomes = SimpleNamespace(sp_clusters=set())

    rtn = sort_by_naming_priority(cur_genomes, cur_genomes, prev_genomes, cur_genomes, {})

    assert rtn == ([], ['S1', 'L1'], ['S2', 'T1'], ['B1'], ['P1'])


def test_set_type_strain(pmc, caplog):
    """Verify type strain clusters are checked against the species they are type material of."""

    cur_genomes = {
        # GTDB epithet differing from SeqCode epithet only by gender is kept
        'T1': StubGenome('s__Halioglobus marinus', seqcode_sp='s__Seongchinamella marina', type_strain=True),
        # placeholder GTDB epithet replaced by SeqCode epithet
        'T2': StubGenome('s__DAHUXN01 sp963811805', seqcode_sp='s__Allingihabitans vicinus', type_strain=True),
        # erroneous suffix of type strain in genus of its SeqCode name is removed
        'T3': StubGenome('s__Limnohabitans lacus_A', seqcode_sp='s__Limnohabitans lacus', type_strain=True),
        # ledger name takes precedence over the SeqCode name
        'T4': StubGenome('s__Methanobrevibacter arboriphilicus', seqcode_sp='s__Methanobinarius arboriphilus',
                         ledger_sp='s__Methanobrevibacter arboriphilicus'),
        # type material without a species name keeps its GTDB name
        'T5': StubGenome('s__Rhodoferax sp000005', type_strain=True)}
    final_taxonomy = taxonomy(cur_genomes)
    unambiguous = {'s__Seongchinamella marina': ('T1', 'TYPE_STRAIN_GENOME'),
                   's__Allingihabitans vicinus': ('T2', 'TYPE_STRAIN_GENOME'),
                   's__Limnohabitans lacus': ('T3', 'TYPE_STRAIN_GENOME'),
                   's__Methanobrevibacter arboriphilicus': ('T4', 'TYPE_STRAIN_GENOME')}

    names = {}
    case_count = defaultdict(int)
    with caplog.at_level(logging.WARNING, logger='rich'):
        for rid in cur_genomes:
            names[rid] = pmc.set_type_strain(rid, final_taxonomy, cur_genomes, unambiguous, case_count)[0]

    assert names == {'T1': 's__Halioglobus marinus',
                     'T2': 's__DAHUXN01 vicinus',
                     'T3': 's__Limnohabitans lacus',
                     'T4': 's__Methanobrevibacter arboriphilicus',
                     'T5': 's__Rhodoferax sp000005'}
    assert dict(case_count) == {'INCONGRUENT_TYPE_STRAIN': 1, 'ERRONEOUS_SUFFIX_TYPE_STRAIN': 1}
    assert 'not designated as an unambiguous' not in caplog.text


def test_set_type_species(pmc, caplog):
    """Verify type species clusters are checked against the species they are type material of."""

    cur_genomes = {
        'S1': StubGenome('s__Acidiparvus lacustris', seqcode_sp='s__Acidiparvus lacustris', type_species=True),
        'S2': StubGenome('s__Acidiparvus lacustris_A', ncbi_sp='s__Polynucleobacter necessarius',
                         seqcode_sp='s__Acidiparvus necessarius', type_species=True)}
    final_taxonomy = taxonomy(cur_genomes)
    unambiguous = {'s__Acidiparvus lacustris': ('S1', 'TYPE_STRAIN_GENOME')}

    case_count = defaultdict(int)
    with caplog.at_level(logging.WARNING, logger='rich'):
        assert pmc.set_type_species('S1', final_taxonomy, cur_genomes, unambiguous, case_count) == (
            's__Acidiparvus lacustris', '')
        assert pmc.set_type_species('S2', final_taxonomy, cur_genomes, unambiguous, case_count)[0] == (
            's__Acidiparvus lacustris_A')

    assert dict(case_count) == {'INCONGRUENT_TYPE_SPECIES_SPECIFIC': 1}
    assert 'Species s__Acidiparvus necessarius represented by type species genome S2 not designated' in caplog.text


def test_resolve_duplicate_type_strain_species(pmc):
    """Verify duplicate type strain names are resolved using the genus of the species genomes are type material of."""

    cur_genomes = {
        # SeqCode type strain in the genus of its SeqCode name keeps the name
        'T1': StubGenome('s__Allorhodoferax lacus', seqcode_sp='s__Allorhodoferax lacus', type_strain=True),
        # NCBI type strain transferred from another genus
        'T2': StubGenome('s__Allorhodoferax lacus', ncbi_sp='s__Rhodoferax lacus', type_strain=True)}
    final_taxonomy = taxonomy(cur_genomes)
    pmc.sp_name_mngr = SimpleNamespace(suffixed_placeholder_sp_epithet=lambda generic, specific: f'{specific}_A')

    resolved = pmc.resolve_duplicate_type_strain_species(final_taxonomy, cur_genomes)

    assert resolved == {'T2'}
    assert final_taxonomy['T1'][-1] == 's__Allorhodoferax lacus'
    assert final_taxonomy['T2'][-1] == 's__Allorhodoferax lacus_A'


def test_resolve_placeholder_specific_epthets(pmc, tmp_path):
    """Verify a previous placeholder specific name is retained as a well-formed binomial, unless the
    name only corrects the Latin part of the previous name to account for a genus transfer."""

    cur_genomes = {
        # different placeholder epithet in the same genus reverts to the previous epithet
        'G1': StubGenome('s__Methanobinarius arboriphilicus_A'),
        # previous placeholder in a different genus is not retained
        'G2': StubGenome('s__Enorma sp000002'),
        # Latin epithet is not replaced by a previous placeholder
        'G3': StubGenome('s__Rhodoferax lacus'),
        # correction made to account for a genus transfer, keeping the suffix, is kept
        'G4': StubGenome('s__Aequitasia blattaphila_A'),
        # correction made to account for a genus transfer that changes the suffix reverts
        'G5': StubGenome('s__Hassallia radiosa_B')}
    prev_genomes = {'P1': StubGenome('s__Methanobinarius arboriphilus_A'),
                    'P2': StubGenome('s__Collinsella sp000002'),
                    'P3': StubGenome('s__Rhodoferax sp000003'),
                    'P4': StubGenome('s__Aequitasia blattaphilus_A'),
                    'P5': StubGenome('s__Hassallia radiosum_A')}
    new_to_prev_rid = {'G1': 'P1', 'G2': 'P2', 'G3': 'P3', 'G4': 'P4', 'G5': 'P5'}
    final_taxonomy = taxonomy(cur_genomes)
    cur_gtdb_sp = {'s__Aequitasia blattaphila_A': ('G4', 'SPECIFIC_EPITHET_GENUS_TRANSFER'),
                   's__Hassallia radiosa_B': ('G5', 'SPECIFIC_EPITHET_GENUS_TRANSFER')}
    pmc.final_name_log = open(tmp_path / 'final_sp_name.log', 'w')

    pmc.resolve_placeholder_specific_epthets(final_taxonomy, prev_genomes, new_to_prev_rid, set(), cur_gtdb_sp)
    pmc.final_name_log.close()

    assert final_taxonomy['G1'][-1] == 's__Methanobinarius arboriphilus_A'
    assert final_taxonomy['G2'][-1] == 's__Enorma sp000002'
    assert final_taxonomy['G3'][-1] == 's__Rhodoferax lacus'
    assert final_taxonomy['G4'][-1] == 's__Aequitasia blattaphila_A'
    assert final_taxonomy['G5'][-1] == 's__Hassallia radiosum_A'
    assert cur_gtdb_sp['s__Methanobinarius arboriphilus_A'] == ('G1', 'RETAIN_PREVIOUS_PLACEHOLDER_NAME')
    assert cur_gtdb_sp['s__Aequitasia blattaphila_A'] == ('G4', 'SPECIFIC_EPITHET_GENUS_TRANSFER')
    assert cur_gtdb_sp['s__Hassallia radiosum_A'] == ('G5', 'RETAIN_PREVIOUS_PLACEHOLDER_NAME')


def test_write_gtdb_synonym_table(pmc, tmp_path):
    """Verify the reference species columns of the synonyms table are written with a header."""

    header = ['Synonym type', 'NCBI species', 'GTDB representative', 'Strain IDs', 'Type sources', 'Priority year',
              'GTDB type species', 'GTDB type strain', 'NCBI assembly type', 'NCBI synonym',
              'Highest-quality synonym genome', 'Synonym strain IDs', 'Synonym type sources',
              'Synonym priority year', 'Synonym GTDB type species', 'Synonym GTDB type strain',
              'Synonym NCBI assembly type', 'ANI', 'AF', 'Warnings']
    row = ['MAJORITY_VOTE_SYNONYM', 's__Pelagibacter ubiqueversans', 'R1'] + [''] * 6 + ['s__Pelagibacter communis', 'G1'] \
        + [''] * 6 + ['99.0', '90.0', '']
    assert len(row) == len(header)
    synonym_file = tmp_path / 'synonyms.tsv'
    synonym_file.write_text('\n'.join([
        '\t'.join(header + ['Reference species', 'Reference source used',
                            'Synonym reference species', 'Synonym reference source used']),
        '\t'.join(row + ['s__Pelagibacter ubiqueversans', 'SeqCode', 's__Pelagibacter communis', 'NCBI'])]) + '\n')
    old_synonym_file = tmp_path / 'synonyms_old.tsv'
    old_synonym_file.write_text('\n'.join(['\t'.join(header), '\t'.join(row[:-1])]) + '\n')

    cur_genomes = SimpleNamespace(full_gid={'R1': 'RS_GCF_000000001.1', 'G1': 'GB_GCA_000000002.1'})
    final_taxonomy = {'R1': BURK + ['g__Pelagibacter', 's__Pelagibacter ubiqueversans']}

    for in_file, ref_columns in [(synonym_file, ['s__Pelagibacter ubiqueversans', 'SeqCode', 's__Pelagibacter communis', 'NCBI']),
                                 (old_synonym_file, ['', '', '', ''])]:
        out_file = tmp_path / 'gtdb_synonyms_final.tsv'
        pmc.write_gtdb_synonym_table(cur_genomes, final_taxonomy, str(in_file), str(out_file))

        out_header, out_row = [line.split('\t') for line in out_file.read_text().splitlines()]
        assert len(out_header) == len(out_row) == 25
        assert out_header[-4:] == ['Reference species', 'Reference source used',
                                   'Synonym reference species', 'Synonym reference source used']
        assert out_row[1] == 's__Pelagibacter ubiqueversans'
        assert out_row[3] == 'RS_GCF_000000001.1'
        assert out_row[-4:] == ref_columns
