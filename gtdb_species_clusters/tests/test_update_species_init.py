# Test naming of species clusters in UpdateSpeciesInit class

import logging
from types import SimpleNamespace

import pytest

from gtdb_species_clusters.genome import Genome
from gtdb_species_clusters.genome_utils import SeqCodeRecord
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters.taxa import Taxa
from gtdb_species_clusters.update_species_init import UpdateSpeciesInit

BURK = 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae'


class StubGenome:
    def __init__(self, gtdb, ncbi, seqcode=None, ledger_sp='', type_species=False, type_strain=False):
        self.gtdb_taxa = Taxa(f'{BURK};{gtdb}')
        self.ncbi_taxa = Taxa(f'{BURK};{ncbi}')
        self.seqcode_taxa = Taxa(f'{BURK};{seqcode}') if seqcode else None
        self.type_strains_ledger_sp = ledger_sp
        self.gtdb_type_species = type_species
        self.type_strain = type_strain or type_species

    @property
    def ref_taxa(self):
        return self.seqcode_taxa if self.seqcode_taxa is not None else self.ncbi_taxa

    @property
    def ref_taxa_source(self):
        return 'SeqCode' if self.seqcode_taxa is not None else 'NCBI'

    def ref_taxa_columns(self):
        seqcode_sp = self.seqcode_taxa.species if self.seqcode_taxa else ''
        seqcode_genus = self.seqcode_taxa.genus if self.seqcode_taxa else ''
        return f'{seqcode_sp}\t{seqcode_genus}\t{self.ref_taxa_source}'

    def type_species(self):
        # mirrors Genome.type_species: ledger species, else SeqCode species, else NCBI species
        if self.type_strains_ledger_sp:
            return self.type_strains_ledger_sp
        return self.ref_taxa.species

    def is_gtdb_type_species(self):
        return self.gtdb_type_species

    def is_effective_type_strain(self):
        return self.type_strain


class StubGenomes(dict):
    def __init__(self, genomes, sp_clusters=()):
        super().__init__(genomes)
        self.sp_clusters = set(sp_clusters)


class StubPriorityManager:
    def __init__(self, genus_years):
        self.genus_years = genus_years

    def genus_priority_year(self, genus):
        return self.genus_years.get(genus, Genome.NO_PRIORITY_YEAR)


@pytest.fixture
def sp_init(tmp_path):
    p = UpdateSpeciesInit(str(tmp_path))
    p.sp_priority_mngr = StubPriorityManager({'g__Limnohabitans': 2010, 'g__Acidiparvus': 2024})
    return p


def run_naming(sp_init, step, cur_genomes, prev_genomes=None, synonyms=None, ledger=None):
    prev_genomes = prev_genomes if prev_genomes is not None else StubGenomes({})
    clusters = {rid: {rid} for rid in cur_genomes}
    rids = sp_init._sort_by_naming_priority(clusters, prev_genomes, cur_genomes, synonyms or {})
    cluster_sp_names = {}
    used_sp_names = {}
    if step == 'binomial':
        step_fn = sp_init.name_binomial_clusters
        step_fn(rids, cluster_sp_names, used_sp_names, clusters, prev_genomes, cur_genomes, synonyms or {})
    else:
        step_fn = getattr(sp_init, f'name_{step}_clusters')
        step_fn(rids, cluster_sp_names, used_sp_names, clusters, prev_genomes, cur_genomes, ledger or {}, synonyms or {})

    return cluster_sp_names


def test_sort_by_naming_priority(sp_init):
    """Verify representatives with a SeqCode name but no NCBI name are not treated as placeholders."""

    cur_genomes = StubGenomes({
        'P1': StubGenome('g__JAAB01;s__JAAB01 sp000001', 'g__;s__'),
        'B1': StubGenome('g__JAAB01;s__JAAB01 sp000002', 'g__;s__', seqcode='g__Acidiparvus;s__Acidiparvus lacustris'),
        'T1': StubGenome('g__;s__', 'g__;s__', seqcode='g__Aquirufa;s__Aquirufa nivalisilvae', type_strain=True),
        'S1': StubGenome('g__;s__', 'g__;s__', seqcode='g__Acidiparvus;s__Acidiparvus lacustris', type_species=True)})
    clusters = {rid: {rid} for rid in cur_genomes}

    rids = sp_init._sort_by_naming_priority(clusters, StubGenomes({}), cur_genomes, {})
    assert rids == ['S1', 'T1', 'B1', 'P1']



def test_sort_by_naming_priority_continuing_clusters(sp_init):
    """Verify clusters containing a previous representative are named before new clusters, with ties by genome ID."""

    cur_genomes = StubGenomes({gid: StubGenome('g__VXYK01;s__VXYK01 sp047516945', 'g__;s__')
                               for gid in ('G9', 'G1', 'G7', 'G3', 'G5')})
    clusters = {'G9': {'G9', 'P2'},     # previous representative P2 replaced by G9
                'G1': {'G1', 'G0'},     # new cluster
                'G7': {'G7', 'P8'},     # previous representative P8 replaced by G7
                'G3': {'G3'},           # new cluster
                'G5': {'G5'}}           # previous representative
    prev_genomes = StubGenomes({}, sp_clusters={'G5', 'P2', 'P8'})

    rids = sp_init._sort_by_naming_priority(clusters, prev_genomes, cur_genomes, {})
    assert rids == ['G5', 'G7', 'G9', 'G1', 'G3']


def test_continuing_cluster_keeps_previous_name(sp_init):
    """Verify a cluster whose previous representative was replaced keeps its previous placeholder name
    when the representative of a new cluster is classified to the same species."""

    gtdbtk_sp = 'g__VXYK01;s__VXYK01 sp047516945'
    cur_genomes = StubGenomes({
        'G060099985': StubGenome(gtdbtk_sp, 'g__;s__'),     # representative of a new cluster
        'G060103125': StubGenome(gtdbtk_sp, 'g__;s__'),
        'G977925895': StubGenome(gtdbtk_sp, 'g__;s__'),     # replaced previous representative G047516945
        'G047516945': StubGenome(gtdbtk_sp, 'g__;s__')})
    clusters = {'G060099985': {'G060099985', 'G060103125'},
                'G977925895': {'G977925895', 'G047516945'}}
    prev_genomes = StubGenomes({}, sp_clusters={'G047516945'})

    rids = sp_init._sort_by_naming_priority(clusters, prev_genomes, cur_genomes, {})
    names = {}
    sp_init.name_placeholder_clusters(rids, names, {}, clusters, prev_genomes, cur_genomes, {})
    sp_init.sp_curation_log.close()

    assert names == {'G977925895': 's__VXYK01 sp047516945',
                     'G060099985': 's__VXYK01 sp060099985'}
    assert cur_genomes['G060103125'].gtdb_taxa.species == 's__VXYK01 sp060099985'


def test_name_type_species_clusters(sp_init, tmp_path, caplog):
    """Verify type species clusters are named from their SeqCode name and conflicts list the source."""

    cur_genomes = StubGenomes({
        # SeqCode type species without an NCBI species and with a placeholder GTDB genus
        'S1': StubGenome('g__JAAB01;s__JAAB01 sp000001', 'g__;s__',
                         seqcode='g__Acidiparvus;s__Acidiparvus lacustris', type_species=True),
        # SeqCode type species in a GTDB genus with another Latin name
        'S2': StubGenome('g__Limnohabitans;s__Limnohabitans sp000002', 'g__;s__',
                         seqcode='g__Aquirufa;s__Aquirufa nivalisilvae', type_species=True),
        # type species of genus without any species name is skipped with a warning
        'S3': StubGenome('g__JAAB03;s__JAAB03 sp000003', 'g__;s__', type_species=True),
        # NCBI type species is unchanged
        'S4': StubGenome('g__Polynucleobacter;s__Polynucleobacter necessarius',
                         'g__Polynucleobacter;s__Polynucleobacter necessarius', type_species=True)})

    with caplog.at_level(logging.WARNING, logger='rich'):
        names = run_naming(sp_init, 'type_species', cur_genomes)
    sp_init.sp_curation_log.close()

    assert names == {'S1': 's__Acidiparvus lacustris',
                     'S2': 's__Limnohabitans nivalisilvae',
                     'S4': 's__Polynucleobacter necessarius'}
    assert cur_genomes['S1'].gtdb_taxa.species == 's__Acidiparvus lacustris'
    assert 'no SeqCode or NCBI species assignment: S3' in caplog.text

    conflicts = (tmp_path / 'conflicting_type_species_of_genus.tsv').read_text().splitlines()
    assert conflicts == ['NCBI genus\tProposed GTDB genera\tReference source used',
                         'g__Aquirufa\tg__Limnohabitans\tSeqCode']

    curation = (tmp_path / 'sp_curation_log.tsv').read_text().splitlines()
    assert curation[0].endswith('\tIssue\tSeqCode species\tSeqCode genus\tReference source used')
    assert len(curation) == 2
    tokens = curation[1].split('\t')
    assert tokens[1] == 'S2'
    assert tokens[6] == 's__Limnohabitans nivalisilvae'
    assert tokens[-3:] == ['s__Aquirufa nivalisilvae', 'g__Aquirufa', 'SeqCode']


def test_name_type_strain_clusters(sp_init):
    """Verify type strain clusters are named from the species they are type material of."""

    cur_genomes = StubGenomes({
        # SeqCode type strain without a GTDB genus takes the SeqCode generic name
        'T1': StubGenome('g__;s__', 'g__;s__', seqcode='g__Aquirufa;s__Aquirufa nivalisilvae', type_strain=True),
        # SeqCode name in another genus with an epithet differing only by gender keeps the GTDB name
        'T2': StubGenome('g__Rhodoferax;s__Rhodoferax planktonica', 'g__;s__',
                         seqcode='g__Limnohabitans;s__Limnohabitans planktonicus', type_strain=True),
        # SeqCode name takes precedence over the NCBI name
        'T3': StubGenome('g__Rhodoferax;s__Rhodoferax sp000003', 'g__Rhodoferax;s__Rhodoferax ferrireducens',
                         seqcode='g__Rhodoferax;s__Rhodoferax antarcticus', type_strain=True),
        # ledger name takes precedence over the SeqCode name
        'T4': StubGenome('g__Rhodoferax;s__Rhodoferax sp000004', 'g__;s__',
                         seqcode='g__Rhodoferax;s__Rhodoferax lacus', ledger_sp='s__Rhodoferax_A lacus',
                         type_strain=True),
        # type material without a species name is skipped
        'T5': StubGenome('g__Rhodoferax;s__Rhodoferax sp000005', 'g__;s__', type_strain=True)})

    names = run_naming(sp_init, 'type_strain', cur_genomes, ledger={'T4': 's__Rhodoferax_A lacus'})

    assert names == {'T1': 's__Aquirufa nivalisilvae',
                     'T2': 's__Rhodoferax planktonica',
                     'T3': 's__Rhodoferax antarcticus',
                     'T4': 's__Rhodoferax_A lacus'}


def test_name_binomial_clusters(sp_init):
    """Verify clusters without type material are named from their SeqCode or NCBI binomial name."""

    cur_genomes = StubGenomes({
        # SeqCode type genome that is not effective type material (e.g. untrustworthy as type)
        'B1': StubGenome('g__Limnohabitans;s__Limnohabitans sp000001', 'g__;s__',
                         seqcode='g__Limnohabitans;s__Limnohabitans lacus'),
        'B2': StubGenome('g__Rhodoferax;s__Rhodoferax sp000002', 'g__Rhodoferax;s__Rhodoferax ferrireducens'),
        'P1': StubGenome('g__Rhodoferax;s__Rhodoferax sp000003', 'g__;s__')})

    names = run_naming(sp_init, 'binomial', cur_genomes)

    assert names == {'B1': 's__Limnohabitans lacus',
                     'B2': 's__Rhodoferax ferrireducens'}


def test_update_log_reference_columns(sp_init, tmp_path):
    """Verify the update log records the SeqCode name and source of the proposed name."""

    cur_genomes = StubGenomes({
        'S1': StubGenome('g__JAAB01;s__JAAB01 sp000001', 'g__;s__',
                         seqcode='g__Acidiparvus;s__Acidiparvus lacustris', type_species=True)})
    prev_genomes = StubGenomes({'S1': cur_genomes['S1']}, sp_clusters={'S1'})

    run_naming(sp_init, 'type_species', cur_genomes, prev_genomes=prev_genomes)
    sp_init.sp_name_log.close()

    rows = (tmp_path / 'sp_name_log.tsv').read_text().splitlines()
    assert rows[0].endswith('\tAction\tSeqCode species\tSeqCode genus\tReference source used')
    assert rows[1].split('\t')[1] == 'S1'
    assert rows[1].split('\t')[-3:] == ['s__Acidiparvus lacustris', 'g__Acidiparvus', 'SeqCode']


def test_seqcode_genus_priority():
    """Verify the SeqCode year of the type species is used for genera without an LPSN entry."""

    records = [SeqCodeRecord('Acidiparvus lacustris', '', 2024, True),
               SeqCodeRecord('Acidiparvus necessarius', '', 2022, False),
               SeqCodeRecord('Aquirufa nivalisilvae', '', 2023, True),
               SeqCodeRecord('Aquirufa antheringensis', '', 2021, True),
               SeqCodeRecord('Limnohabitans lacus', '', 2025, True),
               SeqCodeRecord('Rhodoferax lacus', '', Genome.NO_PRIORITY_YEAR, True)]

    mngr = SpeciesPriorityManager.__new__(SpeciesPriorityManager)
    mngr._fout_ambiguous_priority = SimpleNamespace(close=lambda: None)
    mngr.lpsn_genus_priority = {'g__Limnohabitans': 2010}
    mngr.seqcode_genus_priority = SpeciesPriorityManager._seqcode_genus_priority(records)

    assert mngr.seqcode_genus_priority == {'g__Acidiparvus': 2024,
                                           'g__Aquirufa': 2021,
                                           'g__Limnohabitans': 2025}
    assert mngr.genus_priority_year('g__Limnohabitans') == 2010
    assert mngr.genus_priority_year('g__Acidiparvus') == 2024
    assert mngr.genus_priority_year('g__Rhodoferax') == Genome.NO_PRIORITY_YEAR


def test_set_gtdbtk_classification(tmp_path):
    """Verify new genomes with a previously unseen SeqCode or NCBI species take that genus and species."""

    prev_genomes = StubGenomes({
        'P1': StubGenome('g__Acidiparvus;s__Acidiparvus lacustris', 'g__;s__',
                         seqcode='g__Acidiparvus;s__Acidiparvus lacustris')})

    cur_genomes = Genomes()
    cur_genomes.genomes = {
        # SeqCode species already seen in previous release: GTDB-Tk classification kept
        'C1': StubGenome('g__;s__', 'g__;s__', seqcode='g__Acidiparvus;s__Acidiparvus lacustris'),
        # unseen SeqCode species
        'C2': StubGenome('g__;s__', 'g__Limnohabitans;s__Limnohabitans planktonicus',
                         seqcode='g__Aquirufa;s__Aquirufa nivalisilvae'),
        # unseen NCBI species
        'C3': StubGenome('g__;s__', 'g__Rhodoferax;s__Rhodoferax ferrireducens'),
        # no species name
        'C4': StubGenome('g__;s__', 'g__;s__')}

    gtdbtk_file = tmp_path / 'gtdbtk_classify.tsv'
    lines = ['user_genome\tclassification']
    for gid in cur_genomes.genomes:
        lines.append(f'{gid}\t{BURK};g__Polynucleobacter;s__Polynucleobacter sp{gid}')
    gtdbtk_file.write_text('\n'.join(lines) + '\n')

    num_updated, num_ref_sp = cur_genomes.set_gtdbtk_classification(str(gtdbtk_file), prev_genomes)

    assert (num_updated, num_ref_sp) == (4, 2)
    assert cur_genomes['C1'].gtdb_taxa.species == 's__Polynucleobacter spC1'
    assert cur_genomes['C2'].gtdb_taxa.genus == 'g__Aquirufa'
    assert cur_genomes['C2'].gtdb_taxa.species == 's__Aquirufa nivalisilvae'
    assert cur_genomes['C3'].gtdb_taxa.species == 's__Rhodoferax ferrireducens'
    assert cur_genomes['C4'].gtdb_taxa.species == 's__Polynucleobacter spC4'
