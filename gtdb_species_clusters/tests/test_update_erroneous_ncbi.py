# Test identification of genomes with erroneous species assignments in UpdateErroneousNCBI class

from types import SimpleNamespace

import pytest

from gtdb_species_clusters.update_erroneous_ncbi import UpdateErroneousNCBI


@pytest.fixture
def erroneous_ncbi(tmp_path):
    return UpdateErroneousNCBI(str(tmp_path))


class StubGenome:
    def __init__(self, ncbi_sp='s__', seqcode_sp=None, ledger_sp='', type_strain=False, untrustworthy_sp=False):
        self.ncbi_taxa = SimpleNamespace(species=ncbi_sp)
        self.ncbi_untrustworthy_sp = untrustworthy_sp
        self.seqcode_taxa = SimpleNamespace(species=seqcode_sp, genus='g__' + seqcode_sp[3:].split()[0]) if seqcode_sp else None
        self.type_strains_ledger_sp = ledger_sp
        self.type_strain = type_strain

    @property
    def ref_taxa(self):
        return self.seqcode_taxa if self.seqcode_taxa is not None else self.ncbi_taxa

    def ref_taxa_columns(self):
        if self.seqcode_taxa is not None:
            return f'{self.seqcode_taxa.species}\t{self.seqcode_taxa.genus}\tSeqCode'
        return '\t\tNCBI'

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

    def trusted_species_names(self):
        # mirrors Genome.trusted_species_names
        names = set()
        if self.seqcode_taxa is not None or not self.ncbi_untrustworthy_sp:
            names.add(self.ref_taxa.species)
        if not self.ncbi_untrustworthy_sp:
            names.add(self.ncbi_taxa.species)
        return names


def test_identify_misclassified_genomes_cluster(erroneous_ncbi, tmp_path):
    """Verify genomes are checked against type strain genomes of their SeqCode or NCBI species."""

    cur_genomes = {
        # SeqCode type genome with no NCBI species anchors its SeqCode species in cluster R1
        'S1': StubGenome(seqcode_sp='s__Acidiparvus lacustris', type_strain=True),
        'G1': StubGenome(ncbi_sp='s__Acidiparvus lacustris'),
        'G2': StubGenome(ncbi_sp='s__Acidiparvus lacustris'),
        # SeqCode type genome also anchors its differing NCBI species, which has no other type strain genome
        'S2': StubGenome(ncbi_sp='s__Methanobrevibacter arboriphilus',
                         seqcode_sp='s__Methanobinarius arboriphilus', type_strain=True),
        'G3': StubGenome(ncbi_sp='s__Methanobrevibacter arboriphilus'),
        # SeqCode type genome does not anchor its differing NCBI species when that species has its own
        # type strain genome, and is not misclassified as it is grouped by its SeqCode species
        'N2': StubGenome(ncbi_sp='s__Halarcobacter ebronensis', type_strain=True),
        'S3': StubGenome(ncbi_sp='s__Halarcobacter ebronensis',
                         seqcode_sp='s__Halarcobacter ibericensis', type_strain=True),
        'G7': StubGenome(ncbi_sp='s__Halarcobacter ebronensis'),
        # ledger genome anchors its ledger species and its NCBI species, which has no other type strain genome
        'L1': StubGenome(ncbi_sp='s__Synechococcus lacustris', ledger_sp='s__Synechococcus_D lacustris'),
        'G4': StubGenome(ncbi_sp='s__Synechococcus lacustris'),
        # NCBI type strain genome
        'N1': StubGenome(ncbi_sp='s__Escherichia coli', type_strain=True),
        'G5': StubGenome(ncbi_sp='s__Escherichia coli'),
        'G6': StubGenome(ncbi_sp='s__Escherichia coli')}

    cur_clusters = {'S1': {'S1', 'G1'},
                    'G2': {'G2', 'G3'},
                    'S2': {'S2'},
                    'N2': {'N2', 'G7'},
                    'S3': {'S3'},
                    'L1': {'L1'},
                    'G4': {'G4'},
                    'N1': {'N1', 'G5'},
                    'G6': {'G6'}}

    untrustworthy_sp_ledger = tmp_path / 'ncbi_untrustworthy_sp.tsv'
    untrustworthy_sp_ledger.write_text('Genome ID\tNCBI species\n'
                                       'G5\ts__Escherichia coli\n')

    misclassified = erroneous_ncbi.identify_misclassified_genomes_cluster(
        cur_genomes, cur_clusters, str(untrustworthy_sp_ledger))
    assert misclassified == {'G2', 'G3', 'G4', 'G6'}

    lines = (tmp_path / 'ncbi_misclassified_sp.gtdb_clustering.tsv').read_text().splitlines()
    assert lines[0] == ('Genome ID\tNCBI species\tGenome cluster\tType species cluster\t'
                        'SeqCode species\tSeqCode genus\tReference source used')
    assert sorted(lines[1:]) == ['G2\ts__Acidiparvus lacustris\tG2\tS1\t\t\tNCBI',
                                 'G3\ts__Methanobrevibacter arboriphilus\tG2\tS2\t\t\tNCBI',
                                 'G4\ts__Synechococcus lacustris\tG4\tL1\t\t\tNCBI',
                                 'G5\ts__Escherichia coli\tN1\tManually marked as being erroneous\t\t\tNCBI',
                                 'G6\ts__Escherichia coli\tG6\tN1\t\t\tNCBI']


def test_seqcode_type_strains_in_different_clusters(erroneous_ncbi, tmp_path, caplog):
    """Verify SeqCode and NCBI type strain genomes of the same species in different clusters are reported."""

    cur_genomes = {
        'N1': StubGenome(ncbi_sp='s__Planktophila lacus', type_strain=True),
        'S1': StubGenome(ncbi_sp='s__Planktophila lacus', seqcode_sp='s__Planktophila lacus', type_strain=True)}
    cur_clusters = {'N1': {'N1'}, 'S1': {'S1'}}

    untrustworthy_sp_ledger = tmp_path / 'ncbi_untrustworthy_sp.tsv'
    untrustworthy_sp_ledger.write_text('Genome ID\tNCBI species\n')

    with caplog.at_level('ERROR'):
        misclassified = erroneous_ncbi.identify_misclassified_genomes_cluster(
            cur_genomes, cur_clusters, str(untrustworthy_sp_ledger))
    assert len(misclassified) == 1
    assert any(r.getMessage().startswith('Species s__Planktophila lacus has multiple type strain genomes in different clusters')
               for r in caplog.records)


def test_untrustworthy_ncbi_species_not_anchored(erroneous_ncbi, tmp_path):
    """Verify a type strains ledger genome does not anchor its untrustworthy NCBI species (e.g., G900100765 is the
    type strain of s__Mooreia alkaloidigena, which NCBI assigns to s__Catalinimonas alkaloidigena)."""

    cur_genomes = {
        'L1': StubGenome(ncbi_sp='s__Catalinimonas alkaloidigena', ledger_sp='s__Mooreia alkaloidigena',
                         untrustworthy_sp=True),
        'G1': StubGenome(ncbi_sp='s__Catalinimonas alkaloidigena')}
    cur_clusters = {'L1': {'L1'}, 'G1': {'G1'}}

    untrustworthy_sp_ledger = tmp_path / 'ncbi_untrustworthy_sp.tsv'
    untrustworthy_sp_ledger.write_text('Genome ID\tNCBI species\n'
                                       'L1\ts__Catalinimonas alkaloidigena\n')

    misclassified = erroneous_ncbi.identify_misclassified_genomes_cluster(
        cur_genomes, cur_clusters, str(untrustworthy_sp_ledger))
    assert misclassified == set()

    lines = (tmp_path / 'ncbi_misclassified_sp.gtdb_clustering.tsv').read_text().splitlines()
    assert lines[1:] == ['L1\ts__Catalinimonas alkaloidigena\tL1\tManually marked as being erroneous\t\t\tNCBI']
