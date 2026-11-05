# Test identification of new taxa in UpdateCurationTrees class

from types import SimpleNamespace

from gtdb_species_clusters.taxa import Taxa
from gtdb_species_clusters.update_curation_trees import UpdateCurationTrees


def stub_genome(ncbi_taxonomy, seqcode_taxonomy=None):
    ncbi_taxa = Taxa(ncbi_taxonomy)
    seqcode_taxa = Taxa(seqcode_taxonomy) if seqcode_taxonomy else None
    return SimpleNamespace(ncbi_taxa=ncbi_taxa,
                           ref_taxa=seqcode_taxa if seqcode_taxa else ncbi_taxa,
                           ref_taxa_source='SeqCode' if seqcode_taxa else 'NCBI')


class StubGenomes(dict):
    def __init__(self, genomes):
        super().__init__(genomes)
        self.sp_clusters = {gid: set() for gid in genomes}


def test_new_ncbi_taxa(tmp_path):
    """Verify new taxa are identified from the SeqCode or NCBI classification of genomes."""

    acido = 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae'
    prev_genomes = StubGenomes({
        # genus known from NCBI in previous release
        'P1': stub_genome(f'{acido};g__Limnohabitans;s__Limnohabitans planktonicus'),
        # SeqCode genus known in previous release, although genome has no NCBI genus
        'P2': stub_genome(f'{acido};g__;s__', f'{acido};g__Acidiparvus;s__Acidiparvus lacustris')})

    cur_genomes = {
        'P1': prev_genomes['P1'],
        'P2': prev_genomes['P2'],
        # new SeqCode genus for genome without NCBI genus
        'C1': stub_genome(f'{acido};g__;s__', f'{acido};g__Aquirufa;s__Aquirufa nivalisilvae'),
        # new NCBI genus
        'C2': stub_genome(f'{acido};g__Rhodoferax;s__Rhodoferax ferrireducens'),
        # NCBI genus that is renamed under the SeqCode is not new
        'C3': stub_genome(f'{acido};g__Polynucleobacter;s__Polynucleobacter necessarius',
                          f'{acido};g__Acidiparvus;s__Acidiparvus necessarius')}

    cur_clusters = {'P1': {'P1', 'C2'}, 'P2': {'P2', 'C3'}, 'C1': {'C1'}}

    p = UpdateCurationTrees(str(tmp_path), 'gtdb_test')
    p.new_ncbi_taxa(prev_genomes, cur_genomes, cur_clusters)

    lines = (tmp_path / 'gtdb_test_new_ncbi_genus.tsv').read_text().splitlines()
    assert lines == ['NCBI domain\tNCBI taxa\tNo. representatives\tRepresentative ID(s)\tReference source used',
                     'd__Bacteria\tg__Aquirufa\t1\tC1\tSeqCode',
                     'd__Bacteria\tg__Rhodoferax\t1\tP1\tNCBI']
    assert (tmp_path / 'gtdb_test_new_ncbi_genus.tree').read_text() == '(C1,P1);\n'

    lines = (tmp_path / 'gtdb_test_new_ncbi_species.tsv').read_text().splitlines()
    assert lines[1:] == ['d__Bacteria\ts__Acidiparvus necessarius\t1\tP2\tSeqCode',
                         'd__Bacteria\ts__Aquirufa nivalisilvae\t1\tC1\tSeqCode',
                         'd__Bacteria\ts__Rhodoferax ferrireducens\t1\tP1\tNCBI']

    lines = (tmp_path / 'gtdb_test_new_ncbi_family.tsv').read_text().splitlines()
    assert lines[1:] == []
