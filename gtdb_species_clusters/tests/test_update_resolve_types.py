# Test parsing of LTP taxonomy strings in ResolveTypes class

import pytest

from gtdb_species_clusters.update_resolve_types import ResolveTypes, LTP_METADATA


@pytest.fixture
def resolve_types(tmp_path):
    return ResolveTypes(1, 'ltp', str(tmp_path))


def test_parse_ltp_taxonomy_str(resolve_types):
    """Verify taxa and species are parsed from standard LTP taxonomy strings."""

    taxa, sp = resolve_types._parse_ltp_taxonomy_str(
        'Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia fergusonii;|')
    assert taxa == ['Bacteria', 'Pseudomonadota', 'Gammaproteobacteria',
                    'Enterobacterales', 'Enterobacteriaceae', 'Escherichia', 'Escherichia fergusonii']
    assert sp == 's__Escherichia fergusonii'

    taxa, sp = resolve_types._parse_ltp_taxonomy_str(
        'Bacteria;Pseudomonadota;Alphaproteobacteria;Hyphomicrobiales;Rhizobiaceae;Rhizobium;Rhizobium leguminosarum;type sp.|USDA 2370')
    assert taxa[-1] == 'Rhizobium leguminosarum'
    assert sp == 's__Rhizobium leguminosarum'


def test_parse_ltp_taxonomy_str_cleaning(resolve_types):
    """Verify whitespace and leading or trailing '|' or '.' characters are removed from taxa."""

    taxa, sp = resolve_types._parse_ltp_taxonomy_str(
        'Bacteria; Pseudomonadota; Betaproteobacteria; Burkholderiales;|Oxalobacteraceae; Oxalobacter.;Oxalobacter paeniformigenes;|')
    assert taxa == ['Bacteria', 'Pseudomonadota', 'Betaproteobacteria',
                    'Burkholderiales', 'Oxalobacteraceae', 'Oxalobacter', 'Oxalobacter paeniformigenes']
    assert sp == 's__Oxalobacter paeniformigenes'

    taxa, sp = resolve_types._parse_ltp_taxonomy_str(
        'Bacteria;Pseudomonadota;Alphaproteobacteria;Rhodospirillales;Rhodospirillaceae;Aerophototrophica.;Aerophototrophica crusticola;type sp.|')
    assert taxa[-2:] == ['Aerophototrophica', 'Aerophototrophica crusticola']
    assert sp == 's__Aerophototrophica crusticola'

    _taxa, sp = resolve_types._parse_ltp_taxonomy_str(
        'Bacteria;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas nicosulfuronedens ;|')
    assert sp == 's__Pseudomonas nicosulfuronedens'


def test_parse_ltp_taxonomy_str_invalid(resolve_types):
    """Verify LTP entries without a taxonomy are reported as invalid."""

    assert resolve_types._parse_ltp_taxonomy_str(';;|') == ([], None)


def test_ltp_missing_taxonomy(resolve_types, tmp_path, caplog):
    """Verify LTP sequences without a taxonomy are written to file instead of the log."""

    ltp_taxonomy_file = tmp_path / 'ltp_taxonomy.tsv'
    ltp_taxonomy_file.write_text(
        'AF530475\tBacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia fergusonii;|\n'
        'MH645804\t;;|\n'
        'KC247325\t;;|\n')

    with caplog.at_level('WARNING'):
        ltp_species = resolve_types.ltp_defined_species(str(ltp_taxonomy_file))
    assert ltp_species == {'s__Escherichia fergusonii'}
    assert not caplog.records

    resolve_types.write_ltp_missing_taxonomy()
    lines = (tmp_path / 'ltp_missing_taxonomy.tsv').read_text().splitlines()
    assert lines == ['LTP accession\tLTP taxonomy string',
                     'KC247325\t;;|',
                     'MH645804\t;;|']


def ltp_hit(species):
    """Create a high-confidence LTP hit to the specified species."""

    return LTP_METADATA(accession='X', taxonomy='', taxa=[], species=species, ssu_len=1500,
                        evalue=0.0, bitscore=2000, aln_len=1500, perc_iden=100.0, perc_aln=100.0)


def test_ltp_conflicting_genomes(resolve_types):
    """Verify only retained type strain genomes are reported as having conflicting LTP assignments."""

    ncbi_sp = 's__Curtobacterium herbarum'
    type_gids = {'G1', 'G2', 'G3', 'G4', 'G5'}
    ltp_metadata = {'G1': [ltp_hit(ncbi_sp)],                                # untrustworthy, LTP agrees
                    'G2': [ltp_hit(ncbi_sp)],                                # retained, LTP agrees
                    'G3': [ltp_hit('s__Curtobacterium flaccumfaciens')],     # retained, LTP conflicts
                    'G4': [ltp_hit('s__Curtobacterium oceanosedimentum')],   # untrustworthy, LTP conflicts
                    'G5': [ltp_hit(ncbi_sp), ltp_hit('s__Other species')]}   # retained, LTP includes NCBI species

    conflicts = resolve_types.ltp_conflicting_genomes(ncbi_sp, type_gids, {'G1': '', 'G4': ''}, ltp_metadata)
    assert conflicts == {'G3': {'s__Curtobacterium flaccumfaciens'}}

    # no conflicts once the conflicting genome is also considered untrustworthy
    conflicts = resolve_types.ltp_conflicting_genomes(ncbi_sp, type_gids, {'G1': '', 'G3': '', 'G4': ''}, ltp_metadata)
    assert conflicts == {}


class StubGenome:
    def __init__(self, gtdb_type_designation_sources='LPSN'):
        self.gtdb_type_designation_sources = gtdb_type_designation_sources


def test_resolve_retained_ltp_conflicts(resolve_types):
    """Verify retained genomes with conflicting LTP assignments are considered untrustworthy."""

    ncbi_sp = 's__Curtobacterium herbarum'
    type_gids = {'G1', 'G2', 'G3', 'G4', 'G5', 'G6'}
    ltp_metadata = {'G1': [ltp_hit(ncbi_sp)],                               # retained, LTP agrees
                    'G2': [ltp_hit('s__Curtobacterium flaccumfaciens')],    # retained, LTP conflicts
                    'G3': [],                                               # retained, no LTP hit
                    'G4': [ltp_hit('s__Curtobacterium oceanosedimentum')],  # already untrustworthy
                    'G5': [ltp_hit('s__Curtobacterium pusillum')],          # retained, but in GTDB type strains ledger
                    'G6': [ltp_hit(ncbi_sp), ltp_hit('s__Other species')]}  # retained, LTP includes NCBI species
    cur_genomes = {gid: StubGenome() for gid in type_gids}
    cur_genomes['G5'] = StubGenome('GTDB curator')

    untrustworthy = resolve_types.resolve_retained_ltp_conflicts(
        ncbi_sp, type_gids, {'G4': 'reason'}, ltp_metadata, cur_genomes)
    assert list(untrustworthy) == ['G2']
    assert untrustworthy['G2'] == ('Conflicting 16S rRNA hits to LTP database of s__Curtobacterium flaccumfaciens, '
                                   'while other type strain genomes have hits to s__Curtobacterium herbarum')

    # no genomes considered untrustworthy unless a retained genome has the expected LTP assignment
    untrustworthy = resolve_types.resolve_retained_ltp_conflicts(
        ncbi_sp, type_gids, {'G1': 'reason', 'G6': 'reason'}, ltp_metadata, cur_genomes)
    assert untrustworthy == {}
