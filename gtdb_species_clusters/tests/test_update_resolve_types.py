# Test parsing of LTP taxonomy strings in ResolveTypes class

import pytest

from gtdb_species_clusters.update_resolve_types import ResolveTypes


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
