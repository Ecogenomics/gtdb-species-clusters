# Test functions in genome_utils module

import pytest

from gtdb_species_clusters.genome_utils import parse_ncbi_env_bioproject_ledger


def write_ledger(path, rows, header='NCBI BioProject\tAssembly type\tSubmitter\tDate Added\t'):
    """Write ledger with Windows line endings and trailing tabs as in the curated ledger."""

    path.write_bytes(('\r\n'.join([header] + rows) + '\r\n').encode())
    return str(path)


def test_parse_ncbi_env_bioproject_ledger(tmp_path):
    """Verify BioProjects are split into those with environmental genomes and those with isolates."""

    ledger = write_ledger(tmp_path / 'ncbi_env_bioprojects.tsv',
                          ['PRJEB31003\tEnvironmental\tDHP\t"Nov. 5, 2020"\t',
                           'PRJNA903559\tEnvironmental\tDHP\t"Oct. 18, 2023"\tMajority of assemblies are MAGs',
                           'PRJEB96880\tIsolates\tDHP\t"Oct. 16, 2026"\t',
                           'PRJNA000001\t isolates \tDHP\t"Oct. 16, 2026"\t',
                           ''])

    env_bioprojects, isolate_bioprojects = parse_ncbi_env_bioproject_ledger(ledger)
    assert env_bioprojects == {'PRJEB31003', 'PRJNA903559'}
    assert isolate_bioprojects == {'PRJEB96880', 'PRJNA000001'}


def test_parse_ncbi_env_bioproject_ledger_invalid_type(tmp_path, caplog):
    """Verify processing stops for an unknown assembly type."""

    ledger = write_ledger(tmp_path / 'ncbi_env_bioprojects.tsv',
                          ['PRJEB31003\tEnvironmental\tDHP\t2020\t',
                           'PRJEB96880\tIsolate\tDHP\t2026\t'])

    with pytest.raises(SystemExit):
        parse_ncbi_env_bioproject_ledger(ledger)
    assert 'Invalid assembly type `Isolate` for BioProject PRJEB96880' in caplog.records[-1].getMessage()


def test_parse_ncbi_env_bioproject_ledger_old_format(tmp_path, caplog):
    """Verify processing stops if the ledger lacks the Assembly type column."""

    ledger = write_ledger(tmp_path / 'ncbi_env_bioprojects.tsv',
                          ['PRJEB31003\tDHP\t"Nov. 5, 2020"\t'],
                          header='NCBI BioProject\tSubmitter\tDate Added\t')

    with pytest.raises(SystemExit):
        parse_ncbi_env_bioproject_ledger(ledger)
    assert 'missing the `Assembly type` column' in caplog.records[-1].getMessage()
