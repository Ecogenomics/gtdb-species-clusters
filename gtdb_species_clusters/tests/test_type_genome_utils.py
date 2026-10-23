# Test functions in type_genome_utils

import logging
from types import SimpleNamespace

from gtdb_species_clusters.type_genome_utils import ClusteredGenome, check_ecoli_shigella_type_strains


def stub_genomes(genome_info):
    """Create genomes with an NCBI species and type strain status."""

    return {gid: SimpleNamespace(ncbi_taxa=SimpleNamespace(species=species),
                                 is_effective_type_strain=lambda is_type=is_type: is_type)
            for gid, (species, is_type) in genome_info.items()}


def clustered(*gids):
    return [ClusteredGenome(ani=97.0, af=80.0, gid=gid) for gid in gids]


GENOMES = stub_genomes({'EC_T': ('s__Escherichia coli', True),
                        'EC_2': ('s__Escherichia coli', False),
                        'SF_T': ('s__Shigella flexneri', True),
                        'SD_T': ('s__Shigella dysenteriae', True),
                        'SD_2': ('s__Shigella dysenteriae', False),
                        'OTHER_T': ('s__Escherichia fergusonii', True)})


def test_ecoli_shigella_single_cluster():
    """Verify all type strain genomes in the E. coli cluster passes, ignoring non-type and other species genomes."""

    clusters = {'EC_2': clustered('EC_T', 'SF_T', 'SD_T'),
                'SD_2': clustered(),
                'OTHER_T': clustered()}
    assert check_ecoli_shigella_type_strains(clusters, GENOMES)


def test_ecoli_shigella_split_clusters(caplog):
    """Verify type strain genomes split across clusters fails and reports each cluster."""

    clusters = {'EC_T': clustered('EC_2'),
                'SF_T': clustered('SD_T', 'SD_2')}
    with caplog.at_level(logging.ERROR, logger='rich'):
        assert not check_ecoli_shigella_type_strains(clusters, GENOMES)
    assert 'cluster of EC_T' in caplog.text
    assert 'cluster of SF_T' in caplog.text


def test_ecoli_shigella_unclustered_type_strain(caplog):
    """Verify a type strain genome that is not in any cluster fails."""

    clusters = {'EC_T': clustered('SF_T')}
    with caplog.at_level(logging.ERROR, logger='rich'):
        assert not check_ecoli_shigella_type_strains(clusters, GENOMES)
    assert 'no species cluster: SD_T' in caplog.text


def test_ecoli_shigella_shigella_representative():
    """Verify a single cluster represented by a Shigella genome fails."""

    clusters = {'SF_T': clustered('EC_T', 'SD_T')}
    assert not check_ecoli_shigella_type_strains(clusters, GENOMES)
