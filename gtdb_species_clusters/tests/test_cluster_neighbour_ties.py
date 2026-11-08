# Test that ties in selecting the closest neighbour of a genome do not depend on the
# order in which ANI results are read (u_cluster_named_reps and u_cluster_de_novo).

import logging
import pickle
from collections import namedtuple
from types import SimpleNamespace

import pytest

from gtdb_species_clusters.type_genome_utils import GenomeRadius, is_closer_neighbour
from gtdb_species_clusters.update_cluster_named_reps import UpdateClusterNamedReps
from gtdb_species_clusters.update_cluster_de_novo import UpdateClusterDeNovo


def make_step(cls, output_dir):
    """Create a step class without checking for the skani binary."""

    step = cls.__new__(cls)
    step.log = logging.getLogger('rich')
    step.ani_sp = 95.0
    step.af_sp = 50.0
    step.max_ani_neighbour = 97.0
    step.max_af_neighbour = 50.0
    step.output_dir = str(output_dir)
    return step


def test_is_closer_neighbour_order():
    # higher ANI wins regardless of AF and genome ID
    assert is_closer_neighbour(96.1, 60.0, 'G1', 96.0, 90.0, 'G9')
    assert not is_closer_neighbour(95.9, 90.0, 'G9', 96.0, 60.0, 'G1')

    # equal ANI: higher AF wins
    assert is_closer_neighbour(96.0, 85.0, 'G1', 96.0, 80.0, 'G9')
    assert not is_closer_neighbour(96.0, 80.0, 'G9', 96.0, 85.0, 'G1')

    # equal ANI and AF: higher genome ID wins
    assert is_closer_neighbour(96.0, 80.0, 'G9', 96.0, 80.0, 'G1')
    assert not is_closer_neighbour(96.0, 80.0, 'G1', 96.0, 80.0, 'G9')
    assert not is_closer_neighbour(96.0, 80.0, 'G1', 96.0, 80.0, 'G1')


def test_is_closer_neighbour_default():
    # a tie in ANI never displaces a default radius without an AF
    assert not is_closer_neighbour(95.0, 90.0, 'G1', 95.0, None, None)
    assert is_closer_neighbour(95.1, 90.0, 'G1', 95.0, None, None)

    # a tie in ANI and AF never displaces a value without a genome
    assert is_closer_neighbour(96.0, 85.0, 'G1', 96.0, 80.0, None)
    assert not is_closer_neighbour(96.0, 80.0, 'G1', 96.0, 80.0, None)


PAIRWISE_ROWS = [
    # rep A: B and C tie in ANI, C has the higher AF
    ('A', 'B', 96.0, 80.0, 70.0),
    ('A', 'C', 96.0, 70.0, 85.0),
    # rep D: E and F tie in ANI and AF
    ('D', 'E', 96.5, 80.0, 80.0),
    ('D', 'F', 96.5, 80.0, 80.0),
    # rep G: a tie with the default radius does not set a neighbour
    ('G', 'H', 95.0, 90.0, 90.0),
    # rep H: a higher ANI below the AF threshold is skipped
    ('H', 'B', 96.0, 40.0, 30.0),
]


@pytest.mark.parametrize('rows', [PAIRWISE_ROWS, PAIRWISE_ROWS[::-1]])
def test_rep_radius_ties(tmp_path, rows):
    rep_ani_file = tmp_path / 'gtdb_rep_pairwise_ani.tsv'
    with open(rep_ani_file, 'w') as f:
        f.write('Representative 1\tRepresentative 2\tANI\tAF12\tAF21\n')
        for r1, r2, ani, af12, af21 in rows:
            f.write(f'{r1}\t{r2}\t{ani}\t{af12}\t{af21}\n')

    step = make_step(UpdateClusterNamedReps, tmp_path)
    rep_radius = step._rep_radius(set('ABCDEFGH'), str(rep_ani_file))

    assert rep_radius['A'] == GenomeRadius(ani=96.0, af=85.0, neighbour_gid='C')
    assert rep_radius['D'] == GenomeRadius(ani=96.5, af=80.0, neighbour_gid='F')
    assert rep_radius['G'] == GenomeRadius(ani=95.0, af=None, neighbour_gid=None)
    assert rep_radius['H'] == GenomeRadius(ani=95.0, af=None, neighbour_gid=None)


@pytest.mark.parametrize('reverse', [False, True])
def test_cluster_named_reps_ties(tmp_path, reverse):
    hits = [('R1', (96.0, 80.0, 70.0)), ('R2', (96.0, 70.0, 80.0)), ('R0', (95.5, 99.0, 99.0))]
    if reverse:
        hits = hits[::-1]
    ani_af = {'Q': dict(hits)}
    rep_radius = {rid: GenomeRadius(ani=95.0, af=None, neighbour_gid=None) for rid in ('R0', 'R1', 'R2')}

    step = make_step(UpdateClusterNamedReps, tmp_path)
    step.ClusteredGenome = namedtuple('ClusteredGenome', 'ani af gid')
    clusters = step.cluster(ani_af, ['Q'], rep_radius)

    assert [g.gid for g in clusters['R2']] == ['Q']
    assert not clusters['R0'] and not clusters['R1']


@pytest.mark.parametrize('reverse', [False, True])
def test_nonrep_radius_ties(tmp_path, reverse):
    hits = [('R1', (96.0, 80.0, 70.0)), ('R2', (96.0, 70.0, 80.0)),
            ('R3', (96.0, 60.0, 60.0)), ('R4', (97.0, 40.0, 40.0))]
    if reverse:
        hits = hits[::-1]
    ani_af = {'Q': dict(hits),
              'T': {'R1': (95.0, 90.0, 90.0)}}
    ani_af_file = tmp_path / 'ani_af.pkl'
    with open(ani_af_file, 'wb') as f:
        pickle.dump(ani_af, f)

    step = make_step(UpdateClusterDeNovo, tmp_path)
    radius = step.nonrep_radius({'Q', 'T', 'U'}, {'R1', 'R2', 'R3', 'R4'}, str(ani_af_file))

    assert radius['Q'] == GenomeRadius(ani=96.0, af=80.0, neighbour_gid='R2')
    # a tie with the default radius does not set a neighbour
    assert radius['T'] == GenomeRadius(ani=95.0, af=50.0, neighbour_gid=None)
    assert radius['U'] == GenomeRadius(ani=95.0, af=50.0, neighbour_gid=None)


def stub_genomes(scores):
    return {gid: SimpleNamespace(score_type_strain=lambda score=score: score) for gid, score in scores.items()}


def test_select_rep_genomes_ties(tmp_path):
    # D1 and D2 become de novo reps; Q ties with both in ANI and AF and joins the higher genome ID
    genomes = stub_genomes({'D1': 100, 'D2': 99, 'Q': 50})
    nonrep_radius = {gid: GenomeRadius(ani=95.0, af=50.0, neighbour_gid=None) for gid in genomes}
    unclustered_ani_af = {'Q': {'D1': (97.0, 80.0, 80.0)},
                          'D2': {'Q': (97.0, 80.0, 80.0)}}

    step = make_step(UpdateClusterDeNovo, tmp_path)
    reps = step.select_rep_genomes(genomes, nonrep_radius, set(genomes), unclustered_ani_af)

    assert reps == {'D1', 'D2'}
    assert nonrep_radius['Q'] == GenomeRadius(ani=97.0, af=80.0, neighbour_gid='D2')


def test_select_rep_genomes_named_tie(tmp_path):
    # Q ties in ANI and AF between its closest named rep (N9) and the de novo rep D1;
    # the tie is not broken in favour of D1, so Q forms its own cluster as before
    genomes = stub_genomes({'D1': 100, 'Q': 50})
    nonrep_radius = {'D1': GenomeRadius(ani=95.0, af=50.0, neighbour_gid=None),
                     'Q': GenomeRadius(ani=96.0, af=80.0, neighbour_gid='N9')}
    unclustered_ani_af = {'Q': {'D1': (96.0, 80.0, 80.0)}}

    step = make_step(UpdateClusterDeNovo, tmp_path)
    reps = step.select_rep_genomes(genomes, nonrep_radius, set(genomes), unclustered_ani_af)

    assert reps == {'D1', 'Q'}
    assert nonrep_radius['Q'] == GenomeRadius(ani=96.0, af=80.0, neighbour_gid='N9')
