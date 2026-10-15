# Test methods in RepActions class

from gtdb_species_clusters.update_rep_actions import RepActions
from gtdb_species_clusters import defaults as Defaults


class StubSkani:
    """Record calls to skani."""

    def __init__(self):
        self.calls = []

    def pairs(self, gid_pairs, genome_files, preset, report_progress=True):
        self.calls.append((list(gid_pairs), genome_files, preset, report_progress))
        return {}


class StubGenomes:
    genomic_files = {'R1': 'r1.fna', 'R2': 'r2.fna', 'G1': 'g1.fna', 'G2': 'g2.fna', 'G3': 'g3.fna'}


def test_calculate_rep_anis(tmp_path):
    """Verify ANI between all representatives and their species genomes is calculated in a single batch."""

    rep_actions = RepActions(1, str(tmp_path))
    rep_actions.skani = StubSkani()

    rep_actions.calculate_rep_anis([('R1', {'G1', 'G2'}), ('R2', ['G3']), ('R3', [])], StubGenomes())

    assert len(rep_actions.skani.calls) == 1
    gid_pairs, genome_files, preset, report_progress = rep_actions.skani.calls[0]
    assert sorted(gid_pairs) == [('G1', 'R1'), ('G2', 'R1'), ('G3', 'R2')]
    assert genome_files is StubGenomes.genomic_files
    assert preset == Defaults.SKANI_PRESET
    assert report_progress
