# Test methods in RepActions class

from types import SimpleNamespace

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


class StubSpeciesGenomes(dict):
    """Genomes with an NCBI species and a fixed ANI score."""

    def __init__(self, genome_species, genome_scores):
        super().__init__()
        self.genomic_files = {}
        for gid, species in genome_species.items():
            self.genomic_files[gid] = f'{gid}.fna'
            self[gid] = SimpleNamespace(ncbi_taxa=SimpleNamespace(species=species),
                                        genomic_file=f'{gid}.fna',
                                        score_ani=lambda _ani, score=genome_scores[gid]: score,
                                        is_gtdb_type_strain=lambda: False)


def test_rep_candidates_exclude_shigella(tmp_path):
    """Verify genomes from Shigella species are never candidate representatives."""

    rep_actions = RepActions(1, str(tmp_path))
    cur_genomes = StubSpeciesGenomes({'G1': 's__Escherichia coli',
                                      'G2': 's__Shigella flexneri',
                                      'G3': 's__',
                                      'G4': 's__Shigella sonnei'},
                                     {'G1': 1, 'G2': 1, 'G3': 1, 'G4': 1})

    assert sorted(rep_actions.rep_candidates(['G1', 'G2', 'G3', 'G4'], cur_genomes)) == ['G1', 'G3']


def test_top_ani_score_skips_shigella(tmp_path):
    """Verify a lower scoring genome is selected over a higher scoring Shigella genome."""

    rep_actions = RepActions(1, str(tmp_path))
    rep_actions.skani = StubSkani()
    cur_genomes = StubSpeciesGenomes({'R1': 's__Escherichia coli',
                                      'G1': 's__Escherichia coli',
                                      'G2': 's__Shigella flexneri'},
                                     {'R1': 50, 'G1': 60, 'G2': 90})

    new_rid, top_score, _ani, _af = rep_actions.top_ani_score('R1', ['G1', 'G2'], cur_genomes)
    assert (new_rid, top_score) == ('G1', 60)

    gid_pairs = rep_actions.skani.calls[0][0]
    assert gid_pairs == [('G1', 'R1')]


def test_top_ani_score_prev_rep_only_shigella(tmp_path):
    """Verify no representative is returned when a species cluster only contains Shigella genomes."""

    rep_actions = RepActions(1, str(tmp_path))
    rep_actions.skani = StubSkani()
    prev_genomes = StubSpeciesGenomes({'R1': 's__Shigella flexneri'}, {'R1': 50})
    cur_genomes = StubSpeciesGenomes({'G1': 's__Shigella flexneri',
                                      'G2': 's__Shigella boydii'},
                                     {'G1': 60, 'G2': 70})

    new_rid, _top_score, _ani, _af = rep_actions.top_ani_score_prev_rep('R1', ['G1', 'G2'], prev_genomes, cur_genomes)
    assert new_rid is None
