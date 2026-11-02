# Test handling of the processing pool in the Skani class

import signal

import pytest

from gtdb_species_clusters import skani as skani_module
from gtdb_species_clusters.skani import Skani, _ignore_sigint


def _sigint_handler(_data):
    """Report the SIGINT handler of a pool worker."""

    return signal.getsignal(signal.SIGINT) is signal.SIG_IGN


def _failing_worker(_data):
    """Worker that fails, as skani does when a genome file is missing."""

    raise RuntimeError('skani failed')


@pytest.fixture
def skani():
    s = Skani(None, cpus=2)
    yield s
    s.stop_pool()


def test_workers_ignore_sigint(skani):
    """Verify pool workers ignore Ctrl-C so only the main process handles it."""

    assert all(skani.pool.map(_sigint_handler, range(4)))
    assert signal.getsignal(signal.SIGINT) is not signal.SIG_IGN  # main process still handles it


def test_pool_stopped_when_pairs_fails(skani, monkeypatch):
    """Verify the pool is terminated if calculating ANI between pairs fails."""

    monkeypatch.setattr(skani_module, 'skani_worker', _failing_worker)

    with pytest.raises(RuntimeError, match='skani failed'):
        skani.pairs([('G1', 'G2')], {'G1': 'g1.fna', 'G2': 'g2.fna'}, None, report_progress=False)

    # pool is terminated so no workers are left behind
    for worker in skani.pool._pool:
        worker.join(timeout=5)
        assert not worker.is_alive()


def test_ignore_sigint_sets_handler():
    """Verify the pool initializer ignores SIGINT and the handler can be restored."""

    original = signal.getsignal(signal.SIGINT)
    try:
        _ignore_sigint()
        assert signal.getsignal(signal.SIGINT) is signal.SIG_IGN
    finally:
        signal.signal(signal.SIGINT, original)
