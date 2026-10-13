###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import re
import sys
import shutil
import signal
import socket
import getpass
import logging
import subprocess
from collections import Counter
from datetime import datetime
from typing import Dict, List, Optional, Tuple

from gtdblib.util.shell.execute import check_dependencies

from gtdb_species_clusters.genome_utils import read_qc_file


CANARY_FILE = 'canary.running'
DONE_FILE = 'done'

BATCH_FILE_RE = re.compile(r'^genomes_(\d+)\.lst$')


def _exit_on_signal(signum, _frame):
    """Convert a termination signal into SystemExit so cleanup code runs."""

    raise SystemExit(128 + signum)


class GTDB_Tk():
    """Perform initial classification of new and updated genomes using GTDB-Tk."""

    def __init__(self, cpus: int, output_dir: str):
        """Initialization.

        :param cpus: Number of CPUs to be used by GTDB-Tk.
        :param output_dir: Output directory.
        """

        check_dependencies(['gtdbtk'])

        self.output_dir = output_dir
        self.log = logging.getLogger('rich')

        self.cpus = cpus
        if self.cpus > 64:
            self.log.error(
                'Testing indicates pplacer will stale if used with more than 64 CPUs.')
            sys.exit(1)

        self.host = socket.gethostname()
        self.pid = os.getpid()

    def _summary_files(self, batch_dir: str) -> List[str]:
        """Get GTDB-Tk summary files in batch output directory.

        GTDB-Tk only writes a summary file for a domain if the batch
        contains genomes from that domain.

        :param batch_dir: GTDB-Tk output directory for a batch of genomes.

        :return: Existing archaeal and bacterial summary files.
        """

        summary_files = []
        for summary_file in ['gtdbtk.ar53.summary.tsv', 'gtdbtk.bac120.summary.tsv']:
            summary_file = os.path.join(batch_dir, summary_file)
            if os.path.exists(summary_file):
                summary_files.append(summary_file)

        return summary_files

    @staticmethod
    def _write_key_values(f, values: Dict[str, str]) -> None:
        """Write key/value pairs as tab-separated lines."""

        for key, value in values.items():
            f.write(f'{key}\t{value}\n')

    @staticmethod
    def _read_key_values(path: str) -> Optional[Dict[str, str]]:
        """Read tab-separated key/value pairs, or None if the file does not exist."""

        try:
            with open(path, encoding='utf-8') as f:
                content = f.read()
        except FileNotFoundError:
            return None

        values = {}
        for line in content.splitlines():
            if '\t' in line:
                key, value = line.split('\t', 1)
                values[key] = value

        return values

    def _batch_status(self, batch_dir: str) -> Tuple[str, Optional[Dict[str, str]]]:
        """Determine status of a batch.

        A batch is done if it has a done file. Batches processed before
        done files were introduced are also considered done if they have
        GTDB-Tk summary files and are not being processed.

        :param batch_dir: GTDB-Tk output directory for a batch of genomes.

        :return: Status ('done', 'running', or 'pending') and contents of the canary file for running batches.
        """

        if os.path.exists(os.path.join(batch_dir, DONE_FILE)):
            return 'done', None

        canary = self._read_key_values(os.path.join(batch_dir, CANARY_FILE))
        if canary is not None:
            return 'running', canary

        # the canary may have been removed after the done file was checked
        if os.path.exists(os.path.join(batch_dir, DONE_FILE)) or self._summary_files(batch_dir):
            return 'done', None

        return 'pending', None

    def _is_stale_canary(self, canary: Dict[str, str]) -> bool:
        """Check if a canary was left by a process on this host that is no longer running."""

        if canary.get('host') != self.host:
            return False

        try:
            pid = int(canary.get('pid', ''))
        except ValueError:
            return False

        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return True
        except PermissionError:
            # process exists, but is owned by another user
            return False

        return False

    def _remove_stale_canary(self, canary_file: str, stale_canary: Dict[str, str]) -> None:
        """Remove a stale canary file.

        The canary is first renamed so only one process can remove it. If the
        renamed file is not the stale canary that was read, another process
        claimed the batch in the meantime and the canary is restored.
        """

        tmp_file = f'{canary_file}.stale.{self.host}.{self.pid}'
        try:
            os.rename(canary_file, tmp_file)
        except FileNotFoundError:
            return

        if self._read_key_values(tmp_file) != stale_canary:
            try:
                os.link(tmp_file, canary_file)
            except FileExistsError:
                pass

        os.remove(tmp_file)

    def _claim_batch(self, canary_file: str, cmd: str) -> bool:
        """Atomically create the canary file for a batch.

        :return: True if the batch was claimed by this process, False if it was already claimed.
        """

        try:
            fd = os.open(canary_file, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            return False

        with os.fdopen(fd, 'w', encoding='utf-8') as f:
            self._write_key_values(f, {
                'host': self.host,
                'pid': self.pid,
                'user': getpass.getuser(),
                'started': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'command': cmd})

        return True

    def _run_gtdbtk(self, cmd: List[str]) -> bool:
        """Run GTDB-Tk, terminating it if this process is interrupted.

        :return: True if GTDB-Tk completed successfully.
        """

        # GTDB-Tk is run in its own process group so it and its
        # child processes (e.g., pplacer) can be terminated together
        proc = subprocess.Popen(cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                encoding='utf-8',
                                start_new_session=True)
        try:
            _stdout, stderr = proc.communicate()
        except BaseException:
            # make sure GTDB-Tk is not left writing to a batch that
            # will be released for processing by another instance
            self.log.warning('Terminating GTDB-Tk.')
            try:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=60)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
            except ProcessLookupError:
                pass
            raise

        if proc.returncode != 0:
            self.log.error(stderr)
            self.log.error(f'GTDB-Tk failed with return code {proc.returncode}.')
            return False

        return True

    def _read_batch_files(self, batch_dir: str) -> List[Tuple[int, str]]:
        """Get genome batch files sorted by batch index."""

        batch_files = []
        for f in os.listdir(batch_dir):
            m = BATCH_FILE_RE.match(f)
            if m:
                batch_files.append((int(m.group(1)), os.path.join(batch_dir, f)))

        return sorted(batch_files)

    def _create_batch_files(self,
                            batch_dir: str,
                            genomic_files: List[Tuple[str, str]],
                            new_updated_gids: set,
                            batch_size: int) -> List[Tuple[int, str]]:
        """Create genome batch files, or read existing batch files.

        Batch files are written to a temporary directory that is renamed once
        complete, so other instances never see a partial set of batch files.

        :return: Batch index and path of each genome batch file.
        """

        if not os.path.exists(batch_dir):
            tmp_dir = f'{batch_dir}.tmp.{self.host}.{self.pid}'
            os.makedirs(tmp_dir)
            for batch_idx, start in enumerate(range(0, len(genomic_files), batch_size)):
                genome_batch_file = os.path.join(
                    tmp_dir, f'genomes_{batch_idx}.lst')

                fout = open(genome_batch_file, 'w')
                for i in range(start, min(start+batch_size, len(genomic_files))):
                    gid, gf = genomic_files[i]
                    fout.write('{}\t{}\n'.format(gf, gid))
                fout.close()

            try:
                os.rename(tmp_dir, batch_dir)
                return self._read_batch_files(batch_dir)
            except OSError:
                # another instance created the batch files first
                shutil.rmtree(tmp_dir)

        self.log.warning(
            f'Using existing genome batch files in {batch_dir}.')
        genome_batch_files = self._read_batch_files(batch_dir)

        # check if there are genomes not already in a batch file. Ideally,
        # this would never happen, but sometimes we process past this step
        # and then identify genomes missing in the database. These need to
        # be put into a batch file for processing.
        missing_gids = set(new_updated_gids)
        last_batch_idx = 0
        for batch_idx, batch_file in genome_batch_files:
            last_batch_idx = max(last_batch_idx, batch_idx)

            with open(batch_file) as f:
                for line in f:
                    tokens = line.strip().split('\t')
                    missing_gids.discard(tokens[1])

        if len(missing_gids) > 0:
            batch_idx = last_batch_idx + 1
            genome_batch_file = os.path.join(
                batch_dir, f'genomes_{batch_idx}.lst')

            tmp_file = os.path.join(
                batch_dir, f'.genomes_{batch_idx}.lst.tmp.{self.host}.{self.pid}')
            fout = open(tmp_file, 'w')
            for gid, gf in genomic_files:
                if gid in missing_gids:
                    fout.write('{}\t{}\n'.format(gf, gid))
            fout.close()
            os.replace(tmp_file, genome_batch_file)

            genome_batch_files.append((batch_idx, genome_batch_file))
            self.log.info('Added the batch file {} with {:,} genomes.'.format(
                genome_batch_file,
                len(missing_gids)))

        return genome_batch_files

    def _process_batch(self, batch_idx: int, genome_batch_file: str) -> str:
        """Process a batch of genomes with GTDB-Tk unless it is done or being processed by another instance.

        :return: 'processed' if this instance processed the batch, otherwise 'skipped'.
        """

        out_dir = os.path.join(self.output_dir, f'gtdbtk_batch{batch_idx}')
        canary_file = os.path.join(out_dir, CANARY_FILE)
        done_file = os.path.join(out_dir, DONE_FILE)

        status, canary = self._batch_status(out_dir)
        if status == 'done':
            self.log.info(
                f'Skipping genome batch {batch_idx} as it has been processed.')
            return 'skipped'

        if status == 'running':
            if not self._is_stale_canary(canary):
                self.log.warning('Skipping genome batch {} as it is being processed on {} (PID {}, started {}). '
                                 'If this process is no longer running, delete {} and rerun u_gtdbtk.'.format(
                                     batch_idx,
                                     canary.get('host', 'unknown host'),
                                     canary.get('pid', '?'),
                                     canary.get('started', '?'),
                                     out_dir))
                return 'skipped'

            self.log.warning('Reprocessing genome batch {} as PID {} on this host is no longer running.'.format(
                batch_idx,
                canary.get('pid')))
            self._remove_stale_canary(canary_file, canary)

        os.makedirs(out_dir, exist_ok=True)
        cmd = 'gtdbtk classify_wf --cpus {} --pplacer_cpus {} --force --batchfile {} --out_dir {}'.format(
            self.cpus,
            self.cpus,
            genome_batch_file,
            out_dir)

        if not self._claim_batch(canary_file, cmd):
            self.log.warning(
                f'Skipping genome batch {batch_idx} as it was claimed by another instance.')
            return 'skipped'

        succeeded = False
        try:
            # another instance may have finished the batch between
            # checking its status and claiming it
            if os.path.exists(done_file):
                self.log.info(
                    f'Skipping genome batch {batch_idx} as it has been processed.')
                return 'skipped'

            if len(os.listdir(out_dir)) > 1:
                self.log.warning(
                    f'Rerunning genome batch {batch_idx} as output directory has incomplete GTDB-Tk results.')

            self.log.info(f'Executing: {cmd}')
            succeeded = self._run_gtdbtk(cmd.split())
            if succeeded:
                with open(done_file, 'w', encoding='utf-8') as f:
                    canary = self._read_key_values(canary_file) or {}
                    self._write_key_values(f, {
                        'host': self.host,
                        'started': canary.get('started', '?'),
                        'finished': datetime.now().strftime('%Y-%m-%d %H:%M:%S')})
        finally:
            # done file must be written before the canary is removed so
            # other instances never see a finished batch as pending
            try:
                os.remove(canary_file)
            except FileNotFoundError:
                pass

        if not succeeded:
            self.log.error(
                f'Failed to process genome batch {batch_idx}. It will be rerun by the next instance of u_gtdbtk.')
            sys.exit(1)

        return 'processed'

    def _report_progress(self, genome_batch_files: List[Tuple[int, str]]) -> bool:
        """Report status of all batches.

        :return: True if all batches are done.
        """

        status_count = Counter()
        running_hosts = Counter()
        for batch_idx, _ in genome_batch_files:
            out_dir = os.path.join(self.output_dir, f'gtdbtk_batch{batch_idx}')
            status, canary = self._batch_status(out_dir)
            status_count[status] += 1
            if status == 'running':
                running_hosts[canary.get('host', 'unknown host')] += 1

        running_str = ''
        if running_hosts:
            running_str = ' ({})'.format(', '.join(
                f'{host}: {count:,}' for host, count in sorted(running_hosts.items())))

        self.log.info('Genome batches: {:,} done, {:,} running{}, {:,} remaining.'.format(
            status_count['done'],
            status_count['running'],
            running_str,
            status_count['pending']))

        return status_count['done'] == len(genome_batch_files)

    def _combine_summary_files(self,
                               genome_batch_files: List[Tuple[int, str]],
                               new_updated_gids: set) -> None:
        """Combine GTDB-Tk summary files across all batches."""

        output_file = os.path.join(self.output_dir, 'gtdbtk_classify.tsv')
        tmp_file = f'{output_file}.tmp.{self.host}.{self.pid}'
        fout = open(tmp_file, 'w')
        bHeader = True
        gtdbtk_processed = set()
        for batch_idx, _ in genome_batch_files:
            batch_dir = os.path.join(self.output_dir, f'gtdbtk_batch{batch_idx}')
            summary_files = self._summary_files(batch_dir)
            if not summary_files:
                fout.close()
                os.remove(tmp_file)
                self.log.error(
                    f'No GTDB-Tk summary files found in {batch_dir}. Delete this directory and rerun u_gtdbtk.')
                sys.exit(1)

            for summary_file in summary_files:
                with open(summary_file, encoding='utf-8') as f:
                    header = f.readline()

                    if bHeader:
                        fout.write(header)
                        bHeader = False

                    for line in f:
                        tokens = line.strip().split('\t')
                        gid = tokens[0]
                        if gid in new_updated_gids:
                            # Ideally, this shouldn't be necessary, but
                            # sometimes we process past this step and then
                            # identify genomes missing in the database. This
                            # can result in GTDB-Tk having been applied to
                            # genomes that looked like they were "new", but
                            # really were just erroneously missing from the
                            # database.
                            fout.write(line)
                            gtdbtk_processed.add(gid)

        fout.close()

        # several instances may finish at the same time, so write
        # the combined file atomically
        os.replace(tmp_file, output_file)

        self.log.info(
            'Identified {:,} genomes as being processed by GTDB-Tk.'.format(len(gtdbtk_processed)))
        skipped_gids = new_updated_gids - gtdbtk_processed
        if len(skipped_gids) > 0:
            self.log.warning('Identified {:,} genomes as being skipped by GTDB-Tk.'.format(
                len(skipped_gids)))

    def run(self,
            genomes_new_updated_file: str,
            qc_passed_file: str,
            batch_size: int) -> None:
        """Perform initial classification of new and updated genomes using GTDB-Tk.

        Genomes are processed in batches. Multiple instances can be run at the same
        time on different servers. Each instance claims a batch by creating a canary
        file in the batch output directory, and creates a done file once the batch
        is processed. Batches with a canary or done file are skipped by other
        instances. Use of batches also minimizes the impact of server crashes.
        Results are combined once all batches are done.

        :param genomes_new_updated_file: File indicating new and updated genomes in the current GTDB release.
        :param qc_passed_file: File indicating genomes that have passed QC.
        :param batch_size: Number of genomes to process in a single batch.
        """

        self.log.info(f'Running on {self.host} (PID {self.pid}).')

        # get list of genomes passing QC
        self.log.info('Reading genomes passing QC.')
        gids_pass_qc = read_qc_file(qc_passed_file)
        self.log.info(f' - identified {len(gids_pass_qc):,} genomes')

        # get path to genomes passing QC
        self.log.info(
            'Reading path to genomic file for new/updated genomes passing QC.')
        genomic_files = []
        new_updated_gids = set()
        total_count = 0
        with open(genomes_new_updated_file, encoding='utf-8') as f:
            header = f.readline().strip().split('\t')

            genomic_file_index = header.index('Genomic file')

            for line in f:
                tokens = line.strip().split('\t')

                gid = tokens[0]
                total_count += 1
                if gid in gids_pass_qc:
                    gf = tokens[genomic_file_index]
                    genomic_files.append((gid, gf))
                    new_updated_gids.add(gid)
        self.log.info(
            f' - identified {len(genomic_files):,} of {total_count:,} genomes as passing QC')

        # create batch files
        batch_dir = os.path.join(self.output_dir, 'genome_batch_files')
        genome_batch_files = self._create_batch_files(batch_dir,
                                                      genomic_files,
                                                      new_updated_gids,
                                                      batch_size)

        # process genomes with GTDB-Tk in batches; termination signals are
        # converted to SystemExit so the canary of the current batch is removed
        prev_handlers = {}
        for signum in [signal.SIGTERM, signal.SIGHUP]:
            if signal.getsignal(signum) == signal.SIG_DFL:
                prev_handlers[signum] = signal.signal(signum, _exit_on_signal)

        processed_count = 0
        try:
            for batch_idx, genome_batch_file in genome_batch_files:
                if self._process_batch(batch_idx, genome_batch_file) == 'processed':
                    processed_count += 1
        finally:
            for signum, handler in prev_handlers.items():
                signal.signal(signum, handler)

        self.log.info(
            f'Processed {processed_count:,} genome batches on {self.host}.')

        # combine summary files once all batches are done
        if not self._report_progress(genome_batch_files):
            self.log.warning(
                'Not all genome batches are done. Rerun u_gtdbtk once they are complete to combine the GTDB-Tk results.')
            return

        self._combine_summary_files(genome_batch_files, new_updated_gids)
