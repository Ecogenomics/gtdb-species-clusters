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
import sys
import ntpath
import logging
import re
import subprocess
import sqlite3
import time
import datetime
import multiprocessing as mp
from itertools import combinations
from collections import defaultdict
from typing import Dict

from gtdblib.util.shell.execute import check_dependencies
from gtdblib.util.bio.accession import canonical_gid


class ANIError(Exception):
    """Exception for failed execution of ANI/AF."""
    pass


def skani(qid: str, rid: str, q_gf: str, r_gf: str, preset: str = '--medium'):
    """CalculateANI between a pair of genomes."""

    # run skani and write results to stdout
    try:
        cmd = ['skani', 'dist',
                '-q', q_gf,
                '-r', r_gf,
                '-o', '/dev/stdout']

        if preset is not None:
            cmd += [preset]

        proc = subprocess.Popen(cmd,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    encoding='utf-8')
        stdout, stderr = proc.communicate()

        if proc.returncode != 0:  # skani returned an error code
            print(stderr)
            raise ANIError(f"skani exited with code {proc.returncode}.")
    except Exception as e:
        print(e)
        raise

    result_lines = stdout.splitlines()
    if len(result_lines) == 2:
        tokens = result_lines[1].split('\t')
        ani = float(tokens[2])
        af_r = float(tokens[3])
        af_q = float(tokens[4])
        ani_af = (qid, rid, ani, af_r, af_q)
    else:
        # genomes too divergent to determine ANI and AF
        # with skani so default to zeros
        ani_af = (qid, rid, 0.0, 0.0, 0.0)

    return ani_af


def skani_worker(data):
    """Worker for running skani in parallel."""

    qid, rid, q_gf, r_gf, preset = data
    return skani(qid, rid, q_gf, r_gf, preset)


class Skani():
    """Calculate ANI between genomes using skani with an optional precomputed cache."""

    @staticmethod
    def symmetric_ani_af(ani_af, gid1, gid2):
        """Calculate symmetric ANI and AF statistics between genomes."""

        if gid1 == gid2:
            return 100.0, 100.0

        if gid1 in ani_af and gid2 in ani_af[gid1]:
            ani, af_r, af_q = ani_af[gid1][gid2]
            return ani, max(af_r, af_q)
        elif gid2 in ani_af and gid1 in ani_af[gid2]:
            ani, af_r, af_q = ani_af[gid2][gid1]
            return ani, max(af_r, af_q)

        return 0.0, 0.0

    @staticmethod
    def result(ani_af, gid1, gid2):
        """Calculate ANI and both AF values between genomes.
        
        AF results are return for gid1 vs gid2 followed by gid2 vs. gid1.
        """

        if gid1 == gid2:
            return 100.0, 100.0, 100.0

        if gid1 in ani_af and gid2 in ani_af[gid1]:
            ani, af_r, af_q = ani_af[gid1][gid2]
            return ani, af_r, af_q
        elif gid2 in ani_af and gid1 in ani_af[gid2]:
            ani, af_r, af_q = ani_af[gid2][gid1]
            return ani, af_q, af_r

        return 0.0, 0.0, 0.0

    def __init__(self, ani_cache_file, cpus, read_cache=True):
        """Initialization."""

        check_dependencies(['skani'])

        # do database insertions in batches to reduce number of transactions
        self.DB_BATCH_SIZE = 100

        self.cpus = cpus
        self.log = logging.getLogger('rich')
        self.ani_cache_file = ani_cache_file

        self._create_db()

        self.ani_af_cache = defaultdict(dict)
        if ani_cache_file and read_cache:
            self.log.info('Reading ANI cache into memory for fast access:')
            self.ani_af_cache = self._read_cache()
            self.log.info(' - done')

        self.pool = mp.Pool(self.cpus)

        self.log.info('Using skani v{}.'.format(self.get_version()))

    def __del__(self):
        """Destructor."""

        self._write_cache()
        
        if self.db_conn is not None:
            self.db_conn.close()

        self.pool.close()

    def get_version(self):
        """Returns the version of skani on the system path.

        Returns
        -------
        str
            The string containing the skani version.
        """
        try:
            proc = subprocess.Popen(['skani', '-V'], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, encoding='utf-8')
            stdout, stderr = proc.communicate()
            version = re.search(r'skani (.+)', stdout)
            return version.group(1)
        except Exception as e:
            print(e)
            return 'unknown'

    def _create_db(self):
        """Read previously calculated ANI values."""

        if self.ani_cache_file:
            self.log.info(f'Connecting to ANI DB: {self.ani_cache_file}')
            try:
                self.db_conn = sqlite3.connect(self.ani_cache_file)
            except Exception as e:
                print(e)
                raise

            self.db_cur = self.db_conn.cursor()
            self.db_cur.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table' AND name = 'ani_table'")

            if not self.db_cur.fetchone():
                self.log.info(f' - creating ANI database table')

                self.db_cur.execute('''CREATE TABLE ani_table
                (query_id TEXT NOT NULL,
                ref_id TEXT NOT NULL,
                ani REAL NOT NULL,
                af_r REAL NOT NULL,
                af_q REAL NOT NULL)''')

                self.db_cur.execute(
                    'CREATE INDEX gid_idx ON ani_table(query_id, ref_id)')

                self.db_conn.commit()
            else:
                self.log.info(
                    f' - database contains {self.cache_size():,} entries')
        else:
            self.log.info('Not using an ANI database.')
            self.db_conn = None

    def _read_cache(self):
        """Read all ANI and AF values in the current cache."""

        ani_af = defaultdict(dict)

        if self.db_conn:
            self.db_cur.execute(
                'SELECT query_id, ref_id, ani, af_r, af_q FROM ani_table')
            for row in self.db_cur:
                qid, rid, ani, af_r, af_q = row
                ani_af[qid][rid] = (ani, af_r, af_q)

        return ani_af

    def _write_cache(self):
        """Write cache to file."""

        if self.db_conn:
            self.db_conn.commit()

    def cache_size(self):
        """Get number of rows in the ANI cache database."""

        if self.db_conn:
            self.db_cur.execute("SELECT max(rowid) from ani_table")
            num_rows = self.db_cur.fetchone()[0]
            if num_rows is None:
                # DB could be empty in which case num_rows will be None
                return 0

            return num_rows

        return 0

    def in_cache(self, qid, rid):
        """Check if ANI comparisons is already in the cache."""

        ani_af = None
        if self.db_conn:
            self.db_cur.execute(
                'SELECT ani, af_r, af_q FROM ani_table WHERE query_id=? AND ref_id=?', (qid, rid,))
            ani_af = self.db_cur.fetchone()

        return ani_af

    def get_genome_id(self, genome_path):
        """Extract genome ID from path to genomic file."""

        genome_id = ntpath.basename(genome_path)

        return canonical_gid(genome_id)

    def dist(self,
                query_paths: Dict[str, str],
                ref_paths: Dict[str, str],
                output_dir: str, 
                preset: str,
                min_af: float = 50,
                min_sketch_ani:float = 85):
        """Calculate ANI/AF between query and reference genomes."""

        os.makedirs(output_dir, exist_ok=True)

        self.log.info(f'Calculating ANI between {len(query_paths):,} query and {len(ref_paths):,} reference genomes:')

        try:
            # create files with path to genomes
            self.log.info(' - create file with path to query genomes')
            qry_path_file = os.path.join(output_dir, 'skani_query_genome_paths.tsv')
            fout = open(qry_path_file, 'w')
            for gp in query_paths.values():
                fout.write(f'{gp}\n')
            fout.close()

            self.log.info(' - create file with path to reference genomes')
            ref_path_file = os.path.join(output_dir, 'skani_ref_genome_paths.tsv')
            fout = open(ref_path_file, 'w')
            for gp in ref_paths.values():
                fout.write(f'{gp}\n')
            fout.close()
            
            # compute ANI/AF between genome pairs in lower triangle
            self.log.info(' - calculating ANI and AF between query and reference genomes')
            results_file = os.path.join(output_dir, 'skani_dist.tsv')
            cmd = ['skani', 'dist',
                    '-t', str(self.cpus), 
                    '--ql', qry_path_file,
                    '--rl', ref_path_file,
                    '-o', results_file,
                    '--min-af', str(min_af),
                    '-s', str(min_sketch_ani)]

            if preset is not None:
                cmd += [preset]

            proc = subprocess.Popen(cmd,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    encoding='utf-8')
            _stdout, stderr = proc.communicate()
            if proc.returncode != 0:  # skani returned an error code
                print(stderr)
                raise ANIError(f"skani dist exited with code {proc.returncode}.")

        except Exception as e:
            print(e)
            raise

        # get map from genome path to genome ID
        gp_to_gid_map = {}
        for gid, gp in qry_path_file.items():
            gp_to_gid_map[gp] = gid

        for gid, gp in ref_paths.items():
            gp_to_gid_map[gp] = gid

        # parse skani results
        self.log.info(' - parsing skani results')
        ani_af = defaultdict(dict)
        with open(results_file) as f:
            header = f.readline().strip().split('\t')

            ref_file_idx = header.index('Ref_file')
            query_file_idx = header.index('Query_file')
            ani_idx = header.index('ANI')
            af_r_idx = header.index('Align_fraction_ref')
            af_q_idx = header.index('Align_fraction_query')

            for line in f:
                tokens = line.strip().split('\t')

                rid = gp_to_gid_map[tokens[ref_file_idx]]
                qid = gp_to_gid_map[tokens[query_file_idx]]
                ani = float(tokens[ani_idx])
                af_r = float(tokens[af_r_idx])
                af_q = float(tokens[af_q_idx])

                ani_af[qid][rid] = (ani, af_r, af_q)

        return ani_af

    def triangle(self, genome_paths: Dict[str, str], 
                        output_dir: str, 
                        preset: str,
                        min_af: float = 50,
                        min_sketch_ani:float = 85):
        """Calculate ANI/AF between genome pairs in lower triangle.
        
        This method is fast, but has large memory requirements if large
        numbers of genomes need to be compared.
        """

        os.makedirs(output_dir, exist_ok=True)

        self.log.info(f'Calculating ANI between {len(genome_paths):,} genomes in lower triangle ({preset}; min_af = {min_af}; s = {min_sketch_ani}):')

        try:
            # create file with path to genomes
            self.log.info(' - create file with path to genomes')
            genome_path_file = os.path.join(output_dir, 'skani_genome_paths.tsv')
            fout = open(genome_path_file, 'w')
            for gp in genome_paths.values():
                fout.write(f'{gp}\n')
            fout.close()
            
            # compute ANI/AF between genome pairs in lower triangle
            self.log.info(' - calculating ANI and AF between genomes')
            results_file = os.path.join(output_dir, 'skani_dist.tsv')
            cmd = ['skani', 'triangle',
                    '-t', str(self.cpus), 
                    '-l', genome_path_file,
                    '-o', results_file,
                    '--sparse',
                    '--min-af', str(min_af),
                    '-s', str(min_sketch_ani)]

            if preset is not None:
                cmd += [preset]

            proc = subprocess.Popen(cmd,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    encoding='utf-8')
            _stdout, stderr = proc.communicate()
            if proc.returncode != 0:  # skani returned an error code
                print(stderr)
                raise ANIError(f"skani triangle exited with code {proc.returncode}.")

        except Exception as e:
            print(e)
            raise

        # get map from genome path to genome ID
        gp_to_gid_map = {}
        for gid, gp in genome_paths.items():
            gp_to_gid_map[gp] = gid

        # parse skani results
        self.log.info(' - parsing skani results')
        ani_af = defaultdict(dict)
        with open(results_file) as f:
            header = f.readline().strip().split('\t')

            ref_file_idx = header.index('Ref_file')
            query_file_idx = header.index('Query_file')
            ani_idx = header.index('ANI')
            af_r_idx = header.index('Align_fraction_ref')
            af_q_idx = header.index('Align_fraction_query')

            for line in f:
                tokens = line.strip().split('\t')

                rid = gp_to_gid_map[tokens[ref_file_idx]]
                qid = gp_to_gid_map[tokens[query_file_idx]]
                ani = float(tokens[ani_idx])
                af_r = float(tokens[af_r_idx])
                af_q = float(tokens[af_q_idx])

                ani_af[qid][rid] = (ani, af_r, af_q)

        return ani_af

    def search(self, 
                query_paths: Dict[str, str], 
                ref_paths: Dict[str, str],
                output_dir: str, 
                preset: str,
                min_af: float = 50,
                min_sketch_ani:float = 85):
        """Calculate skani between query and reference genomes.
        
        Since skani is symmetric each pair of genomes is only processed once. That is,
        genome A and genome B will be processed with skani as (A,B) or (B,A), but not
        both since they will give the same result. Representative genomes are first
        sketched and then the query genomes searched against this DB. This is slower
        than using "dist" but keeps memory usage much more reasonable.
        """

        os.makedirs(output_dir, exist_ok=True)

        log_str = f'Calculating ANI between {len(query_paths):,} query and {len(ref_paths):,} reference genomes'
        log_str += f' ({preset}; min_af = {min_af}; s = {min_sketch_ani}):'
        self.log.info(log_str)

        try:
            # create file with path to reference genomes
            self.log.info(' - create file with path to reference genomes')
            rep_path_file = os.path.join(output_dir, 'skani_ref_genome_paths.tsv')
            fout = open(rep_path_file, 'w')
            for gp in ref_paths.values():
                fout.write(f'{gp}\n')
            fout.close()

            # create file with path to query genomes
            self.log.info(' - creating file with path to query genomes')
            query_path_file = os.path.join(output_dir, 'skani_query_genome_paths.tsv')
            fout = open(query_path_file, 'w')
            for gid, gp in query_paths.items():
                if gid in ref_paths and ref_paths[gid] != gp:
                    self.log.error(f'Genome ID {gid} has different genomic FASTA files for query and reference.')
                    sys.exit(1)
                fout.write(f'{gp}\n')
            fout.close()
            
            # sketch the reference genomes
            self.log.info(' - sketching reference genomes')
            sketch_dir = os.path.join(output_dir, 'skani_ref_sketches')
            cmd = ['skani', 'sketch',
                    '-t', str(self.cpus), 
                    '-l', rep_path_file,
                    '-o', sketch_dir]

            if preset is not None:
                cmd += [preset]

            proc = subprocess.Popen(cmd,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        encoding='utf-8')
            _stdout, stderr = proc.communicate()
            if proc.returncode != 0:  # skani returned an error code
                print(stderr)
                raise ANIError(f"skani sketch exited with code {proc.returncode}.")
            
            # search query genomes against reference genome sketches
            self.log.info(' - calculating ANI and AF between reference and query genomes')
            results_file = os.path.join(output_dir, 'skani_query_vs_ref.tsv')
            cmd = ['skani', 'search',
                    '-t', str(self.cpus), 
                    '-d', sketch_dir,
                    '--ql', query_path_file,
                    '-o', results_file,
                    '--min-af', str(min_af),
                    '-s', str(min_sketch_ani)]
            proc = subprocess.Popen(cmd,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        encoding='utf-8')
            _stdout, stderr = proc.communicate()
            if proc.returncode != 0:  # skani returned an error code
                print(stderr)
                raise ANIError(f"skani search exited with code {proc.returncode}.")
        except Exception as e:
            print(e)
            raise

        # get map from genome path to genome ID
        gp_to_gid_map = {}
        for gid, gp in query_paths.items():
            gp_to_gid_map[gp] = gid
        for gid, gp in ref_paths.items():
            gp_to_gid_map[gp] = gid

        # parse skani results
        self.log.info(' - parsing skani results')
        ani_af = defaultdict(dict)
        with open(results_file) as f:
            header = f.readline().strip().split('\t')

            ref_file_idx = header.index('Ref_file')
            query_file_idx = header.index('Query_file')
            ani_idx = header.index('ANI')
            af_r_idx = header.index('Align_fraction_ref')
            af_q_idx = header.index('Align_fraction_query')

            for line in f:
                tokens = line.strip().split('\t')

                rid = gp_to_gid_map[tokens[ref_file_idx]]
                qid = gp_to_gid_map[tokens[query_file_idx]]
                ani = float(tokens[ani_idx])
                af_r = float(tokens[af_r_idx])
                af_q = float(tokens[af_q_idx])

                ani_af[qid][rid] = (ani, af_r, af_q)

        return ani_af

    def pairwise(self, gids, genome_files, preset, report_progress=True):
        """Calculate skani between all genome pairs in parallel.
        
        Since skani is symmetric each pair of genomes is only processed once. That is,
        genome A and genome B will be processed with skani as (A,B) or (B,A), but not
        both since they will give the same result.
        """

        gid_pairs = []
        for gid1, gid2 in combinations(gids, 2):
            gid_pairs.append((gid1, gid2))

        return self.pairs(gid_pairs, genome_files, preset, report_progress)

    def pairs(self, gid_pairs, genome_files, preset, report_progress=True):
        """Calculate skani between specified genome pairs in parallel.
        
        Genome pairs are specified as (query id, reference id).        
        """

        # determine pairs in cache and new pairs that must be processed
        items = []
        processed_pairs = 0
        ani_af = defaultdict(dict)
        for qid, rid in gid_pairs:
            ani, af_r, af_q = self.ani_af_cache.get(qid, {}).get(rid, (None, None, None))
            if ani is not None:
                ani_af[qid][rid] = (ani, af_r, af_q)
                processed_pairs += 1
            else:
                items.append((qid, rid, genome_files[qid], genome_files[rid], preset))

        if report_progress:
            self.log.info(f' - identified {processed_pairs:,} pairs in cache')

        # process pairs in parallel using a processing pool, making sure to
        # add results to the ANI cache and SQL database
        db_ani_af = []
        processed_pairs = 0
        start_time = time.time()
        for qid, rid, ani, af_r, af_q in self.pool.imap_unordered(skani_worker, items, chunksize=1):
            ani_af[qid][rid] = (ani, af_r, af_q)

            # since skani is symmetric add ANI / AF results to cache in both
            # directions to avoid unnecessary calculation in the future
            self.ani_af_cache[qid][rid] = (ani, af_r, af_q)
            self.ani_af_cache[rid][qid] = (ani, af_q, af_r)

            if self.db_conn:
                # place results into database in batches to keep the number of transactions sensible
                db_ani_af.append((qid, rid, ani, af_r, af_q))
                db_ani_af.append((rid, qid, ani, af_q, af_r))
                if len(db_ani_af) == self.DB_BATCH_SIZE:
                    self.db_cur.executemany('INSERT INTO ani_table (query_id, ref_id, ani, af_r, af_q) VALUES (?, ?, ?, ?, ?)',
                                            db_ani_af)
                    self.db_conn.commit()
                    db_ani_af = []

            if report_progress:
                processed_pairs += 1
                elapsed_time = time.time() - start_time
                remaining_time = (elapsed_time/processed_pairs) * \
                    float(len(items) - processed_pairs)

                status = '- processing {:,} of {:,} ({:.2f}%) genome pairs [elapsed: {}, remaining: {}; pairs/min: {:,}]'.format(
                    processed_pairs,
                    len(items),
                    float(processed_pairs*100)/len(items),
                    str(datetime.timedelta(seconds=round(elapsed_time))),
                    datetime.timedelta(seconds=round(remaining_time)),
                    int(processed_pairs/(elapsed_time/60.0)))
                sys.stdout.write('\r\033[K')  # clear line
                sys.stdout.write(f'{status}')
                sys.stdout.flush()

        if report_progress:
            sys.stdout.write('\n')

        if db_ani_af and self.db_conn:
            # place final set of results into database
            self.db_cur.executemany('INSERT INTO ani_table (query_id, ref_id, ani, af_r, af_q) VALUES (?, ?, ?, ?, ?)',
                                    db_ani_af)
            self.db_conn.commit()

        return ani_af

    def calculate(self, gid1, gid2, genome_file1, genome_file2):
        """Calculate symmetric ANI and AF between two genomes."""

        ani_af = skani(gid1, gid2, genome_file1, genome_file2)

        _qid, _rid, ani, af_r, af_q = ani_af
        self.ani_af_cache[gid1][gid2] = (ani, af_r, af_q)
        self.ani_af_cache[gid2][gid1] = (ani, af_q, af_r)

        db_ani_af = []
        db_ani_af.append((gid1, gid2, ani, af_r, af_q))
        db_ani_af.append((gid2, gid1, ani, af_q, af_r))
        self.db_cur.executemany('INSERT INTO ani_table (query_id, ref_id, ani, af_r, af_q) VALUES (?, ?, ?, ?, ?)',
                                db_ani_af)
        self.db_conn.commit()

        return ani, max(af_r, af_q)