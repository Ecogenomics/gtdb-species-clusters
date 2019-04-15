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
import csv
import uuid
import time
import tempfile
import argparse
import operator
import ntpath
import logging
import multiprocessing as mp
from itertools import combinations, permutations
from collections import defaultdict

from biolib.external.execute import check_dependencies, run

class ANI_Cache(object):
    """Calculate average nucleotide identity between genomes using a precomputed cache where possible."""

    def __init__(self, ani_cache_file, cpus):
        """Initialization."""
        
        check_dependencies(['fastANI'])
        
        self.cpus = cpus

        self.logger = logging.getLogger('timestamp')
        
        self.ani_cache_file = ani_cache_file
        self._read_cache()
        
    def __del__(self):
        """Destructor."""
        
        self._write_cache()
        
    def _read_cache(self):
        """Read previously calculated ANI values."""
        
        self.ani_cache = defaultdict(lambda: {})

        if self.ani_cache_file:
            if os.path.exists(self.ani_cache_file):
                cache_size = 0
                for line in open(self.ani_cache_file):
                    line_split = line.strip().split('\t')
                    gid1 = line_split[0]
                    gid2 = line_split[1]
                    ani = float(line_split[2])
                    af = float(line_split[3])
                    self.ani_cache[gid1][gid2] = (ani, af)
                    cache_size += 1
                    
                self.logger.info('Read ANI cache with %d entries.' % cache_size)
            else:
                self.logger.warning('ANI cache file does not exist: %s' % self.ani_cache_file)
            
    def _write_cache(self):
        """Write cache to file."""
        
        if self.ani_cache_file:
            fout = open(self.ani_cache_file, 'w')
            cache_size = 0
            for gid1 in self.ani_cache:
                for gid2 in self.ani_cache[gid1]:
                    ani, af = self.ani_cache[gid1][gid2]
                    fout.write('%s\t%s\t%f\t%f\n' % (gid1, gid2, ani, af))
                    cache_size += 1
            fout.close()
            
            self.logger.info('Wrote ANI cache with %d entries.' % cache_size)

    def _get_genome_id(self, genome_path):
        """Extract genome ID from path to genomic file."""
        
        genome_id = ntpath.basename(genome_path)
        if genome_id.startswith('GCA_') or genome_id.startswith('GCF_'):
            genome_id = '_'.join(genome_id.split('_')[0:2])
            if genome_id.startswith('GCA_'):
                genome_id = 'GB_' + genome_id
            else:
                genome_id = 'RS_' + genome_id
        else:
            genome_id = '_'.join(genome_id.split('_')[0:2])
            
        return genome_id
        
    def fastani(self, qid, rid, genomic_files):
        """Calculate ANI between a pair of genomes using FastANI."""
        
        # check cache
        if qid in self.ani_cache:
            if rid in self.ani_cache[qid]:
                ani, af = self.ani_cache[qid][rid]
                ani_af = (qid, rid, ani, af)
                return ani_af
        
        # create file pointing to representative genome files
        tmp_fastani_file = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
        cmd = 'fastANI -q %s -r %s -o %s 2> /dev/null' % (
                    genomic_files[qid], 
                    genomic_files[rid], 
                    tmp_fastani_file)
            
        run(cmd)

        if os.path.exists(tmp_fastani_file) and os.stat(tmp_fastani_file).st_size > 0:
            for line in open(tmp_fastani_file):
                line_split = line.strip().split()

                query_genome = self._get_genome_id(line_split[0])
                ref_genome = self._get_genome_id(line_split[1])

                ani = float(line_split[2])
                af = float(line_split[3])/int(line_split[4])
                ani_af = (query_genome, ref_genome, ani, af)
        else:
            ani_af = (qid, rid, 0.0, 0.0)

        if os.path.exists(tmp_fastani_file):
            os.remove(tmp_fastani_file)

        return ani_af

    def __fastani_worker(self, genomic_files, queue_in, queue_out):
        """Process each data item in parallel."""

        while True:
            qid, rid = queue_in.get(block=True, timeout=None)
            if qid == None:
                break

            ani_af = self.fastani(qid, rid, genomic_files)

            queue_out.put(ani_af)

    def __fastani_writer(self, all_ani_af, num_pairs, report_progress, queue_writer):
        """Store or write results of worker threads in a single thread."""
        
        full_results = {}
        processed = 0
        while True:
            ani_af = queue_writer.get(block=True, timeout=None)
            if ani_af == None:
                for qid in full_results:
                    all_ani_af[qid] = full_results[qid]
                break

            qid, rid, ani, af = ani_af
            if qid not in full_results:
                full_results[qid] = {}
                
            full_results[qid][rid] = (ani, af)
            
            if report_progress:
                processed += 1
                statusStr = '-> Processing %d of %d (%.2f%%) genome pairs.'.ljust(86) % (
                                    processed, 
                                    num_pairs, 
                                    float(processed*100)/num_pairs)
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()
                
        if report_progress:
            sys.stdout.write('\n')
            
    def fastani_pairwise(self, gids, genome_files):
        """Calculate FastANI between all genome pairs in parallel."""
        
        if not gids:
            return {}
            
        if len(gids) == 2: # skip overhead of setting up queues and processes
            d = defaultdict(lambda: {})
            gids = list(gids)
            
            qid, rid, ani, af = self.fastani(gids[0], gids[1], genome_files)
            d[qid][rid] = (ani, af)
            self.ani_cache[qid][rid] = (ani, af)
            
            qid, rid, ani, af = self.fastani(gids[1], gids[0], genome_files)
            d[qid][rid] = (ani, af)
            self.ani_cache[qid][rid] = (ani, af)

            return d
        
        ani_af = mp.Manager().dict()
        
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()
        
        for gid1, gid2 in permutations(gids, 2):
            worker_queue.put((gid1, gid2))

        for _ in range(self.cpus):
            worker_queue.put((None, None))

        try:
            workerProc = [mp.Process(target = self.__fastani_worker, args = (genome_files,
                                                                                worker_queue, 
                                                                                writer_queue)) for _ in range(self.cpus)]
            writeProc = mp.Process(target = self.__fastani_writer, args = (ani_af, 0, False, writer_queue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writer_queue.put(None)
            writeProc.join()
        except:
            for p in workerProc:
                p.terminate()
            writeProc.terminate()
        
        ani_af = dict(ani_af)
        for qid in ani_af:
            for rid in ani_af[qid]:
                self.ani_cache[qid][rid] = ani_af[qid][rid]
        
        return ani_af
        
    def fastani_pairs(self, gid_pairs, genome_files, report_progress=True):
        """Calculate FastANI between specified genome pairs in parallel."""
        
        if not gid_pairs:
            return {}
            
        if len(gid_pairs) <= 6: # skip overhead of setting up queues and processes
            d = defaultdict(lambda: {})
            for idx in range(0, len(gid_pairs)):
                qid, rid, ani, af = self.fastani(gid_pairs[idx][0], gid_pairs[idx][1], genome_files)
                d[qid][rid] = (ani, af)
                self.ani_cache[qid][rid] = (ani, af)
            return d
        
        ani_af = mp.Manager().dict()
        
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()
        
        for gid1, gid2 in gid_pairs:
            worker_queue.put((gid1, gid2))

        for _ in range(self.cpus):
            worker_queue.put((None, None))

        try:
            workerProc = [mp.Process(target = self.__fastani_worker, args = (genome_files,
                                                                                worker_queue, 
                                                                                writer_queue)) for _ in range(self.cpus)]
            writeProc = mp.Process(target = self.__fastani_writer, args = (ani_af, len(gid_pairs), report_progress, writer_queue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writer_queue.put(None)
            writeProc.join()
        except:
            for p in workerProc:
                p.terminate()
            writeProc.terminate()
        
        ani_af = dict(ani_af)
        for qid in ani_af:
            for rid in ani_af[qid]:
                self.ani_cache[qid][rid] = ani_af[qid][rid]
        
        return ani_af
