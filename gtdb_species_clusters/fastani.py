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
import uuid
import tempfile
import ntpath
import logging
import re
import subprocess
import multiprocessing as mp
from itertools import combinations, permutations
from collections import defaultdict

from biolib.external.execute import check_dependencies, run

from gtdb_species_clusters.genome_utils import canonical_gid
from gtdb_species_clusters.type_genome_utils import symmetric_ani


class FastANI(object):
    """Calculate average nucleotide identity between genomes using a precomputed cache where possible."""

    def __init__(self, ani_cache_file, cpus):
        """Initialization."""
        
        check_dependencies(['fastANI'])
        
        self.cpus = cpus

        self.logger = logging.getLogger('timestamp')
        
        self.ani_cache_file = ani_cache_file
        self._read_cache()
        
        self.logger.info('Using FastANI v{}.'.format(self._get_version()))

    def __del__(self):
        """Destructor."""
        
        self.write_cache()
        
    def _get_version(self):
        """Returns the version of FastANI on the system path.
        
        Returns
        -------
        str
            The string containing the FastANI version.
        """
        try:
            proc = subprocess.Popen(['fastANI', '-v'], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, encoding='utf-8')
            stdout, stderr = proc.communicate()
            if stderr.startswith('Unknown option:'):
                return 'unknown (<1.3)'
                
            version = re.search(r'version (.+)', stderr)
            return version.group(1)
        except Exception as e:
            print(e)
            return 'unknown'
        
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
                    
                self.logger.info(f'Read ANI cache with {cache_size:,} entries.')
            else:
                self.logger.warning(f'ANI cache file does not exist: {self.ani_cache_file}')
            
    def write_cache(self):
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
            
            self.logger.info(f'Wrote ANI cache with {cache_size:,} entries.')

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
            
        return canonical_gid(genome_id)
        
    def fastani(self, qid, rid, q_gf, r_gf):
        """CalculateANI between a pair of genomes."""

        # check cache
        if qid in self.ani_cache:
            if rid in self.ani_cache[qid]:
                ani, af = self.ani_cache[qid][rid]
                ani_af = (qid, rid, ani, af)
                return ani_af
        
        # create file pointing to representative genome files
        tmp_fastani_file = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
        cmd = 'fastANI -q %s -r %s -o %s 2> /dev/null' % (
                    q_gf, 
                    r_gf, 
                    tmp_fastani_file)
            
        run(cmd)

        if os.path.exists(tmp_fastani_file) and os.stat(tmp_fastani_file).st_size > 0:
            for line in open(tmp_fastani_file):
                line_split = line.strip().split()
                ani = float(line_split[2])
                af = float(line_split[3])/int(line_split[4])
                ani_af = (qid, rid, ani, af)
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

            ani_af = self.fastani(qid, rid, 
                                    genomic_files[qid], 
                                    genomic_files[rid])

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
                statusStr = '-> Processing {:,} of {:,} ({:.2f}%) genome pairs.'.format(
                                    processed, 
                                    num_pairs, 
                                    float(processed*100)/num_pairs).ljust(86)
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()
                
        if report_progress:
            sys.stdout.write('\n')
            
    def pairwise(self, gids, genome_files, check_cache=False):
        """Calculate FastANI between all genome pairs in parallel."""
        
        if not gids:
            return {}
            
       # check if all pairs are in cache
        if check_cache:
            ani_af = defaultdict(lambda: {})
            
            in_cache = True
            for qid, rid in permutations(gids, 2):
                if qid in self.ani_cache:
                    if rid in self.ani_cache:
                        ani_af[qid][rid] = self.ani_cache[qid][rid]
                    else:
                        in_cache = False
                        break
                else:
                    in_cache = False
                    break
                    
            if in_cache:
                return ani_af
            
        # calculate required ANI pairs
        if len(gids) == 2: # skip overhead of setting up queues and processes
            d = {}
            gids = list(gids)
            d[gids[0]] = {}
            d[gids[1]] = {}
            
            qid, rid, ani, af = self.fastani(gids[0], gids[1],
                                                genome_files[gids[0]],
                                                genome_files[gids[1]])
            d[qid][rid] = (ani, af)
            self.ani_cache[qid][rid] = (ani, af)
            
            qid, rid, ani, af = self.fastani(gids[1], gids[0],
                                                genome_files[gids[1]],
                                                genome_files[gids[0]])
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
        
    def pairs(self, gid_pairs, genome_files, report_progress=True, check_cache=False):
        """Calculate FastANI between specified genome pairs in parallel."""
        
        if not gid_pairs:
            return {}
            
        # check if all pairs are in cache
        if check_cache:
            ani_af = defaultdict(lambda: {})
            
            in_cache = True
            for qid, rid in gid_pairs:
                if qid in self.ani_cache:
                    if rid in self.ani_cache:
                        ani_af[qid][rid] = self.ani_cache[qid][rid]
                    else:
                        in_cache = False
                        break
                else:
                    in_cache = False
                    break
                    
            if in_cache:
                return ani_af

        # calculate required ANI pairs
        if len(gid_pairs) <= 6: # skip overhead of setting up queues and processes
            d = defaultdict(lambda: {})
            for idx in range(0, len(gid_pairs)):
                qid = gid_pairs[idx][0]
                rid = gid_pairs[idx][1]
                qid, rid, ani, af = self.fastani(qid, rid, 
                                                    genome_files[qid],
                                                    genome_files[rid])
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

    def symmetric_ani_cached(self, gid1, gid2, genome_file1, genome_file2):
        """Calculate symmetric ANI and AF between two genomes."""
        
        ani_af12 = self.fastani(gid1, gid2, genome_file1, genome_file2)
        ani_af21 = self.fastani(gid2, gid1, genome_file2, genome_file1)

        self.ani_cache[gid1][gid2] = ani_af12[2:]
        self.ani_cache[gid2][gid1] = ani_af21[2:]
        
        return symmetric_ani(self.ani_cache, gid1, gid2)
