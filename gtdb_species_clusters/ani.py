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
from itertools import combinations
from collections import defaultdict

from biolib.external.execute import check_dependencies

class ANI(object):
    """Calculate average nucleotide identity between genomes."""

    def __init__(self, cpus):
        """Initialization."""
        
        check_dependencies(['fastANI'])
        
        self.cpus = cpus

        self.logger = logging.getLogger('timestamp')

    def fastani_pairwise(self, gids, genome_files):
        """Calculate FastANI between genomes and representatives in parallel."""
        
        ani_af = mp.Manager().dict()
        
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()
        
        for gid1, gid2 in combinations(gids, 2):
            worker_queue.put((gid1, [gid2]))

        for _ in range(self.cpus):
            worker_queue.put((None, None))

        try:
            workerProc = [mp.Process(target = self.__fastani_worker, args = (genome_files,
                                                                                worker_queue, 
                                                                                writer_queue)) for _ in range(self.cpus)]
            writeProc = mp.Process(target = self.__fastani_writer, args = (ani_af, writer_queue))

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
        
        return dict(ani_af)

    def fastani_reps(self, gids, rep_ids, genome_files):
        """Calculate FastANI between genomes and representatives in parallel."""
        
        ani_af = mp.Manager().dict()
        
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()
        
        rep_ids = list(rep_ids)
        for query_gid in gids:
            # process representatives in batches of 100 to keep
            # memory requirements in check
            for start_pos in range(0, len(rep_ids), 100):
                end_pos = min(start_pos + 100, len(rep_ids))
                worker_queue.put((query_gid, rep_ids[start_pos:end_pos]))

        for _ in range(self.cpus):
            worker_queue.put((None, None))

        try:
            workerProc = [mp.Process(target = self.__fastani_worker, args = (genome_files,
                                                                                worker_queue, 
                                                                                writer_queue)) for _ in range(self.cpus)]
            writeProc = mp.Process(target = self.__fastani_writer, args = (ani_af, writer_queue))

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
        
        return dict(ani_af)
        
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
        
    def _fastani(self, query_gid, rep_ids, genomic_files):
        """Calculate ANI between genomes and representatives genomes using FastANI."""
        
        # create file pointing to representative genome files
        tmp_fastani_file = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
        if len(rep_ids) > 1:
            tmp_rep_file = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
            fout = open(tmp_rep_file, 'w')
            for gid in rep_ids:
                fout.write('%s\n' % genomic_files[gid])
            fout.close()
            
            # ANI must be calculated between each genome and each representative
            cmd = 'fastANI --minFrag -1 -q %s --rl %s -o %s 2> /dev/null' % (
                        genomic_files[query_gid], 
                        tmp_rep_file, 
                        tmp_fastani_file)
                        
            os.remove(tmp_rep_file)
        else:
            cmd = 'fastANI --minFrag -1 -q %s -r %s -o %s 2> /dev/null' % (
                        genomic_files[query_gid], 
                        genomic_files[rep_ids[0]], 
                        tmp_fastani_file)
            
        os.system(cmd)

        ani_af = {}
        if os.path.exists(tmp_fastani_file) and os.stat(tmp_fastani_file).st_size > 0:
            for line in open(tmp_fastani_file):
                line_split = line.strip().split()

                query_genome = self._get_genome_id(line_split[0])
                ref_genome = self._get_genome_id(line_split[1])
                if query_genome not in ani_af:
                    ani_af[query_genome] = {}
                    
                ani = float(line_split[2])
                af = float(line_split[3])/int(line_split[4])
                ani_af[query_genome][ref_genome] = (ani, af)
        else:
            self.logger.warning('FastANI results file not created: %s' % query_gid)
            self.logger.warning('Faile command: %s' % cmd)
            return ani_af
        
        os.remove(tmp_fastani_file)

        return ani_af

    def __fastani_worker(self, genomic_files, queue_in, queue_out):
        """Process each data item in parallel."""

        while True:
            gid, rep_ids = queue_in.get(block=True, timeout=None)
            if gid == None:
                break

            ani_af = self._fastani(gid, rep_ids, genomic_files)

            queue_out.put(ani_af)

    def __fastani_writer(self, all_ani_af, queue_writer):
        """Store or write results of worker threads in a single thread."""
        
        full_results = {}
        while True:
            ani_af = queue_writer.get(block=True, timeout=None)
            if ani_af == None:
                for qid in full_results:
                    all_ani_af[qid] = full_results[qid]
                    
                break

            for qid in ani_af:
                if qid not in full_results:
                    full_results[qid] = ani_af[qid]
                else:
                    full_results[qid].update(ani_af[qid])
