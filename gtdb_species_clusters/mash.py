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

from biolib.external.execute import check_dependencies, run

from gtdb_species_clusters.genome_utils import read_genome_path, canonical_gid

class Mash(object):
    """Calculate Mash distance between genomes."""

    def __init__(self, cpus):
        """Initialization."""
        
        check_dependencies(['mash'])
        
        self.cpus = cpus

        self.logger = logging.getLogger('timestamp')
        
    def _mash_genome_id(self, mash_genome_id):
        """Extract canonical GTDB genome ID from Mash results."""
            
        # get filename and remove information past genome accession
        # (e.g., GCA_002498385.1_ASM249838v1_genomic.fna => GCA_002498385.1)
        mash_genome_id = ntpath.basename(mash_genome_id)
        mash_genome_id = mash_genome_id[0:mash_genome_id.find('_', 4)]
        gid = canonical_gid(mash_genome_id)
            
        return gid

    def sketch(self, gids, genome_files, genome_list_file, sketch_file):
        """Create Mash sketch for genomes."""
        
        # create Mash sketch for potential representative genomes
        if not os.path.exists(sketch_file):
            fout = open(genome_list_file, 'w')
            for gid in gids:
                fout.write(genome_files[gid] + '\n')
            fout.close()

            self.logger.info(f'Creating Mash sketch for {len(gids):,} genomes.')
            cmd = 'mash sketch -l -p %d -k 16 -s 5000 -o %s %s 2> /dev/null' % (self.cpus, 
                                                                                sketch_file, 
                                                                                genome_list_file)
            run(cmd)
        else:
            self.logger.warning('Using previously generated sketch file.')
            
    def dist_pairwise(self, min_dist, sketch_file, dist_file):
        """Calculate pairwise Mash distance between genomes."""

        if not os.path.exists(dist_file):
            self.logger.info('Calculating pairwise Mash distances between genomes (d = %.2f).' % min_dist)
            cmd = 'mash dist -p %d -d %f -v %f %s %s > %s 2> /dev/null' % (self.cpus,
                                                                            min_dist,
                                                                            1e-5,
                                                                            sketch_file, 
                                                                            sketch_file, 
                                                                            dist_file)
            run(cmd)
        else:
            self.logger.warning('Using previously generated pairwise distance file.')
            
    def dist(self, min_dist, ref_sketch_file, query_sketch_file, dist_file):
        """Calculate Mash distance between reference and query genomes."""

        if not os.path.exists(dist_file):
            self.logger.info('Calculating Mash distances between reference and query genomes (d = %.2f).' % min_dist)
            cmd = 'mash dist -p %d -d %f -v %f %s %s > %s 2> /dev/null' % (self.cpus,
                                                                            min_dist,
                                                                            1e-5,
                                                                            ref_sketch_file, 
                                                                            query_sketch_file, 
                                                                            dist_file)
            run(cmd)
        else:
            self.logger.warning('Using previously generated pairwise distance file.')
            
    def read_ani(self, dist_file):
        """Read ANI estimates."""

        mash_ani = defaultdict(lambda: {})
        for line in open(dist_file):
            line_split = line.strip().split('\t')
            
            rid = self._mash_genome_id(line_split[0])
            qid = self._mash_genome_id(line_split[1])
            mash_ani[qid][rid] = 100 - 100*float(line_split[2])

        return mash_ani
