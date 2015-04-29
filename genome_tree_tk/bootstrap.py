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
import logging
import random
import multiprocessing as mp

import biolib.seq_io as seq_io
from biolib.external.fasttree import FastTree
from biolib.common import make_sure_path_exists
from biolib.parallel import Parallel

from genome_tree_tk.tree_support import TreeSupport

from genome_tree.markers.rerootTree import RerootTree


class Bootstrap(object):
    """Assess robustness of genome tree by bootstrapping multiple sequence alignment."""

    def __init__(self, cpus):
        """Initialization."""
        self.logger = logging.getLogger()
        
        self.cpus = cpus
    def _producer(self, replicated_num):
        """Infer tree from bootstrapped multiple sequence alignment."""

        output_msa = os.path.join(self.output_dir, 'bootstrap.msa.' + str(replicated_num) + '.fna')
        self._create_bootstrap_alignment(self.msa, output_msa)

        fast_tree = FastTree(multithreaded=False)
        output_tree = os.path.join(self.output_dir, 'bootstrap.tree.' + str(replicated_num) + '.tre')
        fast_tree_output = os.path.join(self.output_dir, 'bootstrap.fasttree.' + str(replicated_num) + '.out')
        fast_tree.run(output_msa, 'prot', 'wag', output_tree, fast_tree_output)

        return True
            
    def _progress(self, processed_items, total_items):
        """Report progress of replicates."""
        
        return 'Processed %d of %d replicates.' % (processed_items, total_items)

    def _create_bootstrap_alignment(self, msa, output_file):
        """Bootstrap multiple sequence alignment.
        
        Parameters
        ----------
        msa : d[seq_id] -> seq
          Full multiple sequence alignment.
        output_file : str
          File to write bootstrapped alignment.
        """
        alignment_len = len(msa[msa.keys()[0]])
        cols = [random.randint(0, alignment_len - 1) for _ in xrange(alignment_len)]

        fout = open(output_file, 'w')
        for seq_id, seq in msa.iteritems():
            fout.write('>' + seq_id + '\n')
            for col in cols:
                fout.write(seq[col])
            fout.write('\n')
        fout.close()
        
    def run(self, input_tree, msa_file, num_replicates, output_dir, output_tree):
        """Bootstrap multiple sequence alignment to assess robustness of genome tree.
        
        Parameters
        ----------
        input_tree : str
          Tree inferred with all data.
        msa_file : str
          File containing multiple sequence alignment for all taxa.
        num_replicates : int
          Number of replicates to perform.
        output_dir : str
          input_tree directory for bootstrap trees.
        output_tree : str
          Output tree with support values.
        """
        
        self.output_dir = output_dir

        # read full multiple sequence alignment
        self.msa = seq_io.read(msa_file)
        
        # calculate replicates
        parallel = Parallel(self.cpus)
        parallel.run(self._producer, None, xrange(num_replicates), self._progress)
                     
        # calculate support values
        rep_tree_files = []
        for rep_index in xrange(num_replicates):
            rep_tree_files.append(os.path.join(output_dir, 'bootstrap.tree.' + str(rep_index) + '.tre'))

        tree_support = TreeSupport()
        tree_support.common_taxa(input_tree, rep_tree_files, output_tree)
