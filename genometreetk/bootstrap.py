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
import logging

import biolib.seq_io as seq_io
from biolib.external.fasttree import FastTree
from biolib.parallel import Parallel
from biolib.bootstrap import bootstrap_alignment, bootstrap_support
from biolib.common import remove_extension, make_sure_path_exists


class Bootstrap(object):
    """Assess robustness of genome tree by bootstrapping multiple sequence alignment."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
          Number of cpus to use.
        """

        self.logger = logging.getLogger()

        self.cpus = cpus

    def _producer(self, replicated_num):
        """Infer tree from bootstrapped multiple sequence alignment.

        Parameters
        ----------
        replicated_num : int
          Unique replicate number.
        """

        output_msa = os.path.join(self.replicate_dir, 'bootstrap_msa.r_' + str(replicated_num) + '.fna')
        bootstrap_alignment(self.msa, output_msa, frac=self.frac)

        fast_tree = FastTree(multithreaded=False)
        output_tree = os.path.join(self.replicate_dir, 'bootstrap_tree.r_' + str(replicated_num) + '.tree')
        fast_tree_output = os.path.join(self.replicate_dir, 'bootstrap_fasttree.r_' + str(replicated_num) + '.out')
        fast_tree.run(output_msa, self.base_type, self.model, output_tree, fast_tree_output)

        return True

    def _progress(self, processed_items, total_items):
        """Report progress of replicates."""

        return '    Processed %d of %d replicates.' % (processed_items, total_items)

    def run(self, input_tree, msa_file, num_replicates, model, base_type, frac, output_dir):
        """Bootstrap multiple sequence alignment.

        Parameters
        ----------
        input_tree : str
          Tree inferred with all data.
        msa_file : str
          File containing multiple sequence alignment for all taxa.
        num_replicates : int
          Number of replicates to perform.
        model : str
          Desired model of evolution.
        base_type : str
          Indicates if bases are nucleotides or amino acids.
        frac : float
          Fraction of alignment to subsample.
        output_dir : str
          Directory for bootstrap trees.
        """

        assert(model in ['wag', 'jtt'])
        assert(base_type in ['nt', 'prot'])

        self.model = model
        self.base_type = base_type
        self.frac = frac

        self.replicate_dir = os.path.join(output_dir, 'replicates')
        make_sure_path_exists(self.replicate_dir)

        # read full multiple sequence alignment
        self.msa = seq_io.read(msa_file)

        # calculate replicates
        self.logger.info('Calculating bootstrap replicates:')
        parallel = Parallel(self.cpus)
        parallel.run(self._producer, None, xrange(num_replicates), self._progress)

        # calculate support values
        rep_tree_files = []
        for rep_index in xrange(num_replicates):
            rep_tree_files.append(os.path.join(self.replicate_dir, 'bootstrap_tree.r_' + str(rep_index) + '.tree'))

        output_tree = os.path.join(output_dir, remove_extension(input_tree) + '.bootstrap.tree')
        bootstrap_support(input_tree, rep_tree_files, output_tree)

        return output_tree
