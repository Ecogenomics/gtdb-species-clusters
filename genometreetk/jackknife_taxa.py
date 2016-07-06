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
import random
from math import floor

import biolib.seq_io as seq_io
from biolib.external.fasttree import FastTree
from biolib.parallel import Parallel
from biolib.bootstrap import bootstrap_support
from biolib.common import remove_extension, make_sure_path_exists

from genometreetk.tree_support import TreeSupport


class JackknifeTaxa(object):
    """Assess robustness by jackknifing taxa in alignment."""

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
        """Infer tree from jackknifed alignments.

        Parameters
        ----------
        replicated_num : int
          Unique replicate number.
        """

        output_msa = os.path.join(self.replicate_dir, 'jk_taxa.msa.' + str(replicated_num) + '.fna')
        self.jackknife_taxa(self.msa, self.perc_taxa_to_keep, self.outgroup_ids, output_msa)

        fast_tree = FastTree(multithreaded=False)
        output_tree = os.path.join(self.replicate_dir, 'jk_taxa.tree.' + str(replicated_num) + '.tre')
        fast_tree_output = os.path.join(self.replicate_dir, 'jk_taxa.fasttree.' + str(replicated_num) + '.out')
        fast_tree.run(output_msa, 'prot', self.model, output_tree, fast_tree_output)

        return True

    def _progress(self, processed_items, total_items):
        """Report progress of replicates."""

        return '    Processed %d of %d replicates.' % (processed_items, total_items)

    def jackknife_taxa(self, msa, perc_taxa_to_keep, outgroup_ids, output_file):
        """Jackknife alignment to a subset of taxa.

        Parameters
        ----------
        msa : d[seq_id] -> seq
          Full multiple sequence alignment.
        perc_taxa_to_keep : float
          Percentage of marker genes to keep in each replicate.
        outgroup_ids : set
          Labels of outgroup taxa.
        output_file : str
          File to write bootstrapped alignment.
        """

        # randomly select ingroup taxa
        ingroup_taxa = set(msa.keys()) - outgroup_ids
        taxa_to_keep = random.sample(ingroup_taxa, int(floor(len(ingroup_taxa) * perc_taxa_to_keep)))

        taxa_to_keep = set(taxa_to_keep).union(outgroup_ids)

        fout = open(output_file, 'w')
        for seq_id, seq in msa.iteritems():
            if seq_id in taxa_to_keep:
                fout.write('>' + seq_id + '\n')
                fout.write(seq + '\n')
        fout.close()

    def run(self, input_tree, msa_file, outgroup_file, perc_taxa_to_keep, num_replicates, model, output_dir):
        """Jackknife taxa.

        Parameters
        ----------
        input_tree : str
          Tree inferred with all data.
        msa_file : str
          File containing multiple sequence alignment for all taxa.
        outgroup_file : str
          File indicating labels of outgroup taxa.
        perc_taxa_to_keep : float
          Percentage of taxa to keep in each replicate.
        num_replicates : int
          Number of replicates to perform.
        model : str
          Desired model of evolution.
        output_dir : str
          input_tree directory for bootstrap trees.
        """

        assert(model in ['wag', 'jtt'])

        self.perc_taxa_to_keep = perc_taxa_to_keep
        self.model = model
        self.replicate_dir = os.path.join(output_dir, 'replicates')
        make_sure_path_exists(self.replicate_dir)
        # read outgroup taxa
        self.outgroup_ids = set()
        if outgroup_file:
            for line in open(outgroup_file):
                self.outgroup_ids.add(line.strip())

        # read full multiple sequence alignment
        self.msa = seq_io.read(msa_file)

        # calculate replicates
        self.logger.info('')
        self.logger.info('  Calculating jackknife taxa replicates:')
        #parallel = Parallel(self.cpus)
        #parallel.run(self._producer, None, xrange(num_replicates), self._progress)

        # calculate support
        rep_tree_files = []
        for rep_index in xrange(num_replicates):
            rep_tree_files.append(os.path.join(self.replicate_dir, 'jk_taxa.tree.' + str(rep_index) + '.tre'))

        tree_support = TreeSupport()
        output_tree = os.path.join(output_dir, remove_extension(input_tree) + '.jk_taxa.tree')
        tree_support.subset_taxa(input_tree, rep_tree_files, output_tree)

        return output_tree