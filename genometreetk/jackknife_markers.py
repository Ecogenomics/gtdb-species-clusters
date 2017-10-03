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
from biolib.common import remove_extension, make_sure_path_exists
from biolib.bootstrap import bootstrap_support


class JackknifeMarkers(object):
    """Assess robustness by jackkifing genes in alignment."""

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

        output_msa = os.path.join(self.replicate_dir, 'jk_markers.msa.' + str(replicated_num) + '.faa')
        self.jackknife_alignment(self.msa, self.perc_markers_to_keep, self.marker_lengths, output_msa)

        fast_tree = FastTree(multithreaded=False)
        output_tree = os.path.join(self.replicate_dir, 'jk_markers.tree.' + str(replicated_num) + '.tre')
        fast_tree_output = os.path.join(self.replicate_dir, 'jk_markers.fasttree.' + str(replicated_num) + '.out')
        fast_tree.run(output_msa, 'prot', self.model, output_tree, fast_tree_output)

        return True

    def _progress(self, processed_items, total_items):
        """Report progress of replicates."""

        return '==> Processed %d of %d replicates.' % (processed_items, total_items)

    def jackknife_alignment(self, msa, perc_markers_to_keep, marker_lengths, output_file):
        """Jackknife alignment to a subset of marker genes.

        The marker_lengths must be specified in the order
        in which genes were concatenated.

        Parameters
        ----------
        msa : d[seq_id] -> seq
          Full multiple sequence alignment.
        perc_markers_to_keep : float
          Percentage of marker genes to keep in each replicate [0, 1].
        marker_lengths : list
          Length of each marker gene.
        output_file : str
          File to write bootstrapped alignment.
        """
        markers_to_keep = random.sample(xrange(0, len(marker_lengths)), int(floor(perc_markers_to_keep * len(marker_lengths))))

        start_pos = [0]
        for index, ml in enumerate(marker_lengths):
            start_pos.append(start_pos[index] + ml)

        mask = [0] * sum(marker_lengths)
        for marker_index in markers_to_keep:
            start = start_pos[marker_index]
            end = start + marker_lengths[marker_index]
            mask[start:end] = [1] * (end - start)

        fout = open(output_file, 'w')
        for seq_id, seq in msa.iteritems():
            fout.write('>' + seq_id + '\n')
            sub_seq = ''.join([base for base, m in zip(seq, mask) if m == 1])
            fout.write(sub_seq + '\n')
        fout.close()

    def run(self, input_tree, 
                    msa_file, 
                    marker_info_file, 
                    mask_file, 
                    perc_markers_to_keep, 
                    num_replicates, 
                    model,
                    jk_dir,
                    output_dir):
        """Jackknife marker genes.

        Marker file should have the format:
          <marker id>\t<marker name>\t<marker desc>\t<length>\n

        Parameters
        ----------
        input_tree : str
          Tree inferred with all data.
        msa_file : str
          File containing multiple sequence alignment for all taxa.
        marker_info_file : str
          File indicating database id, HMM name, description and length of each marker in the alignment.
        mask_file : str
          File indicating masking of multiple sequence alignment.
        perc_markers_to_keep : float [0, 1]
          Percentage of marker genes to keep in each replicate.
        num_replicates : int
          Number of replicates to perform.
        model : str
          Desired model of evolution.
        output_dir : str
          Output directory for jackkife trees.
        """

        assert(model in ['wag', 'jtt'])

        self.model = model
        self.perc_markers_to_keep = perc_markers_to_keep
        
        
        # determine length of each marker gene in alignment
        rep_tree_files = []
        if not jk_dir:
            self.replicate_dir = os.path.join(output_dir, 'replicates')
            make_sure_path_exists(self.replicate_dir)
            
            marker_lengths = []
            total_len = 0
            with open(marker_info_file) as f:
                f.readline()
                for line in f:
                    line_split = line.split('\t')
                    ml = int(line_split[3])
                    marker_lengths.append(ml)
                    total_len += ml
                    
            self.logger.info('Concatenated length of markers: %d' % total_len)
                    
            # read mask
            mask = open(mask_file).readline().strip()
            start = 0
            self.marker_lengths = []
            total_mask_len = 0
            for ml in marker_lengths:
                end = start + ml
                zeros = mask[start:end].count('0')
                start = end
                
                self.marker_lengths.append(ml - zeros)
                total_mask_len += ml - zeros
                
            self.logger.info('Concatenated length of filtered MSA: %d' % total_mask_len)

            # read full multiple sequence alignment
            self.msa = seq_io.read(msa_file)
            
            if len(self.msa.values()[0]) != total_mask_len:
                self.logger.error('Length of MSA does not meet length of mask.')
                sys.exit()

            # calculate replicates
            self.logger.info('Calculating jackknife marker replicates:')
            parallel = Parallel(self.cpus)
            parallel.run(self._producer, None, xrange(num_replicates), self._progress)

            # calculate support
            self.logger.info('Calculating support for %d replicates.' % num_replicates)
            for rep_index in xrange(num_replicates):
                rep_tree_files.append(os.path.join(self.replicate_dir, 'jk_markers.tree.' + str(rep_index) + '.tre'))
        else:
            for f in os.listdir(jk_dir):
                if f.endswith('.tree') or f.endswith('.tre'):
                    rep_tree_files.append(os.path.join(jk_dir, f))
            self.logger.info('Read %d jackknife replicates.' % len(rep_tree_files))

        output_tree = os.path.join(output_dir, remove_extension(input_tree) + '.jk_markers.tree')
        bootstrap_support(input_tree, rep_tree_files, output_tree)

        return output_tree
