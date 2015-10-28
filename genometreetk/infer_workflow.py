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

from biolib.common import make_sure_path_exists
from biolib.misc.time_keeper import TimeKeeper
from biolib.external.hmmer import HmmModelParser
from biolib.external.fasttree import FastTree

from genometreetk.common import (read_genome_id_file,
                                    read_genome_dir_file,
                                    read_marker_id_file,
                                    create_concatenated_alignment)
from genometreetk.markers.align_markers import AlignMarkers


class InferWorkflow(object):
    """Infer phylogenetic tree from concatenated marker genes."""

    def __init__(self, genome_dir_file, pfam_model_file, tigrfams_model_dir, cpus):
        """Initialization.

        Parameters
        ----------
        genome_dir_file : str
            File specifying directory for each genome.
        pfam_model_file : str
            File containing Pfam HMMs.
        tigrfams_model_dir : str
            Directory containing TIGRFAMs HMMs.
        cpus : int
            Number of cpus to use.
        """

        self.logger = logging.getLogger()

        self.genome_dir_file = genome_dir_file
        self.pfam_model_file = pfam_model_file
        self.tigrfams_model_dir = tigrfams_model_dir

        self.cpus = cpus

    def _fetch_marker_models(self, marker_genes, hmm_model_out, hmm_info_out, output_model_dir):
        """Save PFAM and TIGRFAM marker genes into a single HMM model file.

        Parameters
        ----------
        marker_genes : iterable
            Marker genes to fetch.
        hmm_model_out : str
            File to containing phylogenetically informative HMMs.
        hmm_info_out : str
            File to contain information about HMMs.
        output_model_dir : str
            Directory to write individual HMM model files.
        """

        marker_id_to_name = {}
        for line in open(self.pfam_model_file):
            if 'NAME' in line:
                name = line.split()[1].rstrip()
            elif 'ACC' in line:
                acc = line.split()[1].rstrip()
                marker_id_to_name[acc] = name

        fout_model = open(hmm_model_out, 'w')
        for marker_id in marker_genes:
            output_model_file = os.path.join(output_model_dir, marker_id + '.hmm')
            if 'PF' in marker_id:
                os.system('hmmfetch ' + self.pfam_model_file + ' ' + marker_id_to_name[marker_id] + ' > ' + output_model_file)
            else:
                input_model_file = os.path.join(self.tigrfams_model_dir, marker_id + '.HMM')
                os.system('hmmfetch ' + input_model_file + ' ' + marker_id + ' > ' + output_model_file)

            # write model to file
            for line in open(output_model_file):
                fout_model.write(line)

        fout_model.close()

        self.logger.info('    HMM models written to: ' + hmm_model_out)

        # read HMM model metadata
        hmm_model_parse = HmmModelParser(hmm_model_out)
        hmm_models = hmm_model_parse.models()

        fout_info = open(hmm_info_out, 'w')
        fout_info.write('Model Accession\tName\tDescription\tLength\n')
        for model in hmm_models.values():
            fout_info.write('%s\t%s\t%s\t%s\n' % (model.acc, model.name, model.desc, model.leng))
        fout_info.close()

        self.logger.info('    HMM information written to: ' + hmm_info_out)

    def run(self, genome_id_file,
                    marker_id_file,
                    model,
                    output_dir):
        """Identify phylogenetic tree.

        Parameters
        ----------
        genome_id_file : str
            File specifying unique ids of genomes to include in tree.
        marker_id_file : str
            File specifying unique ids of marker genes  to use for inference.
        model : str ['wag' or 'jtt']
            Model of evolution to use.
        output_dir : str
            Directory to store results.
        """

        time_keeper = TimeKeeper()

        output_alignment_dir = os.path.join(output_dir, 'alignments')
        make_sure_path_exists(output_alignment_dir)

        output_model_dir = os.path.join(output_dir, 'hmm_models')
        make_sure_path_exists(output_model_dir)

        # read directory for each genome
        genome_dirs = read_genome_dir_file(self.genome_dir_file)

        # read genomes within the ingroup
        ncbi_genome_ids, user_genome_ids = read_genome_id_file(genome_id_file)
        genome_ids = ncbi_genome_ids.union(user_genome_ids)
        self.logger.info('Inferring tree for %d genomes.' % len(genome_ids))
        self.logger.info('NCBI genomes: %d' % len(ncbi_genome_ids))
        self.logger.info('User genomes: %d' % len(user_genome_ids))

        # get marker genes
        self.logger.info('Reading marker genes.')
        marker_genes = read_marker_id_file(marker_id_file)
        self.logger.info('Read %d marker genes.' % len(marker_genes))

        # gather all single-copy HMMs into a single model file
        hmm_model_out = os.path.join(output_dir, 'phylo.hmm')
        hmm_info_out = os.path.join(output_dir, 'phylo.tsv')
        self.logger.info('Generating marker gene HMM model files.')
        self._fetch_marker_models(marker_genes, hmm_model_out, hmm_info_out, output_model_dir)

        # align gene sequences
        align_markers = AlignMarkers(self.cpus)
        align_markers.run(genome_ids, genome_dirs, marker_genes, True, output_alignment_dir, output_model_dir)

        # create concatenated alignment file
        self.logger.info('Concatenating alignments.')
        concatenated_alignment_file = os.path.join(output_dir, 'concatenated_alignment.faa')
        marker_file = os.path.join(output_dir, 'concatenated_markers.tsv')
        create_concatenated_alignment(genome_ids, marker_genes, output_alignment_dir, concatenated_alignment_file, marker_file)

        # create concatenated genome tree
        self.logger.info('Inferring concatenated genome tree.')
        concatenated_tree = os.path.join(output_dir, 'concatenated.tree')
        concatenated_tree_log = os.path.join(output_dir, 'concatenated.tree.log')
        log_file = os.path.join(output_dir, 'fasttree.log')
        fast_tree = FastTree(multithreaded=True)
        fast_tree.run(concatenated_alignment_file, 'prot', model, concatenated_tree, concatenated_tree_log, log_file)

        # generate summary report
        report_out = os.path.join(output_dir, 'infer_workflow.log')
        fout = open(report_out, 'w')
        fout.write('[infer]\n')
        fout.write('Genome Id file: %s\n' % genome_id_file)
        fout.write('Marker Id file: %s\n' % marker_id_file)
        fout.write('Model of evolution: %s\n' % model)
        fout.write(time_keeper.get_time_stamp())
        fout.close()
