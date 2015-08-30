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

from biolib.misc.time_keeper import TimeKeeper
from biolib.common import make_sure_path_exists
from biolib.external.hmmer import HmmModelParser

from genome_tree_tk.markers.infer_markers import InferMarkers
from genome_tree_tk.markers.paralog_test import ParalogTest
from genome_tree_tk.markers.lgt_test import LgtTest

import pickle


class MarkerWorkflow(object):
    """Determine phylogenetically informative marker genes."""

    def __init__(self, genome_quality_file, img_genome_dir, pfam_model_file, tigrfams_model_dir, cpus):
        """Initialization.

        Parameters
        ----------
        genome_quality_file : str
            File specifying completeness and contamination of genomes.
        img_genome_dir : str
            Directory with genomes in individual directories.
        pfam_model_file : str
            File containing Pfam HMMs.
        tigrfams_model_dir : str
            Directory containing TIGRFAMs HMMs.
        cpus : int
            Number of cpus to use.
        """

        self.logger = logging.getLogger()

        self.genome_quality_file = genome_quality_file
        self.img_genome_dir = img_genome_dir
        self.pfam_model_file = pfam_model_file
        self.tigrfams_model_dir = tigrfams_model_dir

        self.cpus = cpus

    def _get_hmms(self, hmms_dir, marker_genes, hmm_model_out, hmm_info_out):
        """Create HMM file with marker genes.

        Parameters
        ----------
        hmms_dir : str
            Directory containing individual HMMs.
        marker_genes : iterable
            Unique ids of marker genes.
        hmm_model_out : str
            File to containing phylogenetically informative HMMs.
        hmm_info_out : str
            File to contain information about HMMs.
        """

        # place all phylogenetically informative marker genes into a single model file
        fout = open(hmm_model_out, 'w')
        for marker_id in marker_genes:
            for line in open(os.path.join(hmms_dir, marker_id + '.hmm')):
                fout.write(line)
        fout.close()

        self.logger.info('    HMM models written to: ' + hmm_model_out)

        # read HMM model metadata
        hmm_model_parse = HmmModelParser(hmm_model_out)
        hmm_models = hmm_model_parse.models()
        fout = open(hmm_info_out, 'w')
        fout.write('Model Accession\tName\tDescription\tLength\n')
        for model in hmm_models.values():
            fout.write('%s\t%s\t%s\t%s\n' % (model.acc, model.name, model.desc, model.leng))
        fout.close()

        self.logger.info('    HMM information written to: ' + hmm_info_out)

    def run(self, ingroup_file,
            trusted_comp, trusted_cont,
            ubiquity, single_copy, redundancy,
            # non_conspecific,
            min_support,
            min_per_taxa,
            min_per_splits,
            perc_markers_to_jackknife,
            restrict_marker_list,
            taxonomy,
            output_dir):
        """Identify phylogenetic marker genes and save to HMM model file.

        Parameters
        ----------
        ingroup_file : str
            File specifying unique ids of ingroup genomes.
        trusted_comp : float
            Minimum completeness to trust genome for marker set inference.
        trusted_cont : float
            Maximum contamination to trust genome for marker set inference.
        ubiquity : float
            Threshold for defining ubiquity marker genes.
        single_copy : float
            Threshold for defining a single-copy marker gene.
        redundancy : float
            Threshold for declaring HMMs redundant
        non_conspecific : float
            Non-conspecific threshold for retaining multi-copy gene trees.
        min_support : float
            Minimum jackknife support of splits to use during LGT filtering [0, 1].
        min_per_taxa : float
            Minimum percentage of taxa required to consider a split during LGT filtering [0, 1].
        min_per_splits : float
            Minimum recovered well-supported splits required to retain marker during LGT filtering [0, 1].
        perc_markers_to_jackknife : float
            Percentage of taxa to keep during marker jackknifing [0, 1].
        restrict_marker_list : str
            Restrict marker set to genes listed within file.
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.
        output_dir : str
            Directory to store results.
        """

        time_keeper = TimeKeeper()

        alignment_dir = os.path.join(output_dir, 'alignments')
        make_sure_path_exists(alignment_dir)

        gene_tree_dir = os.path.join(output_dir, 'gene_trees')
        make_sure_path_exists(gene_tree_dir)

        model_dir = os.path.join(output_dir, 'hmm_models')
        make_sure_path_exists(model_dir)

        # identify markers in marker set
        valid_marker_genes = set()
        if restrict_marker_list != None:
            for line in open(restrict_marker_list):
                line_split = line.split()

                marker = line_split[0]
                if marker.startswith('pfam'):
                    marker = marker.replace('pfam', 'PF')

                valid_marker_genes.add(marker)

        # identify genes suitable for phylogenetic inference
        self.logger.info('')
        self.logger.info('  Identifying genes suitable for phylogenetic inference.')
        infer_markers = InferMarkers(self.genome_quality_file, self.img_genome_dir,
                                    self.pfam_model_file, self.tigrfams_model_dir, self.cpus)
        results = infer_markers.identify_marker_genes(ingroup_file,
                                                        trusted_comp, trusted_cont,
                                                        ubiquity, single_copy, redundancy,
                                                        valid_marker_genes,
                                                        alignment_dir, model_dir)
        num_ingroup_genomes, num_img_genome, num_user_genomes, trusted_genome_ids, marker_gene_stats, marker_genes = results
        self.logger.info('    Identifying %s single-copy, ubiquitous marker genes.' % len(marker_genes))

        # gather all single-copy HMMs into a single model file
        hmm_model_out = os.path.join(output_dir, 'single_copy_ubiquitous.hmm')
        hmm_info_out = os.path.join(output_dir, 'single_copy_ubiquitous.tsv')
        self.logger.info('')
        self.logger.info('  Gathering single-copy, ubiquitous HMMs.')
        self._get_hmms(model_dir, marker_genes, hmm_model_out, hmm_info_out)

        # infer gene trees
        self.logger.info('')
        self.logger.info('  Inferring gene trees.')
        infer_markers.infer_gene_trees(alignment_dir, gene_tree_dir, '.aln.masked.faa')

        # test gene trees for paralogs
        if False:
            self.logger.info('')
            self.logger.info('  Testing paralogs in gene trees.')
            paralog_test = ParalogTest()
            retained_paralog_markers = paralog_test.run(marker_genes,
                                                        hmm_model_out,
                                                        gene_tree_dir,
                                                        non_conspecific,
                                                        taxonomy,
                                                        '.tree',
                                                        output_dir)

            # gather all single-copy HMMs into a single model file
            paralog_model_out = os.path.join(output_dir, 'single_copy_ubiquitous.paralog_filtered.hmm')
            paralog_info_out = os.path.join(output_dir, 'paralog.paralog_filtered.tsv')
            self.logger.info('')
            self.logger.info('  Gathering single-copy, ubiquitous HMMs after paralog filtering.')
            self._get_hmms(model_dir, marker_genes, paralog_model_out, paralog_info_out)

        # identify gene trees which fail to reproduce the majority
        # of well-supported splits in a jackknifed genome tree
        # pickle.dump(trusted_genome_ids, open('trusted_genome_ids.txt', "wb"))
        # pickle.dump(retained_paralog_markers, open('retained_markers.txt', "wb"))

        # trusted_genome_ids = pickle.load(open('trusted_genome_ids.txt', 'rb'))
        # retained_paralog_markers = pickle.load(open('retained_markers.txt', 'rb'))

        lgt_test = LgtTest(self.cpus)
        results = lgt_test.run(trusted_genome_ids,
                               marker_genes,
                               hmm_model_out,
                               min_support,
                               min_per_taxa,
                               min_per_splits,
                               perc_markers_to_jackknife,
                               gene_tree_dir,
                               alignment_dir,
                               output_dir)
        retained_markers, distances, num_internal_nodes, well_supported_nodes, well_supported_internal_nodes = results

        # gather phylogenetically informative HMMs into a single model file
        phylo_hmm_model_out = os.path.join(output_dir, 'phylo.hmm')
        phylo_hmm_info_out = os.path.join(output_dir, 'phylo.tsv')
        self.logger.info('')
        self.logger.info('  Gathering phylogenetically informative HMMs.')
        self._get_hmms(model_dir, retained_markers, phylo_hmm_model_out, phylo_hmm_info_out)

        # write out metadata regarding potential marker genes
        marker_info_out = os.path.join(output_dir, 'marker_info.tsv')
        fout = open(marker_info_out, 'w')
        fout.write('Model accession\tName\tDescription\tLength\tUbiquity\tSingle copy\tRecovered splits (%)\tCompatible splits (%)\tNormalized compatible split length\tManhattan\tEuclidean\tRetained\tFiltering criteria\n')

        hmm_model_parse = HmmModelParser(hmm_model_out)
        hmm_models = hmm_model_parse.models()
        for model_acc, model_info in hmm_models.iteritems():
            filtering_criteria = ''
            if model_acc not in retained_markers:
                filtering_criteria = 'removed by LGT filter'

            fout.write('%s\t%s\t%s\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\t%.2f\t%.3f\t%s\t%s\n' % (model_acc,
                                                                               model_info.name,
                                                                               model_info.desc,
                                                                               model_info.leng,
                                                                               marker_gene_stats[model_acc][0],
                                                                               marker_gene_stats[model_acc][1],
                                                                               distances[model_acc][0],
                                                                               distances[model_acc][1],
                                                                               distances[model_acc][2],
                                                                               distances[model_acc][3],
                                                                               distances[model_acc][4],
                                                                               model_acc in retained_markers,
                                                                               filtering_criteria))
        fout.close()

        self.logger.info('')
        self.logger.info('  Information about marker genes written to: %s.' % marker_info_out)

        # generate summary report
        report_out = os.path.join(output_dir, 'marker_workflow.log')
        fout = open(report_out, 'w')
        fout.write('[markers]\n')
        fout.write('Number of ingroup genomes: %d\n' % num_ingroup_genomes)
        fout.write('  Number of IMG genomes: %d\n' % num_img_genome)
        fout.write('  Number of user genomes: %d\n' % num_user_genomes)
        fout.write('Number of trusted ingroup genomes: %d\n' % len(trusted_genome_ids))
        fout.write('')
        fout.write('Ubiquity threshold: %f\n' % ubiquity)
        fout.write('Single-copy threshold: %f\n' % single_copy)
        # fout.write('Non-conspecific threshold: %f\n' % non_conspecific)
        fout.write('')
        fout.write('LGT test minimum support: %f\n' % min_support)
        fout.write('LGT test minimum percent taxa for internal split: %f\n' % min_per_taxa)
        fout.write('LGT test minimum percent recovered splits to retain gene tree: %f\n' % min_per_splits)
        fout.write('LGT test percentage of markers to jackknife: %f\n' % perc_markers_to_jackknife)
        fout.write('')
        fout.write('Number of internal nodes: %d\n' % num_internal_nodes)
        fout.write('Number of well-supported nodes: %d\n' % well_supported_nodes)
        fout.write('Number of well-supported, internal nodes: %d\n' % well_supported_internal_nodes)
        fout.write('')
        fout.write('Initial number of marker genes: %d\n' % len(marker_genes))
        # fout.write('Filtered by paralog test: %d\n' % (len(marker_genes) - len(retained_paralog_markers)))
        fout.write('Filtered by LGT test: %d\n' % (len(marker_genes) - len(retained_markers)))
        fout.write('Final marker genes: %d\n' % len(retained_markers))
        fout.write(time_keeper.get_time_stamp())
        fout.close()

        return phylo_hmm_model_out
