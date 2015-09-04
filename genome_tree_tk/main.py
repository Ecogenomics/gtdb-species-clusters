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

from biolib.misc.time_keeper import TimeKeeper
from biolib.common import check_file_exists, make_sure_path_exists
from biolib.taxonomy import Taxonomy
from biolib.external.execute import check_dependencies

from genome_tree_tk.trusted_genome_workflow import TrustedGenomeWorkflow
from genome_tree_tk.dereplication_workflow import DereplicationWorkflow
from genome_tree_tk.marker_workflow import MarkerWorkflow
from genome_tree_tk.infer_workflow import InferWorkflow
from genome_tree_tk.bootstrap import Bootstrap
from genome_tree_tk.jackknife_markers import JackknifeMarkers
from genome_tree_tk.jackknife_taxa import JackknifeTaxa
from genome_tree_tk.combine_support import CombineSupport

from genome_tree_tk.reroot_tree import RerootTree


class OptionsParser():
    def __init__(self):
        """Initialization"""
        self.logger = logging.getLogger()
        self.time_keeper = TimeKeeper()

    def _read_config_file(self):
        """Read configuration info."""

        cfg_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'genome_tree_tk.cfg')

        d = {}
        for line in open(cfg_file):
            key, value = line.split('=')
            d[key] = value.strip()

        return d

    def trusted(self, options):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [GenomeTreeTk - trusted] Determining trusted genomes.')
        self.logger.info('*******************************************************************************')

        config_data = self._read_config_file()

        trusted_genome_workflow = TrustedGenomeWorkflow(config_data['assembly_metadata_file'],
                                                        config_data['checkm_stats_file'],
                                                        config_data['taxonomy_file'])

        trusted_genome_workflow.run(options.trusted_comp,
                                        options.trusted_cont,
                                        options.max_contigs,
                                        options.min_N50,
                                        options.allow_partial_taxonomy,
                                        options.trusted_genomes_file)

        self.logger.info('')
        self.logger.info('  Trusted genome list written to: %s' % options.trusted_genomes_file)

        self.time_keeper.print_time_stamp()

    def dereplicate(self, options):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [GenomeTreeTk - dereplicate] Dereplicate genomes based on taxonomy.')
        self.logger.info('*******************************************************************************')

        if options.trusted_genomes_file:
            check_file_exists(options.trusted_genomes_file)

        config_data = self._read_config_file()
        dereplication_workflow = DereplicationWorkflow(config_data['assembly_metadata_file'],
                                                        config_data['taxonomy_file'],
                                                        config_data['type_strain_file'])

        dereplication_workflow.run(options.max_species,
                                   options.trusted_genomes_file,
                                   options.derep_genome_file)

        self.logger.info('')
        self.logger.info('  Dereplicated genome list written to: %s' % options.derep_genome_file)

        self.time_keeper.print_time_stamp()

    def markers(self, options):
        self.logger.info('*******************************************************************************')
        self.logger.info(' [GenomeTreeTk - markers] Determining marker genes.')
        self.logger.info('*******************************************************************************')

        make_sure_path_exists(options.output_dir)

        config_data = self._read_config_file()
        taxonomy = Taxonomy().read(config_data['taxonomy_file'])

        marker_workflow = MarkerWorkflow(config_data['genome_dir_file'],
                                            config_data['pfam_model_file'],
                                            config_data['tigrfams_model_dir'],
                                            options.cpus)

        phylo_hmm_out = marker_workflow.run(options.ingroup_file,
                                                options.ubiquity,
                                                options.single_copy,
                                                options.redundancy,
                                                options.min_support,
                                                options.min_per_taxa,
                                                options.perc_markers,
                                                options.restict_marker_list,
                                                taxonomy,
                                                options.output_dir)

        self.logger.info('')
        self.logger.info('  Marker genes written to: %s' % phylo_hmm_out)

        self.time_keeper.print_time_stamp()

    def infer(self, options):
        self.logger.info('*******************************************************************************')
        self.logger.info(' [GenomeTreeTk - infer] Inferring genome tree.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.genome_id_file)
        check_file_exists(options.marker_id_file)
        make_sure_path_exists(options.output_dir)

        config_data = self._read_config_file()
        infer_workflow = InferWorkflow(config_data['genome_dir_file'],
                                            config_data['pfam_model_file'],
                                            config_data['tigrfams_model_dir'],
                                            options.cpus)

        infer_workflow.run(options.genome_id_file,
                                options.marker_id_file,
                                options.model,
                                options.output_dir)

        self.logger.info('')
        self.logger.info('  Results written to: %s' % options.output_dir)

        self.time_keeper.print_time_stamp()

    def bootstrap(self, options):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [GenomeTreeTk - bootstrap] Bootstrap multiple sequence alignment.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.input_tree)
        check_file_exists(options.msa_file)
        make_sure_path_exists(options.output_dir)

        bootstrap = Bootstrap(options.cpus)
        output_tree = bootstrap.run(options.input_tree,
                                    options.msa_file,
                                    options.num_replicates,
                                    options.model,
                                    options.gamma,
                                    options.output_dir)

        self.logger.info('')
        self.logger.info('  Bootstrapped tree written to: %s' % output_tree)

        self.time_keeper.print_time_stamp()

    def jk_markers(self, options):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [GenomeTreeTk - jk_markers] Jackknife marker genes.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.input_tree)
        check_file_exists(options.msa_file)
        make_sure_path_exists(options.output_dir)

        jackknife_markers = JackknifeMarkers(options.cpus)
        output_tree = jackknife_markers.run(options.input_tree,
                                                options.msa_file,
                                                options.gene_length_file,
                                                options.perc_markers,
                                                options.num_replicates,
                                                options.model,
                                                options.gamma,
                                                options.output_dir)

        self.logger.info('')
        self.logger.info('  Jackknifed marker tree written to: %s' % output_tree)

        self.time_keeper.print_time_stamp()

    def jk_taxa(self, options):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [GenomeTreeTk - jk_taxa] Jackknife ingroup taxa.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.input_tree)
        check_file_exists(options.msa_file)
        make_sure_path_exists(options.output_dir)

        jackknife_taxa = JackknifeTaxa(options.cpus)
        output_tree = jackknife_taxa.run(options.input_tree,
                                            options.msa_file,
                                            options.outgroup_ids,
                                            options.perc_taxa,
                                            options.num_replicates,
                                            options.model,
                                            options.gamma,
                                            options.output_dir)

        self.logger.info('')
        self.logger.info('  Jackknifed taxa tree written to: %s' % output_tree)

        self.time_keeper.print_time_stamp()

    def combine(self, options):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [GenomeTreeTk - combine] Combine all support values into a single tree.')
        self.logger.info('*******************************************************************************')

        combineSupport = CombineSupport()
        combineSupport.run(options.support_type,
                            options.bootstrap_tree,
                            options.jk_marker_tree,
                            options.jk_taxa_tree,
                            options.output_tree)

        self.time_keeper.print_time_stamp()

    def support_wf(self, options):
        """"Perform entire tree support workflow."""

        self.bootstrap(options)
        self.jk_markers(options)
        self.jk_taxa(options)
        self.combine(options)

    def midpoint(self, options):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [GenomeTreeTk - midpoint] Rerooting tree at midpoint.')
        self.logger.info('*******************************************************************************')

        reroot = RerootTree()
        reroot.midpoint(options.input_tree, options.output_tree)

        self.time_keeper.print_time_stamp()

    def outgroup(self, options):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [GenomeTreeTk - outgroup] Rerooting tree with outgroup.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.outgroup_file)

        outgroup = set()
        for outgroup_id in open(options.outgroup_file):
            outgroup_id = outgroup_id.split('\t')[0].strip()
            outgroup.add(outgroup_id)

        reroot = RerootTree()
        reroot.root_with_outgroup(options.input_tree, options.output_tree, outgroup)

        self.time_keeper.print_time_stamp()

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        logging.basicConfig(format='', level=logging.INFO)

        check_dependencies(('FastTree', 'hmmsearch'))

        if options.subparser_name == 'trusted':
            self.trusted(options)
        elif options.subparser_name == 'dereplicate':
            self.dereplicate(options)
        elif options.subparser_name == 'markers':
            self.markers(options)
        elif options.subparser_name == 'infer':
            self.infer(options)
        elif options.subparser_name == 'bootstrap':
            self.bootstrap(options)
        elif options.subparser_name == 'jk_markers':
            self.jk_markers(options)
        elif options.subparser_name == 'jk_taxa':
            self.jk_taxa(options)
        elif options.subparser_name == 'combine':
            self.combine(options)
        elif options.subparser_name == 'midpoint':
            self.midpoint(options)
        elif options.subparser_name == 'outgroup':
            self.outgroup(options)
        else:
            self.logger.error('  [Error] Unknown RefineM command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
