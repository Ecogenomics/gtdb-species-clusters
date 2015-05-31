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

import sys
import logging

from biolib.misc.time_keeper import TimeKeeper
from biolib.external.execute import check_dependencies
from biolib.common import check_file_exists, make_sure_path_exists

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

    def bootstrap(self, options):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [GenomeTreeTK - bootstrap] Bootstrap multiple sequence alignment.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.input_tree)
        check_file_exists(options.msa_file)
        make_sure_path_exists(options.output_dir)

        bootstrap = Bootstrap(options.cpus)
        output_tree = bootstrap.run(options.input_tree,
                                    options.msa_file,
                                    options.num_replicates,
                                    options.model,
                                    options.output_dir)

        self.logger.info('')
        self.logger.info('  Bootstrapped tree written to: %s' % output_tree)

        self.time_keeper.print_time_stamp()

    def jk_markers(self, options):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [GenomeTreeTK - jk_markers] Jackknife marker genes.')
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
                                                options.output_dir)

        self.logger.info('')
        self.logger.info('  Jackknifed marker tree written to: %s' % output_tree)

        self.time_keeper.print_time_stamp()

    def jk_taxa(self, options):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [GenomeTreeTK - jk_taxa] Jackknife ingroup taxa.')
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
                                            options.output_dir)

        self.logger.info('')
        self.logger.info('  Jackknifed taxa tree written to: %s' % output_tree)

        self.time_keeper.print_time_stamp()

    def combine(self, options):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [GenomeTreeTK - combine] Combine all support values into a single tree.')
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
        self.logger.info(' [GenomeTreeTK - midpoint] Rerooting tree at midpoint.')
        self.logger.info('*******************************************************************************')

        reroot = RerootTree()
        reroot.midpoint(options.input_tree, options.output_tree)

        self.time_keeper.print_time_stamp()

    def outgroup(self, options):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [GenomeTreeTK - outgroup] Rerooting tree with outgroup.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.outgroup_file)

        outgroup = set()
        for outgroup_id in open(options.outgroup_file):
            outgroup.add(outgroup_id)

        reroot = RerootTree()
        reroot.midpoint(options.input_tree, options.output_tree, outgroup)

        self.time_keeper.print_time_stamp()

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        logging.basicConfig(format='', level=logging.INFO)

        # check_dependencies(('FastTree'))

        if(options.subparser_name == 'bootstrap'):
            self.bootstrap(options)
        elif(options.subparser_name == 'jk_markers'):
            self.jk_markers(options)
        elif(options.subparser_name == 'jk_taxa'):
            self.jk_taxa(options)
        elif(options.subparser_name == 'combine'):
            self.combine(options)
        elif(options.subparser_name == 'midpoint'):
            self.midpoint(options)
        elif(options.subparser_name == 'outgroup'):
            self.outgroup(options)
        else:
            self.logger.error('  [Error] Unknown RefineM command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
