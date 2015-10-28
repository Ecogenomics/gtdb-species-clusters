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

import logging

import dendropy


class RerootTree(object):
    """Reroot tree."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

    def root_with_outgroup(self, input_tree, output_tree, outgroup):
        """Reroot the tree using the given outgroup.

        Parameters
        ----------
        input_tree : str
          File containing Newick tree to rerooted.
        output_tree : str
          Name of file for rerooted tree.
        outgroup : iterable
          Labels of taxa in outgroup.
        """

        tree = dendropy.Tree.get_from_path(input_tree, schema='newick', rooting='force-rooted', preserve_underscores=True)
        for taxa in outgroup:
            node = tree.find_node_with_taxon_label(taxa)
            if not node:
                self.logger.error('  [Error] Missing taxa %d. Tree not rerooted.' % taxa)
                return

        mrca = tree.mrca(taxon_labels=outgroup)

        self.logger.info('')
        if len(mrca.leaf_nodes()) != len(outgroup):
            self.logger.warning('  [Warning] Outgroup is not monophyletic. Tree was rerooted at the MCRA of the outgroup.')
            self.logger.warning('  [Warning] The outgroup consisted of %d taxa, while the MCRA has %d leaf nodes.' % (len(outgroup), len(mrca.leaf_nodes())))
        else:
            self.logger.info('  Outgroup is monophyletic. Tree rerooted.')

        tree.reroot_at_edge(mrca.edge, length1=0.5 * mrca.edge_length, length2=0.5 * mrca.edge_length)
        tree.write_to_path(output_tree, schema='newick', suppress_rooting=True, unquoted_underscores=True)

    def midpoint(self, input_tree, output_tree):
        """Reroot tree bat midpoint.

        Parameters
        ----------
        input_tree : str
          File containing Newick tree to rerooted.
        output_tree : str
          Name of file for rerooted tree.
        """

        tree = dendropy.Tree.get_from_path(input_tree, schema='newick', rooting='force-rooted', preserve_underscores=True)
        tree.reroot_at_midpoint()
        tree.write_to_path(output_tree, schema='newick', suppress_rooting=True, unquoted_underscores=True)
