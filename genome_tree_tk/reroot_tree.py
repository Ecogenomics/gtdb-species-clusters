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

    def _is_monophyletic(self, tree, group):
        """Determine if taxa are monophyletic.

        Parameters
        ----------
        tree : dendropy.Tree
          Tree representation.
        group : iterable
          Labels of taxa.

        Returns
        -------
        boolean
          True if taxa re monophyletic, else False.
        """

        tree.update_bipartitions()

        groupSplit = tree.taxon_set.get_taxa_bitmask(labels=group)
        return groupSplit in tree.split_edges

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

        tree = dendropy.Tree.get_from_path(input_tree, schema='newick', rooting='force-unrooted', preserve_underscores=True)

        if not self._isMonophyletic(tree, outgroup):
            self.logger.info('  [Warning] Outgroup is not monophyletic.')
            return
        else:
            self.logger.info('  Outgroup is monophyletic. Tree rerooted.')
            mrca = tree.mrca(taxon_labels=outgroup)
            tree.reroot_at_edge(mrca.edge, length1=0.5 * mrca.edge_length, length2=0.5 * mrca.edge_length, update_bipartitions=True)

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

        tree = dendropy.Tree.get_from_path(input_tree, schema='newick', rooting='force-unrooted', preserve_underscores=True)
        tree.reroot_at_midpoint(update_bipartitions=True)
        tree.write_to_path(output_tree, schema='newick', suppress_rooting=True, unquoted_underscores=True)
