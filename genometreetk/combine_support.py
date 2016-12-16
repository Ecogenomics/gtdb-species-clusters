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


class CombineSupport(object):
    """Combine all support values into a single tree."""

    def __init__(self):
        """Initialization."""
        
        self.logger = logging.getLogger()

    def _collect_support_values(self, tree):
        """Get support value for each node in preorder traversal.

        Parameters
        ----------
        tree : dendropy.Tree
          Tree to obtain support values from.
        """

        preorder_support = []
        for node in tree.preorder_node_iter():
            if node.is_internal():
                preorder_support.append(node.label)

        return preorder_support

    def run(self, support_type, bootstrap_tree, jk_marker_tree, jk_taxa_tree, output_tree):
        """Create new tree indicating combined support values.

        Tree can either be decorated with the average support value
        or the minimum support value as determine by support_type.

        Parameters
        ----------
        support_type : str => 'average' or 'minimum'
            Type of support value to calculate.
        bootstrap_tree : str
          Tree with bootstrap support values.
        jk_marker_tree : str
          Tree with jackknife marker support values.
        jk_taxa_tree : str
          Tree with jackknife taxa support values.
        output_file : str
          File to write tree with combined support values.
        """

        assert(support_type in ['average', 'minimum'])

        tree = dendropy.Tree.get_from_path(bootstrap_tree, schema='newick', rooting='force-rooted', preserve_underscores=True)
        bootstrap_support = self._collect_support_values(tree)

        tree = dendropy.Tree.get_from_path(jk_marker_tree, schema='newick', rooting='force-rooted', preserve_underscores=True)
        jk_marker_support = self._collect_support_values(tree)

        tree = dendropy.Tree.get_from_path(jk_taxa_tree, schema='newick', rooting='force-rooted', preserve_underscores=True)
        jk_taxa_support = self._collect_support_values(tree)

        internalNodeNum = 0
        for node in tree.preorder_node_iter():
            if node.is_internal():
                if support_type == 'average':
                    support = (int(bootstrap_support[internalNodeNum]) + int(jk_marker_support[internalNodeNum]) + int(jk_taxa_support[internalNodeNum])) / 3.0
                elif support_type == 'minimum':
                    support = min(int(bootstrap_support[internalNodeNum]), int(jk_marker_support[internalNodeNum]), int(jk_taxa_support[internalNodeNum]))

                node.label = '%s' % str(int(support + 0.5))
                internalNodeNum += 1

        tree.write_to_path(output_tree, schema='newick', suppress_rooting=True, unquoted_underscores=True)
