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

import dendropy

from biolib.newick import parse_label, create_label


class RerootTree(object):
    """Reroot tree."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()
        
    def _reroot(self, tree, outgroup_node, max_support=100):
        """Reroot tree taking proper care of bootstrap values."""

        # determine support values for each bipartition
        tree.encode_bipartitions()
        support_values = {}
        for nd in tree:
            support, taxon, aux_info = parse_label(nd.label)
            if nd.is_leaf():
                support_values[nd.bipartition] = max_support
            else:
                if support is not None:
                     support_values[nd.bipartition] = int(support)
                else:
                    support_values[nd.bipartition] = None
     
        # move support values for desired re-rooting
        new_root = outgroup_node.parent_node
        tree.reseed_at(new_root)
        tree.encode_bipartitions()
        for nd in tree:
            _, taxon, aux_info = parse_label(nd.label)
            nd.label = create_label(support_values.get(nd.bipartition, "not_specified"), taxon, aux_info)
        tree.seed_node.edge.length = None
     
        # do a hard re-rooting of the tree
        # (this invalidates the previous bipartitions, so must be done seperately)
        tree.is_rooted = True
        tree.reroot_at_edge(outgroup_node.edge,
                                    length1=0.5 * outgroup_node.edge_length,
                                    length2=0.5 * outgroup_node.edge_length)

        # determine bootstrap for new node
        for child in tree.seed_node.child_node_iter():
            if outgroup_node.is_leaf():
                if not child.is_leaf():
                    support, taxon, aux_info = parse_label(child.label)
                    child.label = create_label(max_support, taxon, aux_info)
            else:
                if child != outgroup_node:
                    support, _taxon, _aux_info = parse_label(outgroup_node.label)
                    _support, taxon, aux_info = parse_label(child.label)
                    child.label = create_label(support, taxon, aux_info)
 
        return tree

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

        outgroup_in_tree = set()
        tree = dendropy.Tree.get_from_path(input_tree, 
                                            schema='newick', 
                                            rooting='force-unrooted', 
                                            preserve_underscores=True)

        outgroup = set(outgroup)
        for n in tree.leaf_node_iter():
            if n.taxon.label in outgroup:
                outgroup_in_tree.add(n.taxon)

        self.logger.info('Identified %d outgroup taxa in the tree.' % len(outgroup_in_tree))

        if len(outgroup_in_tree) == 0:
            self.logger.warning('No outgroup taxa identified in the tree.')
            self.logger.warning('Tree was not rerooted.')
            sys.exit(0)

        mrca = tree.mrca(taxa=outgroup_in_tree)

        if len(mrca.leaf_nodes()) != len(outgroup_in_tree):
            self.logger.info('Outgroup is not monophyletic. Tree will be rerooted at the MRCA of the outgroup.')
            self.logger.info('The outgroup consisted of %d taxa, while the MRCA has %d leaf nodes.' % (len(outgroup_in_tree), len(mrca.leaf_nodes())))
            if len(mrca.leaf_nodes()) == len(tree.leaf_nodes()):
                self.logger.warning('The MRCA spans all taxa in the tree.')
                self.logger.warning('This indicating the selected outgroup is likely polyphyletic in the current tree.')
                self.logger.warning('Polyphyletic outgroups are not suitable for rooting. Try another outgroup.')
        else:
            self.logger.info('Outgroup is monophyletic.')

        if mrca.edge_length is None:
            self.logger.info('Tree appears to already be rooted on this outgroup.')
        else:
            self.logger.info('Rerooting tree.')
            self._reroot(tree, mrca)
            tree.write_to_path(output_tree, 
                                schema='newick', 
                                suppress_rooting=True, 
                                unquoted_underscores=True)
            self.logger.info('Rerooted tree written to: %s' % output_tree)

    def midpoint(self, input_tree, output_tree):
        """Reroot tree bat midpoint.

        Parameters
        ----------
        input_tree : str
          File containing Newick tree to rerooted.
        output_tree : str
          Name of file for rerooted tree.
        """

        tree = dendropy.Tree.get_from_path(input_tree, 
                                            schema='newick', 
                                            rooting='force-rooted', 
                                            preserve_underscores=True)
        tree.reroot_at_midpoint()
        tree.write_to_path(output_tree, 
                            schema='newick', 
                            suppress_rooting=True, 
                            unquoted_underscores=True)
