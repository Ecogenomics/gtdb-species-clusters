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

from math import floor

import dendropy


class TreeSupport():
    """Calculate support values for clades."""

    def __init__(self):
        """Initialize."""
        self.logger = logging.getLogger()

    def common_taxa(self, input_tree, replicate_trees, output_tree):
        """ Calculate support for tree with replicates covering the same taxon set.

        Parameters
        ----------
        input_tree : str
          Tree inferred from complete data.
        replicate_trees : iterable
          Files containing replicate trees.
        output_tree: str
          Name of output tree with support values.
        """

        tree = dendropy.Tree.get_from_path(input_tree, schema='newick', rooting='force-unrooted', preserve_underscores=True)
        tree.encode_bipartitions()

        print 'read input tree'

        rep_trees = dendropy.TreeList(taxon_namespace=tree.taxon_namespace)
        for i, rep_tree_file in enumerate(replicate_trees):
            print 'tree', i  # ***
            rep_tree = dendropy.Tree.get_from_path(rep_tree_file,
                                                         schema='newick',
                                                         rooting='force-unrooted',
                                                         preserve_underscores=True,
                                                         taxon_namespace=tree.taxon_namespace)
            rep_tree.bipartitions = True  # this is a hack because the call to frequency_of_bipartition expects this variable to be defined
            rep_tree.encode_bipartitions()
            rep_trees.append(rep_tree)

        for i, node in enumerate(tree.internal_nodes()):
            bootstrap = int(rep_trees.frequency_of_bipartition(bipartition=node.bipartition, is_bipartitions_updated=True) * 100)

            if node.label:
                node.label = str(bootstrap) + ':' + node.label
            else:
                node.label = str(bootstrap)

        tree.write_to_path(output_tree, schema='newick', suppress_rooting=True, unquoted_underscores=True)

    def subset_taxa(self, input_tree, replicate_trees, output_tree):
        """Calculate support for tree with replicates containing a subset of taxa.

        Parameters
        ----------
        input_tree : str
          Tree inferred from complete data.
        replicate_trees : iterable
          Files containing replicate trees.
        output_tree: str
          Name of output tree with support values.
        """

        tree = dendropy.Tree.get_from_path(input_tree, schema='newick', rooting='force-unrooted', preserve_underscores=True)
        for node in tree.internal_nodes():
            node.label = 0
            node.nontrivial_splits = 0

        for rep_tree_file in replicate_trees:
            rep_tree = dendropy.Tree.get_from_path(rep_tree_file, schema='newick', rooting='force-unrooted', preserve_underscores=True)

            rep_tree_list = dendropy.TreeList([rep_tree])

            rep_tree_taxa_set = set([x.taxon.label for x in rep_tree.leaf_nodes()])

            for node in tree.internal_nodes():
                taxa_labels = set([x.taxon.label for x in node.leaf_nodes()]).intersection(rep_tree_taxa_set)

                split = rep_tree.taxon_namespace.taxa_bitmask(labels=taxa_labels)
                normalized_split = dendropy.Bipartition.normalize_bitmask(
                                    bitmask=split,
                                    fill_bitmask=rep_tree.taxon_namespace.all_taxa_bitmask(),
                                    lowest_relevant_bit=1)

                if len(taxa_labels) > 1:
                    # tabulate results for non-trivial splits
                    node.label += int(rep_tree_list.frequency_of_bipartition(split_bitmask=normalized_split))
                    node.nontrivial_splits += 1

        for node in tree.internal_nodes():
            if node.nontrivial_splits > 0:
                node.label = str(int(floor(node.label * 100.0 / node.nontrivial_splits)))
            else:
                node.label = 'NA'

        tree.write_to_path(output_tree, schema='newick', suppress_rooting=True, unquoted_underscores=True)
