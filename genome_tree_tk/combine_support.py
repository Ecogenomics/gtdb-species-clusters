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

from genome_tree.timeKeeper import TimeKeeper
from genome_tree.defaultValues import DefaultValues
from genome_tree.fasttree import FastTree
from genome_tree.common import readTreeModel, makeSurePathExists, checkDirExists

import dendropy


class CombineSupport(object):
    """Combine all support values into a single tree."""

    def __init__(self, outputPrefix, outputDir, supportType):
        """Setup directories for combined command."""
        self.logger = logging.getLogger()

        if not outputPrefix[-1] == '.' and not outputPrefix[-1] == '_':
            outputPrefix += '.'

        self.bootstrapDir = os.path.join(outputDir, DefaultValues.BOOTSTRAP_DIR)
        checkDirExists(self.bootstrapDir)
        self.bootstrapTreeExternal = os.path.join(self.bootstrapDir, outputPrefix + DefaultValues.BOOTSTRAP_TREE_EXTERNAL)
        self.bootstrapTreeAll = os.path.join(self.bootstrapDir, outputPrefix + DefaultValues.BOOTSTRAP_TREE_ALL)

        self.jkMarkersDir = os.path.join(outputDir, DefaultValues.JK_MARKERS_DIR)
        checkDirExists(self.jkMarkersDir)
        self.jkMarkerTreeExternal = os.path.join(self.jkMarkersDir, outputPrefix + DefaultValues.JK_MARKERS_TREE_EXTERNAL)
        self.jkMarkerTreeAll = os.path.join(self.jkMarkersDir, outputPrefix + DefaultValues.JK_MARKERS_TREE_ALL)
        self.jkTaxaDir = os.path.join(outputDir, DefaultValues.JK_TAXA_DIR)
        checkDirExists(self.jkTaxaDir)
        self.jkTaxaTreeExternal = os.path.join(self.jkTaxaDir, outputPrefix + DefaultValues.JK_TAXA_TREE_EXTERNAL)
        self.jkTaxaTreeAll = os.path.join(self.jkTaxaDir, outputPrefix + DefaultValues.JK_TAXA_TREE_ALL)

        self.combinedTreeExternal = os.path.join(outputDir, outputPrefix + supportType + '.' + DefaultValues.COMBINED_TREE_EXTERNAL)
        self.combinedTreeAll = os.path.join(outputDir, outputPrefix + supportType + '.' + DefaultValues.COMBINED_TREE_ALL)

    def collectSupportValues(self, tree):
        """Get support value for each node in pre-order traversal."""
        preOrderSupport = []
        for node in tree.preorder_node_iter():
            if node.is_internal():
                preOrderSupport.append(node.label)

        return preOrderSupport

    def __combine(self, bootstrapTree, jkMarkerTree, jkTaxaTree, combinedTree, supportType):
        """Create new tree indicating average of support values."""

        tree = dendropy.Tree.get_from_path(bootstrapTree, schema='newick', as_rooted=True, preserve_underscores=True)
        bootstrapSupport = self.collectSupportValues(tree)

        tree = dendropy.Tree.get_from_path(jkMarkerTree, schema='newick', as_rooted=True, preserve_underscores=True)
        jkMarkerSupport = self.collectSupportValues(tree)

        tree = dendropy.Tree.get_from_path(jkTaxaTree, schema='newick', as_rooted=True, preserve_underscores=True)
        jkTaxaSupport = self.collectSupportValues(tree)

        internalNodeNum = 0
        for node in tree.preorder_node_iter():
            if node.is_internal():
                if supportType == 'average':
                    support = (int(bootstrapSupport[internalNodeNum]) + int(jkMarkerSupport[internalNodeNum]) + int(jkTaxaSupport[internalNodeNum])) / 3.0
                elif supportType == 'minimum':
                    support = min(int(bootstrapSupport[internalNodeNum]), int(jkMarkerSupport[internalNodeNum]), int(jkTaxaSupport[internalNodeNum]))

                node.label = '%s' % str(int(support + 0.5))
                internalNodeNum += 1

        tree.write_to_path(combinedTree, schema='newick', suppress_rooting=True, unquoted_underscores=True)

    def run(self, treesToInfer, supportType):
        """Create new trees indicating average of support values.."""

        timeKeeper = TimeKeeper()

        if treesToInfer == 'both' or treesToInfer == 'external':
            self.logger.info('  Calculating %s support for external tree.' % supportType)
            self.__combine(self.bootstrapTreeExternal, self.jkMarkerTreeExternal, self.jkTaxaTreeExternal, self.combinedTreeExternal, supportType)

        if treesToInfer == 'both' or treesToInfer == 'all':
            self.logger.info('  Calculating %s support for all tree.' % supportType)
            self.__combine(self.bootstrapTreeAll, self.jkMarkerTreeAll, self.jkTaxaTreeAll, self.combinedTreeAll, supportType)

        timeKeeper.printTimeStamp()
