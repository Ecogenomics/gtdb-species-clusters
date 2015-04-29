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
import random
import multiprocessing as mp

from math import floor

from genome_tree.fasttree import FastTree
from genome_tree_tk.tree_support import TreeSupport
from genome_tree.common import readTreeModel, makeSurePathExists

from checkm.defaultValues import DefaultValues as DefaultValuesCheckM
from checkm.util.seqUtils import readFasta


class JackknifeTaxa(object):
    """Assess robustness of genome tree by jackknifing marker genes in multiple sequence alignment."""

    def __init__(self, outputPrefix, outputDir):
        """Setup directories for jk_taxa command."""
        self.logger = logging.getLogger()

        self.reportOut = os.path.join(outputDir, DefaultValues.REPORT)
        self.treeModel = readTreeModel(self.reportOut)
        self.validOutgroup = os.path.join(outputDir, DefaultValues.INFER_DIR, DefaultValues.INFER_VALID_OUTGROUP)

        self.msaExternal = os.path.join(outputDir, DefaultValues.INFER_DIR, DefaultValues.INFER_EXTERNAL_MSA)
        self.treeExternal = os.path.join(outputDir, DefaultValues.INFER_EXTERNAL_REROOT_TREE)

        self.msaAll = os.path.join(outputDir, DefaultValues.INFER_DIR, DefaultValues.INFER_ALL_MSA)
        self.treeAll = os.path.join(outputDir, DefaultValues.INFER_ALL_REROOT_TREE)

        self.jkTaxaDir = os.path.join(outputDir, DefaultValues.JK_TAXA_DIR)
        makeSurePathExists(self.jkTaxaDir)

        self.jkTaxaTreesExternalDir = os.path.join(self.jkTaxaDir, DefaultValues.JK_TAXA_TREES_EXTERNAL_DIR)
        makeSurePathExists(self.jkTaxaTreesExternalDir)

        self.jkTaxaTreesAllDir = os.path.join(self.jkTaxaDir, DefaultValues.JK_TAXA_TREES_ALL_DIR)
        makeSurePathExists(self.jkTaxaTreesAllDir)

        if not outputPrefix[-1] == '.' and not outputPrefix[-1] == '_':
            outputPrefix += '.'
        self.jkTaxaTreeExternal = os.path.join(self.jkTaxaDir, outputPrefix + DefaultValues.JK_TAXA_TREE_EXTERNAL)
        self.jkTaxaTreeAll = os.path.join(self.jkTaxaDir, outputPrefix + DefaultValues.JK_TAXA_TREE_ALL)

    def __processingThread(self, fullMSA, jkTaxaTreesDir, percentTaxaToKeep, outgroupIds, queueIn, queueOut):
        """Infer tree from jackknifed alignments."""

        while True:
            repIndex = queueIn.get(block=True, timeout=None)
            if repIndex == None:
                break

            outputMSA = os.path.join(jkTaxaTreesDir, 'jk_taxa.msa.' + str(repIndex) + '.fna')
            self.__createAlignmentForTaxa(fullMSA, percentTaxaToKeep, outgroupIds, outputMSA)

            fastTree = FastTree(bMultithreaded=False)
            outputTree = os.path.join(jkTaxaTreesDir, 'jk_taxa.tree.' + str(repIndex) + '.tre')
            fastTreeOutput = os.path.join(jkTaxaTreesDir, 'jk_taxa.fasttree.' + str(repIndex) + '.out')
            fastTree.run(outputMSA, self.treeModel, outputTree, fastTreeOutput)

            queueOut.put(repIndex)

    def __reportingThread(self, numReplicates, queueIn):
        """Report number of processed replicates."""

        numProcessed = 0

        statusStr = '    Finished processing %d of %d (%.2f%%) replicates.' % (numProcessed, numReplicates, float(numProcessed) * 100 / numReplicates)
        sys.stderr.write('%s\r' % statusStr)
        sys.stderr.flush()

        while True:
            repIndex = queueIn.get(block=True, timeout=None)
            if repIndex == None:
                break

            numProcessed += 1
            if self.logger.getEffectiveLevel() <= logging.INFO:
                statusStr = '    Finished processing %d of %d (%.2f%%) replicates.' % (numProcessed, numReplicates, float(numProcessed) * 100 / numReplicates)
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')

    def __createAlignmentForTaxa(self, fullMSA, percentTaxaToKeep, outgroupIds, outputFile):
        """Create multiple sequence alignment for jackknifed taxa."""

        # randomly select ingroup taxa
        numIngroupTaxa = len(fullMSA) - len(outgroupIds)
        ingroupTaxaToKeep = random.sample(fullMSA.keys(), int(floor(numIngroupTaxa * percentTaxaToKeep)))

        taxaToKeep = set(ingroupTaxaToKeep).union(outgroupIds)

        fout = open(outputFile, 'w')
        for seqId, seq in fullMSA.iteritems():
            if seqId in taxaToKeep:
                fout.write('>' + seqId + '\n')
                fout.write(seq + '\n')
        fout.close()

    def __jackknife(self, msaFile, jkTaxaTreesDir, inputTree, jkTaxaTree, percentTaxaToKeep, numReplicates, threads):
        """Jackknife taxa to assess robustness of tree.

        Parameters
        ----------
        msaFile : str
          File containing multiple sequence alignment for all taxa.
        jkTaxaTreesDir : str
          Output directory for jackknifed trees
        inputTree : str
          Tree inferred with all data.
        jkTaxaTree : str
          Output tree with jackknife taxa support values.
        percentTaxaToKeep : flow
          Percentage of taxa to keep in each jackknifed tree.
        numReplicates : int
          Number of jackknife replicates to perform.
        threads : int
          Number of processors to use.
        """

        # read full MSA
        fullMSA = readFasta(msaFile)

        # read outgroup taxa
        outgroupIds = set()
        for line in open(self.validOutgroup):
            outgroupIds.add(line.strip())

        # infer genome trees with jackknifed taxa
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for repIndex in xrange(numReplicates):
            workerQueue.put(repIndex)

        for _ in range(threads):
            workerQueue.put(None)

        try:
            calcProc = [mp.Process(target=self.__processingThread, args=(fullMSA, jkTaxaTreesDir, percentTaxaToKeep, outgroupIds, workerQueue, writerQueue)) for _ in range(threads)]
            writeProc = mp.Process(target=self.__reportingThread, args=(numReplicates, writerQueue))

            writeProc.start()

            for p in calcProc:
                p.start()

            for p in calcProc:
                p.join()

            writerQueue.put(None)
            writeProc.join()
        except:
            # make sure all processes are terminated
            for p in calcProc:
                p.terminate()

            writeProc.terminate()

        # calculate support
        self.logger.info('')
        self.logger.info('  Calculating support.')

        repTreeFiles = []
        for repIndex in xrange(numReplicates):
            repTreeFiles.append(os.path.join(jkTaxaTreesDir, 'jk_taxa.tree.' + str(repIndex) + '.tre'))

        treeSupport = TreeSupport()
        treeSupport.supportOnSubsetTaxa(inputTree, repTreeFiles, jkTaxaTree)

    def run(self, percentTaxaToKeep, numReplicates, treesToInfer, threads):
        """Calculate jackknife for the external and 'all' tree."""

        timeKeeper = TimeKeeper()

        if treesToInfer == 'both' or treesToInfer == 'external':
            self.logger.info('  Inferring jackknife taxa support for external genome trees.')
            self.__jackknife(self.msaExternal, self.jkTaxaTreesExternalDir, self.treeExternal, self.jkTaxaTreeExternal, percentTaxaToKeep, numReplicates, threads)
            timeKeeper.printTimeStamp()

        if treesToInfer == 'both' or treesToInfer == 'all':
            self.logger.info('')
            self.logger.info('  Inferring jackknife taxa support for complete genome trees.')
            self.__jackknife(self.msaAll, self.jkTaxaTreesAllDir, self.treeAll, self.jkTaxaTreeAll, percentTaxaToKeep, numReplicates, threads)
            timeKeeper.printTimeStamp()

        # generate summary report
        fout = open(self.reportOut, 'a')
        fout.write('\n\n')
        fout.write('[jk_taxa]\n')
        fout.write('Percentage of taxa to keep: %f\n' % percentTaxaToKeep)
        fout.write('Number of replicates: %d\n' % numReplicates)
        fout.write(timeKeeper.getTimeStamp())
        fout.close()

        timeKeeper.printTimeStamp()
