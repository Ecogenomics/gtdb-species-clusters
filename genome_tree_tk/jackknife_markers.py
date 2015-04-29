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

from genome_tree.timeKeeper import TimeKeeper
from genome_tree.defaultValues import DefaultValues
from genome_tree.fasttree import FastTree
from genome_tree.treeSupport import TreeSupport
from genome_tree.common import readTreeModel, makeSurePathExists

from genome_tree.markers.rerootTree import RerootTree

from checkm.defaultValues import DefaultValues as DefaultValuesCheckM
from checkm.util.seqUtils import readFasta


class JackknifeMarkers(object):
    """Assess robustness of genome tree by jackkifing marker genes in multiple sequence alignment."""

    def __init__(self, outputPrefix, outputDir):
        """Setup directories for jk_markers command."""
        self.logger = logging.getLogger()

        self.reportOut = os.path.join(outputDir, DefaultValues.REPORT)
        self.treeModel = readTreeModel(self.reportOut)
        self.markerGeneStats = os.path.join(outputDir, DefaultValues.MARKER_GENE_STATS)
        self.validOutgroup = os.path.join(outputDir, DefaultValues.INFER_DIR, DefaultValues.INFER_VALID_OUTGROUP)

        self.msaExternal = os.path.join(outputDir, DefaultValues.INFER_DIR, DefaultValues.INFER_EXTERNAL_MSA)
        self.treeExternal = os.path.join(outputDir, DefaultValues.INFER_EXTERNAL_REROOT_TREE)

        self.msaAll = os.path.join(outputDir, DefaultValues.INFER_DIR, DefaultValues.INFER_ALL_MSA)
        self.treeAll = os.path.join(outputDir, DefaultValues.INFER_ALL_REROOT_TREE)

        self.jkMarkersDir = os.path.join(outputDir, DefaultValues.JK_MARKERS_DIR)
        makeSurePathExists(self.jkMarkersDir)

        self.jkMarkersTreesExternalDir = os.path.join(self.jkMarkersDir, DefaultValues.JK_MARKERS_TREES_EXTERNAL_DIR)
        makeSurePathExists(self.jkMarkersTreesExternalDir)

        self.jkMarkersTreesAllDir = os.path.join(self.jkMarkersDir, DefaultValues.JK_MARKERS_TREES_ALL_DIR)
        makeSurePathExists(self.jkMarkersTreesAllDir)

        if not outputPrefix[-1] == '.' and not outputPrefix[-1] == '_':
            outputPrefix += '.'
        self.jkMarkerTreeExternal = os.path.join(self.jkMarkersDir, outputPrefix + DefaultValues.JK_MARKERS_TREE_EXTERNAL)
        self.jkMarkerTreeAll = os.path.join(self.jkMarkersDir, outputPrefix + DefaultValues.JK_MARKERS_TREE_ALL)

    def __processingThread(self, fullMSA, jkMarkersTreesDir, percentMarkersToKeep, markerLengths, queueIn, queueOut):
        """Infer tree from jackknifed alignments."""

        while True:
            repIndex = queueIn.get(block=True, timeout=None)
            if repIndex == None:
                break

            outputMSA = os.path.join(jkMarkersTreesDir, 'jk_markers.msa.' + str(repIndex) + '.fna')
            self.__createJackknifedAlignment(fullMSA, percentMarkersToKeep, markerLengths, outputMSA)

            fastTree = FastTree(bMultithreaded=False)
            outputTree = os.path.join(jkMarkersTreesDir, 'jk_markers.tree.' + str(repIndex) + '.tre')
            fastTreeOutput = os.path.join(jkMarkersTreesDir, 'jk_markers.fasttree.' + str(repIndex) + '.out')
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

    def __createJackknifedAlignment(self, fullMSA, percentMarkersToKeep, markerLengths, outputFile):
        """Jackknife multiple sequence alignment to a subset of marker genes."""
        markerToKeep = random.sample(xrange(0, len(markerLengths)), int(floor(percentMarkersToKeep * len(markerLengths))))

        markerStartPos = [0]
        for index, ml in enumerate(markerLengths):
            markerStartPos.append(markerStartPos[index] + ml)

        mask = [0] * sum(markerLengths)
        for markerIndex in markerToKeep:
            start = markerStartPos[markerIndex]
            end = start + markerLengths[markerIndex]
            mask[start:end] = [1] * (end - start)

        fout = open(outputFile, 'w')
        for seqId, seq in fullMSA.iteritems():
            fout.write('>' + seqId + '\n')
            subSeq = "".join(nt if m == 1 else '' for nt, m in zip(seq, mask))
            fout.write(subSeq + '\n')
        fout.close()

    def __jackknife(self, msaFile, jkMarkersTreesDir, inputTree, jkMarkerTree, percentMarkersToKeep, numReplicates, threads):
        """Jackknife marker genes to assess robustness of tree."""

        # determine length of each marker gene in complete multiple sequence alignment
        markerLengths = []
        for line in open(self.markerGeneStats):
            lineSplit = line.split('\t')
            markerLengths.append(int(lineSplit[1]))

        # read full MSA
        fullMSA = readFasta(msaFile)

        # infer genome trees with jackknifed marker sets
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for repIndex in xrange(numReplicates):
            workerQueue.put(repIndex)

        for _ in range(threads):
            workerQueue.put(None)

        try:
            calcProc = [mp.Process(target=self.__processingThread, args=(fullMSA, jkMarkersTreesDir, percentMarkersToKeep, markerLengths, workerQueue, writerQueue)) for _ in range(threads)]
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
            repTreeFiles.append(os.path.join(jkMarkersTreesDir, 'jk_markers.tree.' + str(repIndex) + '.tre'))

        treeSupport = TreeSupport()
        treeSupport.supportOnCommonTaxa(inputTree, repTreeFiles, jkMarkerTree)

    def run(self, percentMarkersToKeep, numReplicates, treesToInfer, threads):
        """Calculate jackknife for the external and 'all' tree."""

        timeKeeper = TimeKeeper()

        if treesToInfer == 'both' or treesToInfer == 'external':
            self.logger.info('  Inferring jackknife marker support for external genome trees.')
            self.__jackknife(self.msaExternal, self.jkMarkersTreesExternalDir, self.treeExternal, self.jkMarkerTreeExternal, percentMarkersToKeep, numReplicates, threads)
            timeKeeper.printTimeStamp()

        if treesToInfer == 'both' or treesToInfer == 'all':
            self.logger.info('')
            self.logger.info('  Inferring jackknife marker support for complete genome trees.')
            self.__jackknife(self.msaAll, self.jkMarkersTreesAllDir, self.treeAll, self.jkMarkerTreeAll, percentMarkersToKeep, numReplicates, threads)
            timeKeeper.printTimeStamp()

        # generate summary report
        fout = open(self.reportOut, 'a')
        fout.write('\n\n')
        fout.write('[jk_markers]\n')
        fout.write('Percentage of markers kept: %f\n' % percentMarkersToKeep)
        fout.write('Number of replicates: %d\n' % numReplicates)
        fout.write(timeKeeper.getTimeStamp())
        fout.close()

        timeKeeper.printTimeStamp()
