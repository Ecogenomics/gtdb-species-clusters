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
from collections import defaultdict

import biolib.seq_io as seq_io
from biolib.external.fasttree import FastTree

from genome_tree_tk.defaultValues import DefaultValues
from genome_tree_tk.jackknife_markers import JackknifeMarkers

import dendropy


class LgtTest(object):
    """Identify gene trees that may have undergone one or more lateral transfer.

    Specifically, the following test is applied:
      1) infer a jackknifed genome tree by randomly subsampling 50% of
         the marker genes under 100 replicates
      2) identify all splits with at least a minimum jackknife support value, and
         where at least a certain percentage of the taxa fall on each side of the split
      3) determine how many of these "well-support, internal" splits are recovered in each gene tree
      4) filter gene trees which do not recover at specific percentage of these splits
    """

    def __init__(self, cpus):
        """Initialize.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """

        self.logger = logging.getLogger()

        self.cpus = cpus

    def _create_concatenated_alignment(self, genome_ids,
                                       marker_genes,
                                       alignment_dir,
                                       concatenated_alignment_file,
                                       marker_file):
        """Create concatenated multiple sequence alignment for all genomes.

        Parameters
        ----------
        genome_ids : iterable
            Genomes of interest.
        marker_genes : iterable
            Unique ids of marker genes.
        alignment_dir : str
            Directory containing multiple sequence alignments.
        concatenated_alignment_file : str
            File to containing concatenated alignment.
        marker_file : str
            File indicating length of each marker in the alignment.
        """

        # Read alignment files. Some genomes may have multiple
        # copies of a marker gene in which case the last one
        # is arbitrarily taken. This is acceptable as all genes
        # are already screen to be conspecific.
        alignments = defaultdict(dict)
        marker_length = {}
        for mg in marker_genes:
            f = mg + '.aln.masked.faa'
            seqs = seq_io.read_fasta(os.path.join(alignment_dir, f))

            for seq_id, seq in seqs.iteritems():
                genome_id = seq_id[0:seq_id.find(DefaultValues.SEQ_CONCAT_CHAR)]

                alignments[mg][genome_id] = seq

                marker_length[mg] = len(seq)

        # create marker file
        fout = open(marker_file, 'w')
        for mg in marker_genes:
            fout.write('%s\t%s\t%s\t%d\n' % (mg, mg, mg, marker_length[mg]))
        fout.close()

        # create concatenated alignment
        concatenated_seqs = {}
        for mg in marker_genes:
            seqs = alignments[mg]

            for genome_id in genome_ids:
                if genome_id in seqs:
                    # append alignment
                    concatenated_seqs[genome_id] = concatenated_seqs.get(genome_id, '') + seqs[genome_id]
                else:
                    # missing gene
                    concatenated_seqs[genome_id] = concatenated_seqs.get(genome_id, '') + '-' * marker_length[mg]

        # save concatenated alignment
        seq_io.write_fasta(concatenated_seqs, concatenated_alignment_file)

    def run(self, genome_ids,
                    marker_genes,
                    hmm_model_file,
                    min_support,
                    min_per_taxa,
                    min_per_splits,
                    perc_markers_to_jackknife,
                    gene_tree_dir,
                    alignment_dir,
                    output_dir):
        """Identify gene trees which do not recover well-support, internal splits in a jackknifed genome tree.

        Parameters
        ----------
        genome_ids : iterable
            Genomes of interest.
        marker_genes : iterable
            Unique ids of marker genes.
        hmm_model_file : str
            File containing HMMs for each marker gene.
        min_support : float
            Minimum jackknife support of splits to use during LGT filtering [0, 1].
        min_per_taxa : float
            Minimum percentage of taxa required to consider a split during LGT filtering [0, 1].
        min_per_splits : float
            Minimum recovered well-supported splits required to retain marker during LGT filtering [0, 1].
        perc_markers_to_jackknife : float
            Percentage of taxa to keep during marker jackknifing [0, 1].
        gene_tree_dir : str
            Directory containing gene trees.
        alignment_dir : str
            Directory containing multiple sequence alignments.
        output_dir : str
            Output directory.
        """

        # create concatenated alignment file
        self.logger.info('')
        self.logger.info('  Concatenating alignments.')
        concatenated_alignment_file = os.path.join(output_dir, 'concatenated_alignment.faa')
        marker_file = os.path.join(output_dir, 'concatenated_markers.tsv')
        self._create_concatenated_alignment(genome_ids, marker_genes, alignment_dir, concatenated_alignment_file, marker_file)

        # create concatenated genome tree
        self.logger.info('  Inferring concatenated genome tree.')
        concatenated_tree = os.path.join(output_dir, 'concatenated.tree')
        concatenated_tree_log = os.path.join(output_dir, 'concatenated.tree.log')
        log_file = os.path.join(output_dir, 'fasttree.log')
        fast_tree = FastTree(True)
        fast_tree.run(concatenated_alignment_file, 'prot', 'wag', concatenated_tree, concatenated_tree_log, log_file)

        # calculate jackknife support values
        self.logger.info('  Calculating jackknife marker support values.')
        output_dir = os.path.join(output_dir, 'jackknife_markers')
        jackknife_markers = JackknifeMarkers(self.cpus)
        jackknife_tree = jackknife_markers.run(concatenated_tree, concatenated_alignment_file, marker_file, perc_markers_to_jackknife, 100, 'wag', output_dir)
        # jackknife_tree = os.path.join(output_dir, 'concatenated.jk_markers.tree')

        # identify well-support, internal splits
        tree = dendropy.Tree.get_from_path(jackknife_tree, schema='newick', rooting='force-unrooted', preserve_underscores=True)
        num_leaves = len(tree.leaf_nodes())

        num_internal_nodes = 0
        num_major_splits = 0
        well_supported_major_splits = 0

        splits = []
        for node in tree.internal_nodes():
            num_internal_nodes += 1

            num_node_leaves = len(node.leaf_nodes())
            if min(num_node_leaves, num_leaves - num_node_leaves) >= max(min_per_taxa * num_leaves, 2):
                num_major_splits += 1

                if int(node.label) > (min_support * 100.0):
                    well_supported_major_splits += 1
                    split = set([x.taxon.label for x in node.leaf_nodes()])
                    splits.append((split, node.edge_length))

        self.logger.info('    # internal nodes: %d' % num_internal_nodes)
        self.logger.info('    # major splits: %d' % num_major_splits)
        self.logger.info('    # well-supported, major splits: %d' % well_supported_major_splits)

        # filter gene trees that do not recover well-support, internal splits
        self.logger.info('')
        self.logger.info('  Filtering gene trees.')

        retained_markers = set()
        distances = {}
        for i, mg in enumerate(sorted(marker_genes)):
            sys.stdout.write('    Processed %d of %d (%.2f) gene trees.\r' % (i + 1, len(marker_genes), (i + 1) * 100.0 / len(marker_genes)))
            sys.stdout.flush()

            # read gene tree
            f = mg + '.tree'
            gene_tree_file = os.path.join(gene_tree_dir, f)
            gene_tree = dendropy.Tree.get_from_path(gene_tree_file, schema='newick', rooting='force-unrooted', preserve_underscores=True)

            # prune gene tree so each genome is present exactly once
            processed_genome_ids = set()
            taxa_to_prune = []
            for node in gene_tree.leaf_nodes():
                genome_id = node.taxon.label.split(DefaultValues.SEQ_CONCAT_CHAR)[0]

                if genome_id in processed_genome_ids or genome_id not in genome_ids:
                    taxa_to_prune.append(node.taxon)

                processed_genome_ids.add(genome_id)

            gene_tree.prune_taxa(taxa_to_prune)

            # rename nodes to contain only genome id
            gene_tree_taxa_set = set()
            for node in gene_tree.leaf_nodes():
                genome_id = node.taxon.label.split(DefaultValues.SEQ_CONCAT_CHAR)[0]
                node.taxon.label = genome_id
                gene_tree_taxa_set.add(genome_id)

            # re-encode the split system over the new taxon namespace
            gene_tree.migrate_taxon_namespace(dendropy.TaxonNamespace(gene_tree_taxa_set))
            gene_tree.encode_bipartitions()
            split_bitmasks = set(b.split_bitmask for b in gene_tree.bipartition_encoding)

            # determine number of splits recovered by or compatible with this gene tree
            recovered_splits = 0
            compatible_splits = 0
            compatible_edge_length = 0
            for split, edge_length in splits:
                common_taxa_labels = split.intersection(gene_tree_taxa_set)

                split = gene_tree.taxon_namespace.taxa_bitmask(labels=common_taxa_labels)
                normalized_split = dendropy.Bipartition.normalize_bitmask(
                                    bitmask=split,
                                    fill_bitmask=gene_tree.taxon_namespace.all_taxa_bitmask(),
                                    lowest_relevant_bit=1)

                if normalized_split in split_bitmasks:
                    recovered_splits += 1

                if gene_tree.is_compatible_with_bipartition(dendropy.Bipartition(bitmask=normalized_split, is_rooted=False)):
                    compatible_splits += 1
                    compatible_edge_length += edge_length

            if compatible_splits >= min_per_splits * len(splits):
                retained_markers.add(mg)

            perc_recovered_splits = recovered_splits * 100.0 / len(splits)
            perc_comp_splits = compatible_splits * 100.0 / len(splits)
            norm_comp_edge_length = float(compatible_edge_length) / sum([s[1] for s in splits])

            # calculate weighted Robinson-Foulds (Manhattan) and Felsenstein's Euclidean
            # distances to the concatenated genome tree
            pruned_tree = tree.clone(depth=2)
            pruned_tree.retain_taxa_with_labels(gene_tree.taxon_namespace.labels())
            pruned_tree.migrate_taxon_namespace(gene_tree.taxon_namespace)
            pruned_tree.encode_bipartitions()

            pruned_tree_edge_len = sum([e.length for e in pruned_tree.edges() if e.length])
            gene_tree_edge_len = sum([e.length for e in gene_tree.edges() if e.length])
            pruned_tree.scale_edges(1.0 / pruned_tree_edge_len)
            gene_tree.scale_edges(1.0 / gene_tree_edge_len)

            manhattan = dendropy.calculate.treecompare.weighted_robinson_foulds_distance(pruned_tree, gene_tree)
            euclidean = dendropy.calculate.treecompare.euclidean_distance(pruned_tree, gene_tree)

            distances[mg] = (perc_recovered_splits, perc_comp_splits, norm_comp_edge_length, manhattan, euclidean)

        sys.stdout.write('\n')
        self.logger.info('    Filtered gene trees: %d' % (len(marker_genes) - len(retained_markers)))
        self.logger.info('    Retained gene trees: %d' % len(retained_markers))

        return retained_markers, distances, num_internal_nodes, num_major_splits, well_supported_major_splits
