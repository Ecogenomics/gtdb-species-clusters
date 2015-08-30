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

import dendropy

from biolib.external.hmmer import HmmModelParser

from genome_tree_tk.reroot_tree import RerootTree


class ParalogTest(object):
    """Identify marker genes where all paralogs are conspecific."""

    def __init__(self):
        """Initialize."""

        self.logger = logging.getLogger()

    def _patristic_dist(self, tree, taxa1, taxa2):
        """Calculate the patristic distance between two taxa.

        Parameters
        ----------
        tree : dendropy.Tree
            Tree to calculate distance over.
        taxa1 : dendropy.Node
            Taxon 1 to calculate distance from.
        taxa2 : dendropy.Node
            Taxon 2 to calculate distance from.
        """

        mrca = tree.mrca(taxon_labels=[taxa1.taxon.label, taxa2.taxon.label])

        if mrca.parent_node == None:
            # MRCA is the root of the tree
            return taxa1.distance_from_root() + taxa2.distance_from_root()
        else:

            dist = taxa1.edge_length
            parent_node = taxa1.parent_node
            while parent_node != mrca:
                dist += parent_node.edge_length
                parent_node = parent_node.parent_node

            dist += taxa2.edge_length
            parent_node = taxa2.parent_node
            while parent_node != mrca:
                dist += parent_node.edge_length
                parent_node = parent_node.parent_node

            return dist

    def run(self, marker_genes, hmm_model_file, gene_tree_dir, accept_per, taxonomy, extension, output_dir):
        """Determine if paralogs are conspecific.

        Parameters
        ----------
        marker_genes : iterable
            Marker genes of interest.
        hmm_model_file : str
            File containing HMMs for each marker gene.
        gene_tree_dir : str
            Directory containing gene trees to process.
        accept_per : float
            Non-conspecific threshold for retaining multi-copy gene trees.
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.
        extension : str
            Extension of gene trees to process.
        """

        # parse HMMER models of interest
        hmm_model_parse = HmmModelParser(hmm_model_file)
        hmm_models = hmm_model_parse.models()

        retained_markers = set()

        paralog_output_file = os.path.join(output_dir, 'paralog_filtering.tsv')
        fout = open(paralog_output_file, 'w')
        fout.write('Model accession\tName\tDescription\tLength\tNon-conspecific genomes\tRetained\n')
        for i, mg in enumerate(marker_genes):
            sys.stdout.write('    Processed %d of %d (%.2f) gene trees.\r' % (i + 1, len(marker_genes), (i + 1) * 100.0 / len(marker_genes)))
            sys.stdout.flush()

            f = mg + extension

            tree = dendropy.Tree.get_from_path(os.path.join(gene_tree_dir, f), schema='newick', rooting='force-unrooted', preserve_underscores=True)

            taxa = tree.leaf_nodes()
            num_taxa = len(taxa)

            # root tree at mid-point
            tree.reroot_at_midpoint(update_bipartitions=True)

            # get species name of each taxa
            leaf_node_to_species_name = {}
            for t in taxa:
                genome_id = t.taxon.label.split('|')[0]
                genus = taxonomy[genome_id][5]
                sp = taxonomy[genome_id][6].lower()

                leaf_node_to_species_name[t.taxon.label] = genus + ' ' + sp

            # find all paralogs
            paralogs = defaultdict(set)
            for i in xrange(0, len(taxa)):
                genome_id = taxa[i].taxon.label.split('|')[0]
                for j in xrange(i + 1, len(taxa)):
                    # genes from the same genome are paralogs, but we filter out
                    # those that are identical (distance of 0 on the tree) to
                    # speed up computation and because these clearly do not
                    # adversely effect phylogenetic inference
                    if genome_id == taxa[j].taxon.label.split('|')[0] and self._patristic_dist(tree, taxa[i], taxa[j]) > 0:
                        paralogs[genome_id].add(taxa[i].taxon.label)
                        paralogs[genome_id].add(taxa[j].taxon.label)

            # check if paralogs are conspecific
            non_conspecific_genomes = []
            for genome_id, taxa_labels in paralogs.iteritems():
                lca_node = tree.mrca(taxon_labels=taxa_labels)

                children = lca_node.leaf_nodes()
                species = set()
                for child in children:
                    child_genome_id = child.taxon.label.split('|')[0]

                    genus = taxonomy[child_genome_id][5]
                    sp = taxonomy[child_genome_id][6].lower()
                    species.add(genus + ' ' + sp)

                if len(species) > 1:
                    non_conspecific_genomes.append(genome_id)

            fout.write('%s\t%s\t%s\t%s\t%d\t%s\n' % (mg,
                                                       hmm_models[mg].name,
                                                       hmm_models[mg].desc,
                                                       hmm_models[mg].leng,
                                                       len(non_conspecific_genomes),
                                                       len(non_conspecific_genomes) <= accept_per * num_taxa))

            if len(non_conspecific_genomes) <= accept_per * num_taxa:
                retained_markers.add(mg)

        sys.stdout.write('\n')

        self.logger.info('')
        self.logger.info('  Paralog filtering results written to: %s' % paralog_output_file)

        self.logger.info('    Filtered gene trees: %d' % (len(marker_genes) - len(retained_markers)))
        self.logger.info('    Retained gene trees: %d' % len(retained_markers))

        return retained_markers
