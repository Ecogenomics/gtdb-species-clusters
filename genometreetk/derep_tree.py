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
import math
import logging

from genometreetk.common import read_gtdb_metadata

from biolib.seq_io import read_seq
from biolib.newick import parse_label

import dendropy


class DereplicateTree(object):
    """Dereplicate tree."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger()
        
    def _derep_msa(self, msa_file, selected_taxa, output_msa):
        """Dereplicate multiple sequence alignment."""
        
        selected_taxa_labels = set()
        for taxon in selected_taxa:
            selected_taxa_labels.add(taxon.label)
        
        fout = open(output_msa, 'w')
        for seq_id, seq, annotation in read_seq(msa_file, keep_annotation=True):
            if seq_id in selected_taxa_labels:
                fout.write('>%s %s\n' % (seq_id, annotation))
                fout.write('%s\n' % seq)
        fout.close()
        
    def _derep_lineage(self, node, num_taxa_to_retain, genome_metadata):
        """Select genomes from lineage."""
        
        # rank all genomes in lineage with GTDB representatives first,
        # followed by genomes of decreasing quality
        rep_list = []
        taxa_list = []
        for leaf in node.leaf_iter():
            comp, cont, rep = genome_metadata[leaf.taxon.label]
            qual = float(comp) - 5*float(cont)
            if rep == 't' or rep == 'True' or rep == 'true':
                rep_list.append((qual, leaf.taxon))
            else:
                taxa_list.append((qual, leaf.taxon))
                
        sorted_list = sorted(rep_list, reverse=True) + sorted(taxa_list, reverse=True)
        
        taxa_to_keep = []
        for _qual, taxon in sorted_list[0:num_taxa_to_retain]:
            taxa_to_keep.append(taxon)
        
        return taxa_to_keep
        
    def _select_taxa(self, tree,
                            node_of_interest, 
                            outgroup_node, 
                            num_taxa_to_retain,
                            keep_unclassified,                            
                            genome_metadata):
        """Select genomes in named lineages on path from ingroup to outgroup."""
        
        # get most recent common ancestor of outgroup and lineage of interest
        outgroup_leaf_taxon = outgroup_node.leaf_iter().next().taxon
        lineage_of_interest_taxon = node_of_interest.leaf_iter().next().taxon
        mrca = tree.mrca(taxa=[outgroup_leaf_taxon, lineage_of_interest_taxon])
        
        # get taxon of lineage of interest
        taxa_of_interest = []
        parent = node_of_interest
        while parent != mrca:
            _support, taxon, _auxiliary_info = parse_label(parent.label)
            if taxon:
                taxa_of_interest.append(taxon)
            parent = parent.parent_node
            
        self.logger.info('Taxonomy for lineage of interest: %s' % ';'.join(taxa_of_interest))

        # select taxa from named lineages by traversing tree
        # in preorder and terminating descent at named taxa
        # not in path to lineage of interest
        selected_taxa = []
        
        stack = []
        for c in mrca.child_node_iter():
            stack.append(c)
            
        while stack:
            cur_node = stack.pop()
            
            if cur_node.is_leaf() and keep_unclassified:
                selected_taxa.append(cur_node.taxon)
            
            _support, taxon, _auxiliary_info = parse_label(cur_node.label)
            
            if taxon and taxon not in taxa_of_interest:
                # select roughly equal taxa from each child lineage to
                # enure we retain the correct depth (and the named node)
                # for this lineage
                derep_taxa = []
                num_children = sum([1 for c in cur_node.child_node_iter()])
                child_taxa_to_sample = int(math.ceil((1.0/num_children)*num_taxa_to_retain))
                for i, c in enumerate(cur_node.child_node_iter()):
                    taxa_to_sample = min(child_taxa_to_sample, num_taxa_to_retain - len(derep_taxa))
                    derep_taxa += self._derep_lineage(c, taxa_to_sample, genome_metadata)
 
                selected_taxa += derep_taxa
                self.logger.info('Selecting %d taxa from %s.' % (len(derep_taxa), taxon))
            elif cur_node == node_of_interest:
                self.logger.info('Retaining all taxa in lineage of interest.')
                for leaf in node_of_interest.leaf_iter():
                    selected_taxa.append(leaf.taxon) 
            else:
                for c in cur_node.child_node_iter():
                    stack.append(c)
                
        return selected_taxa

    def run(self, 
            input_tree, 
            lineage_of_interest, 
            outgroup,
            gtdb_metadata,
            num_taxa_to_retain,
            msa_file,
            keep_unclassified,
            output_dir):
        """Dereplicate tree.

        Parameters
        ----------
        input_tree : str
            Tree to dereplicate
        lineage_of_interest : str
            Named lineage where all taxa should be retain.
        outgroup : str
            Named lineage to use as outgroup.
        gtdb_metadata : str
            File containing metadata for taxa in tree.
        num_taxa_to_retain: int
            Taxa to retain in dereplicated lineages.
        msa_file : str
            Multiple sequence alignment to dereplicate along with tree.
        keep_unclassified : boolean
            Keep all taxa in unclassified lineages.
        output_dir:
            Output dir.
        """
        
        # read GTDB metadata
        self.logger.info('Reading metadata.')
        genome_metadata = read_gtdb_metadata(gtdb_metadata, ['checkm_completeness',
                                                              'checkm_contamination',
                                                              'gtdb_representative'])
        
        # read tree
        self.logger.info('Reading tree.')
        tree = dendropy.Tree.get_from_path(input_tree, 
                                            schema='newick', 
                                            rooting='force-rooted', 
                                            preserve_underscores=True)

        # locate node of interest and outgroup node
        self.logger.info('Identifying lineage of interest and outgroup.')
        node_of_interest = None
        outgroup_node = None
        for node in tree.preorder_node_iter():
            _support, taxon_str, _auxiliary_info = parse_label(node.label)
            
            if not taxon_str:
                continue
                
            for taxon in [t.strip() for t in taxon_str.split(';')]:
                if taxon == lineage_of_interest:
                    node_of_interest = node
                elif taxon == outgroup:
                    outgroup_node = node
                
        if not node_of_interest:
            self.logger.error('Could not find specified lineage of interest: %s' % lineage_of_interest)
            sys.exit()
            
        if not outgroup_node:
            self.logger.error('Could not find outgroup: %s' % outgroup)
            sys.exit()
                       
        # select taxa to retain
        self.logger.info('Selecting taxa to retain.')
        selected_taxa = self._select_taxa(tree,
                                            node_of_interest, 
                                            outgroup_node, 
                                            num_taxa_to_retain, 
                                            keep_unclassified,
                                            genome_metadata)

        self.logger.info('Retaining %d taxa.' % len(selected_taxa))
        
        # prune tree
        self.logger.info('Pruning tree.')
        tree.retain_taxa(selected_taxa)
        
        # dereplicate MSA if requested
        if msa_file:
            self.logger.info('Dereplicating MSA.')
            msa_name, msa_ext = os.path.splitext(os.path.basename(msa_file))
            output_msa = os.path.join(output_dir, msa_name + '.derep' + msa_ext)
            self._derep_msa(msa_file, selected_taxa, output_msa)
        
        # write out results
        tree_name, tree_ext = os.path.splitext(os.path.basename(input_tree))
        output_tree = os.path.join(output_dir, tree_name + '.derep' + tree_ext)
        tree.write_to_path(output_tree, 
                            schema='newick', 
                            suppress_rooting=True, 
                            unquoted_underscores=True)
