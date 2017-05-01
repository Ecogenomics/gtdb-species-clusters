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
from collections import defaultdict

from math import floor

import dendropy

from biolib.common import is_float
from biolib.taxonomy import Taxonomy


class PhylogeneticDiversity():
    """Calculate phylogenetic diversity."""

    def __init__(self):
        """Initialize."""
        self.logger = logging.getLogger()
        
    def _read_taxa_list(self, taxa_list):
        """Read taxa from file."""
        
        taxa = set()
        for line in open(taxa_list):
            taxa.add(line.strip().split('\t')[0])
            
        return taxa
        
    def _total_pd(self, tree):
        """Calculate PD over entire tree."""
        
        total_pd = 0
        total_taxa = 0
        for node in tree.preorder_node_iter():
            if node.parent_node is not None:
                total_pd += node.edge.length
                if node.is_leaf():
                    total_taxa += 1
                    
        return total_pd, total_taxa
        
    def _read_reps(self, rep_list):
        """Read genomes assigned to a representative."""
        
        genome_reps = {}
        
        if rep_list:
            for line in open(rep_list):
                line_split = line.strip().split('\t')
                
                rep_id = line_split[0]
                
                if len(line_split) >= 2:
                    for genome_id in map(str.strip, line_split[1].split(',')):
                        genome_reps[genome_id] = rep_id
                
        return genome_reps
        
    def _include_reps(self, ingroup_taxa, outgroup_taxa, genome_reps, is_ingroup):
        """Expand set of taxa to include representatives of taxa."""
        
        rep_is_ingroup = []
        rep_is_outgroup = []
        if is_ingroup:
            taxa_with_reps = set(ingroup_taxa)
            for genome_id in ingroup_taxa:
                if genome_id in genome_reps:
                    rep_id = genome_reps[genome_id]
                    if rep_id == genome_id:
                        continue
                        
                    taxa_with_reps.add(rep_id)
                    if rep_id in ingroup_taxa:
                        rep_is_ingroup.append(rep_id)
                    else:
                        rep_is_outgroup.append(rep_id)
        else:
            taxa_with_reps = set(outgroup_taxa)
            for genome_id in genome_reps:
                if genome_id in ingroup_taxa:
                    continue
                    
                rep_id = genome_reps[genome_id]
                if rep_id == genome_id:
                    continue
                    
                if rep_id in ingroup_taxa:
                    taxa_with_reps.add(rep_id)
                    rep_is_ingroup.append(rep_id)
                elif rep_id in outgroup_taxa:
                    taxa_with_reps.add(rep_id)
                    rep_is_outgroup.append(rep_id)
                else:
                    # this should never happend
                    self.logger.warning('There is an outgroup taxa (%s) with a representative (%s) not present in the tree.' % (genome_id, rep_id))
                        
        return taxa_with_reps, rep_is_ingroup, rep_is_outgroup
        
    def _taxon_pd(self, tree, ingroup, out_taxa_with_reps, genome_reps):
        """Calculate phylogenetic gain of each ingroup taxon relative to outgroup."""

        pg_taxon = {}
        for taxon in ingroup:
            rep_id = genome_reps.get(taxon, None)
            if rep_id and rep_id not in ingroup:
                pg_taxon[taxon] = [0, rep_id + ' (assigned to outgroup representative)']
            
        for leaf in tree.leaf_node_iter():
            if leaf.taxon.label in ingroup:
                # find first internal node containing an outgroup taxon
                pg = 0
                parent = leaf
                outgroup_taxon = 'None'
                while parent:
                    # check for outgroup taxon
                    stop = False
                    for tip in parent.leaf_iter():
                        if tip.taxon.label in out_taxa_with_reps:
                            outgroup_taxon = tip.taxon.label
                            if outgroup_taxon in ingroup:
                                outgroup_taxon += ' (one or more outgroup taxa are assigned to this ingroup taxon)'
                            stop = True
                            
                    if stop:
                        break
                    
                    pg += parent.edge.length
                    parent = parent.parent_node
                    
                pg_taxon[leaf.taxon.label] = [pg, outgroup_taxon]
                
                # propagate information to genomes represented by this genome_id
                for genome_id, rep_id in genome_reps.iteritems():
                    if genome_id in ingroup and rep_id == leaf.taxon.label:
                        pg_taxon[genome_id] = [pg, outgroup_taxon]

        return pg_taxon

    def pd(self, tree, taxa_list, rep_list, per_taxa_pg_file):
        """Calculate phylogenetic diversity of extant taxa."""

        genome_reps = self._read_reps(rep_list)
        
        self.logger.info('Reading tree.')
        tree = dendropy.Tree.get_from_path(tree,
                                            schema='newick',
                                            rooting='force-rooted',
                                            preserve_underscores=True)
                                                                            
        # get total branch length of tree
        self.logger.info('Calculating total PD.')
        total_pd, total_taxa = self._total_pd(tree)
        self.logger.info('Total PD for %d taxa = %.2f' % (total_taxa, total_pd))

        # get PD of ingroup taxa
        self.logger.info('Calculating total PD for specified taxa.')
        ingroup = self._read_taxa_list(taxa_list)
        
        in_taxa = set()
        for leaf in tree.leaf_node_iter():
            if leaf.taxon.label in ingroup:
                in_taxa.add(leaf.taxon.label)
            
        #tree_in = dendropy.Tree(tree)
        #tree_in.retain_taxa_with_labels(ingroup)
        #in_pd, in_taxa = self._total_pd(tree_in)
        
        self.logger.info('Specified ingroup taxa: %d' % len(ingroup))
        self.logger.info('Ingroup taxa as representatives or singletons in tree: %d' % len(in_taxa))

        # calculate PD for ingroup with additional genomes assigned to a representative
        ingroup_with_reps, rep_is_ingroup, rep_is_outgroup = self._include_reps(ingroup, None, genome_reps, True)
        tree_in_with_reps = dendropy.Tree(tree)
        tree_in_with_reps.retain_taxa_with_labels(ingroup_with_reps)
        in_pd_with_reps, in_taxa_with_reps = self._total_pd(tree_in_with_reps)
        
        self.logger.info('Ingroup taxa represented by another ingroup taxa: %d' % len(rep_is_ingroup))
        self.logger.info('Ingroup taxa represented by an outgroup taxa: %d' % len(rep_is_outgroup))
        self.logger.info('Unique outgroup taxa representing an ingroup taxa: %d' % len(set(rep_is_outgroup)))
        ingroup_taxa = len(in_taxa) + len(rep_is_ingroup) + len(rep_is_outgroup)
        ingroup_taxa_derep = len(in_taxa) + len(set(rep_is_outgroup))
        
        # get PD of outgroup taxa
        outgroup = set()
        for leaf in tree.leaf_node_iter():
            if leaf.taxon.label not in ingroup:
                outgroup.add(leaf.taxon.label)
                
        self.logger.info('Outgroup taxa as representatives or singletons in tree: %d' % len(outgroup))
        
        # calculate PD for outgroup with additional genomes assigned to a representative
        #tree_out = dendropy.Tree(tree)
        #tree_out.retain_taxa_with_labels(outgroup)
        #out_pd, out_taxa = self._total_pd(tree_out)

        outgroup_with_reps, rep_is_ingroup, rep_is_outgroup = self._include_reps(ingroup, outgroup, genome_reps, False)
        tree_out_with_reps = dendropy.Tree(tree)
        tree_out_with_reps.retain_taxa_with_labels(outgroup_with_reps)
        out_pd_with_reps, out_taxa_with_reps = self._total_pd(tree_out_with_reps)
                
        self.logger.info('Outgroup taxa represented by another outgroup taxa: %d' % len(rep_is_outgroup))
        self.logger.info('Outgroup taxa represented by an ingroup taxa: %d' % len(rep_is_ingroup))
        self.logger.info('Unique ingroup taxa representing an outgroup taxa: %d' % len(set(rep_is_ingroup)))
        outgroup_taxa = len(outgroup) + len(rep_is_outgroup) + len(rep_is_ingroup)
        outgroup_taxa_derep = len(outgroup) + len(set(rep_is_ingroup))
        
        # calculate PG of each ingroup taxon relative to outgroup in requested
        if per_taxa_pg_file:
            self.logger.info('Calculating PG of each ingroup taxon relative to outgroup.')
            pg_taxon = self._taxon_pd(tree, ingroup, outgroup.union(rep_is_ingroup), genome_reps)
            
            fout = open(per_taxa_pg_file, 'w')
            fout.write('Taxon\tPG\tPercent PG\tFirst outgroup taxon\n')
            for taxon, pg_stats in pg_taxon.iteritems():
                pg, outgroup_taxon = pg_stats
                fout.write('%s\t%f\t%f\t%s\n' % (taxon, pg, pg * 100.0 / total_pd, outgroup_taxon))
            fout.close()
                   
        return total_pd, ingroup_taxa, ingroup_taxa_derep, in_pd_with_reps, outgroup_taxa, outgroup_taxa_derep, out_pd_with_reps
        
    def _clade_pd(self, tree, ingroup_count, outgroup_count):
        """Calculate PD for named clades."""
        
        pd = {}
        for node in tree.preorder_node_iter():
            if not node.label:
                continue

            taxon = None
            if ':' in node.label:
                _support, taxon = node.label.split(':')
            else:
                if not is_float(node.label):
                    taxon = node.label

            if taxon:
                taxon_pd = 0
                taxon_count = 0
                in_taxon_pd = 0
                in_taxon_count = 0
                in_taxon_derep = 0
                out_taxon_pd = 0
                out_taxon_count = 0
                out_taxon_derep = 0
                for nn in node.postorder_iter():
                    if nn == node:
                        continue
                        
                    # check if group contains taxa from
                    # the ingroup and/or outgroup
                    ingroup_leaves = False
                    outgroup_leaves = False
                    for leaf in nn.leaf_iter():
                        genome_id = leaf.taxon.label
                        if genome_id in ingroup_count:
                            ingroup_leaves = True
                        
                        if genome_id in outgroup_count:
                            outgroup_leaves = True
                            
                    if ingroup_leaves:
                        in_taxon_pd += nn.edge.length

                    if outgroup_leaves:
                        out_taxon_pd += nn.edge.length
                        
                    if nn.is_leaf():
                        genome_id = nn.taxon.label
                        
                        in_taxon_count += ingroup_count.get(genome_id, 0)
                        if genome_id in ingroup_count:
                            in_taxon_derep += 1
                            
                        out_taxon_count += outgroup_count.get(genome_id, 0)
                        if genome_id in outgroup_count:
                            out_taxon_derep += 1
                            
                    taxon_pd += nn.edge.length

                if taxon == 'd__Archaea':
                    print taxon_count, in_taxon_count, out_taxon_count, in_taxon_pd, out_taxon_pd
                pd[taxon] = [taxon_pd, in_taxon_pd, in_taxon_count, in_taxon_derep, out_taxon_pd, out_taxon_count, out_taxon_derep]
                
        return pd
        
    def pd_clade(self, decorated_tree, output_file, taxa_list, rep_list):
        """Calculate phylogenetic diversity of named groups."""
        
        # calculate PD for entire tree
        self.logger.info('Reading tree.')
        tree = dendropy.Tree.get_from_path(decorated_tree,
                                            schema='newick',
                                            rooting='force-rooted',
                                            preserve_underscores=True)
 
        self.logger.info('Calculating total PD.')
        total_pd, total_taxa = self._total_pd(tree)
        
        # get ingroup and outgroup taxa
        ingroup = self._read_taxa_list(taxa_list)
        self.logger.info('Specified ingroup taxa: %d' % len(ingroup))
        
        # get number of ingroup and outgroup genomes 
        # representated by each leaf node
        genome_reps = self._read_reps(rep_list)
        
        ingroup_count = defaultdict(int)
        outgroup_count = defaultdict(int)
        for genome_id in genome_reps:
            rep_id = genome_reps[genome_id]
                
            if genome_id in ingroup:
                ingroup_count[rep_id] += 1
            else:
                outgroup_count[rep_id] += 1
                
        # add count for singletons
        for leaf in tree.leaf_node_iter():
            genome_id = leaf.taxon.label
            
            if genome_id in ingroup and genome_id not in ingroup_count:
                ingroup_count[genome_id] = 1 # ingroup singleton
            elif genome_id not in ingroup and genome_id not in  outgroup_count:
                outgroup_count[genome_id] = 1 # outgroup singleton

        # PD for named groups
        self.logger.info('Calculating PD for named clades.')
        pd_clade = self._clade_pd(tree, ingroup_count, outgroup_count)
        
        print 'ingroup_count, outgroup_count', len(ingroup_count), len(outgroup_count)

        # report results
        fout = open(output_file, 'w')
        fout.write('Clade')
        fout.write('\tTaxa\tTaxa (derep)\tPD\tPercent PD')
        fout.write('\tOut Taxa\tOut Taxa (derep)\tOut PD\tOut Percent PD')
        fout.write('\tIn Taxa\tIn Taxa (derep)\tIn PD\tIn Percent PD')
        fout.write('\tIn PG\tIn Percent PG\n')
        
        ordered_taxa = Taxonomy().sort_taxa(pd_clade.keys())
        for taxon in ordered_taxa:
            taxon_pd, in_taxon_pd, in_taxon_count, in_taxon_derep, out_taxon_pd, out_taxon_count, out_taxon_derep = pd_clade[taxon]
            taxon_count = in_taxon_count + out_taxon_count
            taxon_derep = in_taxon_derep + out_taxon_derep
            in_taxon_pg = taxon_pd - out_taxon_pd
            
            taxon_pd = max(taxon_pd, 1e-9) # make sure PD is never exactly zero to avoid division errors
            
            row = taxon
            row += '\t%d\t%d\t%.2f\t%.2f' % (taxon_count, taxon_derep, taxon_pd, taxon_pd * 100 / total_pd)
            row += '\t%d\t%d\t%.2f\t%.2f' % (out_taxon_count, out_taxon_derep, out_taxon_pd, out_taxon_pd * 100 / taxon_pd)
            row += '\t%d\t%d\t%.2f\t%.2f' % (in_taxon_count, in_taxon_derep, in_taxon_pd, in_taxon_pd * 100 / taxon_pd)
            row += '\t%.2f\t%.2f' % (in_taxon_pg, in_taxon_pg * 100 / taxon_pd)
            fout.write(row + '\n')
        fout.close()
