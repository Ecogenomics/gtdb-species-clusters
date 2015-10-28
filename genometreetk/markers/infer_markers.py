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
import shutil
import logging
import ntpath
from collections import defaultdict

from genometreetk.default_values import DefaultValues
from genometreetk.markers.align_markers import AlignMarkers
from genometreetk.common import read_genome_id_file, read_genome_dir_file

from biolib.external.fasttree import FastTree

import dendropy

import pickle  # ***


class InferMarkers(object):
    """Identify ubiquitous, single-copy marker gene within a set of genomes.

    Currently, this class is specifically designed to work with reference
    genomes from NCBI. In particular, it assumes each genome is contain in
    a separate directory and that Pfam and TIGRFAMs annotations are
    available in a specific format.
    """

    def __init__(self, genome_dir_file, pfam_model_file, tigrfams_model_dir, cpus):
        """Initialization.

        Parameters
        ----------
        genome_dir_file : str
            File specifying directory for each genome.
        pfam_model_file : str
            File containing Pfam HMMs.
        tigrfams_model_dir : str
            Directory containing TIGRFAMs HMMs.
        cpus : int
            Number of cpus to use.
        """

        self.logger = logging.getLogger()

        self.genome_dir_file = genome_dir_file
        self.pfam_model_file = pfam_model_file
        self.tigrfams_model_dir = tigrfams_model_dir

        self.cpus = cpus

        self.protein_file_ext = DefaultValues.PROTEIN_FILE_EXTENSION
        self.pfam_extension = DefaultValues.PFAM_EXTENSION
        self.tigr_extension = DefaultValues.TIGR_EXTENSION

    def _read_gene_hits(self, table, genome_ids, genome_dirs, extension):
        """Read gene annotations from top hit file.

        Parameters
        ----------
        table : defaultdict(lambda: defaultdict(set))
            Table to populate with gene information.
        genome_ids : iterable
            Genomes of interest.
        genome_dirs : d[assembly_accession] -> directory
            Path to files for individual genomes.
        extension : str
            Extension of file containing top hits to each gene.
        """

        for genome_id in genome_ids:
            genome_dir = genome_dirs[genome_id]

            assembly = genome_dir[genome_dir.rfind('/') + 1:]
            tophit_file = os.path.join(genome_dir, assembly + extension)

            with open(tophit_file) as f:
                f.readline()

                for line in f:
                    line_split = line.split('\t')

                    gene_id = line_split[0]
                    hits = line_split[1].split(';')
                    for hit in hits:
                        protein_family, _evalue, _bitscore = hit.split(',')
                        table[protein_family][genome_id].add(gene_id)

    def _gene_count_table(self, genome_ids, genome_dirs):
        """Get Pfam and TIGRFAMs annotations for genomes.

        Parameters
        ----------
        genome_ids : iterable
            Genomes of interest.
        genome_dirs : d[assembly_accession] -> directory
            Path to files for individual genomes.

        Returns
        -------
        d[assembly_accession][genome_id] -> set([gene_id_1, ..., gene_id_N])
            Gene location of protein families within each genome.
        """

        table = defaultdict(lambda: defaultdict(set))
        self._read_gene_hits(table, genome_ids, genome_dirs, self.pfam_extension)
        self._read_gene_hits(table, genome_ids, genome_dirs, self.tigr_extension)

        return table

    def _marker_genes(self, genome_ids, gene_count_table, ubiquity_threshold, single_copy_threshold, output_file):
        """Identify genes meeting ubiquity and single-copy thresholds.

        Parameters
        ----------
        genome_ids : iterable
            Genomes of interest.
        gene_count_table : d[family_id][genome_id] -> set([gene_id_1, ..., gene_id_N])
            Gene location of protein families within each genome.
        ubiquity_threshold : float
            Threshold for defining a ubiquitous marker genes [0, 1].
        single_copy_threshold : float
            Threshold for defining a single-copy marker gene [0, 1].
        output_file : str
            Output file indicating ubiquity and single-copy values for all genes.

        Returns
        -------
        dict: d[protein family] -> (ubiquity, single copy)
            Marker genes satisfying selection criteria.
        """

        if ubiquity_threshold > 1 or single_copy_threshold > 1:
            self.logger.error('Ubiquity or single-copy threshold is invalid: %f, %f' % (ubiquity_threshold, single_copy_threshold))
            sys.exit(0)

        fout = open(output_file, 'w')
        fout.write('Model accession\tUbiquity\tSingle copy\n')

        # find genes meeting ubiquity and single-copy thresholds
        markers = {}
        for protein_family, genes_in_genomes in gene_count_table.iteritems():
            ubiquity = 0
            single_copy = 0

            for genome_id in genome_ids:
                count = len(genes_in_genomes.get(genome_id, []))

                if count > 0:
                    ubiquity += 1

                if count == 1:
                    single_copy += 1

            u = ubiquity * 100.0 / len(genome_ids)
            s = single_copy * 100.0 / ubiquity
            fout.write('%s\t%.1f\t%.1f\n' % (protein_family, u, s))

            if ubiquity >= (ubiquity_threshold * len(genome_ids)) and single_copy >= (single_copy_threshold * ubiquity):
                markers[protein_family] = (u, s)

        fout.close()

        return markers

    def _identify_redundant_hmms(self, marker_genes, gene_count_table, redundancy, output_file):
        """Identify HMMs that consistently hit the same gene.

        This function identifies redundant HMMs both between and within
        the set of Pfam and TIGRFAMs HMMs. Redundancy between these two
        protein families is common. Redundancy within the TIGRFAMs HMMs
        is also common as there are a number of lineage-specific HMMs
        (e.g., archaeal and universal). Preference is given to TIGRFAM
        HMMs over Pfam HMMs as the former often models full genes instead
        of just domains.

        Parameters
        ----------
        marker_genes : iterable
            Marker genes to process for redundancy.
        gene_count_table : d[family_id][genome_id] -> set([gene_id_1, ..., gene_id_N])
            Gene location of protein families within each genome.
        redundancy : float
            Threshold for declaring HMMs redundant.
        output_file : str
            Output file to contain list of HMMs deemed to be redundant.

        Returns
        -------
        set
            Marker genes identified as being redundant.
        """

        if redundancy < 1:
            self.logger.error('Redundancy threshold is invalid: %f.' % redundancy)
            sys.exit(0)

        fout = open(output_file, 'w')
        fout.write('Kept marker\tRedundant marker\n')

        marker_gene_list = list(marker_genes)

        # count number of times HMMs hit the same gene
        redundancy_count = defaultdict(lambda: defaultdict(int))
        for i in xrange(0, len(marker_gene_list)):
            marker_gene_i = marker_gene_list[i]
            genes_in_genomes_i = gene_count_table[marker_gene_i]

            for j in xrange(i + 1, len(marker_genes)):
                marker_gene_j = marker_gene_list[j]
                genes_in_genomes_j = gene_count_table[marker_gene_j]

                for genome_id in genes_in_genomes_i:
                    if genome_id in genes_in_genomes_j:
                        if genes_in_genomes_i[genome_id].intersection(genes_in_genomes_j[genome_id]):
                            redundancy_count[marker_gene_i][marker_gene_j] += 1

        # Identify HMMs consistently hitting the same gene across genomes.
        #
        # Note that the following sets of redundant families is at least possible:
        #  X,Y
        #  X,Z
        #
        # It is unclear what to do in this case since, in general, we would expect
        # Y to also be redundant with Z. The following always resolves each redundant
        # pair in order giving preference to TIGRFAMs and HMMs with lower numbers. This
        # will NOT always result in the largest possible set of HMMs (i.e., Y, Z may
        # be removed when one could just remove X), but this seems fair since it is unclear
        # how to resolve such situations.
        hmms_to_remove = set()
        for marker_gene_i in redundancy_count:
            for marker_gene_j, count in redundancy_count[marker_gene_i].iteritems():
                if count > redundancy:
                    if marker_gene_i in hmms_to_remove or marker_gene_j in hmms_to_remove:
                        # marker gene from this redundant pair is already marked for removal
                        continue

                    # preferentially discard PFAM models
                    if 'PF' in marker_gene_i and not 'PF' in marker_gene_j:
                        hmms_to_remove.add(marker_gene_i)
                        fout.write('%s\t%s\n' % (marker_gene_j, marker_gene_i))
                    elif not 'PF' in marker_gene_i and 'PF' in marker_gene_j:
                        hmms_to_remove.add(marker_gene_j)
                        fout.write('%s\t%s\n' % (marker_gene_i, marker_gene_j))
                    elif 'PF' in marker_gene_i and 'PF' in marker_gene_j:
                        # take Pfam model with lowest number as these tend
                        # to encode better known protein families
                        pfam_num_i = marker_gene_i.replace('PF', '')
                        pfam_num_i = int(pfam_num_i[0:pfam_num_i.find('.')])
                        pfam_num_j = marker_gene_j.replace('PF', '')
                        pfam_num_j = int(pfam_num_j[0:pfam_num_j.find('.')])

                        if pfam_num_i > pfam_num_j:
                            hmms_to_remove.add(marker_gene_i)
                            fout.write('%s\t%s\n' % (marker_gene_j, marker_gene_i))
                        else:
                            hmms_to_remove.add(marker_gene_j)
                            fout.write('%s\t%s\n' % (marker_gene_i, marker_gene_j))
                    else:
                        # take TIGRFAMs model with lowest number as these
                        # tend to be more universal and/or to encode better
                        # known protein families
                        tigr_num_i = int(marker_gene_i.replace('TIGR', ''))
                        tigr_num_j = int(marker_gene_j.replace('TIGR', ''))
                        if tigr_num_i > tigr_num_j:
                            hmms_to_remove.add(marker_gene_i)
                            fout.write('%s\t%s\n' % (marker_gene_j, marker_gene_i))
                        else:
                            hmms_to_remove.add(marker_gene_j)
                            fout.write('%s\t%s\n' % (marker_gene_i, marker_gene_j))

        fout.close()

        return hmms_to_remove

    def _fetch_marker_models(self, marker_genes, output_model_dir):
        """Save PFAM and TIGRFAM marker genes into individual model files.

        Parameters
        ----------
        marker_genes : iterable
            Marker genes to process.
        output_model_dir : str
            Directory to store HMM models.
        """

        marker_id_to_name = {}
        for line in open(self.pfam_model_file):
            if 'NAME' in line:
                name = line.split()[1].rstrip()
            elif 'ACC' in line:
                acc = line.split()[1].rstrip()
                marker_id_to_name[acc] = name

        for marker_id in marker_genes:
            if 'PF' in marker_id:
                os.system('hmmfetch ' + self.pfam_model_file + ' ' + marker_id_to_name[marker_id] + ' > ' + os.path.join(output_model_dir, marker_id + '.hmm'))
            else:
                model_file = os.path.join(self.tigrfams_model_dir, marker_id + '.HMM')
                os.system('hmmfetch ' + model_file + ' ' + marker_id + ' > ' + os.path.join(output_model_dir, marker_id + '.hmm'))

    def identify_marker_genes(self, ingroup_file,
                            ubiquity_threshold, single_copy_threshold, redundancy,
                            valid_marker_genes,
                            output_msa_dir, output_model_dir):
        """Identify ubiquitous, single-copy marker genes.

        Parameters
        ----------
        ingroup_file : str
            File specifying unique ids of ingroup genomes.
        ubiquity_threshold : float
            Threshold for defining ubiquity marker genes.
        single_copy_threshold : float
            Threshold for defining a single-copy marker gene.
        redundancy : float
            Threshold for declaring HMMs redundant.
        valid_marker_genes : iterable
            Restrict marker set to genes within this set.
        output_msa_dir : str
            Directory to store multiple sequence alignment of marker genes.
        output_model_dir : str
            Directory to store HMMs of marker genes.
        """

        # read directory for each genome
        genome_dirs = read_genome_dir_file(self.genome_dir_file)

        # read genomes within the ingroup
        ncbi_genome_ids, user_genome_ids = read_genome_id_file(ingroup_file)
        genome_ids = ncbi_genome_ids.union(user_genome_ids)
        self.logger.info('Ingroup genomes: %d' % len(genome_ids))
        self.logger.info('NCBI genomes: %d' % len(ncbi_genome_ids))
        self.logger.info('User genomes: %d' % len(user_genome_ids))

        # identify marker genes
        self.logger.info('Identifying marker genes.')
        gene_stats_file = os.path.join(output_model_dir, '..', 'gene_stats.all.tsv')
        gene_count_table = self._gene_count_table(genome_ids, genome_dirs)
        marker_gene_stats = self._marker_genes(genome_ids, gene_count_table, ubiquity_threshold, single_copy_threshold, gene_stats_file)

        # with open('tmp_marker_gene_list', 'wb') as f:
        #    pickle.dump(marker_gene_stats, f)

        # with open('tmp_marker_gene_list', 'rb') as f:
        #    marker_gene_stats = pickle.load(f)

        marker_genes = set(marker_gene_stats.keys())
        if valid_marker_genes:
            self.logger.info('Restricting %d identified markers to specified set of valid markers.' % len(marker_genes))
            marker_genes = marker_genes.intersection(valid_marker_genes)
        self.logger.info('Identified ubiquitous, single-copy marker genes: %d' % len(marker_genes))

        redundancy_out_file = os.path.join(output_model_dir, '..', 'redundant_markers.tsv')
        redundancy = redundancy * len(genome_ids)
        redundant_hmms = self._identify_redundant_hmms(marker_genes, gene_count_table, redundancy, redundancy_out_file)
        marker_genes = marker_genes - redundant_hmms
        self.logger.info('Marker genes identified as redundant between TIGRFAM and Pfam: %d' % len(redundant_hmms))
        self.logger.info('Remaining ubiquitous, single-copy marker genes: %d' % len(marker_genes))

        # get HMM for each marker gene
        self.logger.info('Fetching HMM for each marker genes.')
        self._fetch_marker_models(marker_genes, output_model_dir)

        # align gene sequences
        align_markers = AlignMarkers(self.cpus)
        align_markers.run(genome_ids, genome_dirs, marker_genes, False, output_msa_dir, output_model_dir)

        return len(genome_ids), len(ncbi_genome_ids), len(user_genome_ids), genome_ids, marker_gene_stats, marker_genes

    def infer_gene_trees(self, msa_dir, output_dir, extension):
        """Infer gene trees.

        Parameters
        ----------
        msa_dir : str
            Directory containing multiple sequence alignment of marker genes.
        output_dir : str
            Directory to store gene trees.
        extension : str
            Extension of multiple sequence alignment files.
        """

        files = os.listdir(msa_dir)
        msa_files = []
        for f in files:
            if f.endswith(extension):
                msa_file = os.path.join(msa_dir, f)
                msa_files.append(msa_file)

                fin = open(msa_file)
                data = fin.readlines()
                fin.close()

                fout = open(msa_file, 'w')
                for line in data:
                    if line[0] != '>':
                        # remove trailing star
                        if line[-1] == '*':
                            line = line[0:-1]
                    fout.write(line)
                fout.close()

        fasttree = FastTree(multithreaded=False)
        fasttree.parallel_run(msa_files, 'prot', 'wag', output_dir, self.cpus)

        # create gene tree without gene ids for visualization in ARB
        for msa_file in msa_files:
            tree_filename = ntpath.basename(msa_file)
            tree_prefix = tree_filename[0:tree_filename.find('.')]

            if tree_prefix.startswith('PF'):
                # patch up output file for Pfam trees
                old_tree_prefix = tree_prefix
                tree_prefix = '.'.join(tree_filename.split('.')[0:2])
                shutil.move(os.path.join(output_dir, old_tree_prefix + '.tree'),
                                os.path.join(output_dir, tree_prefix + '.tree'))

            gene_tree_file = os.path.join(output_dir, tree_prefix + '.tree')
            gene_tree = dendropy.Tree.get_from_path(gene_tree_file, schema='newick', rooting='force-unrooted', preserve_underscores=True)

            # rename nodes to contain only genome id
            for node in gene_tree.leaf_nodes():
                genome_id = node.taxon.label.split(DefaultValues.SEQ_CONCAT_CHAR)[0]
                node.taxon.label = genome_id

            output_tree_file = os.path.join(output_dir, tree_prefix + '.genome_ids.tree')
            gene_tree.write_to_path(output_tree_file, schema='newick', suppress_rooting=True, unquoted_underscores=True)
