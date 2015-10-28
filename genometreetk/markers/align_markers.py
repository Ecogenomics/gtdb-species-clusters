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
import multiprocessing as mp
import logging
from collections import defaultdict

import biolib.seq_io as seq_io
from biolib.external.hmmer import HMMER

from genometreetk.default_values import DefaultValues


class AlignMarkers(object):
    """Align genes to HMM."""

    def __init__(self, cpus):
        """Initialize.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """

        self.logger = logging.getLogger()

        self.cpus = cpus

        self.protein_file_ext = DefaultValues.PROTEIN_FILE_EXTENSION
        self.pfam_extension = DefaultValues.PFAM_EXTENSION
        self.tigr_extension = DefaultValues.TIGR_EXTENSION

    def _genes_in_genomes(self, genome_ids, genome_dirs):
        """Get genes within genomes.

        Parameters
        ----------
        genome_ids : iterable
            Genomes of interest.
        genome_dirs : d[assembly_accession] -> directory
            Path to files for individual genomes.

        Returns
        -------
        d[genome_id][family_id] -> [(gene_id_1, bitscore), ..., (gene_id_N, bitscore)]
            Genes within each genome.
        """

        genes_in_genome = {}
        for genome_id in genome_ids:
            genome_dir = genome_dirs[genome_id]

            marker_id_to_gene_id = defaultdict(list)

            assembly = genome_dir[genome_dir.rfind('/') + 1:]
            tophit_file = os.path.join(genome_dir, assembly + self.pfam_extension)
            with open(tophit_file) as f:
                f.readline()
                for line in f:
                    line_split = line.split('\t')

                    gene_id = line_split[0]
                    hits = line_split[1].split(';')
                    for hit in hits:
                        pfam_id, _evalue, bitscore = hit.split(',')
                        marker_id_to_gene_id[pfam_id].append((gene_id, float(bitscore)))

            tophit_file = os.path.join(genome_dir, assembly + self.tigr_extension)
            with open(tophit_file) as f:
                f.readline()
                for line in f:
                    line_split = line.split('\t')

                    gene_id = line_split[0]
                    hits = line_split[1].split(';')
                    for hit in hits:
                        tigrfam_id, _evalue, bitscore = hit.split(',')
                        marker_id_to_gene_id[tigrfam_id].append((gene_id, float(bitscore)))

            genes_in_genome[genome_id] = marker_id_to_gene_id

        return genes_in_genome

    def _run_hmm_align(self, genome_ids,
                                genome_dirs,
                                genes_in_genomes,
                                ignore_multi_copy,
                                output_msa_dir,
                                output_model_dir,
                                queue_in,
                                queue_out):
        """Run each marker gene in a separate thread.

        Only the gene with the highest bitscore is used for genomes with
        multiple hits to a given protein family.

        Parameters
        ----------
        genome_ids : iterable
            Genomes of interest.
        genome_dirs : d[assembly_accession] -> directory
            Path to files for individual genomes.
        genes_in_genomes : d[genome_id][family_id] -> [(gene_id_1, bitscore), ..., (gene_id_N, bitscore)]
            Genes within each genome.
        ignore_multi_copy : bool
            Flag indicating if genes with multiple hits should be ignored (True) or the gene with the highest bitscore taken (False).
        output_msa_dir : str
            Output directory for multiple sequence alignment.
        output_model_dir : str
            Output directory for HMMs.
        queue_in : Queue
            Input queue for parallel processing.
        queue_out : Queue
            Output queue for parallel processing.
        """

        while True:
            marker_id = queue_in.get(block=True, timeout=None)
            if marker_id == None:
                break

            marker_seq_file = os.path.join(output_msa_dir, marker_id + '.faa')
            fout = open(marker_seq_file, 'w')
            for genome_id in genome_ids:
                genome_dir = genome_dirs[genome_id]

                assembly = genome_dir[genome_dir.rfind('/') + 1:]
                genes_file = os.path.join(genome_dir, assembly + self.protein_file_ext)
                seqs = seq_io.read_fasta(genes_file)

                hits = genes_in_genomes[genome_id].get(marker_id, None)
                if not hits or (ignore_multi_copy and len(hits) > 1):
                    continue

                # get gene with highest bitscore
                hits.sort(key=lambda x: x[1], reverse=True)
                gene_id, _bitscore = hits[0]

                fout.write('>' + genome_id + DefaultValues.SEQ_CONCAT_CHAR + gene_id + '\n')
                fout.write(seqs[gene_id] + '\n')
            fout.close()

            hmmer = HMMER('align')
            hmmer.align(os.path.join(output_model_dir, marker_id + '.hmm'), marker_seq_file, os.path.join(output_msa_dir, marker_id + '.aln.faa'), trim=False, outputFormat='Pfam')
            self._mask_alignment(os.path.join(output_msa_dir, marker_id + '.aln.faa'), os.path.join(output_msa_dir, marker_id + '.aln.masked.faa'))

            queue_out.put(marker_id)

    def _report_threads(self, num_genes, writer_queue):
        """Report progress of parallel processing.

        Parameters
        ----------
        num_genes : int
            Number of genes being processed.
        writer_queue : Queue
            Output queue for parallel processing.
        """

        num_processed_genes = 0
        while True:
            marker_id = writer_queue.get(block=True, timeout=None)
            if marker_id == None:
                break

            num_processed_genes += 1
            statusStr = '==> Finished processing %d of %d (%.2f%%) marker genes.' % (num_processed_genes, num_genes, float(num_processed_genes) * 100 / num_genes)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def _mask_alignment(self, input_file, output_file):
        """Read HMMER alignment in STOCKHOLM format and output masked alignment in FASTA format.

        Parameters
        ----------
        input_file : str
            Input sequence file in STOCKHOLM format.
        output_file : str
            Output sequence file in FASTA format.
        """

        # read STOCKHOLM alignment
        seqs = {}
        for line in open(input_file):
            line = line.rstrip()
            if line == '' or line[0] == '#' or line == '//':
                if 'GC RF' in line:
                    mask = line.split('GC RF')[1].strip()
                continue
            else:
                line_split = line.split()
                seqs[line_split[0]] = line_split[1].upper().replace('.', '-').strip()

        # output masked sequences in FASTA format
        fout = open(output_file, 'w')
        for seq_id, seq in seqs.iteritems():
            fout.write('>' + seq_id + '\n')

            masked_seq = ''.join([seq[i] for i in xrange(0, len(seq)) if mask[i] == 'x'])
            fout.write(masked_seq + '\n')
        fout.close()

    def run(self, genome_ids, genome_dirs, marker_genes, ignore_multi_copy, output_msa_dir, output_model_dir):
        """Perform multithreaded alignment of marker genes using HMM align.

        Parameters
        ----------
        genome_ids : iterable
            Genomes of interest.
        genome_dirs : d[assembly_accession] -> directory
            Path to files for individual genomes.
        marker_genes : iterable
            Unique ids of marker genes to align.
        ignore_multi_copy : bool
            Flag indicating if genes with multiple hits should be ignored (True) or the gene with the highest bitscore taken (False).
        output_msa_dir : str
            Output directory for multiple sequence alignments.
        output_model_dir : str
            Output directory for HMMs.
        """

        # get mapping of marker ids to gene ids for each genome
        self.logger.info('Determining genes in genomes of interest.')
        genes_in_genomes = self._genes_in_genomes(genome_ids, genome_dirs)

        # align marker genes
        self.logger.info('Aligning marker genes:')
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()

        for _, marker_id in enumerate(marker_genes):
            worker_queue.put(marker_id)

        for _ in range(self.cpus):
            worker_queue.put(None)

        try:
            calc_proc = [mp.Process(target=self._run_hmm_align, args=(genome_ids,
                                                                      genome_dirs,
                                                                      genes_in_genomes,
                                                                      ignore_multi_copy,
                                                                      output_msa_dir,
                                                                      output_model_dir,
                                                                      worker_queue,
                                                                      writer_queue)) for _ in range(self.cpus)]
            write_proc = mp.Process(target=self._report_threads, args=(len(marker_genes), writer_queue))

            write_proc.start()

            for p in calc_proc:
                p.start()

            for p in calc_proc:
                p.join()

            writer_queue.put(None)
            write_proc.join()
        except:
            for p in calc_proc:
                p.terminate()

            write_proc.terminate()
