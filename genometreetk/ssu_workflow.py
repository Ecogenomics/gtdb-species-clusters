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

import biolib.seq_io as seq_io
from biolib.misc.time_keeper import TimeKeeper
from biolib.external.fasttree import FastTree

import genometreetk.ncbi as ncbi
from genometreetk.common import (read_gtdb_genome_quality,
                                    read_genome_dir_file)


class SSU_Workflow(object):
    """Infer SSU trees."""

    def __init__(self, gtdb_metadata_file, genome_dir_file):
        """Initialization.

        Parameters
        ----------
        gtdb_metadata_file : str
            File specifying GTDB metadata for each genome.
        genome_dir_file : str
            File specifying directory for each genome.
        """

        self.logger = logging.getLogger()

        self.gtdb_metadata_file = gtdb_metadata_file
        self.genome_dir_file = genome_dir_file

    def _get_ssu_seqs(self, genomes_to_consider, output_dir):
        """Get 16S sequences from genomes.

        Parameters
        ----------
        genomes_to_consider : iterable
            IDs of genomes to obtain 16S sequences.
        output_dir : str
            Desired output directory.

        Returns
        -------
        str
            Name of file with 16S sequences in FASTA format.
        """


        genome_dirs = read_genome_dir_file(self.genome_dir_file)

        ssu_output_file = os.path.join(output_dir, 'ssu_reps.fna')
        fout = open(ssu_output_file, 'w')

        no_identified_ssu = 0
        for genome_id in genomes_to_consider:
            ssu_hmm_file = os.path.join(genome_dirs[genome_id], 'ssu_gg_2013_08', 'ssu.hmm_summary.tsv')

            if not os.path.exists(ssu_hmm_file):
                no_identified_ssu += 1
                self.logger.warning('RefSeq representative %s has no identified 16S sequence.' % genome_id)
                continue


            with open(ssu_hmm_file) as f:
                headers = f.readline().rstrip().split('\t')
                hmm_model_index = headers.index('HMM')
                length_index = headers.index('SSU gene length')

                longest_query_len = 0
                best_hmm_model = None
                seq_info = None
                seq_id = None
                for line in f:
                    line_split = line.rstrip().split('\t')
                    query_len = int(line_split[length_index])
                    hmm_model = line_split[hmm_model_index]
                    if query_len > longest_query_len or (best_hmm_model == 'euk' and hmm_model != 'euk'):
                        # preferentially select the longest sequence, but aim
                        # for a SSU that had the best match to either the archaeal
                        # or bacterial HMM models as we know these genomes are prokaryotic
                        # (sequences with a best hit to a Eukaryotic model will appear as
                        # Euks when aligned with ssu-align!)
                        longest_query_len = query_len
                        best_hmm_model = hmm_model
                        seq_id = line_split[0]
                        seq_info = line.strip()

                if longest_query_len:
                    ssu_seq_file = os.path.join(genome_dirs[genome_id], 'ssu_gg_2013_08', 'ssu.fna')
                    seqs = seq_io.read_fasta(ssu_seq_file)
                    fout.write('>' + genome_id + ' ' + seq_info + '\n')
                    fout.write(seqs[seq_id] + '\n')
        fout.close()

        self.logger.info('Identified %d RefSeq representatives with no 16S sequence.' % no_identified_ssu)

        return ssu_output_file

    def _trim_seqs(self, input_msa, output_msa, remove_identical=False, min_per_taxa=0.7, min_bp=1000):
        """Trim ends of sequences.

        input_msa : str
            File with MSA to trim.
        output_msa : str
            New file with trimmed MSA.
        remove_identical : boolean
            Flag indicating if identical sequence should be removed.
        min_per_taxa : float [0, 1.0]
            Minimum required taxa to retain leading and trailing columns.
        min_bp : int
            Minimum required length to retain sequence.
        """

        # read seqs
        seqs = seq_io.read_fasta(input_msa)

        # filter identical seqs
        identical_seqs = set()
        if remove_identical:
            self.logger.info('Filtering identical sequences.')

            seq_ids = seqs.keys()
            for i in xrange(0, len(seq_ids)):
                seq_id_I = seq_ids[i]

                if seq_id_I in identical_seqs:
                    continue

                for j in xrange(i + 1, len(seqIds)):
                    seq_id_J = seq_ids[j]
                    if seqs[seq_id_I] == seqs[seq_id_J]:
                      self.logger.info('Seq %s and %s are identical.' % (seq_id_I, seq_id_J))
                      identical_seqs.add(seq_id_J)

            self.logger.info('Identified %d of %d sequences as identical.' % (len(identical_seqs), len(seqs)))

        # trim start and end columns to consensus alignment
        first_char = []
        last_char = []
        for seq_id, seq in seqs.iteritems():
            if seq_id in identical_seqs:
                continue

            for i, ch in enumerate(seq):
                if ch != '.' and ch != '-':
                    first_char.append(i)
                    break

            for i in xrange(len(seq) - 1, -1, -1):
                if seq[i] != '.' and seq[i] != '-':
                    last_char.append(i)
                    break

        first_char.sort()
        last_char.sort(reverse=True)

        trim_index = int((len(first_char) * min_per_taxa) + 1)

        start = first_char[trim_index]
        end = last_char[trim_index]

        self.logger.info('Trimming seqs from %d to %d leaving a %dbp length alignment.' % (start, end, end - start + 1))

        short_seq_file = output_msa + '.short'
        fout = open(output_msa, 'w')
        fout_short = open(short_seq_file, 'w')
        num_filtered_seq = 0
        for seq_id, seq in seqs.iteritems():
            if seq_id in identical_seqs:
                continue

            valid_bp = 0
            for i in xrange(start, min(len(seq), end + 1)):
                ch = seq[i]
                if ch != '.' and ch != '-':
                    valid_bp += 1

            if valid_bp >= min_bp:
                fout.write('>' + seq_id + '\n')
                fout.write(seq[start:end + 1] + '\n')
            else:
                self.logger.info('Filtering seq %s with %d of %d (%.1f%%) aligned bases.' % (seq_id, valid_bp, (end - start + 1), valid_bp * 100.0 / (end - start + 1)))
                num_filtered_seq += 1
                fout_short.write('>' + seq_id + '\n')
                fout_short.write(seq[start:end + 1] + '\n')

        fout.close()
        fout_short.close()

        self.logger.info('Filtered %d of %d sequences due to length.' % (num_filtered_seq, len(seqs) - len(identical_seqs)))
        self.logger.info('Short sequence written to: %s' % short_seq_file)

    def ncbi_rep_tree(self, min_rep_quality, output_dir):
        """Infer 16S tree for NCBI representative and reference genomes.

        Parameters
        ----------
        min_rep_quality : float [0, 100]
            Minimum genome quality for a genome to be a representative.
        output_dir : str
            Directory to store results
        """

        accession_to_taxid, _complete_genomes, representative_genomes = ncbi.read_refseq_metadata(self.gtdb_metadata_file , keep_db_prefix=True)
        self.logger.info('Identified %d RefSeq genomes.' % len(accession_to_taxid))
        self.logger.info('Identified %d representative or reference genomes.' % len(representative_genomes))

        # access genome quality
        genome_quality = read_gtdb_genome_quality(self.gtdb_metadata_file)
        missing_quality = set(accession_to_taxid.keys()) - set(genome_quality.keys())
        if missing_quality:
            self.logger.error('There are %d genomes without metadata information.' % len(missing_quality))
            sys.exit()

        new_genomes_to_consider = []
        poor_quality_ref_genomes = 0
        for genome_id in list(representative_genomes):
            comp, cont, qual = genome_quality.get(genome_id, [-1, -1, -1])
            if qual >= min_rep_quality:
                new_genomes_to_consider.append(genome_id)
            else:
                # check if genome is marked as a representative at NCBI
                if genome_id in representative_genomes:
                    poor_quality_ref_genomes += 1
                    self.logger.warning('Filtered RefSeq representative %s with comp = %.2f, cont = %.2f' % (genome_id, comp, cont))

        genomes_to_consider = new_genomes_to_consider
        self.logger.info('Filtered %d poor quality RefSeq representative genomes.' % poor_quality_ref_genomes)
        self.logger.info('Considering %d genomes after filtering at a quality of %.2f.' % (len(genomes_to_consider), min_rep_quality))

        # get SSU sequences for genomes
        ssu_output_file = self._get_ssu_seqs(genomes_to_consider, output_dir)

        # align sequences
        ssu_align_dir = os.path.join(output_dir, 'ssu_align')
        os.system('ssu-align --dna %s %s' % (ssu_output_file, ssu_align_dir))
        os.system('ssu-mask --afa %s' % ssu_align_dir)

        # trim sequences
        for domain in ['archaea', 'bacteria']:
            input_msa = os.path.join(ssu_align_dir, 'ssu_align.' + domain + '.mask.afa')
            trimmed_msa = os.path.join(output_dir, domain + '.trimmed.fna')
            self._trim_seqs(input_msa, trimmed_msa)

            # infer tree
            output_tree = os.path.join(output_dir, domain + '.ssu_refseq_reps.tree')
            os.system('FastTreeMP -nosupport -nt -gamma %s > %s' % (trimmed_msa, output_tree))
