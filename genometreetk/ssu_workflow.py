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
from biolib.external.blast import Blast
from biolib.taxonomy import Taxonomy

import genometreetk.ncbi as ncbi
from genometreetk.common import (read_gtdb_metadata,
                                    read_genome_dir_file,
                                    read_gtdb_taxonomy)


class SSU_Workflow(object):
    """Infer SSU trees."""

    def __init__(self, gtdb_metadata_file, genome_dir_file, cpus):
        """Initialization.

        Parameters
        ----------
        gtdb_metadata_file : str
            File specifying GTDB metadata for each genome.
        genome_dir_file : str
            File specifying directory for each genome.
        cpus : int
            Maximum number of CPUs to use.
        """

        self.logger = logging.getLogger()

        self.gtdb_metadata_file = gtdb_metadata_file
        self.genome_dir_file = genome_dir_file
        self.cpus = cpus

    def _get_ssu_seqs(self,
                      min_ssu_length,
                      min_ssu_contig,
                      gtdb_taxonomy,
                      genomes_to_consider,
                      output_dir):
        """Get 16S sequences from genomes.

        Parameters
        ----------
        min_ssu_length : int
            Minimum required length of 16S sequences.
        min_ssu_contig : int
            Minimum required length of contig containing 16S sequence.
        gtdb_taxonomy : list
            GTDB taxonomy for each genome assembly.
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

        ssu_output_file = os.path.join(output_dir, 'gtdb_ssu.fna')

        fout = open(ssu_output_file, 'w')
        no_identified_ssu = 0
        for genome_id in genomes_to_consider:
            ssu_hmm_file = os.path.join(genome_dirs[genome_id], 'ssu_gg_2013_08', 'ssu.hmm_summary.tsv')

            if not os.path.exists(ssu_hmm_file):
                no_identified_ssu += 1
                # self.logger.warning('Genome %s has no identified 16S sequence.' % genome_id)
                continue


            with open(ssu_hmm_file) as f:
                headers = f.readline().rstrip().split('\t')
                hmm_model_index = headers.index('HMM')
                ssu_length_index = headers.index('SSU gene length')
                contig_length_index = headers.index('Sequence length')

                longest_contig_len = 0
                best_hmm_model = None
                seq_info = None
                seq_id = None
                for line in f:
                    line_split = line.rstrip().split('\t')

                    ssu_query_len = int(line_split[ssu_length_index])
                    if ssu_query_len < min_ssu_length:
                        continue

                    hmm_model = line_split[hmm_model_index]
                    contig_length = int(line_split[contig_length_index])
                    if contig_length < min_ssu_contig:
                        continue

                    if contig_length > longest_contig_len or (best_hmm_model == 'euk' and hmm_model != 'euk'):
                        # preferentially select the longest contig, but aim
                        # for a SSU that had the best match to either the archaeal
                        # or bacterial HMM models as we know these genomes are prokaryotic
                        # (sequences with a best hit to a Eukaryotic model will appear as
                        # Euks when aligned with ssu-align!)
                        longest_contig_len = contig_length
                        best_hmm_model = hmm_model
                        seq_id = line_split[0]
                        seq_info = line.strip()

                if longest_contig_len:
                    ssu_seq_file = os.path.join(genome_dirs[genome_id], 'ssu_gg_2013_08', 'ssu.fna')
                    seqs = seq_io.read_fasta(ssu_seq_file)
                    fout.write('>' + genome_id + ' ' + ';'.join(gtdb_taxonomy[genome_id]) + '\n')
                    fout.write(seqs[seq_id] + '\n')
        fout.close()

        self.logger.info('Identified %d genomes with no 16S sequence meeting filtering criteria.' % no_identified_ssu)

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

    def _tax_filter(self, ssu_output_file, taxonomy, output_dir):
        """Identify sequence to filter based on taxonomy of best BLAST hit.

        """

        blast = Blast(self.cpus)
        self.logger.info('Creating BLASTN database.')
        blast.create_blastn_db(ssu_output_file)

        self.logger.info('Performing reciprocal BLAST to identify sequences with incongruent taxonomies.')
        blast_table = os.path.join(output_dir, 'blastn.tsv')
        blast.blastn(ssu_output_file,
                     ssu_output_file,
                     blast_table,
                     evalue=1e-10,
                     max_matches=2,
                     output_fmt='custom',
                     task='megablast')

        filter = set()
        order_index = Taxonomy.rank_labels.index('order')
        for hit in blast.read_hit(blast_table, table_fmt='custom'):
            if hit.query_id == hit.subject_id:
                # ignore self hits
                continue

            if hit.perc_identity >= 97 and hit.alignment_len > 800:
                # there is a close hit in the database so verify it has
                # the expected taxonomic order

                order_of_query = taxonomy[hit.query_id][order_index][3:].strip()
                if order_of_query != Taxonomy.rank_prefixes[order_index]:
                    order_of_subject = taxonomy[hit.subject_id][order_index][3:].strip()
                    if order_of_query and order_of_subject and order_of_query != order_of_subject:
                        filter.add(hit.query_id)

        return filter

    def run(self, min_ssu_length,
                    min_ssu_contig,
                    min_quality,
                    max_contigs,
                    min_N50,
                    tax_filter,
                    ncbi_rep_only,
                    user_genomes,
                    genome_list,
                    output_dir):
        """Infer 16S tree spanning select GTDB genomes.

        Parameters
        ----------
        min_ssu_length : int
            Minimum required length of 16S sequences.
        min_ssu_contig : int
            Minimum required length of contig containing 16S sequence.
        min_quality : float [0, 100]
            Minimum genome quality for a genome to be include in tree.
        max_contigs : int
            Maximum number of contigs to include genome.
        min_N50 : int
            Minimum N50 to include genome.
        tax_filter : boolean
            Filter sequences based on incongruent taxonomy classification.
        ncbi_rep_only : boolean
            Restrict tree to NCBI representative and reference genomes.
        user_genomes : boolean
            Include User genomes in tree.
        genome_list : str
            Explict list of genomes to use (ignores --ncbi_rep_only and --user_genomes).
        output_dir : str
            Directory to store results
        """
        
        if genome_list and (ncbi_rep_only or user_genomes):
            self.logger.error("The 'genome_list' flag cannot be used with the 'ncbi_rep_only' and 'user_genomes' flags.")
            sys.exit(-1)

        genome_quality = read_gtdb_metadata(self.gtdb_metadata_file, ['checkm_completeness',
                                                                      'checkm_contamination',
                                                                      'scaffold_count',
                                                                      'n50_contigs'])

        gtdb_taxonomy = read_gtdb_taxonomy(self.gtdb_metadata_file)

        num_user_genomes = 0
        num_ncbi_genomes = 0
        for genome_id in genome_quality:
            if genome_id.startswith('U_'):
                num_user_genomes += 1
            elif genome_id.startswith('RS_') or genome_id.startswith('GB_'):
                num_ncbi_genomes += 1
            else:
                self.logger.warning('Unrecognized genome prefix: %s' % genome_id)
        self.logger.info('Initially considering %d genomes (%d NCBI, %d User).' % (len(genome_quality), num_ncbi_genomes, num_user_genomes))
        
        # filter genomes based on quality and database source
        genomes_in_list = set()
        if genome_list:
            for line in open(genome_list):
                genomes_in_list.add(line.rstrip().split('\t')[0])
            self.logger.info('Restricted genomes to the %d in the genome list.' % len(genomes_in_list))
        else:
            if ncbi_rep_only:
                _accession_to_taxid, _complete_genomes, ncbi_rep_genomes = ncbi.read_refseq_metadata(self.gtdb_metadata_file , keep_db_prefix=True)
                self.logger.info('Considering only NCBI representative or reference genomes.')
                self.logger.info('Identified %d RefSeq genomes.' % len(accession_to_taxid))
                self.logger.info('Identified %d representative or reference genomes.' % len(ncbi_rep_genomes))
                
            if not user_genomes:
                self.logger.info('Filtering User genomes.')

        self.logger.info('Filtering genomes based on specified critieria.')
        new_genomes_to_consider = []
        filtered_genomes = 0
        gq = 0
        sc = 0
        n50 = 0
        for genome_id in genome_quality:                
            if genomes_in_list: 
                if genome_id not in genomes_in_list:
                    continue
            else:
                if not user_genomes and genome_id.startswith('U_'):
                    filtered_genomes += 1
                    continue

                if ncbi_rep_only and genome_id not in ncbi_rep_genomes:
                    filtered_genomes += 1
                    continue
   
            comp, cont, scaffold_count, n50_contigs = genome_quality.get(genome_id, [-1, -1, -1])
            q = float(comp) - float(cont)
            if q < min_quality or int(scaffold_count) > max_contigs or int(n50_contigs) < min_N50:
                if q < min_quality:
                    gq += 1
                    
                if int(scaffold_count) > max_contigs:
                    sc += 1
                    
                if int(n50_contigs) < min_N50:
                    n50 += 1
            
                filtered_genomes += 1
                continue
                
            new_genomes_to_consider.append(genome_id)

        genomes_to_consider = new_genomes_to_consider
        self.logger.info('Filtered %d genomes (%d on genome quality, %d on number of contigs, %d on N50).' % (filtered_genomes, gq, sc, n50))
        self.logger.info('Considering %d genomes after filtering.' % len(genomes_to_consider))

        # get SSU sequences for genomes
        ssu_output_file = self._get_ssu_seqs(min_ssu_length,
                                             min_ssu_contig,
                                             gtdb_taxonomy,
                                             genomes_to_consider,
                                             output_dir)

        # identify erroneous SSU sequences
        if tax_filter:
            self.logger.info('Filtering sequences with incongruent taxonomy strings.')
            filter = self._tax_filter(ssu_output_file, gtdb_taxonomy, output_dir)

            self.logger.info('Filtered %d sequences.' % len(filter))
            print filter

            if len(filter) > 0:
                ssu_filtered_output = os.path.join(output_dir, 'gtdb_ssu.filtered.fna')
                fout = open(ssu_filtered_output, 'w')
                for seq_id, seq, annotation in seq_io.read_seq(ssu_output_file, keep_annotation=True):
                    if seq_id not in filter:
                        fout.write('>' + seq_id + ' ' + annotation + '\n')
                        fout.write(seq + '\n')
                fout.close()

                ssu_output_file = ssu_filtered_output

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
            output_tree = os.path.join(output_dir, domain + '.tree')
            os.system('FastTreeMP -nosupport -nt -gamma %s > %s' % (trimmed_msa, output_tree))
