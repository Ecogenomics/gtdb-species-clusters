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
from biolib.common import remove_extension
from biolib.external.fasttree import FastTree
from biolib.external.blast import Blast
from biolib.taxonomy import Taxonomy

import genometreetk.ncbi as ncbi
from genometreetk.common import (read_gtdb_metadata,
                                    read_genome_dir_file,
                                    read_gtdb_taxonomy)


class RNA_Workflow(object):
    """Infer RNA gene trees."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Maximum number of CPUs to use.
        """

        self.logger = logging.getLogger()
        self.cpus = cpus

    def _get_rna_seqs(self,
                        rna_name,
                        rna_file,
                        min_rna_length,
                        min_scaffold_length,
                        gtdb_taxonomy,
                        genomes_to_consider,
                        output_dir):
        """Get 16S sequences from genomes.
        
        Assumes the FASTA header has the format:
            <genome_id>~<seq_id> <gtdb_taxonomy> <rna_gene_length> <contig_length>

        Parameters
        ----------
        rna_name : str
            Name of rRNA gene.
        rna_file : str
            File with rRNA gene sequences in FASTA format.
        min_rna_length : int
            Minimum required length of rRNA gene sequences.
        min_scaffold_length : int
            Minimum required length of scaffold containing rRNE gene sequence.
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

        rna_output_file = os.path.join(output_dir, '%s.fna' % rna_name)
        fout = open(rna_output_file, 'w')
        filtered_rna_len = 0
        filtered_scaffold_len = 0
        genomes_with_seq = set()
        for seq_id, seq, annotation in seq_io.read_seq(rna_file, keep_annotation=True):
            genome_id, contig_id = seq_id.split('~', 1)
            if genome_id not in genomes_to_consider:
                continue
                
            if len(seq) < min_rna_length:
                filtered_rna_len += 1
                continue
               
            scaffold_len = int(annotation.split(' ')[2].strip())
            if scaffold_len < min_scaffold_length:
                filtered_scaffold_len += 1
                continue
                
            genomes_with_seq.add(genome_id)
            
            fout.write('>%s %s\n' % (genome_id, annotation))
            fout.write(seq + '\n')
            
        fout.close()

        self.logger.info('Filtered %d rRNA gene sequence based on rRNA gene length <%d.' % (filtered_rna_len, min_rna_length))
        self.logger.info('Filtered %d rRNA gene sequence based on scaffold length <%d.' % (filtered_scaffold_len, min_scaffold_length))
        self.logger.info('Identified %d genomes with a valid rRNA gene sequence.' % len(genomes_with_seq))

        return rna_output_file

    def _trim_seqs(self, input_msa, output_msa, remove_identical=False, min_per_taxa=0.5, min_bp=1000):
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

        trim_index = int((len(seqs) * min_per_taxa))

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
        
        extant_taxa = Taxonomy().extant_taxa(taxonomy)

        tax_filter_dir = os.path.join(output_dir, 'tax_filter')
        if not os.path.exists(tax_filter_dir):
            os.makedirs(tax_filter_dir)
        
        blast = Blast(self.cpus)
        self.logger.info('Creating BLASTN database.')
        blast.create_blastn_db(ssu_output_file)

        self.logger.info('Performing reciprocal BLAST to identify sequences with incongruent taxonomies.')
        
        blast_table = os.path.join(tax_filter_dir, 'blastn.tsv')
        blast.blastn(ssu_output_file,
                     ssu_output_file,
                     blast_table,
                     evalue=1e-10,
                     max_matches=2,
                     output_fmt='custom',
                     task='blastn')

        filter = set()
        fout = open(os.path.join(tax_filter_dir, 'filtered_seqs.tsv'), 'w')
        fout.write('Query ID\tQuery Taxonomy\tSubject ID\tSubject Taxonomy\tPerc. Identity\tAlign. Length\tMismatch Rank\tNo. Query Genomes\tNo. Subject Genomes\n')
        for hit in blast.read_hit(blast_table, table_fmt='custom'):
            if hit.query_id == hit.subject_id:
                # ignore self hits
                continue

            if hit.alignment_len > 800:
                query_genome_id = hit.query_id.split('~', 1)[0]
                subject_genome_id = hit.subject_id.split('~', 1)[0]
                    
                # require a (very lenient) percent identity of threshold from Yarza et al., 2014
                
                if hit.perc_identity >= 75: # phylum
                    rank_index = Taxonomy.rank_labels.index('phylum')
                    query_taxa = taxonomy[query_genome_id][rank_index][3:].strip()
                    subject_taxa = taxonomy[subject_genome_id][rank_index][3:].strip()
                    if query_taxa and subject_taxa and query_taxa != subject_taxa:
                        filter.add(hit.query_id)
                        fout.write('%s\t%s\t%s\t%s\t%.2f\t%d\t%s\t%d\t%d\n' % (hit.query_id, 
                                                                ';'.join(taxonomy[query_genome_id]),
                                                                hit.subject_id,
                                                                ';'.join(taxonomy[subject_genome_id]),
                                                                hit.perc_identity,
                                                                hit.alignment_len,
                                                                'Phylum',
                                                                len(extant_taxa['p__' + query_taxa]),
                                                                len(extant_taxa['p__' + subject_taxa])))
                                                                
                if hit.perc_identity >= 78.5: # class
                    rank_index = Taxonomy.rank_labels.index('class')
                    query_taxa = taxonomy[query_genome_id][rank_index][3:].strip()
                    subject_taxa = taxonomy[subject_genome_id][rank_index][3:].strip()
                    if query_taxa and subject_taxa and query_taxa != subject_taxa:
                        filter.add(hit.query_id)
                        fout.write('%s\t%s\t%s\t%s\t%.2f\t%d\t%s\t%d\t%d\n' % (hit.query_id, 
                                                                ';'.join(taxonomy[query_genome_id]),
                                                                hit.subject_id,
                                                                ';'.join(taxonomy[subject_genome_id]),
                                                                hit.perc_identity,
                                                                hit.alignment_len,
                                                                'Class',
                                                                len(extant_taxa['c__' + query_taxa]),
                                                                len(extant_taxa['c__' + subject_taxa])))
                
                if hit.query_id not in filter and hit.perc_identity >= 82: # order
                    rank_index = Taxonomy.rank_labels.index('order')
                    query_taxa = taxonomy[query_genome_id][rank_index][3:].strip()
                    subject_taxa = taxonomy[subject_genome_id][rank_index][3:].strip()
                    if query_taxa and subject_taxa and query_taxa != subject_taxa:
                        filter.add(hit.query_id)
                        fout.write('%s\t%s\t%s\t%s\t%.2f\t%d\t%s\t%d\t%d\n' % (hit.query_id, 
                                                                ';'.join(taxonomy[query_genome_id]),
                                                                hit.subject_id,
                                                                ';'.join(taxonomy[subject_genome_id]),
                                                                hit.perc_identity,
                                                                hit.alignment_len,
                                                                'Order',
                                                                len(extant_taxa['o__' + query_taxa]),
                                                                len(extant_taxa['o__' + subject_taxa])))
                                                                
                if hit.query_id not in filter and hit.perc_identity >= 86.5: # family
                    rank_index = Taxonomy.rank_labels.index('family')
                    query_taxa = taxonomy[query_genome_id][rank_index][3:].strip()
                    subject_taxa = taxonomy[subject_genome_id][rank_index][3:].strip()
                    if query_taxa and subject_taxa and query_taxa != subject_taxa:
                        filter.add(hit.query_id)
                        fout.write('%s\t%s\t%s\t%s\t%.2f\t%d\t%s\t%d\t%d\n' % (hit.query_id, 
                                                                ';'.join(taxonomy[query_genome_id]),
                                                                hit.subject_id,
                                                                ';'.join(taxonomy[subject_genome_id]),
                                                                hit.perc_identity,
                                                                hit.alignment_len,
                                                                'Family',
                                                                len(extant_taxa['f__' + query_taxa]),
                                                                len(extant_taxa['f__' + subject_taxa])))
                                                                
                if hit.query_id not in filter and hit.perc_identity >= 94.5: # genus
                    rank_index = Taxonomy.rank_labels.index('genus')
                    query_taxa = taxonomy[query_genome_id][rank_index][3:].strip()
                    subject_taxa = taxonomy[subject_genome_id][rank_index][3:].strip()
                    if query_taxa and subject_taxa and query_taxa != subject_taxa:
                        filter.add(hit.query_id)
                        fout.write('%s\t%s\t%s\t%s\t%.2f\t%d\t%s\t%d\t%d\n' % (hit.query_id, 
                                                                ';'.join(taxonomy[query_genome_id]),
                                                                hit.subject_id,
                                                                ';'.join(taxonomy[subject_genome_id]),
                                                                hit.perc_identity,
                                                                hit.alignment_len,
                                                                'Genus',
                                                                len(extant_taxa['g__' + query_taxa]),
                                                                len(extant_taxa['g__' + subject_taxa])))

                if hit.query_id not in filter and hit.perc_identity >= 99: # species
                    rank_index = Taxonomy.rank_labels.index('species')
                    query_taxa = taxonomy[query_genome_id][rank_index][3:].strip()
                    subject_taxa = taxonomy[subject_genome_id][rank_index][3:].strip()
                    if query_taxa and subject_taxa and query_taxa != subject_taxa:
                        filter.add(hit.query_id)
                        fout.write('%s\t%s\t%s\t%s\t%.2f\t%d\t%s\t%d\t%d\n' % (hit.query_id, 
                                                                ';'.join(taxonomy[query_genome_id]),
                                                                hit.subject_id,
                                                                ';'.join(taxonomy[subject_genome_id]),
                                                                hit.perc_identity,
                                                                hit.alignment_len,
                                                                'Species',
                                                                len(extant_taxa['s__' + query_taxa]),
                                                                len(extant_taxa['s__' + subject_taxa])))

        fout.close()
        
        return filter

    def run(self, rna_name,
                    gtdb_metadata_file, 
                    rna_file, 
                    min_rna_length,
                    min_scaffold_length,
                    min_quality,
                    max_contigs,
                    min_N50,
                    tax_filter,
                    genome_list,
                    output_dir,
                    align_method='ssu_align'):
        """Infer rRNA gene tree spanning select GTDB genomes.

        Parameters
        ----------
        rna_name : str
            Name of rRNA gene.
        gtdb_metadata_file : str
            File specifying GTDB metadata for each genome.
        rna_file : str
            File with rRNA gene sequences in FASTA format.
        min_rna_length : int
            Minimum required length of rRNA gene sequences.
        min_scaffold_length : int
            Minimum required length of scaffold containing rRNA gene sequence.
        min_quality : float [0, 100]
            Minimum genome quality for a genome to be include in tree.
        max_contigs : int
            Maximum number of contigs to include genome.
        min_N50 : int
            Minimum N50 to include genome.
        tax_filter : boolean
            Filter sequences based on incongruent taxonomy classification.
        genome_list : str
            Explict list of genomes to use (ignores --ncbi_rep_only and --user_genomes).
        output_dir : str
            Directory to store results
        """
                
        if rna_name not in ['ssu', 'lsu']:
            self.logger.error('Unrecognized rRNA gene type: %s' % rna_name)
            sys.exit(-1)

        genome_metadata = read_gtdb_metadata(gtdb_metadata_file, ['checkm_completeness',
                                                                      'checkm_contamination',
                                                                      'scaffold_count',
                                                                      'n50_scaffolds',
                                                                      'organism_name',
                                                                      'gtdb_representative'])

        gtdb_taxonomy = read_gtdb_taxonomy(gtdb_metadata_file)

        user_genomes = set()
        uba_genomes = set()
        ncbi_genomes = set()
        rep_genomes = set()
        for genome_id in genome_metadata:
            org_name = str(genome_metadata[genome_id][4])
            if genome_id.startswith('U_'):
                if '(UBA' in org_name:
                    uba_genomes.add(genome_id)
                else:
                    user_genomes.add(genome_id)
            elif genome_id.startswith('RS_') or genome_id.startswith('GB_'):
                ncbi_genomes.add(genome_id)
            else:
                self.logger.warning('Unrecognized genome prefix: %s' % genome_id)
                
            rep = genome_metadata[genome_id][5] == 't'
            if rep:
                rep_genomes.add(genome_id)
                
        self.logger.info('Initially considering %d genomes (%d NCBI, %d UBA, %d User).' % (len(genome_metadata), 
                                                                                            len(ncbi_genomes), 
                                                                                            len(uba_genomes), 
                                                                                            len(user_genomes)))
        self.logger.info('Identified %d representative genomes.' % len(rep_genomes))
        
        # get genomes specified in genome list by user
        genomes_to_consider = set()
        if genome_list:
            for line in open(genome_list):
                genomes_to_consider.add(line.rstrip().split('\t')[0])
            self.logger.info('Restricting genomes to the %d in the genome list.' % len(genomes_to_consider))
        else:
            # filter genomes based on quality and database source
            self.logger.info('Filtering genomes based on specified critieria.')
            self.logger.info('Filtering on minimum quality <%d.' % min_quality)
            self.logger.info('Filtering on number of contigs >%d.' % max_contigs)
            self.logger.info('Filtering on scaffold N50 <%d.' % min_N50)
            
            new_genomes_to_consider = []
            filtered_genomes = 0
            gt = 0
            gq = 0
            sc = 0
            n50 = 0
            for genome_id in genome_metadata:                
                if genome_id not in rep_genomes:
                    gt += 1
                    filtered_genomes += 1
                    continue
                    
                if genome_id not in ncbi_genomes and genome_id not in uba_genomes:
                    gt += 1
                    filtered_genomes += 1
                    continue
                    
                comp, cont, scaffold_count, n50_contigs, _org_name, _rep = genome_metadata[genome_id]
                q = float(comp) - 5*float(cont)
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
            self.logger.info('Filtered %d genomes (%d on genome type, %d on genome quality, %d on number of contigs, %d on N50).' % (filtered_genomes, gt, gq, sc, n50))
            self.logger.info('Considering %d genomes after filtering.' % len(genomes_to_consider))
            
        # limit taxonomy to genomes being considered
        cur_gtdb_taxonomy = {}
        for gid in genomes_to_consider:
            cur_gtdb_taxonomy[gid] = gtdb_taxonomy[gid]

        # get rRNA gene sequences for each genome
        rna_output_file = self._get_rna_seqs(rna_name,
                                                rna_file,
                                                min_rna_length,
                                                min_scaffold_length,
                                                cur_gtdb_taxonomy,
                                                genomes_to_consider,
                                                output_dir)

        # identify erroneous rRNA gene sequences
        if tax_filter:
            self.logger.info('Filtering sequences with incongruent taxonomy strings.')
            filter = self._tax_filter(rna_output_file, cur_gtdb_taxonomy, output_dir)

            self.logger.info('Filtered %d sequences.' % len(filter))
            if len(filter) > 0:
                rna_filtered_output = os.path.join(output_dir, 'gtdb_%s.tax_filter.fna' % rna_name)
                fout = open(rna_filtered_output, 'w')
                for seq_id, seq, annotation in seq_io.read_seq(rna_output_file, keep_annotation=True):
                    if seq_id not in filter:
                        fout.write('>' + seq_id + ' ' + annotation + '\n')
                        fout.write(seq + '\n')
                fout.close()

                rna_output_file = rna_filtered_output

        # align sequences with ssu-align or mothur
        if rna_name == 'ssu':
            if align_method == 'ssu_align':
                self.logger.info('Aligning sequences with ssu-align.')
                align_dir = os.path.join(output_dir, '%s_align' % rna_name)
                os.system('ssu-align --dna %s %s' % (rna_output_file, align_dir))
                os.system('ssu-mask --afa %s' % align_dir)
            elif align_method == 'mothur':
                self.logger.info('Aligning sequences with mothur.')
                align_dir = os.path.join(output_dir, 'mothur')
                if not os.path.exists(align_dir):
                    os.makedirs(align_dir)

                mothur_cmd = 'mothur "#set.dir(output=%s, blastdir=/srv/sw/Mothur/1.39.5)' % align_dir
                mothur_cmd += '; align.seqs(candidate=%s, template=/srv/db/mothur/silva_128/silva.seed_v128.align, search=blast, flip=t, processors=%d)' % (rna_output_file, self.cpus)
                input_prefix = remove_extension(rna_output_file)
                align_file = os.path.join(align_dir, input_prefix + '.align')
                mothur_cmd += '; filter.seqs(fasta=%s, hard=/srv/db/mothur/silva_128/Lane1349.silva.filter, processors=%d);"' % (align_file, self.cpus)
                os.system(mothur_cmd)
                input_msa = os.path.join(align_dir, input_prefix + '.filter.fasta')
        elif rna_name == 'lsu':
            self.logger.info('Aligning sequences with ssu-align.')
            align_dir = os.path.join(output_dir, '%s_align' % rna_name)
            if not os.path.exists(align_dir):
                os.makedirs(align_dir)
                
            os.system('esl-sfetch --index %s' % rna_output_file)
                  
            # search fo sequences using domain-specific LSU HMMs
            for domain in ['archaea', 'bacteria', 'eukaryote']:
                self.logger.info('Matching LSU rRNA genes to %s-specific HMM.' % domain)
                table_out = os.path.join(align_dir, 'cmsearch.%s.%s.tblout' % (rna_name, domain))
                cm_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'cm_files')            
                cm_file = os.path.join(cm_dir, 'lsu_%s.cm' % domain)
                log_file = os.path.join(align_dir, 'cmsearch.%s.%s.out' % (rna_name, domain))
                os.system('cmsearch --hmmonly --cpu %d --noali --tblout %s %s %s > %s' % (self.cpus, table_out, cm_file, rna_output_file, log_file))
                
            # identify top hits for each domain
            self.logger.info('Identifying best domain-specific HMM for each LSU rRNA gene.')
            top_hits = {}
            for domain in ['archaea', 'bacteria', 'eukaryote']:
                table_out = os.path.join(align_dir, 'cmsearch.%s.%s.tblout' % (rna_name, domain))
                for line in open(table_out):
                    if line[0] == '#':
                        continue
                        
                    line_split = line.split()
                    seq_id = line_split[0]
                    start_seq = int(line_split[7])
                    end_seq = int(line_split[8])
                    bitscore = float(line_split[14])

                    prev_bitscore = top_hits.get(seq_id, [None, 0, 0, 0, 0])[4]
                    if bitscore > prev_bitscore:
                        top_hits[seq_id] = [domain, seq_id, start_seq, end_seq, bitscore]
              
            # create MSA for each bacteria and archaea
            for domain in ['archaea', 'bacteria']:
                # creat file of top hits
                top_hits_out = os.path.join(align_dir, 'top_hits.%s.%s.tsv' % (rna_name, domain))
                fout = open(top_hits_out, 'w')
                num_hits = 0
                for top_domain, seq_id, start_seq, end_seq, bitscore in top_hits.values():
                    if top_domain == domain:
                        fout.write('%s\t%d\t%d\%f\n' % (seq_id, start_seq, end_seq, bitscore))
                        num_hits += 1
                fout.close()
                
                # align top hits
                self.logger.info('Creating MSA for %s LSU rRNA genes (%d sequences).' % (domain, num_hits))
                
                if num_hits > 0:
                    seq_file = os.path.join(align_dir, 'cmsearch.%s.%s.fna' % (rna_name, domain))
                    os.system("grep -v '^#' %s | awk '{print $1, $2, $3, $1}' | esl-sfetch -Cf %s - > %s" % (top_hits_out, rna_output_file, seq_file))
                    
                    align_file = os.path.join(align_dir, 'cmalign.%s.%s.stk' % (rna_name, domain))
                    os.system('cmalign --dnaout --outformat Pfam %s %s > %s' % (cm_file, seq_file, align_file))
                    
                    masked_file = os.path.join(align_dir, 'cmalign.%s.%s.mask.afa' % (rna_name, domain))
                    os.system('esl-alimask -p --outformat AFA %s > %s' % (align_file, masked_file))
            
        # trim sequences and infer tree
        if align_method == 'ssu_align':
            for domain in ['archaea', 'bacteria']:
                if rna_name == 'ssu':
                    input_msa = os.path.join(align_dir, 'ssu_align.' + domain + '.mask.afa')
                elif rna_name == 'lsu':
                    input_msa = os.path.join(align_dir, 'cmalign.%s.%s.mask.afa' % (rna_name, domain))
                    
                if not os.path.exists(input_msa):
                    continue
                    
                trimmed_msa = os.path.join(output_dir, domain + '.trimmed.fna')
                self._trim_seqs(input_msa, trimmed_msa)

                # infer tree
                self.logger.info('Inferring tree for %s genes.' % domain)
                output_tree = os.path.join(output_dir, domain + '.tree')
                os.system('FastTreeMP -nosupport -nt -gamma %s > %s' % (trimmed_msa, output_tree))
        elif align_method == 'mothur':
            trimmed_msa = os.path.join(output_dir, input_prefix + '.trimmed.fna')
            self._trim_seqs(input_msa, trimmed_msa)

            # infer tree
            self.logger.info('Inferring tree for %s genes.')
            output_tree = os.path.join(output_dir, input_prefix + '.tree')
            os.system('FastTreeMP -nosupport -nt -gamma %s > %s' % (trimmed_msa, output_tree))
            
    def combine(self, ssu_msa, ssu_tree, lsu_msa, lsu_tree, output_dir):
        """Infer 16S + 23S tree spanning GTDB genomes."""
        
        # identify common 16S and 23S sequences
        ssu_seqs = {}
        for seq_id, seq, annotation in seq_io.read_seq(ssu_msa, keep_annotation=True):
            genome_id = seq_id.split('~')[0]
            ssu_seqs[genome_id] = [seq, annotation]
        self.logger.info('Read %d SSU rRNA sequences.' % len(ssu_seqs))
            
        lsu_seqs = {}
        for seq_id, seq, annotation in seq_io.read_seq(lsu_msa, keep_annotation=True):
            genome_id = seq_id.split('~')[0]
            lsu_seqs[genome_id] = [seq, annotation]
        self.logger.info('Read %d LSU rRNA sequences.' % len(lsu_seqs))
              
        common_seqs = set(ssu_seqs.keys()).intersection(lsu_seqs.keys())
        self.logger.info('Identified %d sequences in common.' % len(common_seqs))
        
        # identify incongruent taxonomic order classifcations between trees
        self.logger.info('Identifying incongruent order-level taxonomic classifications between trees.')
        ssu_taxonomy = Taxonomy().read_from_tree(ssu_tree)
        lsu_taxonomy = Taxonomy().read_from_tree(lsu_tree)

        order_index = Taxonomy.rank_labels.index('order')

        seqs_to_filter = set()
        for seq_id in common_seqs:
            ssu_order = ssu_taxonomy.get(seq_id)[order_index][3:]
            lsu_order = lsu_taxonomy.get(seq_id)[order_index][3:]
            
            # remove designator of paraphyletic orders
            # (since in the concatenated tree this may be resolved)
            ssu_order = ssu_order.split('_')[0]
            lsu_order = lsu_order.split('_')[0]
            
            if ssu_order != lsu_order:
                seqs_to_filter.add(seq_id)
        
        self.logger.info('Identified %d sequences with incongruent classifcations.' % len(seqs_to_filter))
        common_seqs.difference_update(seqs_to_filter)
        
        # write out MSA
        concatenated_msa = os.path.join(output_dir, 'ssu_lsu_concatenated.fna')
        fout = open(concatenated_msa, 'w')
        for seq_id in common_seqs:
            fout.write('>%s %s %s\n' % (seq_id, ssu_seqs[seq_id][1], lsu_seqs[seq_id][1]))
            fout.write('%s%s\n' % (ssu_seqs[seq_id][0], lsu_seqs[seq_id][0]))
        fout.close()
        
        # infer tree
        output_tree = os.path.join(output_dir, 'ssu_lsu_concatenated.tree')
        os.system('FastTreeMP -nosupport -nt -gamma %s > %s' % (concatenated_msa, output_tree))
