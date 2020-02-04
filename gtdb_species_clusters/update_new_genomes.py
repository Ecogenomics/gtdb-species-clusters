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
import argparse
import logging
from collections import defaultdict

from gtdb_species_clusters.genome_utils import canonical_gid


class NewGenomes(object):
    """Identify new or modified genomes."""

    def __init__(self, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
    def same_genome_accn(self, accn1, accn2, identical_accns):
        """Check if NCBI genome accessions are the same."""
        
        if accn1 == accn2:
            return True
        
        if canonical_gid(accn1) != canonical_gid(accn2):
            self.logger.error('Genomes have different canonical genome IDs: {}, {}, {}, {}'.format(
                                accn1,
                                canonical_gid(accn1),
                                accn2,
                                canonical_gid(accn2)))
            sys.exit(-1)
        
        accn1 = accn1.replace('RS_', '').replace('GB_', '')
        accn2 = accn2.replace('RS_', '').replace('GB_', '')

        if identical_accns.get(accn1, None) == accn2:
            return True

        return False
            
    def run(self, 
                prev_gtdb_metadata_file,
                cur_gtdb_metadata_file,
                cur_genome_paths,
                ncbi_assembly_summary_genbank):
        """Identify new or modified genomes."""
        
        self.logger.info('Reading previous GTDB genomes.')
        prev_accns = {}
        gtdb_taxonomy = {}
        gtdb_rep = {}
        ncbi_genome_category = {}
        with open(prev_gtdb_metadata_file, encoding='utf-8') as f:
            header = f.readline().strip().split('\t')
            
            gtdb_index = header.index('gtdb_taxonomy')
            gtdb_rep_index = header.index('gtdb_representative')
            ncbi_genome_cat_index = header.index('ncbi_genome_category')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gid = line_split[0]
                if gid.startswith('U'): # only concerned with genomes at NCBI
                    continue
                        
                prev_accns[canonical_gid(gid)] = gid
                gtdb_taxonomy[canonical_gid(gid)] = line_split[gtdb_index]
                gtdb_rep[canonical_gid(gid)] = line_split[gtdb_rep_index]
                ncbi_genome_category[canonical_gid(gid)] = line_split[ncbi_genome_cat_index]
                
        self.logger.info(f' ... identified {len(prev_accns):,} genomes.')
                            
        # get genomes in current release
        self.logger.info('Reading current GTDB genomes.')
        cur_accns = {}
        with open(cur_gtdb_metadata_file, encoding='utf-8') as f:
                f.readline()
                for line in f:
                    line_split = line.strip().split('\t')
                    
                    gid = line_split[0]
                    if gid.startswith('U'): # only concerned with genomes at NCBI
                        continue
                        
                    cur_accns[canonical_gid(gid)] = gid
        self.logger.info(f' ... identified {len(cur_accns):,} genomes.')
        
        # get equivalent GenBank and RefSeq genome assemblies
        self.logger.info('Determining identical GenBank and RefSeq accessions.')
        identical_accns = {}
        with open(ncbi_assembly_summary_genbank, encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    if 'assembly_accession' in line:
                        header = line.strip().split('\t')
                        
                        gb_accn_index = header.index('# assembly_accession')
                        rs_accn_index = header.index('gbrs_paired_asm')
                        paired_asm_index = header.index('paired_asm_comp')
                        
                        for line in f:
                            line_split = line.strip().split('\t')
                            
                            paired_asm = line_split[paired_asm_index]
                            if paired_asm.lower() == 'identical':
                                gb_accn = line_split[gb_accn_index]
                                rs_accn = line_split[rs_accn_index]
                                identical_accns[gb_accn] = rs_accn
                                identical_accns[rs_accn] = gb_accn
        
        # identify new and modified genome IDs
        self.logger.info('Identifying new or modified genome IDs.')
        new_gids = set()
        updated_gids = set()
        for cur_gid in cur_accns:
            if cur_gid in prev_accns:
                if not self.same_genome_accn(cur_accns[cur_gid], 
                                                prev_accns[cur_gid], 
                                                identical_accns):
                    updated_gids.add(cur_gid)
            else:
                # genome not present in previous GTDB release
                new_gids.add(cur_gid)
            
        lost_gids = set(prev_accns) - set(cur_accns)
        num_lost_gtdb_reps = sum([1 for gid in lost_gids if gtdb_rep[gid] == 't'])
        self.logger.info(f' ... identified {len(new_gids):,} new, {len(updated_gids):,} updated, and {len(lost_gids):,} lost genomes.')
        self.logger.info(f' ... {num_lost_gtdb_reps:,} lost genomes were GTDB representatives.')

        # get path to current GTDB genome directories
        self.logger.info('Identifying path to genomic files for current GTDB genomes.')
        cur_genome_files = {}
        skipped_genomes = 0
        with open(cur_genome_paths) as f:
            for line in f:
                line_split = line.strip().split('\t')
                accn = line_split[0]
                genome_path = line_split[1]
                gid = line_split[2]
                assert(canonical_gid(accn) == gid)
                
                if gid not in cur_accns:
                    # a genome may not be part of the GTDB release
                    # (e.g., genome has no NCBI taxonomy information)
                    skipped_genomes += 1
                    continue

                assembly_id = os.path.basename(os.path.normpath(genome_path))
                genomic_file = os.path.join(genome_path, assembly_id + '_genomic.fna')
                cur_genome_files[gid] = genomic_file
                #***if not os.path.exists(genomic_file):
                #***    self.logger.warning('Genomic file not found: {}'.format(genomic_file))
                
        self.logger.info(f' ... identified genomic file for {len(cur_genome_files):,} genomes.')
        self.logger.info(f' ... skipped {skipped_genomes:,} genomes without GTDB metadata.')
        
        # write out new or modified genome IDs
        self.logger.info('Writing out and verifying path to new and updated genomic FASTA files.')
        output_file = os.path.join(self.output_dir, 'genomes_new_updated.tsv')
        fout = open(output_file, 'w')
        fout.write('Genome ID\tNCBI accession\tStatus\tGenomic file\n')
        for type_str, gids in [('NEW', new_gids), ('UPDATED', updated_gids)]:
            for gid in gids:
                genomic_file = cur_genome_files[gid]
                fout.write('{}\t{}\t{}\t{}\n'.format(gid, 
                                                    cur_accns[gid], 
                                                    type_str, 
                                                    genomic_file))
        fout.close()
        
        # write out lost genomes
        output_file = os.path.join(self.output_dir, 'genomes_lost.tsv')
        fout = open(output_file, 'w')
        fout.write('Genome ID\tGTDB taxonomy\tGTDB representative\tNCBI genome category\n')
        lost_reps = 0
        for gid in lost_gids:
            fout.write('{}\t{}\t{}\t{}\n'.format(
                        gid,
                        gtdb_taxonomy[gid],
                        gtdb_rep[gid],
                        ncbi_genome_category[gid]))
        fout.close()
        
