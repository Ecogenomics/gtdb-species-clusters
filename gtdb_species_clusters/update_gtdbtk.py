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

from biolib.external.execute import check_dependencies

from gtdb_species_clusters.genome_utils import canonical_gid, read_qc_file


class GTDB_Tk(object):
    """Perform initial classification of new and updated genomes using GTDB-Tk."""

    def __init__(self, cpus, output_dir):
        """Initialization."""
        
        check_dependencies(['gtdbtk'])
        
        self.cpus = cpus
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
            
    def run(self, 
                genomes_new_updated_file,
                qc_passed_file,
                batch_size):
        """Perform initial classification of new and updated genomes using GTDB-Tk."""
        
        # get list of genomes passing QC
        self.logger.info('Reading genomes passing QC.')
        gids_pass_qc = read_qc_file(qc_passed_file)
        self.logger.info(f' ... identified {len(gids_pass_qc):,} genomes.')
        
        # get path to genomes passing QC
        self.logger.info('Reading path to genomic file for new/updated genomes passing QC.')
        genomic_files = []
        total_count = 0
        with open(genomes_new_updated_file, encoding='utf-8') as f:
            header = f.readline().strip().split('\t')
            
            genomic_file_index = header.index('Genomic file')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gid = line_split[0]
                total_count += 1
                if gid in gids_pass_qc:
                    gf = line_split[genomic_file_index]
                    genomic_files.append((gid, gf))
        self.logger.info(f' ... identified {len(genomic_files):,} of {total_count:,} genomes as passing QC.')
        
        # process genomes with GTDB-Tk in batches
        for batch_idx, start in enumerate(range(0, len(genomic_files), batch_size)):
            batch_dir = os.path.join(self.output_dir, 'batch_{}'.format(batch_idx))
            if os.path.exists(batch_dir):
                self.logger.warning(f'Skipping {batch_dir} as directory already exists.')
                continue
                
            os.makedirs(batch_dir)
                
            genome_list_file = os.path.join(batch_dir, 'genomes.lst')
            fout = open(genome_list_file, 'w')
            for i in range(start, start+batch_size):
                if i < len(genomic_files):
                    gid, gf = genomic_files[i]
                    fout.write('{}\t{}\n'.format(gf, gid))
            fout.close()
            
            cmd = 'gtdbtk classify_wf --cpus {} --force --batchfile {} --out_dir {}'.format(
                    self.cpus,
                    genome_list_file,
                    batch_dir)
            print(cmd)
            os.system(cmd)
            
        # combine summary files
        fout = open(os.path.join(self.output_dir, 'gtdbtk_classify.tsv'), 'w')
        bHeader = True
        for batch_dir in os.listdir(self.output_dir):
            if not batch_dir.startswith('batch_'):
                continue
                
            batch_dir = os.path.join(self.output_dir, batch_dir)
            ar_summary = os.path.join(batch_dir, 'gtdbtk.ar122.summary.tsv')
            bac_summary = os.path.join(batch_dir, 'gtdbtk.bac120.summary.tsv')
            
            for summary_file in [ar_summary, bac_summary]:
                with open(summary_file, encoding='utf-8') as f:
                    header = f.readline()
                    
                    if bHeader:
                        fout.write(header)
                        bHeader = False
                        
                    for line in f:
                        fout.write(line)
        
        fout.close()
        
