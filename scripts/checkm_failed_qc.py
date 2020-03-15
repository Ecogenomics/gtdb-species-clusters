#!/usr/bin/env python

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

__prog_name__ = 'checkm_failed_qc.py'
__prog_desc__ = 'run CheckM lineage-specific markers over genomes failing QC using NCBI gene calling'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2018'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse
import tempfile
import ntpath
import shutil


class RunCheckm(object):
    """Tun CheckM lineage-specific markers over genomes failing QC using NCBI gene calling."""

    def __init__(self):
        """Initialization."""
        pass

    def run(self, qc_failed_file, gtdb_genome_files, cpus, output_dir):
        """Applying CheckM to genomes."""

        protein_file_dir = os.path.join(output_dir, 'protein_files')
        if not os.path.exists(protein_file_dir):
            os.makedirs(protein_file_dir)
        
        # get path to NCBI protein files
        protein_files = {}
        for line in open(gtdb_genome_files):
            line_split = line.strip().split('\t')
            gid = line_split[0]
            genome_path = line_split[1]
            accession = os.path.basename(os.path.normpath(genome_path))
            protein_file = os.path.join(genome_path, accession + '_protein.faa')
            protein_files[gid] = protein_file
            
        # link to files for failed genomes
        with open(qc_failed_file) as f:
            header = f.readline().strip().split('\t')
            
            acc_index = header.index('Accession')
            
            for line in f:
                line_split = line.strip().split('\t')
                gid = line_split[acc_index]
            
                protein_file = protein_files[gid]
                if os.path.exists(protein_file):
                    cmd = 'ln -s %s %s' % (protein_file, os.path.join(protein_file_dir, gid + '_ncbi_proteins.faa'))
                    os.system(cmd)
                else:
                    print('Missing protein file for %s.' % gid)

        # run CheckM
        checkm_output_dir = os.path.join(output_dir, 'checkm_lineage_wf')
        os.system('checkm lineage_wf --genes --pplacer_threads %d -x faa -t %d %s %s' % (cpus, cpus, protein_file_dir, checkm_output_dir))

        tree_qa_file = os.path.join(checkm_output_dir, 'tree_qa.o2.tsv')
        os.system('checkm tree_qa -o 2 --tab_table -f %s %s' % (tree_qa_file, checkm_output_dir))

        qa_file = os.path.join(checkm_output_dir, 'qa.tsv')
        os.system('checkm qa -t %d --tab_table -f %s %s %s' % (cpus, qa_file, os.path.join(checkm_output_dir, 'lineage.ms'), checkm_output_dir))

        profile_file = os.path.join(checkm_output_dir, 'profile.tsv')
        os.system('checkm join_tables -f %s %s %s' % (profile_file, qa_file, tree_qa_file))

if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('qc_failed_file', help='file indicating genomes that failed QC')
    parser.add_argument('gtdb_genome_files', help='file indicating path to GTDB genomes')
    parser.add_argument('output_dir', help='output directory')
    parser.add_argument('-c', '--cpus', help='number of processors to use', type=int, default=40)

    args = parser.parse_args()

    try:
        runCheckm = RunCheckm()
        runCheckm.run(args.qc_failed_file, args.gtdb_genome_files, args.cpus, args.output_dir)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
