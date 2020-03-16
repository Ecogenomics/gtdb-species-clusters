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

__prog_name__ = 'checkm_compare.py'
__prog_desc__ = 'compare CheckM estimates'

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


class Compare(object):
    """Compare CheckM estimates."""

    def __init__(self):
        """Initialization."""
        pass

    def run(self, qc_failed_file, checkm_qa_files, output_file):
        """compare CheckM estimates."""
        
        orig_estimates = {}
        with open(qc_failed_file) as f:
            header = f.readline().strip().split('\t')
            
            acc_index = header.index('Accession')
            comp_index = header.index('Completeness (%)')
            cont_index = header.index('Contamination (%)')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gid = line_split[acc_index]
                comp = float(line_split[comp_index])
                cont = float(line_split[cont_index])
                
                orig_estimates[gid] = (comp, cont)
                
        new_estimates = {}
        with open(checkm_qa_files) as f:
            header = f.readline().strip().split('\t')
            
            comp_index = header.index('Completeness')
            cont_index = header.index('Contamination')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gid = line_split[0].replace('_ncbi_proteins', '')
                comp = float(line_split[comp_index])
                cont = float(line_split[cont_index])
                
                new_estimates[gid] = (comp, cont)
        
        fout = open(output_file, 'w')
        fout.write('Accession\tOriginal completeness\tNew completeness\tOriginal contamination\tNew contamination\n')
        for gid in new_estimates:
            orig_comp, orig_cont = orig_estimates[gid]
            new_comp, new_cont = new_estimates[gid]
            
            orig_quality = orig_comp - 5*orig_cont
            if orig_quality >= 50:
                continue
                
            new_quality = new_comp - 5*new_cont
            if new_quality < 50:
                continue
                
            if (new_comp - orig_comp > 5
                or new_cont - orig_cont < -1):
                print(gid, orig_comp, new_comp, orig_cont, new_cont)
                fout.write('%s\t%.2f\t%.2f\t%.2f\t%.2f\n' % (gid, orig_comp, new_comp, orig_cont, new_cont))
        fout.close()

if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('qc_failed_file', help='file indicating genomes that failed QC')
    parser.add_argument('checkm_qa_files', help='file with alternative CheckM estimates')
    parser.add_argument('output_file', help='output directory')

    args = parser.parse_args()

    try:
        p = Compare()
        p.run(args.qc_failed_file, args.checkm_qa_files, args.output_file)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
