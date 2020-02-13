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
import csv
import sys
from copy import deepcopy
from collections import defaultdict, namedtuple

import biolib.seq_io as seq_io
from biolib.taxonomy import Taxonomy


def canonical_gid(gid):
    """Get canonical form of NCBI genome accession.
    
    Example:
        G005435135 -> G005435135
        GCF_005435135.1 -> G005435135
        GCF_005435135.1_ASM543513v1_genomic -> G005435135
        RS_GCF_005435135.1 -> G005435135
        GB_GCA_005435135.1 -> G005435135
    """
    
    if gid.startswith('U'):
        return gid
        
    gid = gid.replace('RS_', '').replace('GB_', '')
    gid = gid.replace('GCA_', 'G').replace('GCF_', 'G')
    if '.' in gid:
        gid = gid[0:gid.find('.')]
    
    return gid
    
    
def select_highest_quality(gids, cur_genomes):
    """Select highest quality genome."""
    
    # sort by decreasing type strain score, followed by 
    # genome ID in order to ensure ties are broken in a
    # deterministic fashion between runs
    q = {k:cur_genomes[k].score_type_strain() for k in gids}
    q_sorted = sorted(q.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
    
    return q_sorted[0][0]


def read_genome_path(genome_path_file):
    """Determine path to genomic FASTA file for each genome."""

    genome_files = {}
    for line in open(genome_path_file):
        line_split = line.strip().split('\t')
        
        gid = line_split[0]
        gid = canonical_gid(gid)
        
        genome_path = line_split[1]
        accession = os.path.basename(os.path.normpath(genome_path))
        
        genome_files[gid] = os.path.join(genome_path, accession + '_genomic.fna')
        
    return genome_files

    
def exclude_from_refseq(genbank_assembly_file):
    """Parse exclude from RefSeq field from NCBI assembly files."""
    
    excluded_from_refseq_note = {}
    with open(genbank_assembly_file, encoding='utf-8') as f:
        for line in f:
            if line[0] == '#':
                if line.startswith('# assembly_accession'):
                    header = line.strip().split('\t')
                    exclude_index = header.index('excluded_from_refseq')
            else:
                line_split = line.strip('\n\r').split('\t')
                gid = canonical_gid(line_split[0])
                excluded_from_refseq_note[gid] = line_split[exclude_index]
    
    return excluded_from_refseq_note
    
    
def read_qc_file(qc_file):
    """Read genomes passing QC from file."""
    
    passed_qc = set()
    with open(qc_file, encoding='utf-8') as f:
        f.readline()
        
        for line in f:
            line_split = line.strip().split('\t')
            gid = canonical_gid(line_split[0])
            passed_qc.add(gid)
            
    return passed_qc
    
    
def read_cur_new_updated(genomes_new_updated_file):
    """Determine new and updated genomes in current GTDB release."""
    
    cur_new = set()
    cur_updated = set()
    with open(genomes_new_updated_file, encoding='utf-8') as f:
        header = f.readline().strip().split('\t')
        
        status_index = header.index('Status')
        
        for line in f:
            line_split = line.strip().split('\t')
            
            gid = line_split[0]
            status = line_split[status_index]
            if status == 'NEW':
                cur_new.add(gid)
            elif status == 'UPDATED':
                cur_updated.add(gid)
                
    return cur_new, cur_updated
    
    
def read_gtdbtk_classifications(gtdbtk_classify_file):
    """Determine classification of genomes according to GTDB-Tk."""
    
    classification = {}
    with open(gtdbtk_classify_file) as f:
        header = f.readline().strip().split('\t')
        
        classification_index = header.index('classification')
        
        for line in f:
            line_split = line.strip().split('\t')
            
            gid = canonical_gid(line_split[0])
            taxa = [t.strip() for t in line_split[classification_index].split(';')]
            
            classification[gid] = taxa
            
    return classification
        
