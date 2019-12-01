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
    
    
def filter_genomes(metadata_file,
                    min_comp,
                    max_cont,
                    min_quality, 
                    max_contigs, 
                    min_N50, 
                    max_ambiguous, 
                    max_gap_length):
    """Indentify genomes passing filtering criteria."""
    
    genome_ids = set()

    csv_reader = csv.reader(open(metadata_file, 'rb'))
    bHeader = True
    for row in csv_reader:
        if bHeader:
            bHeader = False

            genome_index = row.index('accession')
            comp_index = row.index('checkm_completeness')
            cont_index = row.index('checkm_contamination')
            contig_count_index = row.index('contig_count')
            n50_scaffolds_index = row.index('n50_scaffolds')
            ambiguous_bases_index = row.index('ambiguous_bases')
            total_gap_length_index = row.index('total_gap_length')
            
        else:
            genome_id = row[genome_index]
            
            comp = float(row[comp_index])
            cont = float(row[cont_index])
            quality = comp - 5*cont
            
            contig_count = int(row[contig_count_index])
            n50_scaffolds = int(row[n50_scaffolds_index])
            ambiguous_bases = int(row[ambiguous_bases_index])
            total_gap_length = int(row[total_gap_length_index])
            
            if comp >= min_comp and cont <= max_cont and quality >= min_quality:
                if contig_count <= max_contigs and n50_scaffolds >= min_N50:
                    if ambiguous_bases <= max_ambiguous and total_gap_length <= max_gap_length:
                        genome_ids.add(genome_id)

    return genome_ids


 
def read_marker_percentages(gtdb_domain_report):
    """Parse percentage of marker genes for each genome."""
    
    marker_perc = {}
    with open(gtdb_domain_report, encoding='utf-8') as f:
        header = f.readline().rstrip().split('\t')
        
        domain_index = header.index('Predicted domain')
        bac_marker_perc_index = header.index('Bacterial Marker Percentage')
        ar_marker_perc_index = header.index('Archaeal Marker Percentage')
        ncbi_taxonomy_index = header.index('NCBI taxonomy')
        gtdb_taxonomy_index = header.index('GTDB taxonomy')
        
        for line in f:
            line_split = line.strip().split('\t')
            
            gid = canonical_gid(line_split[0])
            domain = line_split[domain_index]
            bac_perc = float(line_split[bac_marker_perc_index])
            ar_perc = float(line_split[ar_marker_perc_index])
            ncbi_domain = [t.strip() for t in line_split[ncbi_taxonomy_index].split(';')][0]
            gtdb_domain = [t.strip() for t in line_split[gtdb_taxonomy_index].split(';')][0]

            marker_perc[gid] = max(bac_perc, ar_perc)
            
            if not gid.startswith('U'):
                if marker_perc[gid] > 10:
                    if ncbi_domain != gtdb_domain and ncbi_domain != 'None':
                        print(f'[WARNING] NCBI and GTDB domains disagree in domain report: {gid}')
                        
                        if ncbi_domain != domain and domain != 'None':
                            print(f' ... NCBI domain {ncbi_domain} also disagrees with predicted domain {domain}.')
                        
                if marker_perc[gid] > 25 and abs(bac_perc - ar_perc) > 5:
                    if domain != gtdb_domain and domain != 'None':
                        print(f'[WARNING] GTDB and predicted domain (Bac = {bac_perc:.1f}%; Ar = {ar_perc:.1f}%) disagree in domain report: {gid}')

    return marker_perc
    
    
def exclude_from_refseq(refseq_assembly_file, genbank_assembly_file):
    """Parse exclude from RefSeq field from NCBI assembly files."""
    
    excluded_from_refseq_note = {}
    for assembly_file in [refseq_assembly_file, genbank_assembly_file]:
        with open(assembly_file, encoding='utf-8') as f:
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
    
    
def read_gtdb_sp_clusters(cluster_file):
    """Read GTDB species cluster file."""
        
    clusters = defaultdict(list)
    species = {}
    with open(cluster_file) as f:
        headers = f.readline().strip().split('\t')
        
        type_sp_index = headers.index('GTDB species')
        type_genome_index = headers.index('Representative genome')
            
        num_clustered_index = headers.index('No. clustered genomes')
        clustered_genomes_index = headers.index('Clustered genomes')
        
        for line in f:
            line_split = line.strip().split('\t')
            
            sp = line_split[type_sp_index]
            rid = canonical_gid(line_split[type_genome_index])
            species[rid] = sp

            num_clustered = int(line_split[num_clustered_index])
            if num_clustered > 0:
                clusters[rid] = [canonical_gid(g.strip()) for g in line_split[clustered_genomes_index].split(',')]
            else:
                clusters[rid] = []
                    
    return clusters, species
    
    
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
        
