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
from collections import defaultdict, namedtuple

import biolib.seq_io as seq_io
from biolib.taxonomy import Taxonomy

from genometreetk.default_values import DefaultValues
from genometreetk.aai import aai_thresholds


# make sure large CSV files can be read
csv.field_size_limit(sys.maxsize)


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

    
def check_domain_assignment(genome_id, gtdb_taxonomy, ncbi_taxonomy, rep_is_bacteria):
    """Check that domain assignment agrees with GTDB and NCBI taxonomy.
    
    Parameters
    ----------
    genome_id : str
        Genome of interest.
    gtdb_taxonomy : list
        GTDB taxonomy.
    ncbi_taxonomy : list
        NCBI taxonomy.
    rep_is_bacteria : boolean
        Flag indicating if genome was classified as a Bacteria (True) or  Archaea (False).
        
    Returns
    -------
    boolean
        True if taxonomies agrees with domain assignment, else False.
    """
    
    gtdb_domain = gtdb_taxonomy[genome_id][0]
    ncbi_domain = ncbi_taxonomy[genome_id][0]
    if gtdb_domain != 'd__' and ncbi_domain != 'd__' and gtdb_domain != ncbi_domain:
        print '[Warning] GTDB and NCBI domain assignments do not agree: %s' % genome_id
        return False
    else:
        for domain in [gtdb_domain, ncbi_domain]:
            if domain == 'd__':
                continue
                
            if rep_is_bacteria and domain != 'd__Bacteria':
                print '[Warning] Taxonomy and predicted bacterial domain assignment do not agree: %s' % genome_id
                print '*%s\t%s' % (genome_id, ';'.join(Taxonomy.rank_prefixes))
                print gtdb_taxonomy[genome_id]
                print ncbi_taxonomy[genome_id]
                return False
            elif not rep_is_bacteria and domain != 'd__Archaea':
                print '[Warning] Taxonomy and predicted archaeal domain assignment do not agree: %s' % genome_id
                print '*%s\t%s' % (genome_id, ';'.join(Taxonomy.rank_prefixes))
                print gtdb_taxonomy[genome_id]
                print ncbi_taxonomy[genome_id]
                return False
                
    return True
    

def predict_bacteria(genome_id, bac_seqs, ar_seqs):
    """Check that domain assignment agrees with GTDB and NCBI taxonomy.
    
    Parameters
    ----------
    genome_id : str
        Genome of interest.
    bac_seqs : dict
        Bacterial sequences for all genomes.
    ar_seqs : dict
        Archaea sequences for all genomes.
    rep_is_bacteria : boolean
        Flag indicating if genome was classified as a Bacteria (True) or  Archaea (False).
        
    Returns
    -------
    boolean
        True if predicted domain is Bacteria, else False for Archaea.
    """
    rep_bac_seq = bac_seqs[genome_id]
    rep_ar_seq = ar_seqs[genome_id]
    
    per_bac_aa = float(len(rep_bac_seq) - rep_bac_seq.count('-')) / len(rep_bac_seq)
    per_ar_aa = float(len(rep_ar_seq) - rep_ar_seq.count('-')) / len(rep_ar_seq)
    
    is_bac = per_bac_aa >= per_ar_aa
    
    return is_bac, per_bac_aa, per_ar_aa
    
    
def reassign_representative(cur_representative_id,
                                cur_aai,
                                new_representative_id,
                                new_aai,
                                trusted_user_genomes):
    """Determines best representative genome.

    Genomes are preferentially assigned to representatives based on
    source repository (RefSeq => GenBank => Trusted User => User) and AAI.
    """
    
    source_order = {'R': 0,  # RefSeq
                    'G': 1,  # GenBank
                    'U': 2}  # User

    if not cur_representative_id:
        # no currently assigned representative
        return new_representative_id, new_aai

    cur_source = source_order[cur_representative_id[0]]
    new_source = source_order[new_representative_id[0]]
    if new_representative_id in trusted_user_genomes:
        new_source = source_order['G']
    
    if new_source < cur_source:
        # give preference to genome source
        return new_representative_id, new_aai
    elif new_source == cur_source:
        # for the same source, find the representative with the highest AAI
        if new_aai > cur_aai:
            return new_representative_id, new_aai

    return cur_representative_id, cur_aai

def assign_rep(rep_id, 
                    genome_id,
                    rep_is_bacteria, 
                    genome_is_bacteria,
                    bac_seqs, 
                    ar_seqs,
                    species,
                    gtdb_taxonomy,
                    genome_aa_count,
                    trusted_user_genomes,
                    aai_threshold,
                    min_matches,
                    assigned_representative, 
                    cur_aai):
    """Detering if genome should be assigned to a representative."""

    # do not cluster genomes from different predicted domains
    if rep_is_bacteria != genome_is_bacteria:
        return assigned_representative, cur_aai
        
    # do not cluster genomes from different named species
    rep_species = species.get(rep_id, None)
    genome_species = species.get(genome_id, None)
    if rep_species and genome_species and rep_species != genome_species:
        return assigned_representative, cur_aai

    # do not cluster NCBI or trusted User genomes with User representatives
    if rep_id.startswith('U_') and rep_id not in trusted_user_genomes: 
        if not genome_id.startswith('U__') or genome_id in trusted_user_genomes:
            return assigned_representative, cur_aai
            
    # do not cluster genomes from different named groups
    skip = False
    for rep_taxon, taxon in zip(gtdb_taxonomy[rep_id], gtdb_taxonomy[genome_id]):
        rep_taxon = rep_taxon[3:]
        taxon = taxon[3:]
        if rep_taxon and taxon and rep_taxon != taxon:
            skip = True
            break
    if skip:
        return assigned_representative, cur_aai
            
    if rep_is_bacteria:
        rep_seq = bac_seqs[rep_id]
        genome_seq = bac_seqs[genome_id]
    else:
        rep_seq = ar_seqs[rep_id]
        genome_seq = ar_seqs[genome_id]
        
    max_mismatches = (1.0 - cur_aai) * genome_aa_count
    aai = aai_thresholds(rep_seq, genome_seq, max_mismatches, min_matches)
    if aai > aai_threshold:  
        assigned_representative, cur_aai = reassign_representative(assigned_representative,
                                                                            cur_aai,
                                                                            rep_id,
                                                                            aai,
                                                                            trusted_user_genomes)
        
    return assigned_representative, cur_aai
    
        
def species_label(gtdb_taxonomy, ncbi_taxonomy, ncbi_organism_name):
    """Determine 'best' species label for each genome.

    Currently, this is just being set to the species label in the
    GTDB taxonomy. In theory, the NCBI taxonomy and organism name
    could also be consulted. However, since the GTDB taxonomy redefines
    some species this might be problematic so isn't currently being
    done.

    Parameters
    ----------
    gtdb_taxonomy : d[assembly_accession] -> [d__, ..., s__]
        GTDB taxonomy of each genome.
    ncbi_taxonomy : d[assembly_accession] -> [d__, ..., s__]
        NCBI taxonomy of each genome.
    ncbi_organism_name : d[assembly_accession] -> name
        NCBI organism name of each genome.

    Return
    ------
    dict : d[assembly_accession] -> species name
        Species name of each genome.
    """

    taxonomy = Taxonomy()

    species = {}
    species_index = Taxonomy.rank_index['s__']
    for genome_id, taxa in gtdb_taxonomy.iteritems():
        sp = taxa[species_index]
        if sp != 's__':
            species[genome_id] = sp

    if False:   # do not consider NCBI information as
                # it may conflict with GTDB information
                # in unwanted ways
        
        for genome_id, taxa in ncbi_taxonomy.iteritems():
            if genome_id in species:
                continue

            sp = taxa[species_index]
            sp = taxonomy.extract_valid_species_name(sp)
            if sp:
                species[genome_id] = sp

        for genome_id, sp in ncbi_organism_name.iteritems():
            if genome_id in species:
                continue

            sp = taxonomy.extract_valid_species_name(sp)
            if sp:
                species[genome_id] = sp

    return species

    
def read_gtdb_metadata(metadata_file, fields):
    """Parse genome quality from GTDB metadata.

    Parameters
    ----------
    metadata_file : str
        Metadata for all genomes in CSV file.
    fields : iterable
        Fields  to read.

    Return
    ------
    dict : d[genome_id] -> namedtuple
        Value for fields indicted by genome IDs.
    """

    gtdb_metadata = namedtuple('gtdb_metadata', ' '.join(fields))
    m = {}

    csv_reader = csv.reader(open(metadata_file, 'rb'))
    bHeader = True
    for row in csv_reader:
        if bHeader:
            headers = row

            genome_index = headers.index('accession')

            indices = []
            for field in fields:
                indices.append(headers.index(field))

            bHeader = False
        else:
            genome_id = row[genome_index]

            values = []
            for i in indices:
                # save values as floats or strings
                try:
                    values.append(float(row[i]))
                except ValueError:
                    values.append(row[i])
            m[genome_id] = gtdb_metadata._make(values)

    return m


def read_gtdb_phylum(metadata_file):
    """Parse GTDB phylum information from GTDB metadata.

    Parameters
    ----------
    metadata_file : str
        Metadata for all genomes.

    Return
    ------
    dict : d[genome_id] -> phyla
    """

    genome_phyla = {}

    csv_reader = csv.reader(open(metadata_file, 'rt'))
    bHeader = True
    for row in csv_reader:
        if bHeader:
            headers = row
            genome_index = headers.index('accession')
            phylum_index = headers.index('gtdb_phylum')
            bHeader = False
        else:
            genome_id = row[genome_index]
            genome_phyla[genome_id] = row[phylum_index]

    return genome_phyla


def read_gtdb_taxonomy(metadata_file):
    """Parse GTDB taxonomy from GTDB metadata.

    Parameters
    ----------
    metadata_file : str
        Metadata for all genomes.

    Return
    ------
    dict : d[genome_id] -> taxonomy list
    """

    taxonomy = {}

    csv_reader = csv.reader(open(metadata_file, 'rt'))
    bHeader = True
    for row in csv_reader:
        if bHeader:
            headers = row
            genome_index = headers.index('accession')
            taxonomy_index = headers.index('gtdb_taxonomy')
            bHeader = False
        else:
            genome_id = row[genome_index]
            taxa_str = row[taxonomy_index].strip()

            if taxa_str:
                taxonomy[genome_id] = map(str.strip, taxa_str.split(';'))
            else:
                taxonomy[genome_id] = list(Taxonomy.rank_prefixes)

    return taxonomy
    
def read_gtdb_representative(metadata_file):
    """Parse GTDB representative from GTDB metadata.

    Parameters
    ----------
    metadata_file : str
        Metadata for all genomes.

    Return
    ------
    dict : d[genome_id] -> True or False
    """

    gtdb_reps = {}

    csv_reader = csv.reader(open(metadata_file, 'rt'))
    bHeader = True
    for row in csv_reader:
        if bHeader:
            headers = row
            genome_index = headers.index('accession')
            gtdb_representative_index = headers.index('gtdb_representative')
            bHeader = False
        else:
            genome_id = row[genome_index]
            is_rep = (row[gtdb_representative_index] == 't')
            gtdb_reps[genome_id] = is_rep

    return gtdb_reps


def read_gtdb_ncbi_taxonomy(metadata_file):
    """Parse NCBI taxonomy from GTDB metadata.

    Parameters
    ----------
    metadata_file : str
        Metadata for all genomes.

    Return
    ------
    dict : d[genome_id] -> taxonomy list
    """

    taxonomy = {}

    csv_reader = csv.reader(open(metadata_file, 'rt'))
    bHeader = True
    for row in csv_reader:
        if bHeader:
            headers = row
            genome_index = headers.index('accession')
            taxonomy_index = headers.index('ncbi_taxonomy')
            bHeader = False
        else:
            genome_id = row[genome_index]
            taxa_str = row[taxonomy_index].strip()

            if taxa_str:
                taxonomy[genome_id] = taxa_str.split(';')
            else:
                taxonomy[genome_id] = list(Taxonomy.rank_prefixes)

    return taxonomy


def read_gtdb_ncbi_organism_name(metadata_file):
    """Parse NCBI organism name from GTDB metadata.

    Parameters
    ----------
    metadata_file : str
        Metadata for all genomes.

    Return
    ------
    dict : d[genome_id] -> organism name
        Organism name of each genome.
    """

    d = {}

    csv_reader = csv.reader(open(metadata_file, 'rt'))
    bHeader = True
    for row in csv_reader:
        if bHeader:
            headers = row
            genome_index = headers.index('accession')
            organism_name_index = headers.index('ncbi_organism_name')
            bHeader = False
        else:
            genome_id = row[genome_index]
            organism_name = row[organism_name_index].strip()

            if organism_name:
                d[genome_id] = organism_name

    return d


def read_gtdb_ncbi_type_strain(metadata_file):
    """Parse NCBI type strain from GTDB metadata.

    Parameters
    ----------
    metadata_file : str
        Metadata for all genomes.

    Return
    ------
    set
        Set of genomes marked as type strains by NCBI.
    """

    type_strains = set()

    csv_reader = csv.reader(open(metadata_file, 'rt'))
    bHeader = True
    for row in csv_reader:
        if bHeader:
            headers = row
            genome_index = headers.index('accession')
            type_strain_index = headers.index('ncbi_type_strain')
            bHeader = False
        else:
            genome_id = row[genome_index]
            if bool(row[type_strain_index]):
                type_strains.add(genome_id)

    return type_strains


def read_tree_model(reportFile):
    for line in open(reportFile):
        if 'Model of evolution:' in line:
            modelStr = line[line.find(':') + 1:].strip()

    return modelStr


def read_genome_dir_file(genome_dir_file):
    """Read genome directories from file.

    Parameters
    ----------
    genome_dir_file : str
        File to read.

    Returns
    -------
    dict : d[genome_id] -> directory
        Directory for each genome.
    """

    # read directory for each genome
    genome_dirs = {}
    for line in open(genome_dir_file):
        line_split = line.split('\t')
        genome_dirs[line_split[0]] = line_split[1].strip()

    return genome_dirs


def read_marker_id_file(marker_id_file):
    """Read marker ids from file.

    Parameters
    ----------
    marker_id_file : str
        File to read.

    Returns
    -------
    set
        Marker ids.
    """

    marker_genes = set()
    for line in open(marker_id_file):
        if line[0] == '#':
            continue

        marker_genes.add(line.split('\t')[0].strip())

    return marker_genes


def read_genome_id_file(genome_id_file):
    """Read genome ids from file.

    Parameters
    ----------
    genome_file : str
        File to read.

    Returns
    -------
    set
        NCBI genome ids.
    set
        User genome ids.
    """

    ncbi_genome_ids = set()
    user_genome_ids = set()
    for line in open(genome_id_file):
        if line[0] == '#':
            continue

        if '\t' in line:
            genome_id = line.split('\t')[0].strip()
        else:
            genome_id = line.split()[0].strip()

        if genome_id.startswith('U_'):
            user_genome_ids.add(genome_id)
        else:
            ncbi_genome_ids.add(genome_id)

    return ncbi_genome_ids, user_genome_ids


def create_concatenated_alignment(genome_ids,
                                   marker_genes,
                                   alignment_dir,
                                   concatenated_alignment_file,
                                   marker_file):
    """Create concatenated multiple sequence alignment for all genomes.

    Parameters
    ----------
    genome_ids : iterable
        Genomes of interest.
    marker_genes : iterable
        Unique ids of marker genes.
    alignment_dir : str
        Directory containing multiple sequence alignments.
    concatenated_alignment_file : str
        File to containing concatenated alignment.
    marker_file : str
        File indicating length of each marker in the alignment.
    """

    # Read alignment files. Some genomes may have multiple
    # copies of a marker gene in which case the last one
    # is arbitrarily taken. This is acceptable as all genes
    # are already screen to be conspecific.
    alignments = defaultdict(dict)
    marker_length = {}
    for mg in marker_genes:
        f = mg + '.aln.masked.faa'
        seqs = seq_io.read_fasta(os.path.join(alignment_dir, f))

        for seq_id, seq in seqs.iteritems():
            genome_id = seq_id[0:seq_id.find(DefaultValues.SEQ_CONCAT_CHAR)]

            alignments[mg][genome_id] = seq

            marker_length[mg] = len(seq)

    # create marker file
    fout = open(marker_file, 'w')
    for mg in marker_genes:
        fout.write('%s\t%s\t%s\t%d\n' % (mg, mg, mg, marker_length[mg]))
    fout.close()

    # create concatenated alignment
    concatenated_seqs = {}
    for mg in marker_genes:
        seqs = alignments[mg]

        for genome_id in genome_ids:
            if genome_id in seqs:
                # append alignment
                concatenated_seqs[genome_id] = concatenated_seqs.get(genome_id, '') + seqs[genome_id]
            else:
                # missing gene
                concatenated_seqs[genome_id] = concatenated_seqs.get(genome_id, '') + '-' * marker_length[mg]

    # save concatenated alignment
    seq_io.write_fasta(concatenated_seqs, concatenated_alignment_file)

