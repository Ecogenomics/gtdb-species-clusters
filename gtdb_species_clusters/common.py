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
    
    
def accession_version(gid):
    """Get version of NCBI genome accession."""
    
    version = 0
    if '.' in gid:
        version = int(gid[gid.find('.')+1:])
    
    return version


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
        

def binomial_species(taxonomy):
    """Get binomial, including Candidatus, species names in NCBI taxonomy."""
    
    binomial_names = defaultdict(set)
    for gid, taxa in taxonomy.items():
        species = taxa[6]
        if species == 's__':
            continue
            
        if any(c.isdigit() for c in species):
            continue
        
        is_candidatus = False
        if 'Candidatus ' in species:
            species = species.replace('Candidatus ', '')
            is_candidatus = True
            
        tokens = species[3:].split()
        if len(tokens) != 2:
            continue
            
        genus, specific = tokens

        if is_candidatus:
            species = species.replace('s__', 's__Candidatus ')
            
        if genus.istitle() and all(c.islower() for c in specific):
            binomial_names[species].add(gid)
        
    return binomial_names


def generic_name(species_name):
    """Get the generic name of species name."""
    
    if species_name == 's__':
        return ''
    
    s = species_name.replace('Candidatus ', '')
    generic, specific = s[3:].split()
    
    return generic
    
    
def specific_epithet(species_name):
    """Get the specific epithet of species name."""
    
    if species_name == 's__':
        return ''
    
    s = species_name.replace('Candidatus ', '')
    generic, specific = s.split()
    
    return specific

    
def genome_species_assignments(ncbi_taxonomy):
    """Get species assignment for each genome."""
    
    ncbi_species = binomial_species(ncbi_taxonomy)
    gid_to_species = {}
    for sp in ncbi_species:
        for gid in ncbi_species[sp]:
            gid_to_species[gid] = sp
            
    return gid_to_species
    
    
def canonical_generic_name(species_name):
    """Get species name with canonicalized generic name."""
    
    generic, specific = species_name[3:].split()
    if '_' in generic:
        generic = generic[0:generic.rfind('_')]
        
    return f's__{generic} {specific}'
        
    
def canonical_species_name(species_name):
    """Get canonical species name from GTDB or NCBI species name."""
    
    if species_name == 's__':
        return species_name
    
    species_name = species_name.replace('Candidatus ', '')
    prefix = species_name[0:3]
    full_name = species_name[3:]
    genus, species = full_name.split(' ')[0:2]
    
    underscore_pos = genus.rfind('_')
    if underscore_pos != -1:
        genus = genus[0:underscore_pos]
        
    underscore_pos = species.rfind('_')
    if underscore_pos != -1:
        species = species[0:underscore_pos]
        
    genus = genus.replace('[','').replace(']','')
    species = species.replace('[','').replace(']','')
        
    return prefix + genus + ' ' + species
    
    
def generic_specific_names(species_name):
    """Get generic and specific names from species name."""
    
    species_name = species_name.replace('s__', '')
    species_name = species_name.replace('Candidatus ', '')
    generic, specific = species_name.split(' ')
    
    return generic, specific
    
    
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
        print('[Warning] GTDB and NCBI domain assignments do not agree: %s' % genome_id)
        return False
    else:
        for domain in [gtdb_domain, ncbi_domain]:
            if domain == 'd__':
                continue
                
            if rep_is_bacteria and domain != 'd__Bacteria':
                print('[Warning] Taxonomy and predicted bacterial domain assignment do not agree: %s' % genome_id)
                print('*%s\t%s' % (genome_id, ';'.join(Taxonomy.rank_prefixes)))
                print(gtdb_taxonomy[genome_id])
                print(ncbi_taxonomy[genome_id])
                return False
            elif not rep_is_bacteria and domain != 'd__Archaea':
                print('[Warning] Taxonomy and predicted archaeal domain assignment do not agree: %s' % genome_id)
                print('*%s\t%s' % (genome_id, ';'.join(Taxonomy.rank_prefixes)))
                print(gtdb_taxonomy[genome_id])
                print(ncbi_taxonomy[genome_id])
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

    with open(metadata_file, encoding='utf-8') as f:
        headers = f.readline().strip().split('\t')

        genome_index = headers.index('accession')

        indices = []
        for field in fields:
            indices.append(headers.index(field))

        for line in f:
            line_split = line.strip().split('\t')
            genome_id = canonical_gid(line_split[genome_index])

            values = []
            for i in indices:
                # save values as floats or strings
                v = line_split[i]
                try:
                    values.append(float(v))
                except ValueError:
                    if v is None or v == '' or v == 'none':
                        values.append(None)
                    elif v == 'f' or v.lower() == 'false':
                        values.append(False)
                    elif v == 't' or v.lower() == 'true':
                        values.append(True)
                    else:
                        values.append(v)
            m[genome_id] = gtdb_metadata._make(values)

    return m


def read_gtdb_accessions(metadata_file):
    """Parse NCBI accession information from GTDB metadata."""
    
    cur_accns = {}
    with open(metadata_file, encoding='utf-8') as f:
        f.readline()
        for line in f:
            line_split = line.strip().split('\t')
            
            gid = line_split[0]
            cur_accns[canonical_gid(gid)] = gid

    return cur_accns
    

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

    with open(metadata_file, encoding='utf-8') as f:
        headers = f.readline().strip().split('\t')
        genome_index = headers.index('accession')
        taxonomy_index = headers.index('gtdb_taxonomy')
        for line in f:
            line_split = line.strip().split('\t')
            genome_id = canonical_gid(line_split[genome_index])
            taxa_str = line_split[taxonomy_index].strip()

            if taxa_str and taxa_str != 'none':
                taxonomy[genome_id] = [t.strip() for t in taxa_str.split(';')]
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

    with open(metadata_file) as f:
        headers = f.readline().strip().split('\t')
        genome_index = headers.index('accession')
        gtdb_representative_index = headers.index('gtdb_representative')

        for line in f:
            line_split = line.strip().split('\t')
            genome_id = line_split[genome_index]
            is_rep = (line_split[gtdb_representative_index] == 't')
            gtdb_reps[genome_id] = is_rep

    return gtdb_reps


def read_gtdb_ncbi_taxonomy(metadata_file, 
                            species_exception_file,
                            genus_exception_file):
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
    with open(metadata_file, encoding='utf-8') as f:
        headers = f.readline().strip().split('\t')
        genome_index = headers.index('accession')
        taxonomy_index = headers.index('ncbi_taxonomy')
        
        for line in f:
            line_split = [token.strip() for token in line.strip().split('\t')]
            gid = canonical_gid(line_split[genome_index])
            taxa_str = line_split[taxonomy_index].strip()
            taxa_str = taxa_str.replace('Candidatus ', '')

            if taxa_str and taxa_str != 'none':
                taxonomy[gid] = [t.strip() for t in taxa_str.split(';')]
            else:
                taxonomy[gid] = list(Taxonomy.rank_prefixes)
    
    ncbi_update_count = 0
    species_updates = {}
    if species_exception_file:
        with open(species_exception_file, encoding='utf-8') as f:
            f.readline()
            for line in f:
                line_split = [token.strip() for token in line.strip().split('\t')]
                gid = canonical_gid(line_split[0])

                sp = line_split[1].replace('Candidatus ', '')
                if gid not in taxonomy:
                    print('Genome in species exception list not defined at NCBI: %s' % gid)
                    sys.exit(-1)
                    
                if not sp.startswith('s__'):
                    sp = 's__' + sp
                    
                taxonomy[gid][6] = sp
                ncbi_update_count += 1
                
                species_updates[gid] = sp
    
    if genus_exception_file:
        with open(genus_exception_file, encoding='utf-8') as f:
            f.readline()
            for line in f:
                line_split = [token.strip() for token in line.strip().split('\t')]
                gid = canonical_gid(line_split[0])
                genus = line_split[1]
                if gid not in taxonomy:
                    print('Genome in genus exception list not defined at NCBI: %s' % gid)
                    sys.exit(-1)
                    
                if genus.startswith('g__'):
                    genus = genus[3:]
                    
                taxonomy[gid][5] = f'g__{genus}'
                
                species = taxonomy[gid][6]
                if species != 's__':
                    generic, specific = generic_specific_names(species)
                    taxonomy[gid][6] = f's__{genus} {specific}'
                ncbi_update_count += 1
            
                # sanity check ledgers
                if gid in species_updates and genus not in species_updates[gid]:
                    self.logger.error(f'Species and genus ledgers have conflicting assignments for {gid}.')
                    sys.exit(-1)

    return taxonomy, ncbi_update_count


def read_ncbi_subsp(metadata_file):
    """Read NCBI subspecies classification for genomes."""
    
    subsp = {}
    with open(metadata_file, encoding='utf-8') as f:
        headers = f.readline().strip().split('\t')
        genome_index = headers.index('accession')
        taxonomy_index = headers.index('ncbi_taxonomy_unfiltered')
        
        for line in f:
            line_split = [token.strip() for token in line.strip().split('\t')]
            gid = canonical_gid(line_split[genome_index])
            taxa_str = line_split[taxonomy_index].strip()
            taxa_str = taxa_str.replace('Candidatus ', '')

            if taxa_str and taxa_str != 'none':
                for taxon in [t.strip() for t in taxa_str.split(';')]:
                    if taxon.startswith('sb__'):
                        subsp[gid] = taxon
                
    return subsp
    

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

    with open(metadata_file) as f:
        headers = f.readline().strip().split('\t')
        genome_index = headers.index('accession')
        organism_name_index = headers.index('ncbi_organism_name')

        for line in f:
            line_split = line.strip().split('\t')
            genome_id = line_split[genome_index]
            organism_name = line_split[organism_name_index].strip()

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

    with open(metadata_file) as f:
        headers = f.readline().strip().split('\t')
        genome_index = headers.index('accession')
        type_strain_index = headers.index('ncbi_type_strain')

        for line in f:
            line_split = line.strip().split('\t')
            genome_id = line_split[genome_index]
            if bool(line_split[type_strain_index]):
                type_strains.add(genome_id)

    return type_strains


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
        
        
def rep_change_gids(rep_change_summary_file, field, value):
    """Get genomes with a specific change."""
    
    gids = {}
    with open(rep_change_summary_file) as f:
        header = f.readline().strip().split('\t')
        
        field_index = header.index(field)
        prev_sp_index = header.index('Previous GTDB species')
        
        for line in f:
            line_split = line.strip().split('\t')
            
            v = line_split[field_index]
            if v == value:
                prev_sp = line_split[prev_sp_index]
                gids[line_split[0]] = prev_sp
                
    return gids
    
    
def find_new_type_strains(rid, 
                            prev_gids_type_strain,
                            cur_gids_type_strain,
                            expanded_sp_clusters):
    """Check if species cluster has new genomes from type strain."""
    
    new_type_strains = set()
    for gid in expanded_sp_clusters[rid]:
        if gid not in prev_gids_type_strain and gid in cur_gids_type_strain:
            new_type_strains.add(gid)
            
    return new_type_strains
    
    
def expand_sp_clusters(sp_clusters, 
                        gtdb_species, 
                        gtdbtk_classifications, 
                        cur_new,
                        cur_updated,
                        gids_pass_qc):
    """Expand species clusters to contain new genomes and verify assignment of updated genomes."""
    
    # create mapping between species and representatives
    sp_rid_map = {sp: rid for rid, sp in gtdb_species.items()}
    
    # create mapping between all genomes and species
    gid_sp_map = {}
    for rid, cids in sp_clusters.items():
        sp = gtdb_species[rid]
        for cid in cids:
            gid_sp_map[cid] = sp
                            
    # expand species clusters
    expanded_sp_clusters = deepcopy(sp_clusters)
    
    sp_assignments = 0
    failed_qc = 0
    prev_genome_count = 0
    for gid, taxa in gtdbtk_classifications.items():
        sp = taxa[6]
        if sp == 's__':
            continue
            
        if gid not in gids_pass_qc:
            # ***HACK: this should not be necessary, except GTDB-Tk was run external
            # of complete workflow for R95
            failed_qc += 1
            continue
            
        if sp not in sp_rid_map:
            self.logger.error(f'GTDB-Tk results indicated a new species for {gid}: {sp}')
            sys.exit(-1)
            
        rid = sp_rid_map[sp]
        sp_assignments += 1
        if gid in cur_new:
            expanded_sp_clusters[rid].append(gid)
        elif gid in cur_updated:
            prev_sp = gid_sp_map[gid]
            if prev_sp != sp:
                self.logger.warning(f'Updated genomes {gid} reassigned from {prev_sp} to {sp}.')
                sys.exit(-1)
                # Really, should handle this case. This will be fine so long as the genomes
                # isn't a species representative. If a species representative has changed to
                # the point where it no longer clusters with its previous genome that requires
                # some real thought.
        else:
            # ***HACK: should be an error except GTDB-Tk was run external to workflow in R95
            #self.logger.error(f"Genome {gid} specified in GTDB-Tk results is neither 'new' or 'updated'")
            #sys.exit(-1)
            prev_genome_count += 1

    # ***HACK: this should not be necessary, except GTDB-Tk was run external
    # of complete workflow for R95
    print('failed_qc', failed_qc)
    print('prev_genome_count', prev_genome_count)
        
    return expanded_sp_clusters, sp_assignments
