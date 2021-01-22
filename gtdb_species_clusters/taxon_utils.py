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
import logging
from copy import deepcopy
from itertools import takewhile
from collections import defaultdict, namedtuple

import biolib.seq_io as seq_io
from biolib.taxonomy import Taxonomy

from gtdb_species_clusters.genome_utils import canonical_gid


def ncbi_to_gtdb_synonyms(ncbi_synonym_file, final_gtdb_taxonomy):
    """Convert synonyms defined in terms of NCBI species to GTDB species."""
    
    # parse NCBI synonyms
    ncbi_synonyms = {}
    with open(ncbi_synonym_file) as f:
        headers = f.readline().strip().split('\t')
        
        ncbi_sp_index = headers.index('NCBI species')
        rid_index = headers.index('GTDB representative')
        ncbi_synonym_index = headers.index('NCBI synonym')
        synonym_gid_index = headers.index('Highest-quality synonym genome')
        
        for line in f:
            tokens = line.strip().split('\t')
            
            ncbi_species = tokens[ncbi_sp_index]
            rid = tokens[rid_index]
            ncbi_synonym = tokens[ncbi_synonym_index]
            synonym_gid = tokens[synonym_gid_index]
            
            ncbi_synonyms[synonym_gid] = (rid, ncbi_species, ncbi_synonym)

    # convert to GTDB synonyms
    gtdb_synonyms = {}
    for synonym_gid, (rid, ncbi_species, ncbi_synonym) in ncbi_synonyms.items():
        if rid not in final_gtdb_taxonomy:
            continue # must be from other domain
            
        gtdb_species = final_gtdb_taxonomy[rid][Taxonomy.SPECIES_INDEX]
        gtdb_specific = specific_epithet(gtdb_species)
        ncbi_specific = specific_epithet(ncbi_species)

        gtdb_synonyms[ncbi_synonym] = gtdb_species

        if not test_same_epithet(gtdb_specific, ncbi_specific):
            logging.getLogger('timestamp').warning('Specific name of species with synonyms differs between GTDB and NCBI for {}: {} {}'.format(
                                                    rid,
                                                    gtdb_species,
                                                    ncbi_species))
    
    return gtdb_synonyms


def longest_common_prefix(*s):
    """Find the longest common prefix.
    
    https://rosettacode.org/wiki/Longest_common_prefix#Python
    """
    
    return ''.join(a for a,b in takewhile(lambda x: x[0] == x[1], zip(min(s), max(s))))


def test_same_epithet(epithet1, epithet2):
    """Test if species epithet are the same, except for changes due to difference in the gender of the genus."""
    
    if epithet1 == epithet2:
        return True
        
    if is_placeholder_sp_epithet(epithet1) or is_placeholder_sp_epithet(epithet2):
        return False
        
    lcp = longest_common_prefix(epithet1, epithet2)
    if ((len(lcp) >= max(len(epithet1), len(epithet2)) - 3)
        and lcp != epithet1 and lcp != epithet2):
        # a small change to the suffix presumably reflecting
        # a change in the gender of the genus. The lcp != epithet1/2
        # tests ensure the lcp is not just a substring of one of the
        # epithets (e.g., humilis and humi)
        return True
        
    return False


def canonical_species(species):
    """Get species name without suffix on generic or specific names.
    
    e.g., s__Alpha_A Beta_B -> s__Alpha Beta
    """
    
    assert species.startswith('s__')
    
    generic, specific = species[3:].split()

    if '_' in generic:
        generic = generic.rsplit('_', 1)[0]
    if '_' in specific:
        specific = specific.rsplit('_', 1)[0]
        
    return 's__{} {}'.format(generic, specific)
    

def canonical_taxon(taxon):
    """Get taxon name without suffix.
    
    Strips final alphabetic suffix from taxon names. Note, 
    the canonical name of a species in this context will be:
    s__Alpha_A Beta_B -> s__Alpha_A Beta
    
    That is, any suffix on the generic name is not removed. If this
    should be removed, used the method canonical_species().
    """

    if taxon.startswith('s__'):
        generic, specific = taxon[3:].split()

        if '_' in specific:
            specific = specific.rsplit('_', 1)[0]
            
        return 's__{} {}'.format(generic, specific)

    rank_prefix = ''
    if taxon[1:3] == '__':
        rank_prefix = taxon[0:3]
        taxon = taxon[3:]
        
    if '_' in taxon:
        taxon = taxon.rsplit('_', 1)[0]
   
    return rank_prefix + taxon


def taxon_suffix(taxon):
    """Return alphabetic suffix of taxon."""
    
    if taxon[1:3] == '__':
        taxon = taxon[3:]
    
    if '_' in taxon:
        suffix = taxon.rsplit('_', 1)[1]
        return suffix
        
    return None


def specific_epithet(species_name):
    """Set specific epithet."""

    if species_name == 's__':
        return ''
    
    generic, specific = species_name.split()
    
    return specific
    

def generic_name(species_name):
    """Set specific epithet."""

    if species_name == 's__':
        return ''
    
    generic, specific = species_name.split()
    
    return generic.replace('s__', '')
    
    
def is_latin_sp_epithet(sp_epithet):
    """Check if specific epithet is Latin (non-placeholder)."""

    return not is_placeholder_sp_epithet(sp_epithet)
    
    
def is_placeholder_sp_epithet(sp_epithet):
    """Check if species epithet is a placeholder."""
    
    assert '__' not in sp_epithet
    
    if sp_epithet == '':
        return True
    
    if any(c.isdigit() for c in sp_epithet):
        return True

    if any(c.isupper() for c in sp_epithet):
        return True
    
    return False


def is_placeholder_taxon(taxon):
    """Check if taxon name is a placeholder."""
    
    assert '__' in taxon # expect taxon name to have rank prefix
    
    test_taxon = taxon[3:].replace('[', '').replace(']', '')
    if test_taxon == '':
        return True
        
    if any(c.isdigit() for c in test_taxon):
        return True

    if any(c.isupper() for c in test_taxon[1:]):
        return True
    
    return False


def is_alphanumeric_taxon(taxon):
    """Check if taxon name is an alphanumeric placeholder.
    
    Example: g__9cdg1, f__UBA123, not not g__Prochlorococcus_A
    """
    
    return is_placeholder_taxon(taxon) and taxon_suffix(taxon) is None
    

def is_alphanumeric_sp_epithet(sp_epithet):
    """Check if specific name is an alphanumeric placeholder.
    
    Example: sp012345678, but not coli or coli_B
    """
    
    return is_placeholder_sp_epithet(sp_epithet) and taxon_suffix(sp_epithet) is None


def is_suffixed_taxon(taxon):
    """Check if taxon name is a suffixed placeholder.
    
    Example: g__Prochlorococcus_A, but not g__9cdg1 or f__UBA123
    """
    
    return is_placeholder_taxon(taxon) and taxon_suffix(taxon) is not None


def is_suffixed_sp_epithet(sp_epithet):
    """Check if specific name is a suffixed placeholder.
    
    Example: coli_B, but not coli or sp012345678
    """
    
    return is_placeholder_sp_epithet(sp_epithet) and taxon_suffix(sp_epithet) is not None


def taxon_type(sp_epithet):
    """Determine if taxon is Latin, suffixed Latin, or alphanumeric."""
    
    if is_alphanumeric_taxon(sp_epithet):
        return 'ALPHANUMERIC'
    elif is_suffixed_taxon(sp_epithet):
        return 'SUFFIXED_LATIN'
    
    return 'LATIN'
    
    
def specific_epithet_type(sp_epithet):
    """Determine if specific name is Latin, suffixed Latin, or alphanumeric."""
    
    if is_alphanumeric_sp_epithet(sp_epithet):
        return 'ALPHANUMERIC'
    elif is_suffixed_sp_epithet(sp_epithet):
        return 'SUFFIXED_LATIN'
    
    return 'LATIN'
    

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


def gtdb_merged_genera(prev_genomes, sp_priority_mngr, output_dir):
    """Determine genera previous merged in GTDB."""
    
    # get all GTDB genera defined in previous release
    prev_gtdb_genera = set()
    genera_with_type_species = set()
    for rid in prev_genomes.sp_clusters:
        canonical_genus = canonical_taxon(prev_genomes[rid].gtdb_taxa.genus)
        prev_gtdb_genera.add(canonical_genus)
        
        ncbi_genus = prev_genomes[rid].ncbi_taxa.genus
        gtdb_genus = prev_genomes[rid].gtdb_taxa.genus
        if ncbi_genus == gtdb_genus and prev_genomes[rid].is_gtdb_type_species():
            genera_with_type_species.add(ncbi_genus)
    
    # determine GTDB synonyms
    synonyms = defaultdict(lambda: defaultdict(list))
    for rid in prev_genomes.sp_clusters:
        ncbi_genus = prev_genomes[rid].ncbi_taxa.genus
        gtdb_genus = prev_genomes[rid].gtdb_taxa.genus
        
        if ncbi_genus in prev_gtdb_genera:
            # can't be a synonym if the genus exists in the GTDB taxonomy
            continue
            
        if is_placeholder_taxon(ncbi_genus):
            # doesn't make sense to consider a GTDB genus a 
            # synonym of a NCBI placeholder name
            continue
            
        synonyms[gtdb_genus][ncbi_genus].append(rid)
        
    # write out GTDB synonyms
    fout = open(os.path.join(output_dir, 'gtdb_prev_genera_synonyms.tsv'), 'w')
    fout.write('GTDB genus\tNCBI/synonym genus\tPriority\tSynonym priority\tPriority conflict')
    fout.write('\tNo. GTDB representatives\tGTDB genus has type species\tContains NCBI type species of genus\n')
    
    merged_genera = {}
    for gtdb_genus in synonyms:
        for ncbi_genus in synonyms[gtdb_genus]:
            merged_genera[gtdb_genus] = ncbi_genus
            
            has_type_species = False
            for rid in synonyms[gtdb_genus][ncbi_genus]:
                if prev_genomes[rid].is_gtdb_type_species():
                    has_type_species = True
                    break
                    
            ncbi_year = sp_priority_mngr.genus_priority_year(ncbi_genus)
            gtdb_year = sp_priority_mngr.genus_priority_year(gtdb_genus)
            
            fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        gtdb_genus,
                        ncbi_genus,
                        gtdb_year,
                        ncbi_year,
                        gtdb_year > ncbi_year,
                        len(synonyms[gtdb_genus][ncbi_genus]),
                        gtdb_genus in genera_with_type_species,
                        has_type_species))

    fout.close()
    
    return merged_genera
        
        
def sort_by_naming_priority(rids_of_interest,
                            sp_clusters,
                            prev_genomes, 
                            cur_genomes, 
                            mc_species):
    """Sort representatives by naming priority."""
    
    # group by naming priority
    manual_curation = []
    type_species = []
    type_strains = []
    binomial = []
    placeholder = []
    for rid in sp_clusters:
        if rid not in rids_of_interest:
            continue
            
        ncbi_genus = cur_genomes[rid].ncbi_taxa.genus
        ncbi_sp = cur_genomes[rid].ncbi_taxa.species
        gtdb_genus = cur_genomes[rid].gtdb_taxa.genus
        gtdb_sp = cur_genomes[rid].gtdb_taxa.species
        
        if rid in mc_species:
            manual_curation.append(rid)
        elif (cur_genomes[rid].is_gtdb_type_species()
                and gtdb_genus == ncbi_genus):
            type_species.append(rid)
        elif cur_genomes[rid].is_effective_type_strain():
            type_strains.append(rid)
        elif not is_placeholder_taxon(ncbi_sp) or not is_placeholder_taxon(gtdb_sp):
            binomial.append(rid)
        else:
            placeholder.append(rid)
            
    assert len(rids_of_interest) == len(manual_curation) + len(type_species) + len(type_strains) + len(binomial) + len(placeholder)
            
    # sort groups so previous GTDB representatives are processed first
    manual_curation_sorted = []
    type_species_sorted = []
    type_strains_sorted = []
    binomial_sorted = []
    placeholder_sorted = []
    for d, d_sorted in [(manual_curation, manual_curation_sorted), 
                        (type_species, type_species_sorted), 
                        (type_strains, type_strains_sorted),
                        (binomial, binomial_sorted),
                        (placeholder, placeholder_sorted)]:
        prev_reps = []
        new_reps = []
        for rid in d:
            if rid in prev_genomes.sp_clusters:
                prev_reps.append(rid)
            else:
                new_reps.append(rid)
            
        d_sorted.extend(prev_reps)
        d_sorted.extend(new_reps)

    return manual_curation_sorted, type_species_sorted, type_strains_sorted, binomial_sorted, placeholder_sorted
        

def genome_species_assignments(ncbi_taxonomy):
    """Get species assignment for each genome."""
    
    ncbi_species = binomial_species(ncbi_taxonomy)
    gid_to_species = {}
    for sp in ncbi_species:
        for gid in ncbi_species[sp]:
            gid_to_species[gid] = sp
            
    return gid_to_species
    
    
def generic_specific_names(species_name):
    """Get generic and specific names from species name."""
    
    species_name = species_name.replace('s__', '')
    species_name = species_name.replace('Candidatus ', '')
    generic, specific = species_name.split(' ')
    
    return generic, specific
    

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
