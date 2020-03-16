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
from itertools import takewhile
from collections import defaultdict, namedtuple

import biolib.seq_io as seq_io
from biolib.taxonomy import Taxonomy

from gtdb_species_clusters.genome_utils import canonical_gid


def longest_common_prefix(*s):
    """Find the longest common prefix.
    
    https://rosettacode.org/wiki/Longest_common_prefix#Python
    """
    
    return ''.join(a for a,b in takewhile(lambda x: x[0] == x[1], zip(min(s), max(s))))


def canonical_taxon(taxon):
    """Get taxon name without suffix."""

    if taxon.startswith('s__'):
        generic, specific = taxon.split()
        if '_' in specific:
            canonical_specific = specific.rsplit('_', 1)[0]
            return '{} {}'.format(generic, canonical_specific)
        else:
            return taxon

    rank_prefix = ''
    if taxon[1:3] == '__':
        rank_prefix = taxon[0:3]
        taxon = taxon[3:]
        
    if '_' in taxon:
        taxon = taxon.rsplit('_', 1)[0]
   
    return rank_prefix + taxon


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
    
    
def parse_synonyms(synonym_file):
    """Parse synonyms."""
    
    synonyms = {}
    with open(synonym_file) as f:
        headers = f.readline().strip().split('\t')
        
        ncbi_sp_index = headers.index('NCBI species')
        ncbi_synonym_index = headers.index('NCBI synonym')
        
        for line in f:
            line_split = line.strip().split('\t')
            
            ncbi_sp = line_split[ncbi_sp_index]
            ncbi_synonym = line_split[ncbi_synonym_index]
            synonyms[ncbi_synonym] = ncbi_sp
            
    return synonyms
    
    
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
