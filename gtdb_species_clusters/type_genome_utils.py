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
import ntpath
from collections import defaultdict, namedtuple

from numpy import (mean as np_mean)

from gtdb_species_clusters.common import read_gtdb_metadata


NCBI_TYPE_SPECIES = set(['assembly from type material', 
                        'assembly from neotype material',
                        'assembly designated as neotype'])
NCBI_PROXYTYPE = set(['assembly from proxytype material'])
NCBI_TYPE_SUBSP = set(['assembly from synonym type material'])

GTDB_TYPE_SPECIES = set(['type strain of species', 'type strain of neotype'])
GTDB_TYPE_SUBSPECIES = set(['type strain of subspecies', 'type strain of heterotypic synonym'])
GTDB_NOT_TYPE_MATERIAL = set(['not type material'])

GenomeRadius = namedtuple('GenomeRadius', 'ani af neighbour_gid')


def ncbi_species(unfiltered_ncbi_taxonomy):
    """Get species designation for unfiltered NCBI taxonomy string."""
    
    ncbi_unfiltered_taxa = [t.strip() for t in unfiltered_ncbi_taxonomy.split(';')]
    ncbi_species = None
    ncbi_subspecies = None
    for taxon in ncbi_unfiltered_taxa:
        if taxon.startswith('s__'):
            ncbi_species = taxon[3:]
        elif taxon.startswith('sb__'):
            ncbi_subspecies = taxon[4:]
            
            # fix odd designation impacting less than a dozen genomes 
            ncbi_subspecies = ncbi_subspecies.replace(' pv. ', ' subsp. ') 

    if ncbi_subspecies:
        if 'subsp.' not in ncbi_subspecies:
            print("NCBI subspecies name without 'subsp.' definition: %s" % ncbi_subspecies)
            
        return ncbi_subspecies
    
    if ncbi_species:
        return ncbi_species
        
    return None
    
    
def check_ncbi_subsp(unfiltered_ncbi_taxonomy):
    """Determine if genome is the 'type strain of species' or 'type strain of subspecies'."""
    
    sp_name = ncbi_species(unfiltered_ncbi_taxonomy)
    if 'subsp.' not in sp_name:
        return False
    else:
        tokens = sp_name.split()
        subsp_index = tokens.index('subsp.')
        if tokens[subsp_index - 1] == tokens[subsp_index + 1]:
            return False

    return True
    

def parse_canonical_sp(sp):
    """Get canonical binomial species name."""
    
    sp = sp.replace('Candidatus ', '')
    sp = ' '.join(sp.split()[0:2]).strip()
    
    return sp
        

def symmetric_ani(ani_af, gid1, gid2):
    """Calculate symmetric ANI statistics between genomes."""
    
    if gid1 == gid2:
        return 100.0, 1.0
    
    if (gid1 not in ani_af
        or gid2 not in ani_af 
        or gid1 not in ani_af[gid2]
        or gid2 not in ani_af[gid1]):
        return 0.0, 0.0
    
    cur_ani, cur_af = ani_af[gid1][gid2]
    rev_ani, rev_af = ani_af[gid2][gid1]
    
    # ANI should be the larger of the two values as this
    # is the most conservative circumscription and reduces the
    # change of creating polyphyletic species clusters
    ani = max(rev_ani, cur_ani)
    
    # AF should be the larger of the two values in order to 
    # accomodate incomplete and contaminated genomes
    af = max(rev_af, cur_af)
    
    return ani, af
    
    
def is_isolate(gid, quality_metadata):
    """Check if genome is an isolate."""
    
    m = quality_metadata[gid]
    if m.ncbi_genome_category:
        if ('metagenome' in m.ncbi_genome_category.lower()
            or 'environmental' in m.ncbi_genome_category.lower()
            or 'single cell' in m.ncbi_genome_category.lower()):
            return False
            
    return True
    
    
def is_type_strain(gid, quality_metadata):
    return quality_metadata[gid].gtdb_type_designation in GTDB_TYPE_SPECIES
    
    
def is_ncbi_type_strain(gid, quality_metadata):
    m = quality_metadata[gid]
    if m.ncbi_type_material_designation:
        if m.ncbi_type_material_designation.lower() in NCBI_TYPE_SPECIES:
            if not check_ncbi_subsp(m.ncbi_taxonomy_unfiltered):
                return True
            else:
                # genome is marked as 'assembled from type material', but
                # is a subspecies according to the NCBI taxonomy so is
                # the type strain of a subspecies
                # (e.g. Alpha beta subsp. gamma)
                return False
                
    return False
    
    
def is_complete_genome(gid, quality_metadata):
    m = quality_metadata[gid]
    return (m.ncbi_assembly_level 
            and m.ncbi_assembly_level.lower() in ['complete genome', 'chromosome']
            and m.ncbi_genome_representation
            and m.ncbi_genome_representation.lower() == 'full'
            and m.scaffold_count == m.ncbi_molecule_count
            and m.ncbi_unspanned_gaps == 0
            and m.ncbi_spanned_gaps <= 10
            and m.ambiguous_bases <= 1e4
            and m.total_gap_length <= 1e4
            and m.ssu_count >= 1)


def quality_score_update(gids, quality_metadata):
    """"Calculate quality score of genomes for updating representatives."""

    score = {}
    for gid in gids:
        metadata = quality_metadata[gid]
        
        # set base quality so genomes have the following priority order:
        #  type strain genome
        #  NCBI assembled from type material
        #  NCBI reference or representative genome
        #  type subspecies genome
        q = 0
        if is_type_strain(gid, quality_metadata):
            q = 1e4
        
        if is_ncbi_type_strain(gid, quality_metadata):
            q = 1e3

        # check if genome appears to complete consist of only an unspanned
        # chromosome and unspanned plasmids and thus should be considered
        # very high quality
        if is_complete_genome(gid, quality_metadata):
            q += 100
            
        q += metadata.checkm_completeness - 5*metadata.checkm_contamination
        q -= 5*float(metadata.contig_count)/100
        q -= 5*float(metadata.ambiguous_bases)/1e5
        
        if not is_isolate(gid, quality_metadata):
            if 'single cell' in metadata.ncbi_genome_category.lower():
                q -= 100 # genome is a SAG
            else: 
                q -= 200 # genome is a MAG
        
        # check for near-complete 16S rRNA gene
        gtdb_domain = metadata.gtdb_taxonomy[0]
        min_ssu_len = 1200
        if gtdb_domain == 'd__Archaea':
            min_ssu_len = 900
            
        if metadata.ssu_length and metadata.ssu_length >= min_ssu_len:
            q += 10
            
        score[gid] = q
        
    return score
    
    
def quality_score(gids, quality_metadata):
    """"Calculate quality score for genomes."""

    score = {}
    for gid in gids:
        metadata = quality_metadata[gid]
    
        # check if genome appears to complete consist of only an unspanned
        # chromosome and unspanned plasmids and thus should be considered
        # very high quality
        if (metadata.ncbi_assembly_level 
                and metadata.ncbi_assembly_level.lower() in ['complete genome', 'chromosome']
                and metadata.ncbi_genome_representation
                and metadata.ncbi_genome_representation.lower() == 'full'
                and metadata.scaffold_count == metadata.ncbi_molecule_count
                and metadata.ncbi_unspanned_gaps == 0
                and metadata.ncbi_spanned_gaps <= 10
                and metadata.ambiguous_bases <= 1e4
                and metadata.total_gap_length <= 1e4
                and metadata.ssu_count >= 1):
            q = 100
        else:
            q = 0
            
        q += metadata.checkm_completeness - 5*metadata.checkm_contamination
        q += 200*(metadata.ncbi_type_material_designation is not None 
                    and metadata.ncbi_type_material_designation.lower() in NCBI_TYPE_SPECIES)
        q += 10*((metadata.ncbi_type_material_designation is not None 
                    and metadata.ncbi_type_material_designation.lower() in NCBI_PROXYTYPE)
                  or (metadata.ncbi_refseq_category is not None 
                    and ('representative' in metadata.ncbi_refseq_category.lower()
                        or 'reference' in metadata.ncbi_refseq_category.lower())))

        q -= 5*float(metadata.contig_count)/100
        q -= 5*float(metadata.ambiguous_bases)/1e5
        
        if metadata.ncbi_genome_category:
            if ('metagenome' in metadata.ncbi_genome_category.lower()
                or 'environmental' in metadata.ncbi_genome_category.lower()): # environmental check added Nov. 4, 2019
                q -= 200
            if 'single cell' in metadata.ncbi_genome_category.lower():
                q -= 100
        
        # check for near-complete 16S rRNA gene
        gtdb_domain = metadata.gtdb_taxonomy[0]
        min_ssu_len = 1200
        if gtdb_domain == 'd__Archaea':
            min_ssu_len = 900
            
        if metadata.ssu_length and metadata.ssu_length >= min_ssu_len:
            q += 10
            
        score[gid] = q
        
    return score


def pass_qc(qc, 
            marker_perc,
            min_comp,
            max_cont,
            min_quality,
            sh_exception,
            min_perc_markers,
            max_contigs,
            min_N50,
            max_ambiguous,
            failed_tests):
    """Check if genome passes QC."""
    
    failed = False
    if qc.checkm_completeness < min_comp:
        failed_tests['comp'] += 1
        failed = True
    
    if qc.checkm_strain_heterogeneity_100 >= sh_exception:
        if qc.checkm_contamination > 20:
            failed_tests['cont'] += 1
            failed = True
        q = qc.checkm_completeness - 5*qc.checkm_contamination*(1.0 - qc.checkm_strain_heterogeneity_100/100.0)
        if q < min_quality:
            failed_tests['qual'] += 1
            failed = True
    else:
        if qc.checkm_contamination > max_cont:
            failed_tests['cont'] += 1
            failed = True
        q = qc.checkm_completeness - 5*qc.checkm_contamination
        if q < min_quality:
            failed_tests['qual'] += 1
            failed = True
            
    if marker_perc < min_perc_markers:
        failed_tests['marker_perc'] += 1
        failed = True
            
    if qc.contig_count > max_contigs:
        failed_tests['contig_count'] += 1
        failed = True
    if qc.n50_contigs < min_N50:
        failed_tests['N50'] += 1
        failed = True
    if qc.ambiguous_bases > max_ambiguous:
        failed_tests['ambig'] += 1
        failed = True
    
    return not failed
    
    
def ncbi_type_strain_of_species(type_metadata):
    """Determine genomes considered type strain of species at NCBI."""
    
    type_gids = set()
    for gid in type_metadata:
        if type_metadata[gid].ncbi_type_material_designation in NCBI_TYPE_SPECIES:
            type_gids.add(gid)
            
    return type_gids
    
    
def gtdb_type_strain_of_species(type_metadata):
    """Determine genomes considered type strain of species by GTDB."""
    
    type_gids = set()
    for gid in type_metadata:
        if type_metadata[gid].gtdb_type_designation in GTDB_TYPE_SPECIES:
            type_gids.add(gid)
            
    return type_gids
            

def write_clusters(clusters, type_radius, species, out_file):
    """Write out clustering information."""

    fout = open(out_file, 'w')
    fout.write('NCBI species\tType genome')
    fout.write('\tClosest species\tClosest type genome\tANI radius\tAF closest')
    fout.write('\tNo. clustered genomes\tMean ANI\tMin ANI\tMean AF\tMin AF\tClustered genomes\n')
    for gid in sorted(clusters, key=lambda x: len(clusters[x]), reverse=True):
        if len(clusters[gid]):
            mean_ani = '%.2f' % np_mean([d.ani for d in clusters[gid]])
            min_ani = '%.2f' % min([d.ani for d in clusters[gid]])
            mean_af = '%.2f' % np_mean([d.af for d in clusters[gid]])
            min_af = '%.2f' % min([d.af for d in clusters[gid]])
        else:
            mean_ani = min_ani = mean_af = min_af = 'N/A'
        fout.write('%s\t%s' % (
                        species.get(gid, 'unclassified'), 
                        gid))
                        
        ani, af, closest_gid = type_radius[gid]
        if not af:
            af = 0
            
        if not closest_gid or closest_gid == 'N/A':
            closest_gid = 'N/A'
            closest_sp = 'N/A'
        else:
            closest_sp = species[closest_gid]
        
        fout.write('\t%s\t%s\t%.2f\t%.2f' % (closest_sp,
                                                closest_gid,
                                                ani,
                                                af))
                        
        fout.write('\t%d\t%s\t%s\t%s\t%s\t%s\n' % (
                        len(clusters[gid]),
                        mean_ani, min_ani,
                        mean_af, min_af,
                        ','.join([d.gid for d in clusters[gid]])))
    fout.close()
    
    
def write_type_radius(type_radius, species, out_file):
    """Write out ANI radius for each type genomes."""

    fout = open(out_file, 'w')
    fout.write('NCBI species\tType genome\tANI\tAF\tClosest species\tClosest type genome\n')
    
    for gid in type_radius:
        ani, af, neighbour_gid = type_radius[gid]
        if not af:
            af = 0
            
        if not neighbour_gid or neighbour_gid == 'N/A':
            neighbour_gid = 'N/A'
            neighbour_sp = 'N/A'
        else:
            neighbour_sp = species[neighbour_gid]
        
        fout.write('%s\t%s\t%.2f\t%.2f\t%s\t%s\n' % (species[gid],
                                                        gid,
                                                        ani,
                                                        af,
                                                        neighbour_sp,
                                                        neighbour_gid))
    fout.close()


def read_quality_metadata(metadata_file):
    """Read statistics needed to determine genome quality."""
    
    return read_gtdb_metadata(metadata_file, ['gtdb_taxonomy',
                                                'checkm_completeness',
                                                'checkm_contamination',
                                                'checkm_strain_heterogeneity_100',
                                                'genome_size',
                                                'contig_count',
                                                'n50_contigs',
                                                'scaffold_count',
                                                'ambiguous_bases',
                                                'total_gap_length',
                                                'ssu_count',
                                                'ssu_length',
                                                'mimag_high_quality',
                                                'gtdb_type_designation',
                                                'ncbi_assembly_level',
                                                'ncbi_genome_representation',
                                                'ncbi_refseq_category',
                                                'ncbi_type_material_designation',
                                                'ncbi_molecule_count',
                                                'ncbi_unspanned_gaps',
                                                'ncbi_spanned_gaps',
                                                'ncbi_genome_category'])
                                                

def read_clusters(cluster_file):
    """Read cluster file."""
        
    clusters = defaultdict(list)
    species = {}
    rep_radius = {}
    with open(cluster_file) as f:
        headers = f.readline().strip().split('\t')
        
        type_sp_index = headers.index('NCBI species')
        type_genome_index = headers.index('Type genome')
            
        num_clustered_index = headers.index('No. clustered genomes')
        clustered_genomes_index = headers.index('Clustered genomes')
        closest_sp_index = headers.index('Closest species')
        closest_gid_index = headers.index('Closest type genome')
        closest_ani_index = headers.index('ANI radius')
        closest_af_index = headers.index('AF closest')
        
        for line in f:
            line_split = line.strip().split('\t')
            
            sp = line_split[type_sp_index]
            rid = line_split[type_genome_index]
            species[rid] = sp

            num_clustered = int(line_split[num_clustered_index])
            if num_clustered > 0:
                clusters[rid] = [g.strip() for g in line_split[clustered_genomes_index].split(',')]
            else:
                clusters[rid] = []
                
            closest_ani = float(line_split[closest_ani_index])
            closest_af = float(line_split[closest_af_index])
            closest_gid = line_split[closest_gid_index]
            rep_radius[rid] = GenomeRadius(ani = closest_ani, 
                                             af = closest_af,
                                             neighbour_gid = closest_gid)
                    
    return clusters, species, rep_radius
    
    

