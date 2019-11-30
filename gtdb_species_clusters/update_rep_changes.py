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
from copy import deepcopy
from collections import defaultdict

from gtdb_species_clusters.common import (specific_epithet,
                                            canonical_gid, 
                                            read_gtdb_accessions,
                                            read_gtdb_sp_clusters,
                                            read_gtdb_ncbi_taxonomy,
                                            read_gtdb_metadata,
                                            read_cur_new_updated,
                                            read_gtdbtk_classifications,
                                            read_qc_file,
                                            read_gtdb_taxonomy,
                                            find_new_type_strains,
                                            expand_sp_clusters)
                                            
from gtdb_species_clusters.type_genome_utils import gtdb_type_strain_of_species

from gtdb_species_clusters.genomes import Genomes

                                                        
class RepChanges(object):
    """Identify species representatives that have changed from previous release."""

    def __init__(self, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
    def run(self, 
            prev_gtdb_metadata_file,
            prev_sp_cluster_file,
            cur_gtdb_metadata_file,
            genomes_new_updated_file,
            qc_passed_file,
            gtdbtk_classify_file,
            species_exception_file,
            genus_exception_file):
        """Identify species representatives that have changed from previous release."""
        
        self.logger.info('Creating previous and current GTDB genome set.')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file)
        self.logger.info(f' ...previous genome set contains {len(prev_genomes):,} genomes.')
        
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file)
        self.logger.info(f' ...current genome set contains {len(cur_genomes):,} genomes.')

        # get previous and current genomes from type strains
        self.logger.info('Determining genomes identified as being assembled from type strain.')
        cur_gids_type_strain = prev_genomes.gtdb_type_strain_genomes()
        prev_gids_type_strain = gtdb_type_strain_genomes()
        self.logger.info(' ...identified {:,} previous and {:,} current genomes from type strain.'.format(
                            len(prev_gids_type_strain),
                            len(cur_gids_type_strain)))
        
        sys.exit(-1)
        
        # get current NCBI taxonomy
        self.logger.info('Reading current NCBI taxonomy from GTDB metadata file.')
        cur_ncbi_taxonomy, _ = read_gtdb_ncbi_taxonomy(cur_gtdb_metadata_file, 
                                                        species_exception_file,
                                                        genus_exception_file)
        
        # get previous GTDB species clusters
        self.logger.info('Reading previous GTDB species clusters.')
        prev_sp_clusters, prev_gtdb_species = read_gtdb_sp_clusters(prev_sp_cluster_file)
        self.logger.info(' ... identified {:,} species clusters spanning {:,} genomes.'.format(
                            len(prev_sp_clusters),
                            sum([len(cids) for cids in prev_sp_clusters.values()])))
                            
        # get previous NCBI taxonomy
        self.logger.info('Reading previous NCBI and GTDB taxonomies from GTDB metadata file.')
        prev_ncbi_taxonomy, _ = read_gtdb_ncbi_taxonomy(prev_gtdb_metadata_file, 
                                                        species_exception_file,
                                                        genus_exception_file)
                                                        
        prev_gtdb_taxonomy = read_gtdb_taxonomy(prev_gtdb_metadata_file)
                                                        
        # get new and updated genomes in current GTDB release
        self.logger.info('Reading new and updated genomes in current GTDB release.')
        cur_new, cur_updated = read_cur_new_updated(genomes_new_updated_file)
        self.logger.info(f' ... identified {len(cur_new):,} new and {len(cur_updated):,} updated genomes.')
        
        # get list of genomes passing QC
        self.logger.info('Reading genomes passing QC.')
        gids_pass_qc = read_qc_file(qc_passed_file)
        self.logger.info(f' ... identified {len(gids_pass_qc):,} genomes.')
        
        # read GTDB-Tk classifications for new and updated genomes
        self.logger.info('Reading GTDB-Tk classifications.')
        gtdbtk_classifications = read_gtdbtk_classifications(gtdbtk_classify_file)
        self.logger.info(f' ... identified {len(gtdbtk_classifications):,} classifications.')
        
        # expand previous GTDB species clusters to contain new genomes, 
        # and verify assignment of updated genomes
        self.logger.info('Expanding previous species clusters to contain new genomes based on GTDB-Tk classifications.')
        expanded_sp_clusters, sp_assignments = expand_sp_clusters(prev_sp_clusters, 
                                                                    prev_gtdb_species, 
                                                                    gtdbtk_classifications,
                                                                    cur_new,
                                                                    cur_updated,
                                                                    gids_pass_qc)
        self.logger.info(f' ... {sp_assignments:,} of {len(cur_new):,} new genomes had a GTDB-Tk species assignment.')
                            
        # determine status of each previous GTDB representative
        self.logger.info('Determining status of each previous GTDB representative.')
        
        fout_summary = open(os.path.join(self.output_dir, 'rep_change_summary.tsv'), 'w')
        fout_summary.write('Genome ID\tPrevious GTDB species\tNo. genomes in cluster')
        fout_summary.write('\tGENOMIC_CHANGE\tNCBI_SPECIES_CHANGE\tTYPE_STRAIN_CHANGE')
        fout_summary.write('\tNew type strains\tRepresentative changed\n')
        
        fout_detailed  = open(os.path.join(self.output_dir, 'rep_change_detailed.tsv'), 'w')
        fout_detailed.write('Genome ID\tPrevious GTDB species\tChange type\tChange\n')
        
        unchanged_genome = set()
        updated_genome = set()
        lost_genome = set()
        user_genome = set()
        unchanged_sp = set()
        reassigned_sp = set()
        unchanged_type_strain = set()
        lost_type_strain = set()
        gain_type_strain = set()
        new_type_strain = set()
        num_rep_changes = 0
        for prev_rid in prev_sp_clusters:
            prev_sp = prev_gtdb_species[prev_rid]
            if prev_rid.startswith('U'):
                user_genome.add(prev_rid)
                fout_summary.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                    prev_rid, prev_sp, 
                                    len(expanded_sp_clusters[prev_rid]),
                                    'N/A','N/A','N/A','N/A'))
                continue
                
            fout_summary.write(f'{prev_rid}\t{prev_sp}\t{len(expanded_sp_clusters[prev_rid])}')
            if prev_rid in cur_accns:
                if prev_rid in cur_updated:
                    updated_genome.add(prev_rid)
                    fout_summary.write('\tUPDATED')
                    fout_detailed.write(f'{prev_rid}\t{prev_sp}\tGENOMIC_CHANGE:UPDATED\tNCBI accession updated from {prev_accns[prev_rid]} to {cur_accns[prev_rid]}\n')
                else:
                    unchanged_genome.add(prev_rid)
                    fout_summary.write('\tUNCHANGED')
                    
                prev_sp = prev_ncbi_taxonomy[prev_rid][6]
                cur_sp = cur_ncbi_taxonomy[prev_rid][6]
                if specific_epithet(prev_sp) == specific_epithet(cur_sp):
                    unchanged_sp.add(prev_rid)
                    fout_summary.write('\tUNCHANGED')
                else:
                    reassigned_sp.add(prev_rid)
                    fout_summary.write('\tREASSIGNED')
                    fout_detailed.write(f'{prev_rid}\t{prev_sp}\tNCBI_SPECIES_CHANGE:REASSIGNED\tSpecies reassigned from {prev_sp} to {cur_sp}\n')

                if prev_rid in prev_gids_type_strain and prev_rid in cur_gids_type_strain:
                    unchanged_type_strain.add(prev_rid)
                    fout_summary.write('\tUNCHANGED')
                elif prev_rid not in prev_gids_type_strain and prev_rid not in cur_gids_type_strain:
                    unchanged_type_strain.add(prev_rid)
                    fout_summary.write('\tUNCHANGED')
                elif prev_rid in prev_gids_type_strain and prev_rid not in cur_gids_type_strain:
                    lost_type_strain.add(prev_rid)
                    fout_summary.write('\tLOST')
                    fout_detailed.write(f'{prev_rid}\t{prev_sp}\tTYPE_STRAIN_CHANGE:LOST\tNo longer considered a genome from type strain\n')
                elif prev_rid not in prev_gids_type_strain and prev_rid in cur_gids_type_strain:
                    gain_type_strain.add(prev_rid)
                    fout_summary.write('\tGAINED')
                    fout_detailed.write(f'{prev_rid}\t{prev_sp}\tTYPE_STRAIN_CHANGE:GAINED\tNow considered a genome from type strain\n')
                else:
                    assert(False)
                    
                new_ts = find_new_type_strains(prev_rid, 
                                                prev_gids_type_strain,
                                                cur_gids_type_strain,
                                                expanded_sp_clusters)
                if new_ts:
                    new_type_strain.add(prev_rid)
                    fout_detailed.write('{}\t{}\tNEW_TYPE_STRAINS::NEW\tSpecies cluster has {:,} new genomes from type strain: {}\n'.format(
                                            prev_rid,
                                            prev_sp,
                                            len(new_ts),
                                            ','.join(new_ts)))
                    
                fout_summary.write(f'\t{len(new_ts)}')
                
                if (prev_rid in unchanged_genome 
                    and prev_rid in unchanged_sp
                    and prev_rid in unchanged_type_strain):
                    fout_summary.write('\tNO')
                else:
                    fout_summary.write('\tYES')
                    num_rep_changes += 1
                
                fout_summary.write('\n')
            else:
                lost_genome.add(prev_rid)
                fout_summary.write('\t{}\t{}\t{}\t{}\t{}\n'.format('LOST','N/A','N/A','N/A', 'YES'))
                fout_detailed.write('\tGENOMIC_CHANGE:LOST\tGenome not present in current GTDB release\n')
                num_rep_changes += 1
                
                
        fout_summary.close()
        fout_detailed.close()
        
        num_rep_changes_perc = num_rep_changes*100.0/len(prev_sp_clusters)
        self.logger.info(f' ... identified {num_rep_changes:,} ({num_rep_changes_perc:.1f}%) species with a change to the representative genome.')

        self.logger.info('Genomic changes:')
        unchanged_perc = len(unchanged_genome)*100.0 / len(prev_sp_clusters)
        updated_perc = len(updated_genome)*100.0 / len(prev_sp_clusters)
        lost_perc = len(lost_genome)*100.0 / len(prev_sp_clusters)
        user_perc = len(user_genome)*100.0 / len(prev_sp_clusters)
        self.logger.info(f'  unchanged_genome: {len(unchanged_genome):,} ({unchanged_perc:.1f}%)')
        self.logger.info(f'  updated_genome: {len(updated_genome):,} ({updated_perc:.1f}%)')
        self.logger.info(f'  lost_genome: {len(lost_genome):,} ({lost_perc:.1f}%)')
        self.logger.info(f'  user_genome: {len(user_genome):,} ({user_perc:.1f}%)')
        
        self.logger.info('NCBI species assignment changes:')
        cur_sp_count = len(unchanged_genome) + len(updated_genome)
        unchanged_sp_perc = len(unchanged_sp)*100.0 / cur_sp_count
        reassigned_sp_perc = len(reassigned_sp)*100.0 / cur_sp_count
        self.logger.info(f'  unchanged_sp: {len(unchanged_sp):,} ({unchanged_sp_perc:.1f}%)')
        self.logger.info(f'  reassigned_sp: {len(reassigned_sp):,} ({reassigned_sp_perc:.1f}%)')
        
        self.logger.info('Genome from type strain changes:')
        prev_ts_count = len(unchanged_type_strain) + len(lost_type_strain)
        unchanged_type_strain_perc = len(unchanged_type_strain)*100.0 / prev_ts_count
        lost_type_strain_perc = len(lost_type_strain)*100.0 / prev_ts_count
        gain_type_strain_perc = len(gain_type_strain)*100.0 / prev_ts_count
        new_type_strain_perc = len(new_type_strain)*100.0 / prev_ts_count
        self.logger.info(f'  unchanged_type_strain: {len(unchanged_type_strain):,} ({unchanged_type_strain_perc:.1f}%)')
        self.logger.info(f'  lost_type_strain: {len(lost_type_strain):,} ({lost_type_strain_perc:.1f}%)')
        self.logger.info(f'  gain_type_strain: {len(gain_type_strain):,} ({gain_type_strain_perc:.1f}%)')
        self.logger.info(f'  new_type_strain: {len(new_type_strain):,} ({new_type_strain_perc:.1f}%)')
