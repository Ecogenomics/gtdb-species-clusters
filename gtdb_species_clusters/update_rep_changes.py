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
from collections import defaultdict

from gtdb_species_clusters.genomes import Genomes


class RepChanges(object):
    """Identify species representatives that have changed from previous release."""

    def __init__(self, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
    def run(self, 
            prev_gtdb_metadata_file,
            cur_gtdb_metadata_file,
            genomes_new_updated_file,
            qc_passed_file,
            gtdbtk_classify_file,
            species_exception_file,
            genus_exception_file):
        """Identify species representatives that have changed from previous release."""
        
        # create previous and current GTDB genome sets
        self.logger.info('Creating previous GTDB genome set.')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                                species_exception_file,
                                                genus_exception_file)
        self.logger.info(f' ...previous genome set contains {len(prev_genomes):,} genomes.')
        self.logger.info(' ...previous genome set has {:,} species clusters spanning {:,} genomes.'.format(
                            len(prev_genomes.sp_clusters),
                            prev_genomes.sp_clusters.total_num_genomes()))

        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                                species_exception_file,
                                                genus_exception_file,
                                                create_sp_clusters=False)
        self.logger.info(f' ...current genome set contains {len(cur_genomes):,} genomes.')

        # get previous and current genomes from type strains
        self.logger.info('Determining genomes identified as being assembled from type strain.')
        prev_type_strain_gids = prev_genomes.gtdb_type_strain_genomes()
        cur_type_strain_gids = cur_genomes.gtdb_type_strain_genomes()
        new_type_strain_gids = cur_type_strain_gids - prev_type_strain_gids
        self.logger.info(' ...identified {:,} previous and {:,} current genomes from type strain.'.format(
                            len(prev_type_strain_gids),
                            len(cur_type_strain_gids)))
        self.logger.info(' ...{:,} type strain genomes are new to the current genome set.'.format(
                            len(new_type_strain_gids)))

        # expand previous GTDB species clusters to contain new genomes, 
        # and verify assignment of updated genomes
        self.logger.info('Expanding previous species clusters to contain new genomes based on GTDB-Tk classifications.')
        prev_genomes.sp_clusters.expand_sp_clusters(genomes_new_updated_file,
                                                    qc_passed_file,
                                                    gtdbtk_classify_file)
        self.logger.info(' ...expanded genome set has {:,} species clusters spanning {:,} genomes.'.format(
                            len(prev_genomes.sp_clusters),
                            prev_genomes.sp_clusters.total_num_genomes()))

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
        for prev_rid, prev_gtdb_sp in prev_genomes.sp_clusters.species():
            fout_summary.write(f'{prev_rid}\t{prev_sp}\t{len(prev_genomes.sp_clusters[prev_rid])}')
            if prev_rid in cur_genomes:
                if prev_rid in prev_genomes.sp_clusters.updated_gids:
                    updated_genome.add(prev_rid)
                    fout_summary.write('\tUPDATED')
                    prev_ncbi_accn = prev_genomes[prev_rid].ncbi_accn
                    cur_ncbi_accn = cur_genomes[prev_rid].ncbi_accn
                    assert(prev_ncbi_accn != cur_ncbi_accn)
                    fout_detailed.write((f'{prev_rid}\t{prev_gtdb_sp}\tGENOMIC_CHANGE:UPDATED\tNCBI accession updated from '
                                            f'{prev_genomes[prev_rid].ncbi_accn} to {cur_genomes[prev_rid].ncbi_accn}\n'))
                else:
                    unchanged_genome.add(prev_rid)
                    fout_summary.write('\tUNCHANGED')
                    
                prev_ncbi_sp = prev_genomes[prev_rid].ncbi_species()
                cur_ncbi_sp = cur_genomes[prev_rid].ncbi_species()
                if prev_genomes[prev_rid].ncbi_specific_epithet() == cur_genomes[prev_rid].ncbi_specific_epithet():
                    unchanged_sp.add(prev_rid)
                    fout_summary.write('\tUNCHANGED')
                else:
                    reassigned_sp.add(prev_rid)
                    fout_summary.write('\tREASSIGNED')
                    fout_detailed.write(f'{prev_rid}\t{prev_gtdb_sp}\tNCBI_SPECIES_CHANGE:REASSIGNED\tNCBI species reassigned from {prev_ncbi_sp} to {cur_ncbi_sp}\n')

                if prev_rid in prev_type_strain_gids and prev_rid in cur_type_strain_gids:
                    unchanged_type_strain.add(prev_rid)
                    fout_summary.write('\tUNCHANGED')
                elif prev_rid not in prev_type_strain_gids and prev_rid not in cur_type_strain_gids:
                    unchanged_type_strain.add(prev_rid)
                    fout_summary.write('\tUNCHANGED')
                elif prev_rid in prev_type_strain_gids and prev_rid not in cur_type_strain_gids:
                    lost_type_strain.add(prev_rid)
                    fout_summary.write('\tLOST')
                    fout_detailed.write(f'{prev_rid}\t{prev_gtdb_sp}\tTYPE_STRAIN_CHANGE:LOST\tNo longer considered a genome from type strain\n')
                elif prev_rid not in prev_type_strain_gids and prev_rid in cur_type_strain_gids:
                    gain_type_strain.add(prev_rid)
                    fout_summary.write('\tGAINED')
                    fout_detailed.write(f'{prev_rid}\t{prev_gtdb_sp}\tTYPE_STRAIN_CHANGE:GAINED\tNow considered a genome from type strain\n')
                else:
                    assert(False)
                    
                new_ts = new_type_strain_gids.intersection(prev_genomes.sp_clusters[prev_rid])

                if new_ts:
                    new_type_strain.add(prev_rid)
                    fout_detailed.write('{}\t{}\tNEW_TYPE_STRAINS::NEW\tSpecies cluster has {:,} new genomes from type strain: {}\n'.format(
                                            prev_rid,
                                            prev_gtdb_sp,
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
        
        num_prev_sp_clusters = len(prev_genomes.sp_clusters)
        num_rep_changes_perc = num_rep_changes*100.0/num_prev_sp_clusters
        self.logger.info(f' ... identified {num_rep_changes:,} ({num_rep_changes_perc:.1f}%) species with a change to the representative genome.')

        self.logger.info('Genomic changes:')
        unchanged_perc = len(unchanged_genome)*100.0 / num_prev_sp_clusters
        updated_perc = len(updated_genome)*100.0 / num_prev_sp_clusters
        lost_perc = len(lost_genome)*100.0 / num_prev_sp_clusters
        user_perc = len(user_genome)*100.0 / num_prev_sp_clusters
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
