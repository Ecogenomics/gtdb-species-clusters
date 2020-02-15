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
import pickle
from collections import defaultdict, Counter

from numpy import (mean as np_mean, std as np_std)

from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_name_manager import SpeciesNameManager
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters.type_genome_utils import read_clusters
from gtdb_species_clusters.taxon_utils import parse_synonyms, is_placeholder_taxon


class UpdateSummaryStats(object):
    """Summary statistics indicating changes to GTDB species clusters."""

    def __init__(self, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
    def _parse_updated_sp_reps(self, updated_sp_rep_file):
        """Determine GTDB species clusters with new representatives."""
        
        self.logger.info('Identifying updated GTDB representatives.')
        
        updated_rids = {}
        num_updated_rids = 0
        num_lost_rids = 0
        with open(updated_sp_rep_file) as f:
            header = f.readline().strip().split('\t')
            
            prev_rid_index = header.index('Previous representative ID')
            new_rid_index = header.index('New representative ID')
            
            for line in f:
                tokens = line.strip().split('\t')
                
                prev_rid = tokens[prev_rid_index]
                new_rid = tokens[new_rid_index]
                if new_rid == 'None':
                    new_rid = None
                    num_lost_rids += 1
                elif prev_rid != new_rid:
                    num_updated_rids += 1
                    
                updated_rids[prev_rid] = new_rid

        self.logger.info(f' ... identified {num_updated_rids:,} updated and {num_lost_rids:,} lost representatives.')
        
        return updated_rids

    def run(self, 
            updated_sp_rep_file,
            gtdb_clusters_file,
            prev_gtdb_metadata_file,
            cur_gtdb_metadata_file,
            uba_genome_paths,
            qc_passed_file,
            gtdbtk_classify_file,
            ncbi_genbank_assembly_file,
            untrustworthy_type_file,
            synonym_file,
            gtdb_type_strains_ledger):
        """Summary statistics indicating changes to GTDB species clusters."""

        # create previous and current GTDB genome sets
        self.logger.info('Creating previous GTDB genome set.')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                uba_genome_file=uba_genome_paths,
                                                ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                                untrustworthy_type_ledger=untrustworthy_type_file)
        self.logger.info(' ... previous genome set has {:,} species clusters spanning {:,} genomes.'.format(
                            len(prev_genomes.sp_clusters),
                            prev_genomes.sp_clusters.total_num_genomes()))

        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                create_sp_clusters=False,
                                                uba_genome_file=uba_genome_paths,
                                                qc_passed_file=qc_passed_file,
                                                ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                                untrustworthy_type_ledger=untrustworthy_type_file,
                                                gtdbtk_classify_file=gtdbtk_classify_file)
        self.logger.info(f' ... current genome set contains {len(cur_genomes):,} genomes.')
        
        # report changes in genome sets
        self.logger.info('Comparing previous and current genome sets.')
        prev_gids = set(prev_genomes)
        new_gids = set(cur_genomes)
        num_same_genomes = len(prev_gids.intersection(new_gids))
        num_lost_genomes = len(prev_gids - new_gids)
        num_new_genomes = len(new_gids - prev_gids)
        self.logger.info(f' ... identified {num_same_genomes:,} genomes as being present in both genome sets.')
        self.logger.info(f' ... identified {num_lost_genomes:,} genomes as being lost from the previous genome set.')
        self.logger.info(f' ... identified {num_new_genomes:,} genomes as being new to the current genome set.')
 
        # get changes to representatives of previous GTDB species clusters
        updated_rids = self._parse_updated_sp_reps(updated_sp_rep_file)
        
        # get new GTDB species clusters
        self.logger.info('Reading current GTDB clusters.')
        new_clusters, _ = read_clusters(gtdb_clusters_file)
        self.logger.info(' ... current genome set has {:,} species clusters spanning {:,} genomes.'.format(
                            len(new_clusters),
                            sum(len(cids) for cids in new_clusters.values())))

        new_rid_map = {}
        for rid, cids in new_clusters.items():
            for cid in cids:
                new_rid_map[cid] = rid
                
        # UBA genome sanity check
        prev_uba_count = 0
        for gid in prev_genomes:
            if gid.startswith('UBA'):
                prev_uba_count += 1

        cur_uba_count = 0
        for gid in cur_genomes:
            if gid.startswith('UBA'):
                cur_uba_count += 1
                
        new_uba_count = 0
        for rid, cids in new_clusters.items():
            for cid in cids:
                if cid.startswith('UBA'):
                    new_uba_count += 1
                    
        self.logger.info(f'Verified all genome / cluster sets contain the same number of UBA genomes: {prev_uba_count:,}')
        assert prev_uba_count == cur_uba_count == new_uba_count

        # tabulate changes in GTDB species clusters
        self.logger.info('Calculating statistics of GTDB species clusters.')
        
        fout = open(os.path.join(self.output_dir, 'gtdb_sp_clusters_change_stats.tsv'), 'w')
        fout.write('Previous representative\tPrevious name\tNew representative\tNew name\tRepresentative status\tName status')
        fout.write('\tNo. previous genomes\tNo. current genomes\tNo. same\tNo. lost\tNo. new\tNo. migrated in\tNo. migrated out\tNote\n')
        
        rep_lost_count = 0
        rep_changed_count = 0
        rep_unchanged_count = 0
        rep_merged_count = 0
        
        name_lost_count = 0
        name_changed_count = 0
        name_unchanged_count = 0
        name_merged_count = 0
        
        prev_cluster_ids = set()
        total_num_same = 0
        total_num_lost = 0
        total_num_new = 0
        total_num_migrated_in = 0
        total_num_migrated_out = 0
        for prev_rid, prev_cids in prev_genomes.sp_clusters.items():
            prev_gtdb_sp = prev_genomes[prev_rid].gtdb_taxa.species
            
            new_rid = updated_rids[prev_rid]
            prev_cluster_ids.add(new_rid)
            note = ''
            if new_rid is None:
                new_rid = 'none'
                new_sp = 'none'
                rep_status = 'LOST'
                name_status = 'LOST' # what does this mean; presumable a species name can be recycled elsewhere!
                
                new_cluster = set()
                
                rep_lost_count += 1
                name_lost_count += 1
            elif new_rid not in new_clusters:
                # representative must have been merged when selecting
                # representatives for NCBI species
                merged_rid = new_rid_map[new_rid]
                merged_sp = cur_genomes[merged_rid].gtdb_taxa.species
                note = 'merged with {} with representative {}'.format(merged_sp, merged_rid)
                
                new_rid = 'none'
                rep_status = 'MERGED'
                name_status = 'MERGED'
                
                new_cluster = set()
                
                rep_merged_count += 1
                name_merged_count += 1
            else:
                new_gtdb_sp = cur_genomes[new_rid].gtdb_taxa.species
                new_cluster = new_clusters[new_rid]
                
                if prev_rid == new_rid:
                    rep_status = 'UNCHANGED'
                    rep_unchanged_count += 1
                else:
                    rep_status = 'CHANGED'
                    rep_changed_count += 1
                    
                
                if prev_gtdb_sp == new_gtdb_sp:
                    name_status = 'UNCHANGED'
                    name_unchanged_count += 1
                else:
                    name_status = 'CHANGED'
                    name_changed_count += 1
            
            fout.write('{}\t{}\t{}\t{}\t{}\t{}'.format(
                        prev_rid,
                        prev_gtdb_sp,
                        new_rid,
                        new_gtdb_sp,
                        rep_status,
                        name_status))
            
            num_same = len(new_cluster.intersection(prev_cids))
            num_lost = len(prev_cids - new_gids)
            num_new = len(new_cluster - prev_gids)
            num_migrated_in = len((new_cluster - prev_cids).intersection(prev_gids))
            num_migrated_out = len((prev_cids - new_cluster).intersection(new_gids))
            assert len(new_cluster) == len(prev_cids) - num_lost + num_new + num_migrated_in - num_migrated_out
            assert len(prev_cids) == num_same + num_lost + num_migrated_out
            
            fout.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        len(prev_cids),
                        len(new_cluster),
                        num_same,
                        num_lost,
                        num_new,
                        num_migrated_in,
                        num_migrated_out,
                        note))
                        
            total_num_same += num_same
            total_num_lost += num_lost
            total_num_new += num_new
            total_num_migrated_in += num_migrated_in
            total_num_migrated_out += num_migrated_out

        # add in new GTDB species clusters
        new_cluster_count = 0
        for new_rid in new_clusters:
            if new_rid in prev_cluster_ids:
                continue
                
            new_gtdb_sp = cur_genomes[new_rid].gtdb_taxa.species
            rep_status = 'NEW'
            name_status = 'NEW'
            new_cluster_count += 1
            
            fout.write('{}\t{}\t{}\t{}\t{}\t{}'.format(
                        'n/a',
                        'n/a',
                        new_rid,
                        new_gtdb_sp,
                        rep_status,
                        name_status))
                        
            num_new = len(new_clusters[new_rid] - prev_gids)
            num_migrated_in = len(new_clusters[new_rid].intersection(prev_gids))
            assert len(new_clusters[new_rid]) == num_new + num_migrated_in
            fout.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        0,
                        len(new_clusters[new_rid]),
                        0,
                        0,
                        num_new,
                        num_migrated_in,
                        0,
                        ''))
                        
            total_num_new += num_new
            total_num_migrated_in += num_migrated_in
            
        # report genome statistics
        num_union = len(new_gids.union(prev_gids))
        assert len(new_gids.union(prev_gids)) == total_num_same + total_num_lost + total_num_new + total_num_migrated_in
        assert total_num_migrated_in == total_num_migrated_out
        self.logger.info(f'There were {len(prev_gids):,} genomes in the previous genome sets.')
        self.logger.info(' ... identified {:,} ({:.2f}%) genomes that were assigned to same species cluster.'.format(
                            total_num_same,
                            total_num_same*100.0/len(prev_gids)))
        self.logger.info(' ... identified {:,} ({:.2f}%) genomes that were lost from the species cluster.'.format(
                            total_num_lost,
                            total_num_lost*100.0/len(prev_gids)))
        self.logger.info(' ... identified {:,} ({:.2f}%) genomes that migrated between species cluster.'.format(
                            total_num_migrated_in,
                            total_num_migrated_in*100.0/len(prev_gids)))
        self.logger.info(' ... identified {:,} new genomes which is a {:.2f}% increase.'.format(
                            total_num_new,
                            len(new_gids)*100.0/len(prev_gids) - 100))

        # report representative statistics
        assert len(new_clusters) == len(prev_genomes.sp_clusters) + new_cluster_count - rep_lost_count  - rep_merged_count 
        self.logger.info(f'There are {len(new_clusters):,} total GTDB species representatives.')
        self.logger.info(' ... identified {:,} ({:.2f}%) unchanged representatives.'.format(
                            rep_unchanged_count,
                            rep_unchanged_count*100.0/len(prev_genomes.sp_clusters)))
        self.logger.info(' ... identified {:,} ({:.2f}%) changed representatives.'.format(
                            rep_changed_count,
                            rep_changed_count*100.0/len(prev_genomes.sp_clusters)))
        self.logger.info(' ... identified {:,} ({:.2f}%) lost representatives.'.format(
                            rep_lost_count,
                            rep_lost_count*100.0/len(prev_genomes.sp_clusters)))
        self.logger.info(' ... identified {:,} ({:.2f}%) merged representatives.'.format(
                            rep_merged_count,
                            rep_merged_count*100.0/len(prev_genomes.sp_clusters)))
        self.logger.info(' ... identified {:,} new representatives which is a {:.2f}% increase.'.format(
                            new_cluster_count,
                            len(new_clusters)*100.0/len(prev_genomes.sp_clusters) - 100))

        self.logger.info(' ... identified {:,} ({:.2f}%) cluster names.'.format(
                            name_unchanged_count,
                            name_unchanged_count*100.0/len(prev_genomes.sp_clusters)))
        self.logger.info(' ... identified {:,} ({:.2f}%) changed cluster names.'.format(
                            name_changed_count,
                            name_changed_count*100.0/len(prev_genomes.sp_clusters)))
        self.logger.info(' ... identified {:,} ({:.2f}%) lost cluster names.'.format(
                            name_lost_count,
                            name_lost_count*100.0/len(prev_genomes.sp_clusters)))
        self.logger.info(' ... identified {:,} ({:.2f}%) merged cluster names.'.format(
                            name_merged_count,
                            name_merged_count*100.0/len(prev_genomes.sp_clusters)))
        self.logger.info(' ... identified {:,} ({:.2f}%) new cluster names.'.format(
                            new_cluster_count,
                            new_cluster_count*100.0/len(prev_genomes.sp_clusters)))