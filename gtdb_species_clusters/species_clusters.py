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
import logging
from collections import defaultdict

from biolib.taxonomy import Taxonomy

from gtdb_species_clusters.genome_utils import (canonical_gid,
                                                read_cur_new_updated,
                                                read_qc_file,
                                                read_gtdbtk_classifications)
                                            
                                            
class SpeciesClusters(object):
    """Genomes comprising GTDB-defined species clusters.
    
        Each species cluster is defined relative to a representative genome,
        contains 1 or more genomes (inclusive of the representative), and has
        a unique species name.
    """
    
    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        self.sp_clusters = defaultdict(set)
        self.species_names = {}
        self.genome_rid = {}
        
        self.new_gids = set()
        self.updated_gids = set()

    def __str__(self):
        """User-friendly string representation."""
        
        return '{{num_sp_clusters:{}}}'.format(len(self.sp_clusters))
        
    def __iter__(self):
        """Iterate over species cluster representative IDs."""
        
        for rid in self.sp_clusters:
            yield rid
            
    def items(self):
        """Iterate over representative IDs and clustered genome IDs."""
        
        for rid, cids in self.sp_clusters.items():
            yield (rid, cids)
    
    def clusters(self):
        """Iterate over representative IDs and clustered genome IDs."""
        
        for rid, cids in self.sp_clusters.items():
            yield (rid, cids)
            
    def species(self):
        """Iterate over representative IDs and species names"""
        
        for rid, sp_name in self.species_names.items():
            yield (rid, sp_name)

    def __getitem__(self, rid):
        """Get species cluster."""
        
        return self.sp_clusters[rid]
        
    def __contains__(self, rid):
        """Check if genome defined a species cluster."""
        
        return rid in self.sp_clusters

    def __len__(self):
        """Get number of species clusters."""
        
        return len(self.sp_clusters)

    def total_num_genomes(self):
        """Get total number of genomes across all species clusters."""
        
        return sum([len(cids) for cids in self.sp_clusters.values()])

    def get_species(self, rid):
        """Get species name for cluster."""
        
        return self.species_names[rid]
        
    def get_representative(self, gid):
        """Get representative of genome."""
        
        return self.genome_rid.get(gid, None)
        
    def update_sp_cluster(self, rid, gid, sp_name):
        """Update species cluster."""
        
        self.sp_clusters[rid].add(gid)
        assert(rid not in self.species_names 
                or self.species_names[rid] == sp_name)
                
        self.species_names[rid] = sp_name
        self.genome_rid[gid] = rid
        self.genome_rid[rid] = rid
        
    def load_from_sp_cluster_file(self, cluster_file):
        """Create species clusters from file."""

        with open(cluster_file) as f:
            headers = f.readline().strip().split('\t')
            
            sp_index = headers.index('GTDB species')
            rid_index = headers.index('Representative genome')
                
            num_clustered_index = headers.index('No. clustered genomes')
            cluster_index = headers.index('Clustered genomes')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                sp = line_split[sp_index]
                rid = canonical_gid(line_split[rid_index])
                self.genome_rid[rid] = rid
                
                self.sp_clusters[rid] = set([rid])
                self.species_names[rid] = sp

                num_clustered = int(line_split[num_clustered_index])
                if num_clustered > 0:
                    cids = [canonical_gid(cid.strip()) for cid in line_split[cluster_index].split(',')]
                    self.sp_clusters[rid].update([cid for cid in cids])
                    for cid in cids:
                        self.genome_rid[cid] = rid
                    
    def create_expanded_clusters(self, 
                                original_sp_clusters,
                                genomes_new_updated_file,
                                qc_passed_file,
                                gtdbtk_classify_file):
        """Expand species clusters to include genome in current GTDB release."""
        
        assert(not self.new_gids and not self.updated_gids)
        
        # read GTDB-Tk classifications for new and updated genomes
        gtdbtk_classifications = read_gtdbtk_classifications(gtdbtk_classify_file)
        self.logger.info(f' ... identified {len(gtdbtk_classifications):,} classifications.')
        
        # get new and updated genomes in current GTDB release
        self.new_gids, self.updated_gids = read_cur_new_updated(genomes_new_updated_file)
        self.logger.info(f' ... identified {len(self.new_gids):,} new and {len(self.updated_gids):,} updated genomes.')
        
        # get list of genomes passing QC
        gids_pass_qc = read_qc_file(qc_passed_file)
        new_pass_qc = len(self.new_gids.intersection(gids_pass_qc))
        updated_pass_qc = len(self.updated_gids.intersection(gids_pass_qc))
        self.logger.info(f' ... identified {new_pass_qc:,} new and {updated_pass_qc:,} updated genomes as passing QC.')

        # create mapping between species and representatives
        orig_sp_rid_map = {sp: rid for rid, sp in original_sp_clusters.species_names.items()}
        
        # create mapping between all genomes and species
        orig_gid_sp_map = {}
        for rid, cids in original_sp_clusters.sp_clusters.items():
            sp = original_sp_clusters.species_names[rid]
            for cid in cids:
                orig_gid_sp_map[cid] = sp
                                
        # expand species clusters
        failed_qc = 0
        new_sp = 0
        prev_genome_count = 0
        for gid, taxa in gtdbtk_classifications.items():
            if gid not in gids_pass_qc:
                # ***HACK: this should not be necessary, except GTDB-Tk was run external
                # of complete workflow for R95
                failed_qc += 1
                continue
                
            sp = taxa[6]
            if sp == 's__':
                new_sp += 1
                continue

            if sp not in orig_sp_rid_map:
                self.logger.error(f'GTDB-Tk results indicated a new species for {gid}: {sp}')
                sys.exit(-1)

            orig_rid = orig_sp_rid_map[sp]
            if gid in self.new_gids:
                self.update_sp_cluster(orig_rid, gid, sp)
            elif gid in self.updated_gids:
                self.update_sp_cluster(orig_rid, gid, sp)
                
                orig_sp = orig_gid_sp_map[gid]
                if orig_sp != sp:
                    self.logger.warning(f'Updated genomes {gid} reassigned from {orig_sp} to {sp}.')
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
        
        self.logger.info(f' ... identified {new_sp:,} genomes not assigned to an existing GTDB species cluster')

        assert len(self.sp_clusters) == len(self.species_names)
        