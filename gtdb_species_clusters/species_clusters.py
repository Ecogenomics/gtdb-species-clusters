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

import sys
import logging
from collections import defaultdict

from gtdblib.util.bio.accession import canonical_gid

from gtdb_species_clusters.genome_utils import (read_cur_new_updated,
                                                read_qc_file,
                                                read_gtdbtk_classifications)


class SpeciesClusters():
    """Genomes comprising GTDB-defined species clusters.

        Each species cluster is defined relative to a representative genome,
        contains 1 or more genomes (inclusive of the representative), and has
        a unique species name.
    """

    def __init__(self):
        """Initialization."""

        self.log = logging.getLogger('rich')

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

        if rid in self.species_names and self.species_names[rid] != sp_name:
            self.log.warning('GTDB representative {} appears to have two names: {} {} {}'.format(
                rid, self.species_names[rid], sp_name, gid))
            # sys.exit(-1)

        self.sp_clusters[rid].add(gid)
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
                    cids = [canonical_gid(cid.strip())
                            for cid in line_split[cluster_index].split(',')]
                    self.sp_clusters[rid].update([cid for cid in cids])
                    for cid in cids:
                        self.genome_rid[cid] = rid

    def create_expanded_clusters(self,
                                 prev_genomes,
                                 genomes_new_updated_file,
                                 qc_passed_file,
                                 gtdbtk_classify_file):
        """Expand species clusters to include genome in current GTDB release."""

        assert (not self.new_gids and not self.updated_gids)

        # read GTDB-Tk classifications for new and updated genomes
        gtdbtk_classifications = read_gtdbtk_classifications(
            gtdbtk_classify_file)
        self.log.info(
            f' - identified {len(gtdbtk_classifications):,} classifications')

        # get new and updated genomes in current GTDB release
        self.new_gids, self.updated_gids = read_cur_new_updated(
            genomes_new_updated_file)
        self.log.info(
            f' - identified {len(self.new_gids):,} new and {len(self.updated_gids):,} updated genomes')

        # get list of genomes passing QC
        gids_pass_qc = read_qc_file(qc_passed_file)
        new_pass_qc = len(self.new_gids.intersection(gids_pass_qc))
        updated_pass_qc = len(self.updated_gids.intersection(gids_pass_qc))
        self.log.info(
            f' - identified {new_pass_qc:,} new and {updated_pass_qc:,} updated genomes as passing QC')

        # create mapping between species and representatives
        original_sp_clusters = prev_genomes.sp_clusters
        orig_sp_rid_map = {sp: rid for rid,
                           sp in original_sp_clusters.species_names.items()}

        # create mapping between all genomes and species
        orig_gid_sp_map = {}
        for rid, cids in original_sp_clusters.sp_clusters.items():
            sp = original_sp_clusters.species_names[rid]
            for cid in cids:
                orig_gid_sp_map[cid] = sp

        # expand species clusters
        new_sp = 0
        for gid, taxa in gtdbtk_classifications.items():
            sp = taxa[6]
            if sp == 's__':
                new_sp += 1
                continue

            if sp not in orig_sp_rid_map:
                self.log.error(
                    f'GTDB-Tk results indicated a new species for {gid}: {sp}')
                # sys.exit(-1)

            orig_rid = orig_sp_rid_map[sp]
            if gid in self.new_gids:
                self.update_sp_cluster(orig_rid, gid, sp)
            elif gid in self.updated_gids:
                self.update_sp_cluster(orig_rid, gid, sp)

                orig_sp = orig_gid_sp_map[gid]
                if orig_sp != sp:
                    if prev_genomes[gid].is_gtdb_sp_rep():
                        # A GTDB species representative can changed to the point where
                        # it no longer clusters with its previous genome assembly and
                        # instead is placed into a different species cluster
                        # (e.g. GCA_005518205.1 was updated to GCA_005518205.2 in R207,
                        # and as a result clustered with the species cluster represented by
                        # GCA_005518235.1).
                        self.log.warning(
                            f'Updated GTDB representative {gid} reassigned from {orig_sp} to {sp} (manual inspection required to ensure this is properly resolved).')

            else:
                self.log.error(
                    f"Genome {gid} specified in GTDB-Tk results is neither 'new' or 'updated'")
                sys.exit(-1)

        self.log.info(
            f' - identified {new_sp:,} genomes not assigned to an existing GTDB species cluster')

        assert len(self.sp_clusters) == len(self.species_names)

    def updated_representatives(self, new_sp_clusters):
        """Determine the updated representative for each GTDB species cluster."""

        # get mapping between genomes and representative for
        # new GTDB species clusters
        new_gid_to_rid = {}
        for rid, cids in new_sp_clusters.items():
            for cid in cids:
                new_gid_to_rid[cid] = rid

        # get representative of species cluster containing
        # each of the previous representatives
        prev_to_new_rids = {}
        for rid in self.sp_clusters:
            prev_to_new_rids[rid] = new_gid_to_rid.get(rid, None)

        return prev_to_new_rids
