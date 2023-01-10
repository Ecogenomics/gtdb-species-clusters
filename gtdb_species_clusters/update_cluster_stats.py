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
import logging
from collections import defaultdict

from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.type_genome_utils import read_clusters


class UpdateClusterStats(object):
    """Summary statistics indicating changes to GTDB species cluster membership."""

    def __init__(self, output_dir):
        """Initialization."""

        self.output_dir = output_dir
        self.log = logging.getLogger('rich')

    def _parse_updated_sp_reps(self, updated_sp_rep_file):
        """Determine GTDB species clusters with new representatives."""

        self.log.info('Identifying updated GTDB representatives.')

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

        self.log.info(
            f' - identified {num_updated_rids:,} updated and {num_lost_rids:,} lost representatives.')

        return updated_rids

    def run(self,
            gtdb_clusters_file,
            prev_gtdb_metadata_file,
            cur_gtdb_metadata_file,
            qc_passed_file,
            ncbi_genbank_assembly_file,
            untrustworthy_type_file,
            gtdb_type_strains_ledger,
            ncbi_env_bioproject_ledger):
        """Summary statistics indicating changes to GTDB species cluster membership."""

        # create previous and current GTDB genome sets
        self.log.info('Creating previous GTDB genome set.')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                             gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                             ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                             untrustworthy_type_ledger=untrustworthy_type_file,
                                             ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)
        self.log.info(' - previous genome set has {:,} species clusters spanning {:,} genomes.'.format(
            len(prev_genomes.sp_clusters),
            prev_genomes.sp_clusters.total_num_genomes()))

        self.log.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                            gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                            create_sp_clusters=False,
                                            qc_passed_file=qc_passed_file,
                                            ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                            untrustworthy_type_ledger=untrustworthy_type_file,
                                            ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)

        # report changes in genome sets
        self.log.info('Comparing previous and current genome sets:')
        prev_gids = set(prev_genomes)
        new_gids = set(cur_genomes)
        num_same_genomes = len(prev_gids.intersection(new_gids))
        num_lost_genomes = len(prev_gids - new_gids)
        num_new_genomes = len(new_gids - prev_gids)
        self.log.info(
            f' - identified {num_same_genomes:,} genomes as being present in both genome sets')
        self.log.info(
            f' - identified {num_lost_genomes:,} genomes as being lost from the previous genome set')
        self.log.info(
            f' - identified {num_new_genomes:,} genomes as being new to the current genome set')

        # get new GTDB species clusters
        self.log.info('Reading current GTDB clusters:')
        new_clusters, _ = read_clusters(gtdb_clusters_file)
        self.log.info(' - current genome set has {:,} species clusters spanning {:,} genomes'.format(
            len(new_clusters),
            sum(len(cids) for cids in new_clusters.values())))

        new_rid_map = {}
        for rid, cids in new_clusters.items():
            for cid in cids:
                new_rid_map[cid] = rid

        # get mapping of previous GTDB representatives to new GTDB species clusters
        self.log.info(
            'Mapping previous GTDB representatives to new representatives.')
        prev_to_new_rid = prev_genomes.sp_clusters.updated_representatives(
            new_clusters)
        self.log.info(
            ' - mapped {:,} previous representatives.'.format(len(prev_to_new_rid)))

        new_to_prev_rids = defaultdict(list)
        for prev_rid, new_rid in prev_to_new_rid.items():
            new_to_prev_rids[new_rid].append(prev_rid)

        # tabulate changes in GTDB species clusters
        self.log.info('Calculating statistics of GTDB species clusters.')

        fout = open(os.path.join(self.output_dir,
                                 'gtdb_sp_clusters_change_stats.tsv'), 'w')
        fout.write(
            'New representative\tPrevious representative(s)\tPrevious name(s)\tRepresentative status')
        fout.write(
            '\tNo. previous genomes\tNo. current genomes\tNo. same\tNo. lost\tNo. new\tNo. migrated in\tNo. migrated out\n')

        rep_lost_count = 0
        rep_changed_count = 0
        rep_unchanged_count = 0
        rep_merger_count = 0

        prev_cluster_ids = set()
        total_num_same = 0
        total_num_lost = 0
        total_num_new = 0
        total_num_migrated_in = 0
        total_num_migrated_out = 0

        for new_rid, prev_rids in new_to_prev_rids.items():
            prev_cluster_ids.add(new_rid)

            prev_gtdb_sp = [
                prev_genomes[prev_rid].gtdb_taxa.species for prev_rid in prev_rids]

            prev_cids = set()
            for prev_rid in prev_rids:
                prev_cids.update(prev_genomes.sp_clusters[prev_rid])

            if new_rid is None:
                new_rid = 'none'
                rep_status = 'LOST'
                new_cluster = set()
                rep_lost_count += len(prev_rids)
            else:
                new_cluster = new_clusters[new_rid]

                if len(prev_rids) == 1:
                    if prev_rids[0] == new_rid:
                        rep_status = 'UNCHANGED'
                        rep_unchanged_count += 1
                    else:
                        rep_status = 'CHANGED'
                        rep_changed_count += 1
                else:
                    rep_status = 'MERGER'
                    rep_merger_count += len(prev_rids)

            fout.write('{}\t{}\t{}\t{}'.format(
                new_rid,
                ', '.join(prev_rids),
                ', '.join(prev_gtdb_sp),
                rep_status))

            num_same = len(new_cluster.intersection(prev_cids))
            num_new = len(new_cluster - prev_gids)
            num_lost = len(prev_cids - new_gids)

            num_migrated_in = len(
                (new_cluster - prev_cids).intersection(prev_gids))
            num_migrated_out = len(
                (prev_cids - new_cluster).intersection(new_gids))

            assert len(new_cluster) == len(prev_cids) - num_lost + \
                num_new + num_migrated_in - num_migrated_out
            assert len(prev_cids) == num_same + num_lost + num_migrated_out

            fout.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                len(prev_cids),
                len(new_cluster),
                num_same,
                num_lost,
                num_new,
                num_migrated_in,
                num_migrated_out))

            total_num_same += num_same
            total_num_lost += num_lost
            total_num_new += num_new
            total_num_migrated_in += num_migrated_in
            total_num_migrated_out += num_migrated_out

        assert len(prev_genomes.sp_clusters) == rep_unchanged_count + \
            rep_changed_count + rep_merger_count + rep_lost_count

        # add in new GTDB species clusters
        new_cluster_count = 0
        for new_rid in new_clusters:
            if new_rid in prev_cluster_ids:
                continue

            new_gtdb_sp = cur_genomes[new_rid].gtdb_taxa.species
            rep_status = 'NEW'
            new_cluster_count += 1

            fout.write('{}\t{}\t{}\t{}'.format(
                new_rid,
                'n/a',
                new_gtdb_sp,
                rep_status))

            num_new = len(new_clusters[new_rid] - prev_gids)
            num_migrated_in = len(
                new_clusters[new_rid].intersection(prev_gids))
            assert len(new_clusters[new_rid]) == num_new + num_migrated_in
            fout.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                0,
                len(new_clusters[new_rid]),
                0,
                0,
                num_new,
                num_migrated_in,
                0))

            total_num_new += num_new
            total_num_migrated_in += num_migrated_in

        assert len(new_gids.union(prev_gids)) == total_num_same + \
            total_num_lost + total_num_new + total_num_migrated_in
        assert total_num_migrated_in == total_num_migrated_out

        # report genome statistics
        assert len(prev_gids) == total_num_same + \
            total_num_lost + total_num_migrated_in
        self.log.info(
            f'There were {len(prev_gids):,} genomes in the previous release:')
        self.log.info(' - identified {:,} ({:.2f}%) genomes that were assigned to same species cluster'.format(
            total_num_same,
            total_num_same*100.0/len(prev_gids)))
        self.log.info(' - identified {:,} ({:.2f}%) genomes that were lost from the species cluster'.format(
            total_num_lost,
            total_num_lost*100.0/len(prev_gids)))
        self.log.info(' - identified {:,} ({:.2f}%) genomes that migrated between species cluster'.format(
            total_num_migrated_in,
            total_num_migrated_in*100.0/len(prev_gids)))
        self.log.info('Identified {:,} new genomes which is a {:.2f}% increase.'.format(
            total_num_new,
            len(new_gids)*100.0/len(prev_gids) - 100))

        # report representative statistics
        assert len(prev_genomes.sp_clusters) == rep_unchanged_count + \
            rep_changed_count + rep_lost_count + rep_merger_count
        self.log.info(
            f'There were {len(prev_genomes.sp_clusters):,} previous GTDB species representatives:')
        self.log.info(' - identified {:,} ({:.2f}%) unchanged representatives'.format(
            rep_unchanged_count,
            rep_unchanged_count*100.0/len(prev_genomes.sp_clusters)))
        self.log.info(' - identified {:,} ({:.2f}%) changed representatives'.format(
            rep_changed_count,
            rep_changed_count*100.0/len(prev_genomes.sp_clusters)))
        self.log.info(' - identified {:,} ({:.2f}%) lost representatives'.format(
            rep_lost_count,
            rep_lost_count*100.0/len(prev_genomes.sp_clusters)))
        self.log.info(' - identified {:,} ({:.2f}%) merged representatives'.format(
            rep_merger_count,
            rep_merger_count*100.0/len(prev_genomes.sp_clusters)))
        self.log.info('Identified {:,} new representatives which is a {:.2f}% increase.'.format(
            new_cluster_count,
            len(new_clusters)*100.0/len(prev_genomes.sp_clusters) - 100))
