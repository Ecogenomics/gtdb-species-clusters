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
import pickle
from collections import namedtuple

from numpy import (mean as np_mean)

from gtdblib.util.shell.execute import check_dependencies

from gtdb_species_clusters.mash import Mash
from gtdb_species_clusters.skani import Skani

from gtdb_species_clusters.genomes import Genomes

from gtdb_species_clusters.type_genome_utils import (GenomeRadius,
                                                     write_rep_radius,
                                                     write_clusters)
from gtdb_species_clusters import defaults as Defaults


class UpdateClusterNamedReps(object):
    """Cluster genomes to selected GTDB representatives."""

    def __init__(self, ani_sp, af_sp, ani_cache_file, cpus, output_dir):
        """Initialization."""

        check_dependencies(['skani', 'mash'])

        self.cpus = cpus
        self.output_dir = output_dir

        self.log = logging.getLogger('rich')

        self.ani_sp = ani_sp
        self.af_sp = af_sp

        self.max_ani_neighbour = Defaults.ANI_SYNONYMS
        self.max_af_neighbour = Defaults.AF_SP
        self.min_mash_ani = Defaults.MASH_MIN_ANI

        self.ClusteredGenome = namedtuple('ClusteredGenome', 'ani af gid')

        self.skani = Skani(ani_cache_file, cpus)

    def _rep_radius(self, rep_gids, rep_ani_file):
        """Calculate circumscription radius for representative genomes."""

        # set radius for all representative genomes to default values
        rep_radius = {}
        for gid in rep_gids:
            rep_radius[gid] = GenomeRadius(ani=self.ani_sp,
                                           af=None,
                                           neighbour_gid=None)

        # determine closest ANI neighbour and restrict ANI radius as necessary
        af_warning_count = 0
        with open(rep_ani_file) as f:
            header = f.readline().strip().split('\t')

            rep_gid1_index = header.index('Representative 1')
            rep_gid2_index = header.index('Representative 2')
            ani_index = header.index('ANI')
            af12_index = header.index('AF12')
            af21_index = header.index('AF21')

            for line in f:
                line_split = line.strip().split('\t')

                rep_gid1 = line_split[rep_gid1_index]
                rep_gid2 = line_split[rep_gid2_index]

                if rep_gid1 not in rep_gids or rep_gid2 not in rep_gids:
                    continue

                ani = float(line_split[ani_index])
                af12 = float(line_split[af12_index])
                af21 = float(line_split[af21_index])
                af = max(af12, af21)

                if ani >= self.max_ani_neighbour and af >= self.max_af_neighbour:
                    # typically, representative genomes should not exceed this ANI and AF
                    # criteria as they should have been declared synonyms in
                    # the u_sel_reps step if they are this similar to each other.
                    # However, a 'fudge factor' is used to allow previous GTDB clusters
                    # to remain as seperate clusters if they exceed these thresholds by
                    # a small margin as this can simply be due to differences in the
                    # version of FastANI / skani used to calculate ANI and AF.
                    self.log.warning('ANI neighbours {} and {} have ANI={:.2f} and AF={:.1f}.'.format(
                        rep_gid1, rep_gid2,
                        ani, af))

                if ani > rep_radius[rep_gid1].ani:
                    if af < self.af_sp:
                        af_warning_count += 1
                        # self.log.warning('ANI for {} and {} is >{:.2f}, but AF <{:.2f} [pair skipped].'.format(
                        #                        rep_gid1,
                        #                        rep_gid2,
                        #                        ani, af))
                        continue

                    rep_radius[rep_gid1] = GenomeRadius(ani=ani,
                                                        af=af,
                                                        neighbour_gid=rep_gid2)

        self.log.info(' - ANI circumscription radius: min={:.2f}, mean={:.2f}, max={:.2f}'.format(
            min([d.ani for d in rep_radius.values()]),
            np_mean([d.ani for d in rep_radius.values()]),
            max([d.ani for d in rep_radius.values()])))

        self.log.warning(' - identified {:,} genome pairs meeting ANI radius criteria, but with an AF <{:.1f}'.format(
            af_warning_count,
            self.af_sp))

        return rep_radius

    def calculate_ani(self, cur_genomes, non_rep_gids, rep_gids):
        """Calculate ANI between non-representative and representative genomes.
        
        Returns a dictionary of dictionaries indicating ANI/AF values between
        non-representative and representative genomes. Since skani is symmetric a pair
        is only recorded once with the non-representative used as the first key, e.g.
        d[non-rep genome ID][rep genome ID] = (ani, af)
        """

        if True:  # ***DEBUGGING
            mash = Mash(self.cpus)

            # create Mash sketch for representative genomes
            rep_genome_list_file = os.path.join(
                self.output_dir, 'gtdb_reps.lst')
            rep_mash_sketch_file = os.path.join(
                self.output_dir, 'gtdb_reps.msh')
            mash.sketch(rep_gids, cur_genomes.genomic_files,
                        rep_genome_list_file, rep_mash_sketch_file)

            # create Mash sketch for non-representative genomes
            nonrep_genome_list_file = os.path.join(
                self.output_dir, 'gtdb_nonreps.lst')
            nonrep_genome_sketch_file = os.path.join(
                self.output_dir, 'gtdb_nonreps.msh')
            mash.sketch(non_rep_gids, cur_genomes.genomic_files,
                        nonrep_genome_list_file, nonrep_genome_sketch_file)

            # get Mash distance between representative and non-representative genomes
            mash_dist_file = os.path.join(
                self.output_dir, 'gtdb_reps_vs_nonreps.dst')
            mash.dist(float(100 - self.min_mash_ani)/100,
                      rep_mash_sketch_file,
                      nonrep_genome_sketch_file,
                      mash_dist_file)

            # read Mash distances
            mash_ani = mash.read_ani(mash_dist_file)

            # get genome pairs above Mash threshold
            mash_ani_pairs = []
            for non_rid in mash_ani:
                for rid in mash_ani[non_rid]:
                    if non_rid != rid and mash_ani[non_rid][rid] >= self.min_mash_ani:
                        # skani is symmetric so only need to store pair once; here
                        # the non-representative genome is put first to make for 
                        # easier downstream processing
                        mash_ani_pairs.append((non_rid, rid))

            self.log.info('Identified {:,} genome pairs with a Mash ANI >= {:.1f}%.'.format(
                len(mash_ani_pairs), self.min_mash_ani))

            # calculate ANI between pairs
            self.log.info(
                'Calculating ANI between {:,} genome pairs:'.format(len(mash_ani_pairs)))
            ani_af = self.skani.pairs(
                mash_ani_pairs, cur_genomes.genomic_files)

            pkl_file = os.path.join(self.output_dir, 'ani_af_nonrep_vs_rep.pkl')
            pickle.dump(ani_af, open(pkl_file, 'wb'))
        else:
            self.log.warning('Using previously calculated results in: {}'.format(
                'ani_af_nonrep_vs_rep.pkl'))
            ani_af = pickle.load(
                open(os.path.join(self.output_dir, 'ani_af_nonrep_vs_rep.pkl'), 'rb'))

        return ani_af

    def cluster(self, ani_af, non_rep_gids, rep_radius):
        """Cluster non-representative to representative genomes using species-specific ANI thresholds."""

        clusters = {}
        for rep_id in rep_radius:
            clusters[rep_id] = []

        num_clustered = 0
        for idx, non_rid in enumerate(non_rep_gids):
            if idx % 100 == 0:
                sys.stdout.write('==> Processed {:,} of {:,} genomes [no. clustered = {:,}].\r'.format(
                    idx+1,
                    len(non_rep_gids),
                    num_clustered))
                sys.stdout.flush()

            if non_rid not in ani_af:
                continue

            closest_rid = None
            closest_ani = 0
            closest_af = 0
            for rid in ani_af.get(non_rid, []):
                ani, af = Skani.symmetric_ani_af(ani_af, non_rid, rid)

                if af >= self.af_sp:
                    if ani > closest_ani or (ani == closest_ani and af > closest_af):
                        closest_rid = rid
                        closest_ani = ani
                        closest_af = af

            if closest_rid:
                if closest_ani > rep_radius[closest_rid].ani:
                    num_clustered += 1
                    clusters[closest_rid].append(self.ClusteredGenome(gid=non_rid,
                                                                      ani=closest_ani,
                                                                      af=closest_af))

        sys.stdout.write('==> Processed {:,} of {:,} genomes [no. clustered = {:,}].\r'.format(
            len(non_rep_gids),
            len(non_rep_gids),
            num_clustered))
        sys.stdout.flush()
        sys.stdout.write('\n')

        num_unclustered = len(non_rep_gids) - num_clustered
        self.log.info('Assigned {:,} genomes to {:,} representatives; {:,} genomes remain unclustered.'.format(
            sum([len(clusters[rid]) for rid in clusters]),
            len(clusters),
            num_unclustered))

        return clusters

    def run(self,
            named_rep_file,
            cur_gtdb_metadata_file,
            cur_genomic_path_file,
            qc_passed_file,
            ncbi_genbank_assembly_file,
            untrustworthy_type_file,
            rep_ani_file,
            gtdb_type_strains_ledger,
            ncbi_env_bioproject_ledger):
        """Cluster genomes to selected GTDB representatives."""

        # create current GTDB genome sets
        self.log.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                            gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                            qc_passed_file=qc_passed_file,
                                            ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                            untrustworthy_type_ledger=untrustworthy_type_file,
                                            ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger,
                                            create_sp_clusters=False)

        # get path to previous and current genomic FASTA files
        self.log.info('Reading path to current genomic FASTA files.')
        cur_genomes.load_genomic_file_paths(cur_genomic_path_file)

        # get representative genomes
        rep_gids = set()
        with open(named_rep_file) as f:
            header = f.readline().strip().split('\t')
            rep_index = header.index('Representative')

            for line in f:
                line_split = line.strip().split('\t')
                gid = line_split[rep_index]
                assert gid in cur_genomes
                rep_gids.add(gid)

        self.log.info(
            'Identified representative genomes for {:,} species.'.format(len(rep_gids)))

        # calculate circumscription radius for representative genomes
        self.log.info(
            'Determining ANI species circumscription for {:,} representative genomes:'.format(len(rep_gids)))
        rep_radius = self._rep_radius(rep_gids, rep_ani_file)
        write_rep_radius(rep_radius, cur_genomes, os.path.join(
            self.output_dir, 'gtdb_rep_ani_radius.tsv'))

        assert set(rep_radius) == rep_gids

        # calculate ANI between representative and non-representative genomes
        non_rep_gids = set(cur_genomes.genomes) - rep_gids
        self.log.info(
            'Calculating ANI between representative and non-representative genomes.')
        ani_af = self.calculate_ani(cur_genomes, non_rep_gids, rep_gids)

        # cluster remaining genomes to representatives
        self.log.info(
            'Clustering {:,} non-representatives to {:,} representatives using species-specific ANI radii.'.format(len(non_rep_gids), len(rep_radius)))
        clusters = self.cluster(ani_af, non_rep_gids, rep_radius)

        # write out clusters
        write_clusters(clusters,
                       rep_radius,
                       cur_genomes,
                       os.path.join(self.output_dir, 'gtdb_named_rep_clusters.tsv'))
