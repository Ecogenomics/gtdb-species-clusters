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
from itertools import permutations
from collections import defaultdict, namedtuple

from biolib.external.execute import check_dependencies

from numpy import (mean as np_mean,
                   std as np_std)

from gtdb_species_clusters.fastani import FastANI
from gtdb_species_clusters.genomes import Genomes


class RepGenomicSimilarity(object):
    """Calculate ANI/AF betwenn GTDB representative genomes with the same genus."""

    def __init__(self,
                 ani_cache_file,
                 cpus,
                 output_dir):
        """Initialization."""

        check_dependencies(['fastANI'])

        self.cpus = cpus
        self.output_dir = output_dir

        self.log = logging.getLogger('timestamp')

        self.fastani = FastANI(ani_cache_file, cpus)

    def run(self, gtdb_metadata_file,
            genomic_path_file):
        """Dereplicate GTDB species clusters using ANI/AF criteria."""

        # create GTDB genome sets
        self.log.info('Creating GTDB genome set.')
        genomes = Genomes()
        genomes.load_from_metadata_file(gtdb_metadata_file)
        genomes.load_genomic_file_paths(genomic_path_file)
        self.log.info(' - genome set has {:,} species clusters spanning {:,} genomes.'.format(
            len(genomes.sp_clusters),
            genomes.sp_clusters.total_num_genomes()))

        # get GTDB representatives from same genus
        self.log.info('Identifying GTDB representatives in the same genus.')
        genus_gids = defaultdict(list)
        num_reps = 0
        for gid in genomes:
            if not genomes[gid].gtdb_is_rep:
                continue

            gtdb_genus = genomes[gid].gtdb_taxa.genus
            genus_gids[gtdb_genus].append(gid)
            num_reps += 1
        self.log.info(
            f' - identified {len(genus_gids):,} genera spanning {num_reps:,} representatives')

        # get all intragenus comparisons
        self.log.info('Determining all intragenus comparisons.')
        gid_pairs = []
        for gids in genus_gids.values():
            if len(gids) < 2:
                continue

            for g1, g2 in permutations(gids, 2):
                gid_pairs.append((g1, g2))
        self.log.info(
            f' - identified {len(gid_pairs):,} intragenus comparisons')

        # calculate FastANI ANI/AF between target genomes
        self.log.info('Calculating ANI between intragenus pairs.')
        ani_af = self.fastani.pairs(
            gid_pairs, genomes.genomic_files, report_progress=True, check_cache=True)
        self.fastani.write_cache(silence=True)

        # write out results
        fout = open(os.path.join(self.output_dir,
                    'intragenus_ani_af_reps.tsv'), 'w')
        fout.write(
            'Query ID\tQuery species\tTarget ID\tTarget species\tANI\tAF\n')
        for qid in ani_af:
            for rid in ani_af:
                ani, af = FastANI.symmetric_ani(ani_af, qid, rid)

                fout.write('{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\n'.format(
                    qid,
                    genomes[qid].gtdb_taxa.species,
                    rid,
                    genomes[rid].gtdb_taxa.species,
                    ani,
                    af))
        fout.close()
