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
from itertools import combinations

from gtdblib.util.shell.execute import check_dependencies

from gtdb_species_clusters.skani import Skani
from gtdb_species_clusters.genomes import Genomes


class ANI_SpeciesPair():
    """Calculate all pairwise ANI/AF values between genomes in two species."""

    def __init__(self, ani_cache_file, cpus, output_dir):
        """Initialization."""

        check_dependencies(['skani'])

        self.cpus = cpus
        self.output_dir = output_dir

        self.log = logging.getLogger('rich')

        self.skani = Skani(ani_cache_file, cpus)

    def run(self, gtdb_metadata_file,
            genome_path_file,
            species1,
            species2):
        """Calculate all pairwise ANI/AF values between genomes in two species."""

        # read GTDB species clusters
        self.log.info('Reading GTDB species clusters.')
        genomes = Genomes()
        genomes.load_from_metadata_file(gtdb_metadata_file)
        genomes.load_genomic_file_paths(genome_path_file)
        self.log.info(' - identified {:,} species clusters spanning {:,} genomes.'.format(
            len(genomes.sp_clusters),
            genomes.sp_clusters.total_num_genomes()))

        # find representatives for species of interest
        rid1 = None
        rid2 = None
        gid_to_sp = {}
        for gid, species in genomes.sp_clusters.species():
            if species == species1:
                rid1 = gid
            elif species == species2:
                rid2 = gid

        if rid1 is None:
            self.log.error(
                f'Unable to find representative genome for {species1}.')
            sys.exit(-1)

        if rid2 is None:
            self.log.error(
                f'Unable to find representative genome for {species2}.')
            sys.exit(-1)

        self.log.info(' - identified {:,} genomes in {}.'.format(
            len(genomes.sp_clusters[rid1]),
            species1))
        self.log.info(' - identified {:,} genomes in {}.'.format(
            len(genomes.sp_clusters[rid2]),
            species2))

        # calculate pairwise ANI between all genomes
        self.log.info(f'Calculating pairwise ANI between all genomes.')
        all_gids = genomes.sp_clusters[rid1].union(genomes.sp_clusters[rid2])
        gid_pairs = [(gid1, gid2) for gid1, gid2 in combinations(all_gids, 2)]
        ani_af = self.skani.pairs(gid_pairs, genomes.genomic_files)

        fout = open(os.path.join(self.output_dir, "ani_sp_pairwise.tsv"), 'w')
        fout.write('Genome ID 1\tSpecies 1\tGenome ID 2\tSpecies 2\tANI\tAF\n')
        for gid1 in ani_af:
            for gid2 in ani_af[gid1]:
                ani, af = ani_af[gid1][gid2]
                fout.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    gid1,
                    genomes[gid1].gtdb_taxa.species,
                    gid2,
                    genomes[gid2].gtdb_taxa.species,
                    ani,
                    af))
        fout.close()
