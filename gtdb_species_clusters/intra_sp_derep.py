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
import operator
from itertools import permutations

from gtdblib.util.shell.execute import check_dependencies

from gtdb_species_clusters.mash import Mash
from gtdb_species_clusters.skani import Skani
from gtdb_species_clusters.genomes import Genomes


class IntraSpeciesDereplication(object):
    """Dereplicate GTDB species clusters using ANI/AF criteria."""

    def __init__(self,
                 derep_ani,
                 derep_af,
                 max_genomes_per_sp,
                 ani_cache_file,
                 cpus,
                 output_dir):
        """Initialization."""

        check_dependencies(['skani', 'mash'])

        self.cpus = cpus
        self.output_dir = output_dir

        self.log = logging.getLogger('rich')

        self.max_genomes_per_sp = max_genomes_per_sp
        self.derep_ani = derep_ani
        self.derep_af = derep_af

        # minimum MASH ANI value for dereplicating within a species
        self.min_mash_intra_sp_ani = derep_ani - 1.0

        self.mash = Mash(self.cpus)
        self.skani = Skani(ani_cache_file, cpus)

    def mash_sp_ani(self, gids, genomes, output_prefix):
        """Calculate pairwise Mash ANI estimates between genomes."""

        INIT_MASH_ANI_FILTER = 95.0

        # create Mash sketch for all genomes
        mash_sketch_file = f'{output_prefix}.msh'
        genome_list_file = f'{output_prefix}.lst'
        self.mash.sketch(gids,
                         genomes.genomic_files,
                         genome_list_file,
                         mash_sketch_file,
                         silence=True)

        # get Mash distances
        mash_dist_file = f'{output_prefix}.dst'
        self.mash.dist_pairwise(float(100 - INIT_MASH_ANI_FILTER)/100,
                                mash_sketch_file,
                                mash_dist_file,
                                silence=True)

        # read Mash distances
        mash_ani = self.mash.read_ani(mash_dist_file)

        count = 0
        for qid in mash_ani:
            for rid in mash_ani[qid]:
                if qid != rid:
                    count += 1

        self.log.info(' - identified {:,} pairs passing Mash filtering of ANI >= {:.1f}%.'.format(
            count,
            INIT_MASH_ANI_FILTER))

        return mash_ani

    def priority_score(self, gid, genomes):
        """Get priority score of genome."""

        score = genomes[gid].score_assembly()
        if genomes[gid].is_gtdb_type_subspecies():
            score += 1e4

        return score

    def order_genomes_by_priority(self, gids, genomes):
        """Order genomes by overall priority. """

        genome_priority = {}
        for gid in gids:
            genome_priority[gid] = self.priority_score(gid, genomes)

        sorted_by_priority = sorted(genome_priority.items(),
                                    key=operator.itemgetter(1),
                                    reverse=True)

        return [d[0] for d in sorted_by_priority]

    def mash_sp_dereplicate(self, mash_ani, sorted_gids, ani_threshold):
        """Dereplicate genomes in species using Mash distances."""

        # perform greedy selection of new representatives
        sp_reps = []
        for gid in sorted_gids:
            clustered = False
            for rep_id in sp_reps:
                if gid in mash_ani:
                    ani = mash_ani[gid].get(rep_id, 0)
                else:
                    ani = 0

                if ani >= ani_threshold:
                    clustered = True
                    break

            if not clustered:
                # genome was not assigned to an existing representative,
                # so make it a new representative genome
                sp_reps.append(gid)

        return sp_reps

    def dereplicate_species(self,
                            species,
                            rid,
                            cids,
                            genomes,
                            mash_out_dir):
        """Dereplicate genomes within a GTDB species."""

        # greedily dereplicate genomes based on genome priority
        sorted_gids = self.order_genomes_by_priority(cids.difference([rid]),
                                                     genomes)
        sorted_gids = [rid] + sorted_gids

        # calculate Mash ANI between genomes
        mash_ani = []
        if len(sorted_gids) > 1:
            # calculate MASH distances between genomes
            out_prefix = os.path.join(
                mash_out_dir, species[3:].lower().replace(' ', '_'))
            mash_ani = self.mash_sp_ani(sorted_gids,
                                        genomes,
                                        out_prefix)

        # perform initial dereplication using Mash for species with excessive
        # numbers of genomes
        if len(sorted_gids) > self.max_genomes_per_sp:
            self.log.info(' - limiting species to <={:,} genomes based on priority and Mash dereplication.'.format(
                self.max_genomes_per_sp))

            prev_mash_rep_gids = None
            for ani_threshold in [99.75, 99.5, 99.25, 99.0, 98.75, 98.5, 98.25, 98.0, 97.75, 97.5, 97.0, 96.5, 96.0, 95.0, None]:
                if ani_threshold is None:
                    self.log.warning(
                        ' - delected {:,} highest priority genomes from final Mash dereplication.' % self.max_genomes_per_sp)
                    sorted_gids = mash_rep_gids[0:self.max_genomes_per_sp]
                    break

                mash_rep_gids = self.mash_sp_dereplicate(mash_ani,
                                                         sorted_gids,
                                                         ani_threshold)

                self.log.info(' - dereplicated {} from {:,} to {:,} genomes at {:.2f}% ANI using Mash.'.format(
                    species,
                    len(cids),
                    len(mash_rep_gids),
                    ani_threshold))

                if len(mash_rep_gids) <= self.max_genomes_per_sp:
                    if not prev_mash_rep_gids:
                        # corner case where dereplication is occurring at 99.75%
                        prev_mash_rep_gids = sorted_gids

                    # select maximum allowed number of genomes by taking all genomes in the
                    # current Mash dereplicated set and then the highest priority genomes in the
                    # previous Mash dereplicated set which have not been selected
                    cur_sel_gids = set(mash_rep_gids)
                    prev_sel_gids = set(prev_mash_rep_gids)
                    num_prev_to_sel = self.max_genomes_per_sp - \
                        len(cur_sel_gids)
                    num_prev_selected = 0
                    sel_sorted_gids = []
                    for gid in sorted_gids:
                        if gid in cur_sel_gids:
                            sel_sorted_gids.append(gid)
                        elif (gid in prev_sel_gids
                                and num_prev_selected < num_prev_to_sel):
                            num_prev_selected += 1
                            sel_sorted_gids.append(gid)

                        if len(sel_sorted_gids) == self.max_genomes_per_sp:
                            break

                    assert len(cur_sel_gids - set(sel_sorted_gids)) == 0
                    assert num_prev_to_sel == num_prev_selected
                    assert len(sel_sorted_gids) == self.max_genomes_per_sp

                    sorted_gids = sel_sorted_gids
                    self.log.info(' - selected {:,} highest priority genomes from Mash dereplication at an ANI = {:.2f}%.'.format(
                        len(sorted_gids),
                        ani_threshold))
                    break

                prev_mash_rep_gids = mash_rep_gids
                prev_ani_threshold = ani_threshold

        # calculate skani ANI/AF between genomes passing Mash filtering
        ani_pairs = set()
        for gid1, gid2 in permutations(sorted_gids, 2):
            if gid1 in mash_ani and gid2 in mash_ani[gid1]:
                if mash_ani[gid1][gid2] >= self.min_mash_intra_sp_ani:
                    # skani ANI is symmetric and the AF is calculated in both
                    # directions so only need to run a given pair once
                    if gid1 < gid2:
                        ani_pairs.append((gid1, gid2))
                    else:
                        ani_pairs.append((gid2, gid1))

        self.log.info(' - calculating skani between {:,} pairs with Mash ANI >= {:.1f}%.'.format(
            len(ani_pairs),
            self.min_mash_intra_sp_ani))
        ani_af = self.skani.pairs(ani_pairs,
                                    genomes.genomic_files,
                                    report_progress=False,
                                    check_cache=True)
        self.skani.write_cache(silence=True)

        # perform greedy dereplication
        sp_reps = []
        for idx, gid in enumerate(sorted_gids):
            # determine if genome clusters with existing representative
            clustered = False
            for rid in sp_reps:
                ani, af = Skani.symmetric_ani(ani_af, gid, rid)

                if ani >= self.derep_ani and af >= self.derep_af:
                    clustered = True
                    break

            if not clustered:
                sp_reps.append(gid)

        self.log.info(' - dereplicated {} from {:,} to {:,} genomes.'.format(
            species,
            len(sorted_gids),
            len(sp_reps)))

        # assign clustered genomes to most similar representative
        subsp_clusters = {}
        for rid in sp_reps:
            subsp_clusters[rid] = [rid]

        non_rep_gids = set(sorted_gids) - set(sp_reps)
        for gid in non_rep_gids:
            closest_rid = None
            max_ani = 0
            max_af = 0
            for rid in sp_reps:
                ani, af = Skani.symmetric_ani(ani_af, gid, rid)
                if ((ani > max_ani and af >= self.derep_af)
                        or (ani == max_ani and af >= max_af and af >= self.derep_af)):
                    max_ani = ani
                    max_af = af
                    closest_rid = rid

            assert closest_rid is not None
            subsp_clusters[closest_rid].append(gid)

        return subsp_clusters

    def derep_sp_clusters(self, genomes):
        """Dereplicate each GTDB species cluster."""

        mash_out_dir = os.path.join(self.output_dir, 'mash')
        if not os.path.exists(mash_out_dir):
            os.makedirs(mash_out_dir)

        derep_genomes = {}
        for rid, cids in genomes.sp_clusters.items():
            species = genomes[rid].gtdb_taxa.species

            self.log.info('Dereplicating {} with {:,} genomes [{:,} of {:,} ({:.2f}%) species].'.format(
                species,
                len(cids),
                len(derep_genomes),
                len(genomes.sp_clusters),
                len(derep_genomes)*100.0/len(genomes.sp_clusters)))

            subsp_clusters = self.dereplicate_species(species,
                                                      rid,
                                                      cids,
                                                      genomes,
                                                      mash_out_dir)

            derep_genomes[species] = subsp_clusters

        return derep_genomes

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

        # dereplicate each species cluster
        self.log.info('Performing dereplication with ANI={:.1f}, AF={:.2f}, Mash ANI={:.2f}, max genomes={:,}.'.format(
            self.derep_ani,
            self.derep_af,
            self.min_mash_intra_sp_ani,
            self.max_genomes_per_sp))
        derep_genomes = self.derep_sp_clusters(genomes)

        # write out `subspecies` clusters
        out_file = os.path.join(self.output_dir, 'subsp_clusters.tsv')
        fout = open(out_file, 'w')
        fout.write(
            'Genome ID\tGTDB Species\tGTDB Taxonomy\tPriority score\tNo. clustered genomes\tNo. clustered genomes\tClustered genomes\n')
        for species, subsp_clusters in derep_genomes.items():
            for rid, cids in subsp_clusters.items():
                assert species == genomes[rid].gtdb_taxa.species
                fout.write('{}\t{}\t{}\t{:.3f}\t{}\t{}\n'.format(
                    rid,
                    genomes[rid].gtdb_taxa.species,
                    genomes[rid].gtdb_taxa,
                    self.priority_score(rid, genomes),
                    len(cids),
                    ','.join(cids)))
