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
from collections import defaultdict, Counter

from biolib.taxonomy import Taxonomy

from gtdb_species_clusters.genome import Genome
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters.taxon_utils import specific_epithet
from gtdb_species_clusters.fastani import FastANI


class NCBI_SpeciesManager():
    """Manage NCBI species names assigned to genomes.

     - establish NCBI species considered to be synonyms under the GTDB
     - identify genomes considered to have erroneous NCBI species assignments
     - classify each NCBI species as unambiguous, ambiguous, or synonym
     - validate NCBI type strain genomes are restrcited to a single GTDB cluster
    """

    # NCBI species classifications
    UNAMBIGUOUS = 'UNAMBIGUOUS'
    AMBIGUOUS = 'AMBIGUOUS'

    # Assignment types for unambiguous species
    MAJORITY_VOTE = 'MAJORITY_VOTE'
    TYPE_STRAIN_GENOME = 'TYPE_STRAIN_GENOME'

    def __init__(self, cur_genomes, cur_clusters, mc_species, output_dir):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        self.output_dir = output_dir
        self.cur_genomes = cur_genomes
        self.cur_clusters = cur_clusters

        # specific names that should effectively be treated as genomes
        # lacking an NCBI species assignment. Ideally, all such cases
        # would be caught ahead of time when initially parsing the NCBI
        # taxonomy so this only includes last minute additions.
        self.forbidden_specific_names = set(['cyanobacterium', 'Family'])

        # identify NCBI species with genomes assembled from type strain of species
        self.ncbi_sp_type_strain_genomes = cur_genomes.ncbi_sp_effective_type_genomes()
        self.logger.info(' - identified effective type strain genomes for {:,} NCBI species.'.format(
            len(self.ncbi_sp_type_strain_genomes)))

        self.validate_type_strain_clustering(mc_species)

    def validate_type_strain_clustering(self, mc_species):
        """Validate that all type strain genomes for an NCBI species occur in a single GTDB cluster."""

        self.logger.info(
            'Verifying that all type strain genomes for a NCBI species occur in a single GTDB cluster.')
        rid_map = {}
        for rid, gids in self.cur_clusters.items():
            rid_map[rid] = rid
            for gid in gids:
                rid_map[gid] = rid

        for ncbi_sp, type_gids in self.ncbi_sp_type_strain_genomes.items():
            gtdb_rids = set([rid_map[gid] for gid in type_gids])

            gtdb_rids = set()
            for gid in type_gids:
                rid = rid_map[gid]
                ncbi_specific = specific_epithet(
                    self.cur_genomes[rid].ncbi_taxa.species)

                if (rid in mc_species
                        and specific_epithet(mc_species[rid]) != ncbi_specific):
                    # skip this genome as it has been manually changed, likely
                    # to resolve this NCBI species having multiple type strain genomes
                    continue

                gtdb_rids.add(rid)

            if len(gtdb_rids) > 1:
                self.logger.error('Type strain genomes from NCBI species {} were assigned to {:,} GTDB species clusters: {}.'.format(
                    ncbi_sp,
                    len(gtdb_rids),
                    [(rid, self.cur_genomes[rid].gtdb_taxa.species) for rid in gtdb_rids]))
                sys.exit(-1)

    def identify_type_strain_synonyms(self, ncbi_misclassified_gids):
        """Identify synonyms arising from multiple type strain genomes residing in the same GTDB cluster.

        A type strain genome residing in a GTDB species cluster with a different NCBI assignment is considered
        a type strain synonym.
        """

        # identify type strain synonyms
        self.logger.info('Identifying type strain synonyms.')
        type_strain_synonyms = defaultdict(list)
        failed_type_strain_priority = 0
        for rid, gids in self.cur_clusters.items():
            if not self.cur_genomes[rid].is_effective_type_strain():
                continue

            rep_ncbi_sp = self.cur_genomes[rid].ncbi_taxa.species

            # find species that are a type strain synonym to the current representative,
            # using the best quality genome for each species to establish
            # synonym statistics such as ANI and AF
            ncbi_sp_gids = set(gids) - ncbi_misclassified_gids
            type_gids = [
                gid for gid in ncbi_sp_gids if self.cur_genomes[gid].is_effective_type_strain()]
            q = {gid: self.cur_genomes[gid].score_type_strain()
                 for gid in type_gids}
            q_sorted = sorted(q.items(), key=lambda kv: (
                kv[1], kv[0]), reverse=True)
            processed_sp = set([rep_ncbi_sp])
            for gid, _quality in q_sorted:
                cur_ncbi_sp = self.cur_genomes[gid].ncbi_taxa.species

                if cur_ncbi_sp in processed_sp:
                    continue

                type_strain_synonyms[rid].append(gid)
                processed_sp.add(cur_ncbi_sp)

        self.logger.info(' - identified {:,} GTDB representatives resulting in {:,} type strain synonyms.'.format(
            len(type_strain_synonyms),
            sum([len(gids) for gids in type_strain_synonyms.values()])))

        if failed_type_strain_priority:
            self.logger.warning(
                f'Identified {failed_type_strain_priority:,} non-type strain representatives that failed to priotize an effective type strain genome.')

        return type_strain_synonyms

    def identify_consensus_synonyms(self, ncbi_misclassified_gids):
        """Identify synonyms arising from all genomes with an NCBI species classification
            which lack a type strain genome being contained in a GTDB cluster defined by
            a type strain genome."""

        # get genomes with NCBI species assignment
        ncbi_species_gids = defaultdict(list)
        ncbi_all_sp_gids = set()
        for gid in self.cur_genomes:
            if gid in ncbi_misclassified_gids:
                continue

            ncbi_species = self.cur_genomes[gid].ncbi_taxa.species
            ncbi_specific = specific_epithet(ncbi_species)
            if ncbi_species != 's__' and ncbi_specific not in self.forbidden_specific_names:
                ncbi_species_gids[ncbi_species].append(gid)
                ncbi_all_sp_gids.add(gid)

        # identify consensus synonyms
        consensus_synonyms = defaultdict(list)

        for ncbi_species, ncbi_sp_gids in ncbi_species_gids.items():
            if ncbi_species in self.ncbi_sp_type_strain_genomes:
                # NCBI species define by type strain genome so should
                # be either the representative of a GTDB species cluster
                # or a type strain synonym
                continue

            # determine if all genomes in NCBI species are contained
            # in a single GTDB species represented by a type strain genome
            for gtdb_rid, cids in self.cur_clusters.items():
                assert gtdb_rid in cids

                if not self.cur_genomes[gtdb_rid].is_effective_type_strain():
                    continue

                ncbi_cur_sp_in_cluster = cids.intersection(ncbi_sp_gids)
                ncbi_cur_sp_in_cluster_perc = len(
                    ncbi_cur_sp_in_cluster)*100.0/len(ncbi_sp_gids)
                if ncbi_cur_sp_in_cluster_perc == 100:
                    # using the best quality genome in NCBI species to establish
                    # synonym statistics such as ANI and AF
                    q = {gid: self.cur_genomes[gid].score_type_strain(
                    ) for gid in ncbi_sp_gids}
                    q_sorted = sorted(q.items(), key=lambda kv: (
                        kv[1], kv[0]), reverse=True)
                    consensus_synonyms[gtdb_rid].append(q_sorted[0][0])
                    break

        self.logger.info(' - identified {:,} GTDB representatives resulting in {:,} majority vote synonyms.'.format(
            len(consensus_synonyms),
            sum([len(gids) for gids in consensus_synonyms.values()])))

        return consensus_synonyms

    def write_synonym_table(self,
                            type_strain_synonyms,
                            consensus_synonyms,
                            ani_af,
                            sp_priority_ledger,
                            genus_priority_ledger,
                            lpsn_gss_file):
        """Create table indicating species names that should be considered synonyms."""

        sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger,
                                                  genus_priority_ledger,
                                                  lpsn_gss_file,
                                                  self.output_dir)

        out_file = os.path.join(self.output_dir, 'synonyms.tsv')
        fout = open(out_file, 'w')
        fout.write(
            'Synonym type\tNCBI species\tGTDB representative\tStrain IDs\tType sources\tPriority year')
        fout.write('\tGTDB type species\tGTDB type strain\tNCBI assembly type')
        fout.write(
            '\tNCBI synonym\tHighest-quality synonym genome\tSynonym strain IDs\tSynonym type sources\tSynonym priority year')
        fout.write(
            '\tSynonym GTDB type species\tSynonym GTDB type strain\tSynonym NCBI assembly type')
        fout.write('\tANI\tAF\tWarnings\n')

        incorrect_priority = 0
        failed_type_strain_priority = 0
        for synonyms, synonym_type in [(type_strain_synonyms, 'TYPE_STRAIN_SYNONYM'),
                                       (consensus_synonyms, 'MAJORITY_VOTE_SYNONYM')]:
            for rid, synonym_ids in synonyms.items():
                for gid in synonym_ids:
                    ani, af = FastANI.symmetric_ani(ani_af, rid, gid)

                    fout.write(synonym_type)
                    fout.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                        self.cur_genomes[rid].ncbi_taxa.species,
                        rid,
                        ','.join(sorted(self.cur_genomes[rid].strain_ids())),
                        ','.join(sorted(self.cur_genomes[rid].gtdb_type_sources())).upper().replace(
                            'STRAININFO', 'StrainInfo'),
                        sp_priority_mngr.species_priority_year(
                            self.cur_genomes, rid),
                        self.cur_genomes[rid].is_gtdb_type_species(),
                        self.cur_genomes[rid].is_gtdb_type_strain(),
                        self.cur_genomes[rid].ncbi_type_material))

                    synonym_priority_year = sp_priority_mngr.species_priority_year(
                        self.cur_genomes, gid)
                    if synonym_priority_year == Genome.NO_PRIORITY_YEAR:
                        synonym_priority_year = 'n/a'

                    fout.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                        self.cur_genomes[gid].ncbi_taxa.species,
                        gid,
                        ','.join(sorted(self.cur_genomes[gid].strain_ids())),
                        ','.join(sorted(self.cur_genomes[gid].gtdb_type_sources())).upper().replace(
                            'STRAININFO', 'StrainInfo'),
                        synonym_priority_year,
                        self.cur_genomes[gid].is_gtdb_type_species(),
                        self.cur_genomes[gid].is_gtdb_type_strain(),
                        self.cur_genomes[gid].ncbi_type_material))
                    fout.write('\t{:.3f}\t{:.4f}'.format(ani, af))

                    if self.cur_genomes[rid].is_effective_type_strain() and self.cur_genomes[gid].is_effective_type_strain():
                        priority_gid, note = sp_priority_mngr.species_priority(
                            self.cur_genomes, rid, gid)
                        if priority_gid != rid:
                            incorrect_priority += 1
                            fout.write('\tIncorrect priority: {}'.format(note))
                    elif not self.cur_genomes[rid].is_gtdb_type_strain() and self.cur_genomes[gid].is_gtdb_type_strain():
                        failed_type_strain_priority += 1
                        fout.write(
                            '\tFailed to prioritize type strain of species')

                    fout.write('\n')

        if incorrect_priority:
            self.logger.warning(
                f' - identified {incorrect_priority:,} synonyms with incorrect priority.')

        if failed_type_strain_priority:
            self.logger.warning(
                f' - identified {failed_type_strain_priority:,} synonyms that failed to priotize the type strain of the species.')

    def parse_synonyms_table(self, synonym_file):
        """Parse synonyms."""

        synonyms = {}
        with open(synonym_file) as f:
            headers = f.readline().strip().split('\t')

            ncbi_sp_index = headers.index('NCBI species')
            ncbi_synonym_index = headers.index('NCBI synonym')

            for line in f:
                line_split = line.strip().split('\t')

                ncbi_sp = line_split[ncbi_sp_index]
                ncbi_synonym = line_split[ncbi_synonym_index]
                synonyms[ncbi_synonym] = ncbi_sp

        return synonyms

    def parse_ncbi_misclassified_table(self, ncbi_misclassified_file):
        """Parse genomes with erroneous NCBI species assignments."""

        misclassified_gids = set()
        with open(ncbi_misclassified_file) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                misclassified_gids.add(tokens[0])

        return misclassified_gids

    def ncbi_sp_gtdb_cluster_table(self, final_taxonomy):
        """Create table indicating GTDB species clusters for each NCBI species."""

        # get map between genomes and representatives
        gid_to_rid = {}
        for rid, cids in self.cur_clusters.items():
            for cid in cids:
                gid_to_rid[cid] = rid

        # get genomes with NCBI species assignment
        ncbi_species_gids = defaultdict(list)
        for rid, cids in self.cur_clusters.items():
            if rid not in final_taxonomy:
                continue  # other domain

            for cid in cids:
                ncbi_species = self.cur_genomes[cid].ncbi_taxa.species
                ncbi_specific = specific_epithet(ncbi_species)
                if ncbi_species != 's__' and ncbi_specific not in self.forbidden_specific_names:
                    ncbi_species_gids[ncbi_species].append(cid)

        # write out table
        fout = open(os.path.join(self.output_dir,
                                 'ncbi_sp_to_gtdb_cluster_map.tsv'), 'w')
        fout.write(
            'NCBI species\tNo. genomes\tNo. GTDB clusters\tHighest prevalence (%)\tGTDB species clusters\n')

        for ncbi_sp, gids in ncbi_species_gids.items():
            gtdb_rep_list = [gid_to_rid[gid] for gid in gids]
            gtdb_rep_count = Counter(gtdb_rep_list)

            gtdb_rep_str = []
            for idx, (rid, count) in enumerate(gtdb_rep_count.most_common()):
                gtdb_rep_str.append('{} ({}): {}'.format(
                    final_taxonomy[rid][Taxonomy.SPECIES_INDEX],
                    rid,
                    count))

                if idx == 0:
                    highest_prevalence = count*100.0/len(gids)

            fout.write('{}\t{}\t{}\t{:.2f}\t{}\n'.format(
                ncbi_sp,
                len(gids),
                len(set(gtdb_rep_list)),
                highest_prevalence,
                '; '.join(gtdb_rep_str)))

        fout.close()

    def classify_ncbi_species(self, ncbi_synonyms, misclassified_gids):
        """Classify each NCBI species as either unambiguously, ambiguously, or synonym."""

        # get genomes with NCBI species assignment
        ncbi_all_species = set()
        ncbi_species_gids = defaultdict(list)
        ncbi_all_species_gids = set()
        for gid in self.cur_genomes:
            if gid in misclassified_gids:
                continue

            ncbi_species = self.cur_genomes[gid].ncbi_taxa.species
            ncbi_specific = specific_epithet(ncbi_species)
            if ncbi_species != 's__' and ncbi_specific not in self.forbidden_specific_names:
                ncbi_all_species.add(ncbi_species)
                ncbi_species_gids[ncbi_species].append(gid)
                ncbi_all_species_gids.add(gid)

        # process NCBI species with type strain genomes followed by
        # remaining NCBI species
        ncbi_type_strain_species = []
        ncbi_nontype_species = []
        for ncbi_species, ncbi_sp_gids in ncbi_species_gids.items():
            has_type_strain_genome = False
            for gid in ncbi_sp_gids:
                if self.cur_genomes[gid].is_effective_type_strain():
                    has_type_strain_genome = True
                    break

            if has_type_strain_genome:
                ncbi_type_strain_species.append(ncbi_species)
            else:
                ncbi_nontype_species.append(ncbi_species)
        self.logger.info(' - identified {:,} NCBI species with effective type strain genome and {:,} NCBI species without a type strain genome.'.format(
            len(ncbi_type_strain_species), len(ncbi_nontype_species)))

        # establish if NCBI species can be unambiguously assigned to a GTDB species cluster
        fout = open(os.path.join(self.output_dir,
                                 'ncbi_sp_classification.tsv'), 'w')
        fout.write('NCBI taxon\tClassification\tAssignment type\tGTDB representative\t% NCBI species in cluster\t% cluster with NCBI species\tNCBI classifications in cluster\tNote\n')
        unambiguous_ncbi_sp = {}
        ambiguous_ncbi_sp = set()
        unambiguous_rids = set()
        for ncbi_species in ncbi_type_strain_species + ncbi_nontype_species:
            if ncbi_species in ncbi_synonyms:
                continue

            ncbi_sp_gids = ncbi_species_gids[ncbi_species]

            # get type strains for NCBI species that are GTDB species representatives
            type_strain_rids = [gid for gid in ncbi_sp_gids
                                if (self.cur_genomes[gid].is_effective_type_strain() and gid in self.cur_clusters)]

            if len(type_strain_rids) == 0:
                highest_sp_in_cluster_perc = 0
                highest_cluster_perc = 0
                highest_gtdb_rid = None
                for gtdb_rid, cids in self.cur_clusters.items():
                    assert gtdb_rid in cids

                    ncbi_cur_sp_in_cluster = cids.intersection(ncbi_sp_gids)
                    ncbi_any_sp_in_cluster = cids.intersection(
                        ncbi_all_species_gids)

                    if len(ncbi_any_sp_in_cluster) > 0:
                        ncbi_cur_sp_in_cluster_perc = len(
                            ncbi_cur_sp_in_cluster)*100.0/len(ncbi_sp_gids)
                        cluster_perc = len(
                            ncbi_cur_sp_in_cluster)*100.0/len(ncbi_any_sp_in_cluster)

                        if (ncbi_cur_sp_in_cluster_perc > highest_sp_in_cluster_perc
                                or (ncbi_cur_sp_in_cluster_perc == highest_sp_in_cluster_perc
                                    and cluster_perc > highest_cluster_perc)):
                            highest_sp_in_cluster_perc = ncbi_cur_sp_in_cluster_perc
                            highest_cluster_perc = cluster_perc
                            highest_gtdb_rid = gtdb_rid

                sp_in_cluster_threshold = 50
                cluster_threshold = 50

                if (highest_gtdb_rid not in unambiguous_rids
                        and highest_sp_in_cluster_perc > sp_in_cluster_threshold
                        and highest_cluster_perc > cluster_threshold):
                    unambiguous_ncbi_sp[ncbi_species] = (
                        highest_gtdb_rid, NCBI_SpeciesManager.MAJORITY_VOTE)
                else:
                    ambiguous_ncbi_sp.add(ncbi_species)

                ncbi_sp_str = []
                for sp, count in Counter([self.cur_genomes[gid].ncbi_taxa.species for gid in self.cur_clusters[highest_gtdb_rid]]).most_common():
                    ncbi_sp_str.append('{} ({:,})'.format(sp, count))
                ncbi_sp_str = '; '.join(ncbi_sp_str)

                fout.write('{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\t{}\t{}\n'.format(
                    ncbi_species,
                    NCBI_SpeciesManager.UNAMBIGUOUS if ncbi_species in unambiguous_ncbi_sp else NCBI_SpeciesManager.AMBIGUOUS,
                    NCBI_SpeciesManager.MAJORITY_VOTE if ncbi_species in unambiguous_ncbi_sp else NCBI_SpeciesManager.AMBIGUOUS,
                    highest_gtdb_rid,
                    highest_sp_in_cluster_perc,
                    highest_cluster_perc,
                    ncbi_sp_str,
                    ''))

            elif len(type_strain_rids) == 1:
                gtdb_rid = type_strain_rids[0]
                unambiguous_ncbi_sp[ncbi_species] = (
                    gtdb_rid, NCBI_SpeciesManager.TYPE_STRAIN_GENOME)
                unambiguous_rids.add(gtdb_rid)

                cids = self.cur_clusters[gtdb_rid]
                ncbi_cur_sp_in_cluster = cids.intersection(ncbi_sp_gids)
                ncbi_any_sp_in_cluster = cids.intersection(
                    ncbi_all_species_gids)
                cur_sp_in_cluster_perc = len(
                    ncbi_cur_sp_in_cluster)*100.0/len(ncbi_sp_gids)
                cluster_perc = len(ncbi_cur_sp_in_cluster) * \
                    100.0/len(ncbi_any_sp_in_cluster)

                ncbi_sp_str = []
                for sp, count in Counter([self.cur_genomes[gid].ncbi_taxa.species for gid in self.cur_clusters[gtdb_rid]]).most_common():
                    ncbi_sp_str.append('{} ({:,})'.format(sp, count))
                ncbi_sp_str = '; '.join(ncbi_sp_str)

                fout.write('{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\t{}\t{}\n'.format(
                    ncbi_species,
                    NCBI_SpeciesManager.UNAMBIGUOUS,
                    NCBI_SpeciesManager.TYPE_STRAIN_GENOME,
                    gtdb_rid,
                    cur_sp_in_cluster_perc,
                    cluster_perc,
                    ncbi_sp_str,
                    ''))
            else:
                self.logger.error('Multiple GTDB species clusters represented by type strain genomes of {}: {}'.format(
                    ncbi_species, type_strain_rids))
                sys.exit(-1)

        # get NCBI synonyms
        for ncbi_species, ncbi_sp_gids in ncbi_species_gids.items():
            if ncbi_species not in ncbi_synonyms:
                continue

            fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                ncbi_species,
                'SYNONYM',
                'SYNONYM',
                'n/a',
                'n/a',
                'n/a',
                'n/a',
                'Synonym of {} under GTDB'.format(ncbi_synonyms[ncbi_species])))

        fout.close()

        # sanity check results
        assert set(unambiguous_ncbi_sp).intersection(ncbi_synonyms) == set()
        assert set(ambiguous_ncbi_sp).intersection(ncbi_synonyms) == set()
        assert set(unambiguous_ncbi_sp).intersection(
            ambiguous_ncbi_sp) == set()
        assert set(ncbi_all_species) - set(unambiguous_ncbi_sp) - \
            set(ambiguous_ncbi_sp) - set(ncbi_synonyms) == set()

        gtdb_rep_ncbi_sp_assignments = {}
        for ncbi_species, (gtdb_rid, _classification) in unambiguous_ncbi_sp.items():
            if gtdb_rid in gtdb_rep_ncbi_sp_assignments:
                self.logger.error('GTDB representative {} assigned to multiple unambiguous NCBI species: {} {}'.format(
                    gtdb_rid,
                    gtdb_rep_ncbi_sp_assignments[gtdb_rid],
                    ncbi_species))
                # ***sys.exit(-1)
            gtdb_rep_ncbi_sp_assignments[gtdb_rid] = ncbi_species

        return ncbi_all_species, unambiguous_ncbi_sp, ambiguous_ncbi_sp

    def classify_ncbi_subspecies(self):
        """Classify each NCBI subspecies as either unambiguously or ambiguously."""

        # identify GTDB species cluster containing multiple type strain of subspecies genomes,
        # and subspecies with multiple type strains in different GTDB clusters
        gtdb_rid_type_subspecies = defaultdict(set)
        ncbi_type_subspecies_rid = defaultdict(set)
        for rid, cids in self.cur_clusters.items():
            for cid in cids:
                if self.cur_genomes[cid].is_ncbi_type_subspecies():
                    ncbi_subspecies = self.cur_genomes[cid].ncbi_unfiltered_taxa.subspecies

                    if ncbi_subspecies:
                        ncbi_type_subspecies_rid[ncbi_subspecies].add(rid)
                        gtdb_rid_type_subspecies[rid].add(ncbi_subspecies)

        # get genomes with NCBI subspecies assignment
        ncbi_all_subspecies = set()
        ncbi_subspecies_gids = defaultdict(list)
        ncbi_all_subspecies_gids = set()
        for gid in self.cur_genomes:
            ncbi_subspecies = self.cur_genomes[gid].ncbi_unfiltered_taxa.subspecies
            if ncbi_subspecies:
                ncbi_all_subspecies.add(ncbi_subspecies)
                ncbi_subspecies_gids[ncbi_subspecies].append(gid)
                ncbi_all_subspecies_gids.add(gid)

        # process NCBI subspecies with type genomes followed by
        # remaining NCBI subspecies
        ncbi_type_subspecies = []
        ncbi_nontype_subspecies = []
        for ncbi_subspecies, ncbi_sp_gids in ncbi_subspecies_gids.items():
            has_type_genome = False
            for gid in ncbi_sp_gids:
                if self.cur_genomes[gid].is_ncbi_type_subspecies():
                    has_type_genome = True
                    break

            if has_type_genome:
                ncbi_type_subspecies.append(ncbi_subspecies)
            else:
                ncbi_nontype_subspecies.append(ncbi_subspecies)
        self.logger.info(' - identified {:,} NCBI species with effective type strain genome and {:,} NCBI species without a type strain genome.'.format(
            len(ncbi_type_subspecies), len(ncbi_nontype_subspecies)))

        # establish if NCBI subspecies can be unambiguously assigned to a GTDB species cluster
        fout = open(os.path.join(self.output_dir,
                                 'ncbi_subsp_classification.tsv'), 'w')
        fout.write(
            'NCBI taxon\tClassification\tAssignment type\tGTDB representative')
        fout.write(
            '\t% NCBI subspecies in cluster\t% cluster with NCBI subspecies\tNCBI classifications in cluster')
        fout.write('\tType strains\tNote\n')
        unambiguous_ncbi_subsp = {}
        ambiguous_ncbi_subsp = set()
        unambiguous_rids = set()
        skipped_subspecies = 0
        for ncbi_subspecies in ncbi_type_subspecies + ncbi_nontype_subspecies:
            # do not consider cases such as Serratia marcescens subsp. marcescens,
            # since this isn't a distinct subspecies of S. marcescens that might
            # be promoted to a different species
            specific = ncbi_subspecies.split()[1]
            subsp = ncbi_subspecies.split()[-1]

            if specific == subsp:
                skipped_subspecies += 1
                continue

            ncbi_subsp_gids = ncbi_subspecies_gids[ncbi_subspecies]

            # get type strains for NCBI species that are GTDB species representatives
            type_rids = ncbi_type_subspecies_rid[ncbi_subspecies]

            if len(type_rids) == 0:
                highest_subsp_in_cluster_perc = 0
                highest_cluster_perc = 0
                highest_gtdb_rid = None
                for gtdb_rid, cids in self.cur_clusters.items():
                    assert gtdb_rid in cids

                    ncbi_cur_subsp_in_cluster = cids.intersection(
                        ncbi_subsp_gids)
                    ncbi_any_subsp_in_cluster = cids.intersection(
                        ncbi_all_subspecies_gids)

                    if len(ncbi_any_subsp_in_cluster) > 0:
                        ncbi_cur_sp_in_cluster_perc = len(
                            ncbi_cur_subsp_in_cluster)*100.0/len(ncbi_subsp_gids)
                        cluster_perc = len(
                            ncbi_cur_subsp_in_cluster)*100.0/len(ncbi_any_subsp_in_cluster)

                        if (ncbi_cur_sp_in_cluster_perc > highest_subsp_in_cluster_perc
                                or (ncbi_cur_sp_in_cluster_perc == highest_subsp_in_cluster_perc
                                    and cluster_perc > highest_cluster_perc)):
                            highest_subsp_in_cluster_perc = ncbi_cur_sp_in_cluster_perc
                            highest_cluster_perc = cluster_perc
                            highest_gtdb_rid = gtdb_rid

                if (len(gtdb_rid_type_subspecies[highest_gtdb_rid]) == 0  # if cluster contains a type subspecies, can't assign via majority vote
                        and highest_gtdb_rid not in unambiguous_rids
                        and highest_subsp_in_cluster_perc > 50
                        and highest_cluster_perc > 50):
                    unambiguous_ncbi_subsp[ncbi_subspecies] = (
                        highest_gtdb_rid, NCBI_SpeciesManager.MAJORITY_VOTE)
                else:
                    ambiguous_ncbi_subsp.add(ncbi_subspecies)

                ncbi_subsp_str = []
                for sp, count in Counter([self.cur_genomes[gid].ncbi_unfiltered_taxa.subspecies for gid in self.cur_clusters[highest_gtdb_rid]]).most_common():
                    ncbi_subsp_str.append('{} ({:,})'.format(sp, count))
                ncbi_subsp_str = '; '.join(ncbi_subsp_str)

                fout.write('{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\t{}\t{}\t{}\n'.format(
                    ncbi_subspecies,
                    NCBI_SpeciesManager.UNAMBIGUOUS if ncbi_subspecies in unambiguous_ncbi_subsp else NCBI_SpeciesManager.AMBIGUOUS,
                    NCBI_SpeciesManager.MAJORITY_VOTE if ncbi_subspecies in unambiguous_ncbi_subsp else NCBI_SpeciesManager.AMBIGUOUS,
                    highest_gtdb_rid,
                    highest_subsp_in_cluster_perc,
                    highest_cluster_perc,
                    ncbi_subsp_str,
                    '', ''))

            elif (len(type_rids) == 1
                  and len(gtdb_rid_type_subspecies[list(type_rids)[0]]) == 1):  # Ideally, should allow clusters with multiple type subspecies and resolve by priority

                gtdb_rid = type_rids.pop()

                unambiguous_ncbi_subsp[ncbi_subspecies] = (
                    gtdb_rid, NCBI_SpeciesManager.TYPE_STRAIN_GENOME)
                unambiguous_rids.add(gtdb_rid)

                cids = self.cur_clusters[gtdb_rid]
                ncbi_cur_subsp_in_cluster = cids.intersection(ncbi_subsp_gids)
                ncbi_any_subsp_in_cluster = cids.intersection(
                    ncbi_all_subspecies_gids)
                cur_sp_in_cluster_perc = len(
                    ncbi_cur_subsp_in_cluster)*100.0/len(ncbi_subsp_gids)
                cluster_perc = len(ncbi_cur_subsp_in_cluster) * \
                    100.0/len(ncbi_any_subsp_in_cluster)

                ncbi_subsp_str = []
                for sp, count in Counter([self.cur_genomes[gid].ncbi_unfiltered_taxa.subspecies for gid in self.cur_clusters[gtdb_rid]]).most_common():
                    ncbi_subsp_str.append('{} ({:,})'.format(sp, count))
                ncbi_subsp_str = '; '.join(ncbi_subsp_str)

                type_strain_gids = [cid for cid in self.cur_clusters[gtdb_rid]
                                    if self.cur_genomes[cid].is_ncbi_type_subspecies()]

                fout.write('{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\t{}\t{}\t{}\n'.format(
                    ncbi_subspecies,
                    NCBI_SpeciesManager.UNAMBIGUOUS,
                    NCBI_SpeciesManager.TYPE_STRAIN_GENOME,
                    gtdb_rid,
                    cur_sp_in_cluster_perc,
                    cluster_perc,
                    ncbi_subsp_str,
                    ','.join(type_strain_gids),
                    ''))
            else:
                if len(type_rids) > 1:
                    warn_str = 'Multiple GTDB species clusters contain type genomes of {}: {}'.format(
                        ncbi_subspecies,
                        type_rids)
                    self.logger.warning(warn_str)
                else:
                    warn_str = 'GTDB species clusters contain multiple subspecies with type strain.'

                ambiguous_ncbi_subsp.add(ncbi_subspecies)
                fout.write('{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\t{}\t{}\t{}\n'.format(
                    ncbi_subspecies,
                    NCBI_SpeciesManager.UNAMBIGUOUS if ncbi_subspecies in unambiguous_ncbi_subsp else NCBI_SpeciesManager.AMBIGUOUS,
                    NCBI_SpeciesManager.MAJORITY_VOTE if ncbi_subspecies in unambiguous_ncbi_subsp else NCBI_SpeciesManager.AMBIGUOUS,
                    'n/a',
                    0,
                    0,
                    'n/a',
                    'n/a',
                    warn_str))

        fout.close()

        # sanity check results
        assert set(unambiguous_ncbi_subsp).intersection(
            ambiguous_ncbi_subsp) == set()
        assert len(ncbi_all_subspecies) == len(unambiguous_ncbi_subsp) + \
            len(ambiguous_ncbi_subsp) + skipped_subspecies

        return ncbi_all_subspecies, unambiguous_ncbi_subsp, ambiguous_ncbi_subsp

    def parse_ncbi_classification_table(self, ncbi_classification_table):
        """Parse NCBI classification table."""

        ncbi_all_sp = set()
        unambiguous_ncbi_sp = {}
        ambiguous_ncbi_sp = set()
        with open(ncbi_classification_table) as f:
            header = f.readline().strip().split('\t')

            ncbi_sp_idx = header.index('NCBI taxon')
            classification_idx = header.index('Classification')
            assignment_idx = header.index('Assignment type')
            rid_idx = header.index('GTDB representative')

            for line in f:
                tokens = line.strip().split('\t')

                ncbi_sp = tokens[ncbi_sp_idx]
                ncbi_all_sp.add(ncbi_sp)

                classification = tokens[classification_idx]
                if classification == NCBI_SpeciesManager.UNAMBIGUOUS:
                    gtdb_rid = tokens[rid_idx]
                    assignment_type = tokens[assignment_idx]
                    unambiguous_ncbi_sp[ncbi_sp] = (gtdb_rid, assignment_type)
                elif classification == NCBI_SpeciesManager.AMBIGUOUS:
                    ambiguous_ncbi_sp.add(ncbi_sp)

        return ncbi_all_sp, unambiguous_ncbi_sp, ambiguous_ncbi_sp
