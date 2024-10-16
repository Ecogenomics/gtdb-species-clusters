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
from collections import defaultdict

from numpy import (mean as np_mean, std as np_std)

from gtdb_species_clusters import defaults as Defaults
from gtdb_species_clusters.skani import Skani
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_clusters import SpeciesClusters
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters.genome_utils import select_highest_quality, parse_ncbi_bioproject


class RepActions():
    """Perform initial actions required for changed representatives."""

    # ANI and AF thresholds for determining that an updated genome is sufficiently
    # similar to previous genome that it should be effectively treated as unchanged
    GENOMIC_UPDATE_ANI = 99.0
    GENOMIC_UPDATE_AF = 0.80

    # increase in ANI score require to select new representative genome
    NEW_REP_QC_THRESHOLD = 10

    def __init__(self, ani_cache_file, cpus, output_dir):
        """Initialization."""

        self.output_dir = output_dir
        self.log = logging.getLogger('rich')

        self.skani = Skani(ani_cache_file, cpus)

        # create log indicate actions taken to GTDB
        # representative genomes
        self.action_log = open(os.path.join(
            self.output_dir, 'action_log.tsv'), 'w')
        self.action_log.write(
            'Genome ID\tPrevious GTDB species\tAction\tParameters\n')

        self.new_reps = {}
        self.sp_priority_mngr = None

    def rep_change_gids(self, rep_change_summary_file, field, value):
        """Get genomes with a specific change."""

        gids = {}
        with open(rep_change_summary_file) as f:
            header = f.readline().strip().split('\t')

            field_index = header.index(field)
            prev_sp_index = header.index('Previous GTDB species')

            for line in f:
                line_split = line.strip().split('\t')

                v = line_split[field_index]
                if v == value:
                    prev_sp = line_split[prev_sp_index]
                    gids[line_split[0]] = prev_sp

        return gids

    def top_ani_score_prev_rep(self,
                               prev_rid,
                               sp_cids,
                               prev_genomes,
                               cur_genomes):
        """Identify genome in cluster with highest balanced ANI score to genomic file of representative in previous GTDB release."""

        # calculate ANI to all genomes in species
        gid_pairs = []
        genome_paths = {}
        for cid in sp_cids:
            gid_pairs.append((f'{prev_rid}-P', f'{cid}-C'))
            genome_paths[f'{prev_rid}-P'] = prev_genomes[prev_rid].genomic_file
            genome_paths[f'{cid}-C'] = cur_genomes[cid].genomic_file

        ani_af = self.skani.pairs(
            gid_pairs, 
            genome_paths, 
            preset = Defaults.SKANI_PRESET,
            report_progress=False)

        # determine highest ANI score
        max_score = -1e6
        max_rid = None
        max_ani = None
        max_af = None
        for cid in sp_cids:
            ani, af = Skani.symmetric_ani_af(
                ani_af, f'{prev_rid}-P', f'{cid}-C')

            cur_score = cur_genomes[cid].score_ani(ani)
            if (cur_score > max_score
                    or (cur_score == max_score and ani > max_ani)):
                max_score = cur_score
                max_rid = cid
                max_ani = ani
                max_af = af

        return max_rid, max_score, max_ani, max_af

    def top_ani_score(self,
                      prev_rid,
                      sp_cids,
                      cur_genomes):
        """Identify genome in cluster with highest balanced ANI score to representative genome."""

        # calculate ANI between representative and genomes in species cluster
        gid_pairs = []
        for cid in sp_cids:
            gid_pairs.append((cid, prev_rid))

        ani_af = self.skani.pairs(gid_pairs,
                                    cur_genomes.genomic_files,
                                    preset = Defaults.SKANI_PRESET,
                                    report_progress=False)

        # find genome with top ANI score
        max_score = -1e6
        max_rid = None
        max_ani = None
        max_af = None
        for cid in sp_cids:
            ani, af = Skani.symmetric_ani_af(ani_af, prev_rid, cid)

            cur_score = cur_genomes[cid].score_ani(ani)

            if self.violates_naming_priority(cur_genomes, prev_rid, cid):
                # skip genomes that result in a naming priority violation
                continue

            if cur_score > max_score:
                max_score = cur_score
                max_rid = cid
                max_ani = ani
                max_af = af

        return max_rid, max_score, max_ani, max_af

    def get_updated_rid(self, prev_rid):
        """Get updated representative."""

        if prev_rid in self.new_reps:
            gid, _action = self.new_reps[prev_rid]
            return gid

        return prev_rid

    def update_rep(self, prev_rid, new_rid, action):
        """Update representative genome for GTDB species cluster."""

        if prev_rid in self.new_reps and self.new_reps[prev_rid][0] != new_rid:
            self.log.warning('Representative {} was reassigned multiple times: {} {}.'.format(
                prev_rid, self.new_reps[prev_rid], (new_rid, action)))
            self.log.warning(
                'Assuming last reassignment of {}: {} has priority.'.format(new_rid, action))

        self.new_reps[prev_rid] = (new_rid, action)

    def genomes_in_current_sp_cluster(self,
                                      prev_rid,
                                      prev_genomes,
                                      new_updated_sp_clusters,
                                      cur_genomes):
        """Get genomes in current species cluster."""

        assert prev_rid in prev_genomes.sp_clusters

        sp_cids = prev_genomes.sp_clusters[prev_rid]
        if prev_rid in new_updated_sp_clusters:
            sp_cids = sp_cids.union(new_updated_sp_clusters[prev_rid])
        sp_cids = sp_cids.intersection(cur_genomes)

        return sp_cids

    def violates_naming_priority(self, cur_genomes, prev_rid, new_rid):
        """Check if proposed change in GTDB representative conflicts with naming priority."""

        if (cur_genomes[prev_rid].is_gtdb_type_strain()
                and cur_genomes[prev_rid].ncbi_taxa.specific_epithet != cur_genomes[new_rid].ncbi_taxa.specific_epithet
                and self.sp_priority_mngr.test_species_priority(cur_genomes, prev_rid, new_rid)):
            # Previous representative has naming priority over proposed new representative
            # which is a violations

            if False:  # *** DEBUGGING
                self.log.warning('Reassignments to type strain genome with lower naming priority is not allowed: {}/{}/{}, {}/{}/{}'.format(
                    prev_rid,
                    cur_genomes[prev_rid].ncbi_taxa.species,
                    cur_genomes[prev_rid].year_of_priority(),
                    new_rid,
                    cur_genomes[new_rid].ncbi_taxa.species,
                    cur_genomes[new_rid].year_of_priority()))

            return True

        return False

    def action_genomic_lost(self,
                            rep_change_summary_file,
                            prev_genomes,
                            cur_genomes,
                            new_updated_sp_clusters):
        """Handle species with lost representative genome."""

        # get genomes with specific changes
        self.log.info(
            'Identifying species with lost representative genome:')
        genomic_lost_rids = self.rep_change_gids(rep_change_summary_file,
                                                 'GENOMIC_CHANGE',
                                                 'LOST')
        self.log.info(f' - identified {len(genomic_lost_rids):,} genomes')

        # calculate ANI between previous and current genomes
        for idx, (prev_rid, prev_gtdb_sp) in enumerate(genomic_lost_rids.items()):
            sp_cids = self.genomes_in_current_sp_cluster(prev_rid,
                                                         prev_genomes,
                                                         new_updated_sp_clusters,
                                                         cur_genomes)

            params = {}
            if sp_cids:
                action = 'GENOMIC_CHANGE:LOST:REPLACED'

                new_rid, _top_score, ani, af = self.top_ani_score_prev_rep(prev_rid,
                                                                           sp_cids,
                                                                           prev_genomes,
                                                                           cur_genomes)
                assert new_rid != prev_rid

                params['new_rid'] = new_rid
                params['ani'] = ani
                params['af'] = af
                params['new_assembly_quality'] = cur_genomes[new_rid].score_assembly()
                params['prev_assembly_quality'] = prev_genomes[prev_rid].score_assembly()

                self.update_rep(prev_rid, new_rid, action)
            else:
                action = 'GENOMIC_CHANGE:LOST:SPECIES_RETIRED'
                self.update_rep(prev_rid, None, action)

            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                prev_rid,
                prev_gtdb_sp,
                action,
                params))

            status_str = '-> processed {:,} of {:,} ({:.2f}%) species'.format(
                idx+1,
                len(genomic_lost_rids),
                float(idx+1)*100/len(genomic_lost_rids))
            sys.stdout.write('\r\033[K')  # clear line
            sys.stdout.write(f'{status_str}')
            sys.stdout.flush()

        sys.stdout.write('\n')

    def action_disband_cluster(self, rep_change_summary_file):
        """Handle species clusters explicitly indicated as needing to be disbanded."""

        # get genomes with specific changes
        self.log.info(
            'Identifying species explicitly marked to be disbanded:')
        disbanded_rids = self.rep_change_gids(rep_change_summary_file,
                                              'DISBANDED_CHECK',
                                              'TRUE')
        self.log.info(f' - identified {len(disbanded_rids):,} genomes')

        # log disbanded genomes
        for prev_rid, prev_gtdb_sp in disbanded_rids.items():
            params = {}
            action = 'EXPLICIT_UPDATE:DISBANDED'
            self.update_rep(prev_rid, None, action)

            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                prev_rid,
                prev_gtdb_sp,
                action,
                params))

        return disbanded_rids

    def action_genomic_update(self,
                              rep_change_summary_file,
                              prev_genomes,
                              cur_genomes,
                              new_updated_sp_clusters,
                              disbanded_rids):
        """Handle representatives with updated genomes."""

        # get genomes with specific changes
        self.log.info(
            'Identifying representatives with updated genomic files:')
        genomic_update_gids = self.rep_change_gids(rep_change_summary_file,
                                                   'GENOMIC_CHANGE',
                                                   'UPDATED')
        self.log.info(
            f' - identified {len(genomic_update_gids):,} genomes')

        # calculate ANI between previous and current genomes
        assembly_score_change = []
        for prev_rid, prev_gtdb_sp in genomic_update_gids.items():
            if prev_rid in disbanded_rids:
                continue

            # check that genome hasn't been lost which should
            # be handled differently
            assert prev_rid in cur_genomes

            ani, af = self.skani.calculate(f'{prev_rid}-P', f'{prev_rid}-C',
                                                        prev_genomes[prev_rid].genomic_file,
                                                        cur_genomes[prev_rid].genomic_file)

            params = {}
            params['ani'] = ani
            params['af'] = af
            params['prev_ncbi_accession'] = prev_genomes[prev_rid].ncbi_accn
            params['cur_ncbi_accession'] = cur_genomes[prev_rid].ncbi_accn
            assert prev_genomes[prev_rid].ncbi_accn != cur_genomes[prev_rid].ncbi_accn

            if ani >= RepActions.GENOMIC_UPDATE_ANI and af >= RepActions.GENOMIC_UPDATE_AF:
                params['prev_assembly_quality'] = prev_genomes[prev_rid].score_assembly()
                params['new_assembly_quality'] = cur_genomes[prev_rid].score_assembly()
                action = 'GENOMIC_CHANGE:UPDATED:MINOR_CHANGE'

                d = (cur_genomes[prev_rid].score_assembly() -
                     prev_genomes[prev_rid].score_assembly())
                assembly_score_change.append(d)
            else:
                sp_cids = self.genomes_in_current_sp_cluster(prev_rid,
                                                             prev_genomes,
                                                             new_updated_sp_clusters,
                                                             cur_genomes)

                if sp_cids:
                    new_rid, _top_score, ani, af = self.top_ani_score_prev_rep(prev_rid,
                                                                               sp_cids,
                                                                               prev_genomes,
                                                                               cur_genomes)

                    if new_rid == prev_rid:
                        params['prev_assembly_quality'] = prev_genomes[prev_rid].score_assembly(
                        )
                        params['new_assembly_quality'] = cur_genomes[prev_rid].score_assembly(
                        )
                        action = 'GENOMIC_CHANGE:UPDATED:RETAINED'
                    else:
                        action = 'GENOMIC_CHANGE:UPDATED:REPLACED'
                        params['new_rid'] = new_rid
                        params['ani'] = ani
                        params['af'] = af
                        params['new_assembly_quality'] = cur_genomes[new_rid].score_assembly()
                        params['prev_assembly_quality'] = prev_genomes[prev_rid].score_assembly(
                        )

                        self.update_rep(prev_rid, new_rid, action)
                else:
                    action = 'GENOMIC_CHANGE:UPDATED:SPECIES_RETIRED'
                    self.update_rep(prev_rid, None, action)

            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                prev_rid,
                prev_gtdb_sp,
                action,
                params))

        self.log.info(' - change in assembly score for updated genomes: {:.2f} +/- {:.2f}'.format(
            np_mean(assembly_score_change),
            np_std(assembly_score_change)))

    def action_type_strain_lost(self,
                                rep_change_summary_file,
                                prev_genomes,
                                cur_genomes,
                                new_updated_sp_clusters,
                                disbanded_rids):
        """Handle representatives which have lost type strain genome status."""

        # get genomes with new NCBI species assignments
        self.log.info(
            'Identifying representatives that lost type strain genome status:')
        ncbi_type_species_lost = self.rep_change_gids(rep_change_summary_file,
                                                      'TYPE_STRAIN_CHANGE',
                                                      'LOST')
        self.log.info(
            f' - identified {len(ncbi_type_species_lost):,} genomes')

        for prev_rid, prev_gtdb_sp in ncbi_type_species_lost.items():
            if prev_rid in disbanded_rids:
                continue

            # check that genome hasn't been lost which should
            # be handled differently
            assert prev_rid in cur_genomes

            sp_cids = self.genomes_in_current_sp_cluster(prev_rid,
                                                         prev_genomes,
                                                         new_updated_sp_clusters,
                                                         cur_genomes)

            prev_rep_score = cur_genomes[prev_rid].score_ani(100)
            new_rid, top_score, ani, af = self.top_ani_score(prev_rid,
                                                             sp_cids,
                                                             cur_genomes)

            params = {}
            params['prev_rid_prev_strain_ids'] = prev_genomes[prev_rid].ncbi_strain_identifiers
            params['prev_rid_cur_strain_ids'] = cur_genomes[prev_rid].ncbi_strain_identifiers
            params['prev_rid_prev_gtdb_type_designation'] = prev_genomes[prev_rid].gtdb_type_designation
            params['prev_rid_cur_gtdb_type_designation'] = cur_genomes[prev_rid].gtdb_type_designation
            params['prev_rid_prev_gtdb_type_designation_sources'] = prev_genomes[prev_rid].gtdb_type_designation_sources
            params['prev_rid_cur_gtdb_type_designation_sources'] = cur_genomes[prev_rid].gtdb_type_designation_sources

            if top_score > prev_rep_score:
                action = 'TYPE_STRAIN_CHANGE:LOST:REPLACED'
                assert prev_rid != new_rid

                params['new_rid'] = new_rid
                params['ani'] = ani
                params['af'] = af
                params['new_assembly_quality'] = cur_genomes[new_rid].score_assembly()
                params['prev_assembly_quality'] = prev_genomes[prev_rid].score_assembly()

                params['new_rid_strain_ids'] = cur_genomes[new_rid].ncbi_strain_identifiers
                params['new_rid_gtdb_type_designation'] = cur_genomes[new_rid].gtdb_type_designation
                params['new_rid_gtdb_type_designation_sources'] = cur_genomes[new_rid].gtdb_type_designation_sources

                self.update_rep(prev_rid, new_rid, action)
            else:
                action = 'TYPE_STRAIN_CHANGE:LOST:RETAINED'

            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                prev_rid,
                prev_gtdb_sp,
                action,
                params))

    def action_domain_change(self,
                             rep_change_summary_file,
                             prev_genomes,
                             cur_genomes,
                             disbanded_rids):
        """Handle representatives which have new domain assignments."""

        # get genomes with new NCBI species assignments
        self.log.info(
            'Identifying representative with new domain assignments:')
        domain_changed = self.rep_change_gids(rep_change_summary_file,
                                              'DOMAIN_CHECK',
                                              'REASSIGNED')
        self.log.info(f' - identified {len(domain_changed):,} genomes')

        for prev_rid, prev_gtdb_sp in domain_changed.items():
            if prev_rid in disbanded_rids:
                continue

            action = 'DOMAIN_CHECK:REASSIGNED'
            params = {}
            params['prev_gtdb_domain'] = prev_genomes[prev_rid].gtdb_taxa.domain
            params['cur_gtdb_domain'] = cur_genomes[prev_rid].gtdb_taxa.domain

            self.update_rep(prev_rid, None, action)
            self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                prev_rid,
                prev_gtdb_sp,
                action,
                params))

    def action_ncbi_anomalous_assemblies(self,
                                         rep_change_summary_file,
                                         prev_genomes,
                                         cur_genomes,
                                         new_updated_sp_clusters,
                                         disbanded_rids):
        """Check if representative considered to have an anomalous assembly at NCBI should be replaced."""

        # get genomes with specific changes
        self.log.info('Identifying genomes marked as problematic at NCBI:')
        ncbi_problematic_rids = self.rep_change_gids(rep_change_summary_file,
                                                     'NCBI_ASSEMBLY_QUALITY',
                                                     'NCBI_ANOMALOUS_ASSEMBLY')
        self.log.info(
            f' - identified {len(ncbi_problematic_rids):,} genomes')

        anis = []
        afs = []
        num_frameshifted_proteins = 0
        num_anomalous = 0
        for idx, (prev_rid, prev_gtdb_sp) in enumerate(ncbi_problematic_rids.items()):
            if prev_rid in disbanded_rids:
                continue

            # check that genome hasn't been lost which should be handled differently
            assert prev_rid in cur_genomes

            sp_cids = self.genomes_in_current_sp_cluster(prev_rid,
                                                         prev_genomes,
                                                         new_updated_sp_clusters,
                                                         cur_genomes)

            status_str = '-> processing {:,} of {:,} ({:.2f}%) species [{}: {:,} genomes].'.format(
                idx+1,
                len(ncbi_problematic_rids),
                float(idx+1)*100/len(ncbi_problematic_rids),
                prev_gtdb_sp,
                len(sp_cids))
            sys.stdout.write('\r\033[K')  # clear line
            sys.stdout.write(f'{status_str}')
            sys.stdout.flush()

            # get latest representative of GTDB species clusters as it may
            # have been updated by a previous update rule
            prev_updated_rid = self.get_updated_rid(prev_rid)
            if prev_updated_rid is None:
                self.log.error(
                    "Returned 'None' as updated representative: {prev_rid}")
                sys.exit(-1)

            prev_rep_score = cur_genomes[prev_rid].score_ani(100)
            new_rid, top_score, ani, af = self.top_ani_score(prev_rid,
                                                             sp_cids,
                                                             cur_genomes)

            if top_score > prev_rep_score and prev_updated_rid != new_rid:
                action = 'NCBI_ANOMALOUS_ASSEMBLY:REPLACED:HIGHER_QS'

                params = {}
                params['new_rid'] = new_rid
                params['ani'] = ani
                params['af'] = af
                params['new_assembly_quality'] = cur_genomes[new_rid].score_assembly()
                params['prev_assembly_quality'] = cur_genomes[prev_updated_rid].score_assembly()
                params['new_gtdb_type_of_species'] = cur_genomes[new_rid].is_gtdb_type_strain()
                params['prev_gtdb_type_of_species'] = cur_genomes[prev_updated_rid].is_gtdb_type_strain()
                params['new_ncbi_type_of_species'] = cur_genomes[new_rid].is_ncbi_type_strain()
                params['prev_ncbi_type_of_species'] = cur_genomes[prev_updated_rid].is_ncbi_type_strain()
                params['new_gtdb_type_of_subsp'] = cur_genomes[new_rid].is_gtdb_type_subspecies()
                params['prev_gtdb_type_of_subsp'] = cur_genomes[prev_updated_rid].is_gtdb_type_subspecies()
                params['new_ncbi_type_of_subsp'] = cur_genomes[new_rid].is_ncbi_type_subspecies()
                params['prev_ncbi_type_of_subsp'] = cur_genomes[prev_updated_rid].is_ncbi_type_subspecies()

                anis.append(ani)
                afs.append(af)

                if cur_genomes[prev_rid].is_ncbi_many_frameshifted_proteins():
                    num_frameshifted_proteins += 1

                if cur_genomes[prev_rid].is_ncbi_anomalous_assembly():
                    num_anomalous += 1

                self.update_rep(prev_rid, new_rid, action)

                self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                    prev_rid,
                    prev_gtdb_sp,
                    action,
                    params))

        sys.stdout.write('\n')

        self.log.info(
            f' - identified {len(anis):,} species with improved representatives')
        self.log.info(
            f'   - {num_frameshifted_proteins:,} marked as having many frameshifted proteins')
        self.log.info(
            f'   - {num_anomalous:,} marked as being an anomalous assembly')
        self.log.info(
            f' - ANI = {np_mean(anis):.2f} +/- {np_std(anis):.2f}%; AF = {np_mean(afs):.2f} +/- {np_std(afs):.2f}%')

    def action_improved_rep(self,
                            cur_genomes,
                            new_updated_sp_clusters,
                            disbanded_rids,
                            ncbi_genbank_assembly_file):
        """Check if representative should be replace with higher quality genome."""

        self.log.info(
            'Identifying improved representatives for GTDB species clusters:')

        ncbi_bioproject = parse_ncbi_bioproject(ncbi_genbank_assembly_file)

        num_gtdb_ncbi_type_sp = 0
        num_gtdb_type_sp = 0
        num_ncbi_type_sp = 0
        num_ncbi_reps = 0
        num_complete = 0
        num_isolate = 0
        anis = []
        afs = []
        improved_reps = {}
        ncbi_bioproject_count = defaultdict(int)
        for idx, (prev_rid, cids) in enumerate(new_updated_sp_clusters.clusters()):
            if prev_rid not in cur_genomes or prev_rid in disbanded_rids:
                # indicates genome has been lost or disbanded
                continue

            prev_gtdb_sp = new_updated_sp_clusters.get_species(prev_rid)
            status_str = '-> processing {:,} of {:,} ({:.2f}%) species [{}: {:,} new/updated genomes].'.format(
                idx+1,
                len(new_updated_sp_clusters),
                float(idx+1)*100/len(new_updated_sp_clusters),
                prev_gtdb_sp,
                len(cids))
            sys.stdout.write('\r\033[K')  # clear line
            sys.stdout.write(f'{status_str}')
            sys.stdout.flush()

            # get latest representative of GTDB species clusters as it may
            # have been updated by a previous update rule
            prev_updated_rid = self.get_updated_rid(prev_rid)

            prev_rep_score = cur_genomes[prev_updated_rid].score_ani(100)
            new_rid, top_score, ani, af = self.top_ani_score(prev_updated_rid,
                                                             cids,
                                                             cur_genomes)

            if top_score > prev_rep_score + RepActions.NEW_REP_QC_THRESHOLD:
                assert prev_updated_rid != new_rid

                action = 'IMPROVED_REP:REPLACED:HIGHER_QS'

                params = {}
                params['new_rid'] = new_rid
                params['ani'] = ani
                params['af'] = af
                params['new_assembly_quality'] = cur_genomes[new_rid].score_assembly()
                params['prev_assembly_quality'] = cur_genomes[prev_updated_rid].score_assembly()
                params['new_gtdb_type_of_species'] = cur_genomes[new_rid].is_gtdb_type_strain()
                params['prev_gtdb_type_of_species'] = cur_genomes[prev_updated_rid].is_gtdb_type_strain()
                params['new_ncbi_type_of_species'] = cur_genomes[new_rid].is_ncbi_type_strain()
                params['prev_ncbi_type_of_species'] = cur_genomes[prev_updated_rid].is_ncbi_type_strain()
                params['new_gtdb_type_of_subsp'] = cur_genomes[new_rid].is_gtdb_type_subspecies()
                params['prev_gtdb_type_of_subsp'] = cur_genomes[prev_updated_rid].is_gtdb_type_subspecies()
                params['new_ncbi_type_of_subsp'] = cur_genomes[new_rid].is_ncbi_type_subspecies()
                params['prev_ncbi_type_of_subsp'] = cur_genomes[prev_updated_rid].is_ncbi_type_subspecies()

                anis.append(ani)
                afs.append(af)

                improvement_list = []
                gtdb_type_sp_improv = cur_genomes[new_rid].is_gtdb_type_strain(
                ) and not cur_genomes[prev_updated_rid].is_gtdb_type_strain()
                ncbi_type_sp_improv = cur_genomes[new_rid].is_ncbi_type_strain(
                ) and not cur_genomes[prev_updated_rid].is_ncbi_type_strain()
                if gtdb_type_sp_improv and ncbi_type_sp_improv:
                    num_gtdb_ncbi_type_sp += 1
                    improvement_list.append(
                        'replaced with genome from type strain of species according to GTDB and NCBI')
                elif gtdb_type_sp_improv:
                    num_gtdb_type_sp += 1
                    improvement_list.append(
                        'replaced with genome from type strain of species according to GTDB')
                elif ncbi_type_sp_improv:
                    num_ncbi_type_sp += 1
                    improvement_list.append(
                        'replaced with genome from type strain of species according to NCBI')

                if cur_genomes[new_rid].is_ncbi_representative() and not cur_genomes[prev_updated_rid].is_ncbi_representative():
                    num_ncbi_reps += 1
                    improvement_list.append(
                        'replaced with genome considered to be NCBI representative of species')

                gtdb_type_subsp_improv = cur_genomes[new_rid].is_gtdb_type_subspecies(
                ) and not cur_genomes[prev_updated_rid].is_gtdb_type_subspecies()
                ncbi_type_subsp_improv = cur_genomes[new_rid].is_ncbi_type_subspecies(
                ) and not cur_genomes[prev_updated_rid].is_ncbi_type_subspecies()
                if gtdb_type_subsp_improv and ncbi_type_subsp_improv:
                    num_gtdb_ncbi_type_sp += 1
                    improvement_list.append(
                        'replaced with genome from type strain of subspecies according to GTDB and NCBI')
                elif gtdb_type_subsp_improv:
                    num_gtdb_type_sp += 1
                    improvement_list.append(
                        'replaced with genome from type strain of subspecies according to GTDB')
                elif ncbi_type_subsp_improv:
                    num_ncbi_type_sp += 1
                    improvement_list.append(
                        'replaced with genome from type strain of subspecies according to NCBI')

                if cur_genomes[new_rid].is_isolate() and not cur_genomes[prev_updated_rid].is_isolate():
                    num_isolate += 1
                    improvement_list.append('MAG/SAG replaced with isolate')
                    if new_rid in ncbi_bioproject:
                        ncbi_bioproject_count[ncbi_bioproject[new_rid]] += 1

                if cur_genomes[new_rid].is_complete_genome() and not cur_genomes[prev_updated_rid].is_complete_genome():
                    num_complete += 1
                    improvement_list.append('replaced with complete genome')

                if len(improvement_list) == 0:
                    improvement_list.append(
                        'replaced with higher quality genome')

                params['improvements'] = '; '.join(improvement_list)

                self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                    prev_rid,
                    prev_gtdb_sp,
                    action,
                    params))

                improved_reps[prev_rid] = (new_rid, action)

        sys.stdout.write('\n')
        self.log.info(
            f' - identified {len(improved_reps):,} species with improved representatives')
        self.log.info(
            f'   - {num_gtdb_ncbi_type_sp:,} replaced with GTDB/NCBI genome from type strain')
        self.log.info(
            f'   - {num_gtdb_type_sp:,} replaced with GTDB genome from type strain')
        self.log.info(
            f'   - {num_ncbi_type_sp:,} replaced with NCBI genome from type strain')
        self.log.info(
            f'   - {num_ncbi_reps:,} replaced with NCBI species representative genome')
        self.log.info(
            f'   - {num_isolate:,} replaced MAG/SAG with isolate')
        self.log.info(
            f'   - {num_complete:,} replaced with complete genome assembly')
        self.log.info(
            f' - ANI = {np_mean(anis):.2f} +/- {np_std(anis):.2f}%; AF = {np_mean(afs):.2f} +/- {np_std(afs):.2f}%')

        # report NCBI BioProject resulting in large numbers of GTDB reps being replaced with presumed isolates
        for bioproject, count in ncbi_bioproject_count.items():
            if count >= 2:
                self.log.warning('BioProject {} responsible for {:,} representatives being replaced with presumed isolates. Should verify these genomes are isolates.'.format(
                    bioproject, count))

        return improved_reps

    def action_naming_priority(self,
                               prev_genomes,
                               cur_genomes,
                               new_updated_sp_clusters,
                               disbanded_rids):
        """Check if representative should be replace with genome with higher nomenclatural priority."""

        self.log.info(
            'Identifying genomes with naming priority in GTDB species clusters:')

        out_file = os.path.join(self.output_dir, 'updated_naming_priority.tsv')
        fout = open(out_file, 'w')
        fout.write('NCBI species\tGTDB species\tRepresentative\tStrain IDs\tRepresentative type sources\tPriority year\tGTDB type species\tGTDB type strain\tNCBI assembly type')
        fout.write('\tNCBI synonym\tGTDB synonym\tSynonym genome\tSynonym strain IDs\tSynonym type sources\tPriority year\tGTDB type species\tGTDB type strain\tSynonym NCBI assembly type')
        fout.write('\tANI\tAF\tPriority note\n')

        num_higher_priority = 0
        assembly_score_change = []
        anis = []
        afs = []
        for prev_rid in prev_genomes.sp_clusters:
            if prev_rid in disbanded_rids:
                continue

            # get type strain genomes in GTDB species cluster, including genomes new to this release
            type_strain_gids = [gid for gid in prev_genomes.sp_clusters[prev_rid]
                                if gid in cur_genomes and cur_genomes[gid].is_effective_type_strain()]
            if prev_rid in new_updated_sp_clusters:
                new_type_strain_gids = [
                    gid for gid in new_updated_sp_clusters[prev_rid] if cur_genomes[gid].is_effective_type_strain()]
                type_strain_gids.extend(new_type_strain_gids)

            if len(type_strain_gids) == 0:
                continue

            # check if representative has already been updated
            updated_rid = self.get_updated_rid(prev_rid)
            if updated_rid is None:
                # indicates genome has explicitly been flagged to be disbanded
                continue

            type_strain_sp = set(
                [cur_genomes[gid].ncbi_taxa.species for gid in type_strain_gids])
            if len(type_strain_sp) == 1 and updated_rid in type_strain_gids:
                continue

            updated_sp = cur_genomes[updated_rid].ncbi_taxa.species
            highest_priority_gid = updated_rid

            if updated_rid not in type_strain_gids:
                highest_priority_gid = None
                if updated_sp in type_strain_sp:
                    sp_gids = [gid for gid in type_strain_gids
                               if cur_genomes[gid].ncbi_taxa.species == updated_sp]
                    hq_gid = select_highest_quality(sp_gids, cur_genomes)
                    highest_priority_gid = hq_gid

                # self.log.warning('Representative is a non-type strain genome even though type strain genomes exist in species cluster: {}: {}, {}: {}'.format(
                #                    prev_rid, cur_genomes[prev_rid].is_effective_type_strain(), updated_rid, cur_genomes[updated_rid].is_effective_type_strain()))
                #self.log.warning('Type strain genomes: {}'.format(','.join(type_strain_gids)))

            # find highest priority genome
            note = ''
            for sp in type_strain_sp:
                if sp == updated_sp:
                    continue

                # get highest quality genome from species
                sp_gids = [gid for gid in type_strain_gids
                           if cur_genomes[gid].ncbi_taxa.species == sp]
                hq_gid = select_highest_quality(sp_gids, cur_genomes)

                if highest_priority_gid is None:
                    highest_priority_gid = hq_gid
                else:
                    highest_priority_gid, note = self.sp_priority_mngr.species_priority(cur_genomes,
                                                                                        highest_priority_gid,
                                                                                        hq_gid)

            # check if representative should be updated
            if highest_priority_gid != updated_rid:
                num_higher_priority += 1

                ani, af = self.skani.calculate(updated_rid,
                                                highest_priority_gid,
                                                cur_genomes[updated_rid].genomic_file,
                                                cur_genomes[highest_priority_gid].genomic_file)

                anis.append(ani)
                afs.append(af)

                d = cur_genomes[highest_priority_gid].score_assembly(
                ) - cur_genomes[updated_rid].score_assembly()
                assembly_score_change.append(d)

                action = 'NOMENCLATURE_PRIORITY:REPLACED'
                params = {}
                params['prev_ncbi_species'] = cur_genomes[updated_rid].ncbi_taxa.species
                params['prev_year_of_priority'] = cur_genomes[updated_rid].year_of_priority()
                params['new_ncbi_species'] = cur_genomes[highest_priority_gid].ncbi_taxa.species
                params['new_year_of_priority'] = cur_genomes[highest_priority_gid].year_of_priority()
                params['new_rid'] = highest_priority_gid
                params['ani'] = ani
                params['af'] = af
                params['priority_note'] = note

                self.update_rep(prev_rid, highest_priority_gid, action)
                self.action_log.write('{}\t{}\t{}\t{}\n'.format(
                    prev_rid,
                    cur_genomes[updated_rid].gtdb_taxa.species,
                    action,
                    params))

                fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                    cur_genomes[highest_priority_gid].ncbi_taxa.species,
                    cur_genomes[highest_priority_gid].gtdb_taxa.species,
                    highest_priority_gid,
                    ','.join(
                        sorted(cur_genomes[highest_priority_gid].strain_ids())),
                    ','.join(sorted(cur_genomes[highest_priority_gid].gtdb_type_sources())).upper(
                    ).replace('STRAININFO', 'StrainInfo'),
                    cur_genomes[highest_priority_gid].year_of_priority(),
                    cur_genomes[highest_priority_gid].is_gtdb_type_species(),
                    cur_genomes[highest_priority_gid].is_gtdb_type_strain(),
                    cur_genomes[highest_priority_gid].ncbi_type_material))
                fout.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                    cur_genomes[updated_rid].ncbi_taxa.species,
                    cur_genomes[updated_rid].gtdb_taxa.species,
                    updated_rid,
                    ','.join(sorted(cur_genomes[updated_rid].strain_ids())),
                    ','.join(sorted(cur_genomes[updated_rid].gtdb_type_sources())).upper().replace(
                        'STRAININFO', 'StrainInfo'),
                    cur_genomes[updated_rid].year_of_priority(),
                    cur_genomes[updated_rid].is_gtdb_type_species(),
                    cur_genomes[updated_rid].is_gtdb_type_strain(),
                    cur_genomes[updated_rid].ncbi_type_material))
                fout.write('\t{:.3f}\t{:.4f}\t{}\n'.format(ani, af, note))

        fout.close()

        self.log.info(
            f' - identified {num_higher_priority:,} species with representative changed to genome with higher nomenclatural priority')
        self.log.info(' - change in assembly score for new representatives: {:.2f} +/- {:.2f}'.format(
            np_mean(assembly_score_change),
            np_std(assembly_score_change)))
        self.log.info(' - ANI: {:.2f} +/- {:.2f}'.format(
            np_mean(anis),
            np_std(anis)))
        self.log.info(' - AF: {:.2f} +/- {:.2f}'.format(
            np_mean(afs),
            np_std(afs)))

    def write_updated_clusters(self,
                               prev_genomes,
                               cur_genomes,
                               new_reps,
                               new_updated_sp_clusters,
                               out_file):
        """Write out updated GTDB species clusters."""

        self.log.info(
            'Writing updated GTDB species clusters to file: {}'.format(out_file))

        fout = open(out_file, 'w')
        fout.write(
            'Representative genome\tGTDB species\tNo. clustered genomes\tClustered genomes\n')

        cur_genome_set = set(cur_genomes)

        num_clusters = 0
        for prev_rid in prev_genomes.sp_clusters:

            new_rid, _action = new_reps.get(prev_rid, [prev_rid, None])
            if new_rid is None:
                continue

            sp_cids = self.genomes_in_current_sp_cluster(prev_rid,
                                                         prev_genomes,
                                                         new_updated_sp_clusters,
                                                         cur_genome_set)

            fout.write('{}\t{}\t{}\t{}\n'.format(
                new_rid,
                prev_genomes.sp_clusters.get_species(prev_rid),
                len(sp_cids),
                ','.join(sp_cids)))
            num_clusters += 1

        fout.close()

        self.log.info(f' - wrote {num_clusters:,} clusters')

    def run(self,
            rep_change_summary_file,
            prev_gtdb_metadata_file,
            prev_genomic_path_file,
            cur_gtdb_metadata_file,
            cur_genomic_path_file,
            genomes_new_updated_file,
            qc_passed_file,
            gtdbtk_classify_file,
            ncbi_genbank_assembly_file,
            untrustworthy_type_file,
            gtdb_type_strains_ledger,
            sp_priority_ledger,
            genus_priority_ledger,
            ncbi_env_bioproject_ledger,
            lpsn_gss_file):
        """Perform initial actions required for changed representatives."""

        # create previous and current GTDB genome sets; the gtdb_type_strains_ledger is not
        # set for the previous genome since we want to compare names without the influence
        # of this file; moreover, the species names for the previous set must match the GTDB-Tk
        # classifications which may not happen if species names are changed by the gtdb_type_strains_ledger 
        self.log.info('Creating previous GTDB genome set:')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                             ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                             untrustworthy_type_ledger=untrustworthy_type_file,
                                             ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)
        self.log.info(' - previous genome set has {:,} species clusters spanning {:,} genomes'.format(
            len(prev_genomes.sp_clusters),
            prev_genomes.sp_clusters.total_num_genomes()))

        self.log.info('Creating current GTDB genome set:')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                            gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                            create_sp_clusters=False,
                                            qc_passed_file=qc_passed_file,
                                            ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                            untrustworthy_type_ledger=untrustworthy_type_file,
                                            ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)

        # get path to previous and current genomic FASTA files
        self.log.info(
            'Reading path to previous and current genomic FASTA files.')
        prev_genomes.load_genomic_file_paths(prev_genomic_path_file)
        cur_genomes.load_genomic_file_paths(cur_genomic_path_file)

        # create expanded previous GTDB species clusters
        new_updated_sp_clusters = SpeciesClusters()

        self.log.info(
            'Creating species clusters of new and updated genomes based on GTDB-Tk classifications:')
        new_updated_sp_clusters.create_expanded_clusters(prev_genomes,
                                                         genomes_new_updated_file,
                                                         qc_passed_file,
                                                         gtdbtk_classify_file)

        self.log.info('Identified {:,} expanded species clusters spanning {:,} genomes.'.format(
            len(new_updated_sp_clusters),
            new_updated_sp_clusters.total_num_genomes()))

        # initialize species priority manager
        self.sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger,
                                                       genus_priority_ledger,
                                                       lpsn_gss_file,
                                                       self.output_dir)

        # take action required for each changed representatives
        self.action_genomic_lost(rep_change_summary_file,
                                 prev_genomes,
                                 cur_genomes,
                                 new_updated_sp_clusters)

        disbanded_rids = self.action_disband_cluster(rep_change_summary_file)

        self.action_genomic_update(rep_change_summary_file,
                                   prev_genomes,
                                   cur_genomes,
                                   new_updated_sp_clusters,
                                   disbanded_rids)

        self.action_type_strain_lost(rep_change_summary_file,
                                     prev_genomes,
                                     cur_genomes,
                                     new_updated_sp_clusters,
                                     disbanded_rids)

        self.action_domain_change(rep_change_summary_file,
                                  prev_genomes,
                                  cur_genomes,
                                  disbanded_rids)

        self.action_ncbi_anomalous_assemblies(rep_change_summary_file,
                                              prev_genomes,
                                              cur_genomes,
                                              new_updated_sp_clusters,
                                              disbanded_rids)

        if True:  # ***DEBUG FLAG
            improved_reps = self.action_improved_rep(cur_genomes,
                                                     new_updated_sp_clusters,
                                                     disbanded_rids,
                                                     ncbi_genbank_assembly_file)

            pickle.dump(improved_reps, open(os.path.join(
                self.output_dir, 'improved_reps.pkl'), 'wb'))
        else:
            self.log.warning(
                'Reading improved_reps for pre-cached file. Generally used only for debugging.')
            improved_reps = pickle.load(
                open(os.path.join(self.output_dir, 'improved_reps.pkl'), 'rb'))

        for prev_rid, (new_rid, action) in improved_reps.items():
            self.update_rep(prev_rid, new_rid, action)

        self.action_naming_priority(prev_genomes,
                                    cur_genomes,
                                    new_updated_sp_clusters,
                                    disbanded_rids)

        # report basic statistics
        num_retired_sp = sum(
            [1 for v in self.new_reps.values() if v[0] is None])
        num_replaced_rids = sum(
            [1 for rid, v in self.new_reps.items() if v[0] is not None and rid in prev_genomes.sp_clusters])
        self.log.info(f'Identified {num_retired_sp:,} retired species.')
        self.log.info(
            f'Identified {num_replaced_rids:,} species with a modified representative genome.')

        fout = open(os.path.join(self.output_dir,
                                 'gtdb_type_sp_untrustworthy_at_ncbi.tsv'), 'w')
        fout.write(
            'Genome ID\tGTDB species\tNCBI species\tSame species\tSame genus\tSame family\tGTDB taxonomy\tNCBI taxonomy\n')
        num_gtdb_type_sp_untrustworthy_at_ncbi = 0
        num_diff_sp = 0
        for prev_rid in prev_genomes.sp_clusters:
            new_rid, action = self.new_reps.get(prev_rid, [prev_rid, None])
            if new_rid is None:
                continue

            if cur_genomes[new_rid].is_gtdb_type_strain() and cur_genomes[new_rid].is_ncbi_untrustworthy_as_type():
                num_gtdb_type_sp_untrustworthy_at_ncbi += 1
                fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    new_rid,
                    cur_genomes[new_rid].gtdb_taxa.species,
                    cur_genomes[new_rid].ncbi_taxa.species,
                    cur_genomes[new_rid].gtdb_taxa.species == cur_genomes[new_rid].ncbi_taxa.species,
                    cur_genomes[new_rid].gtdb_taxa.genus == cur_genomes[new_rid].ncbi_taxa.genus,
                    cur_genomes[new_rid].gtdb_taxa.family == cur_genomes[new_rid].ncbi_taxa.family,
                    cur_genomes[new_rid].gtdb_taxa,
                    cur_genomes[new_rid].ncbi_taxa))

                if cur_genomes[new_rid].gtdb_taxa.species != cur_genomes[new_rid].ncbi_taxa.species:
                    num_diff_sp += 1

        fout.close()
        self.log.info(
            f'Identified {num_gtdb_type_sp_untrustworthy_at_ncbi:,} type strain representatives considered untrustworthy as type at NCBI.')
        if num_diff_sp > 0:
            self.log.warning(
                f' - {num_diff_sp:,} genomes have incongruent GTDB and NCBI species assignments')

        self.action_log.close()

        # write out representatives for existing species clusters
        fout = open(os.path.join(self.output_dir,
                                 'updated_species_reps.tsv'), 'w')
        fout.write(
            'Previous representative ID\tNew representative ID\tAction\tRepresentative status\n')
        for rid in prev_genomes.sp_clusters:
            if rid in self.new_reps:
                new_rid, action = self.new_reps[rid]
                if new_rid is not None:
                    fout.write(f'{rid}\t{new_rid}\t{action}\tREPLACED\n')
                elif rid in disbanded_rids:
                    fout.write(f'{rid}\t{new_rid}\t{action}\tDISBANDED\n')
                else:
                    fout.write(f'{rid}\t{new_rid}\t{action}\tLOST\n')
            else:
                fout.write(f'{rid}\t{rid}\tNONE\tUNCHANGED\n')

        fout.close()

        # write out updated species clusters
        out_file = os.path.join(self.output_dir, 'updated_sp_clusters.tsv')
        self.write_updated_clusters(prev_genomes,
                                    cur_genomes,
                                    self.new_reps,
                                    new_updated_sp_clusters,
                                    out_file)
