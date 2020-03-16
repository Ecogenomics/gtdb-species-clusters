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

from biolib.taxonomy import Taxonomy

from gtdb_species_clusters.fastani import FastANI
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_clusters import SpeciesClusters
from gtdb_species_clusters.species_name_manager import SpeciesNameManager
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters.genome_utils import canonical_gid
from gtdb_species_clusters.type_genome_utils import (symmetric_ani, 
                                                        read_clusters, 
                                                        write_clusters, 
                                                        write_rep_radius)
from gtdb_species_clusters.taxon_utils import (generic_name,
                                                specific_epithet,
                                                canonical_taxon,
                                                parse_synonyms,
                                                longest_common_prefix,
                                                is_placeholder_taxon,
                                                is_placeholder_sp_epithet)


class UpdateSpeciesNames(object):
    """Update names of GTDB species clusters."""

    def __init__(self, ani_cache_file, cpus, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
        self.fastani = FastANI(ani_cache_file, cpus)
        
        self.sp_name_log = open(os.path.join(self.output_dir, 'sp_name_log.tsv'), 'w')
        self.sp_name_log.write('GTDB domain\tGenome ID\tPrevious GTDB species\tNew GTDB species\tAction\n')
        
        self.sp_curation_log = open(os.path.join(self.output_dir, 'sp_curation_log.tsv'), 'w')
        self.sp_curation_log.write('GTDB domain\tGenome ID\tPrevious NCBI species\tCurrent NCBI species')
        self.sp_curation_log.write('\tPrevious GTDB genus\tPrevious GTDB species\tProposed GTDB species')
        self.sp_curation_log.write('\tMandatory fix\tGTDB type strain ledger\tIssue\n')
        
    def update_log(self, gid, cur_genomes, prev_genomes, sp, action):
        """Add entry to update log."""

        prev_gtdb_sp = 'n/a'
        if gid in prev_genomes:
            prev_gtdb_sp = prev_genomes[gid].gtdb_taxa.species

        self.sp_name_log.write('{}\t{}\t{}\t{}\t{}\n'.format(
            cur_genomes[gid].gtdb_taxa.domain,
            gid,
            prev_gtdb_sp,
            sp,
            action))
            
    def curation_log(self, 
                        gid, 
                        cur_genomes, 
                        prev_genomes, 
                        proposed_gtdb_sp, 
                        mandatory, 
                        issue, 
                        in_gtdb_type_strain_ledger):
        """Add entry to curation log."""
        
        prev_ncbi_sp = 'n/a'
        if gid in prev_genomes:
            prev_ncbi_sp = prev_genomes[gid].ncbi_taxa.species
            
        prev_gtdb_sp = 'n/a'
        prev_gtdb_genus = 'n/a'
        if gid in cur_genomes:
            prev_gtdb_sp = cur_genomes[gid].gtdb_taxa.species
            prev_gtdb_genus = cur_genomes[gid].gtdb_taxa.genus
            
        cur_ncbi_sp = cur_genomes[gid].ncbi_taxa.species

        self.sp_curation_log.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                    cur_genomes[gid].gtdb_taxa.domain,
                                    gid,
                                    prev_ncbi_sp,
                                    cur_ncbi_sp,
                                    prev_gtdb_genus,
                                    prev_gtdb_sp,
                                    proposed_gtdb_sp,
                                    mandatory,
                                    in_gtdb_type_strain_ledger,
                                    issue))
                                    
                                    
    def _test_same_epithet(self, epithet1, epithet2):
        """Test if species epithet are the same, except for changes due to difference in the gender of the genus."""
        
        if epithet1 == epithet2:
            return True
            
        lcp = longest_common_prefix(epithet1, epithet2)
        if len(lcp) >= min(len(epithet1), len(epithet2)) - 3:
            # a small change to the suffix presumably reflecting
            # a change in the gender of the genus
            return True
            
        return False
                                    
    def _determine_sp_epithet(self, rid, cur_genomes, synonyms):
        """Determine most appropriate epithet for species name.
        
        This method assumes the epithet of the NCBI species name
        should be favoured; presumable because the genome is 
        assembled from the type strain of the species.
        """
        
        ncbi_sp = cur_genomes[rid].ncbi_taxa.species
        if ncbi_sp in synonyms:
            ncbi_sp = synonyms[ncbi_sp]
        ncbi_genus = 'g__' + generic_name(ncbi_sp)
        
        gtdb_sp = cur_genomes[rid].gtdb_taxa.species
        if gtdb_sp in synonyms:
            gtdb_sp = synonyms[gtdb_sp]
        gtdb_genus = 'g__' + generic_name(gtdb_sp)
        
        gtdb_sp_epithet = specific_epithet(gtdb_sp)
        ncbi_sp_epithet = specific_epithet(ncbi_sp)
        if gtdb_sp_epithet != ncbi_sp_epithet:
            if gtdb_genus == ncbi_genus or is_placeholder_taxon(gtdb_genus):
                # take NCBI specific name since there is no concern about
                # a suffix change due to a change in the gender of the genus
                sp_epithet = specific_epithet(ncbi_sp)
            else:
                # genera do not agree so the suffix of the
                # specific name may have been changed to 
                # reflect the gender of the GTDB genus
                if self._test_same_epithet(gtdb_sp_epithet, ncbi_sp_epithet):
                    # should keep GTDB specific name since it is only
                    # a small change to the suffix presumably reflecting
                    # a change in the gender of the genus
                    sp_epithet = gtdb_sp_epithet
                else:
                    # looks like a real conflict so favour the NCBI
                    # specific name since this is a type strain genome;
                    # there may be an issue with the gender of the GTDB
                    # genus, but this will be resolved later since the
                    # GTDB genus names are ultimately set by manual curation
                    sp_epithet = ncbi_sp_epithet
        else:
            sp_epithet = gtdb_sp_epithet
            
        return sp_epithet
        
    def name_type_species_clusters(self,
                                        rids_by_naming_priority,
                                        cluster_sp_names,
                                        used_sp_names,
                                        clusters,
                                        prev_genomes, 
                                        cur_genomes, 
                                        gtdb_type_strain_ledger,
                                        synonyms):
        """Assign names to clusters represented by a type species of genus genome."""

        num_updates = 0
        num_conflicting_type_sp = 0
        for rid in rids_by_naming_priority:
            if not cur_genomes[rid].is_gtdb_type_species():
                continue

            ncbi_sp = cur_genomes[rid].ncbi_taxa.species
            if ncbi_sp == 's__':
                self.logger.error('NCBI type species genome has no NCBI species assignment: {}'.format(rid))
                sys.exit(-1)
                
            if ncbi_sp in synonyms:
                ncbi_sp = synonyms[ncbi_sp]
                
            if rid in gtdb_type_strain_ledger:
                proposed_gtdb_sp = gtdb_type_strain_ledger[rid]
            else:
                # update GTDB assignment only if it is unassigned or a placeholder name
                gtdb_genus = cur_genomes[rid].gtdb_taxa.genus
                gtdb_sp = cur_genomes[rid].gtdb_taxa.species
                if gtdb_sp in synonyms:
                    gtdb_sp = synonyms[gtdb_sp]
                    
                genus = gtdb_genus
                if is_placeholder_taxon(genus):
                    genus = 'g__' + generic_name(ncbi_sp)
                    
                sp_epithet = self._determine_sp_epithet(rid, cur_genomes, synonyms)

                proposed_gtdb_sp = '{} {}'.format(genus.replace('g__', 's__'), sp_epithet)
                
                assert genus != 'g__' and sp_epithet != ''

            ncbi_genus = 'g__' + generic_name(ncbi_sp)
            proposed_gtdb_genus = 'g__' + generic_name(proposed_gtdb_sp)
            if proposed_gtdb_genus != ncbi_genus:
                num_conflicting_type_sp += 1
                ncbi_year = self.sp_priority_mngr.genus_priority(ncbi_genus)
                gtdb_year = self.sp_priority_mngr.genus_priority(proposed_gtdb_genus)

                if gtdb_year < ncbi_year:
                    self.curation_log(rid, 
                                        cur_genomes, 
                                        prev_genomes, 
                                        proposed_gtdb_sp,
                                        False,
                                        "GTDB genus assignment does not reflect type species of genus",
                                        rid in gtdb_type_strain_ledger)
                else:
                    self.curation_log(rid, 
                                        cur_genomes, 
                                        prev_genomes, 
                                        proposed_gtdb_sp,
                                        True,
                                        "GTDB genus assignment does not reflect type species of genus and may violate naming priority ({}: {}, {}: {})".format(
                                        ncbi_genus, ncbi_year, proposed_gtdb_genus, gtdb_year),
                                        rid in gtdb_type_strain_ledger)

            if proposed_gtdb_sp != cur_genomes[rid].gtdb_taxa.species:
                if rid in prev_genomes:
                    num_updates += 1
                    self.update_log(rid, cur_genomes, prev_genomes, proposed_gtdb_sp, "Updated to reflect genome assembled from type species of genus.")
                for cid in clusters[rid]:
                    #***cur_genomes[cid].gtdb_taxa.genus = 'g__' + generic_name(proposed_gtdb_sp)
                    cur_genomes[cid].gtdb_taxa.species = proposed_gtdb_sp
                    
            if proposed_gtdb_sp in used_sp_names:
                self.logger.error('NCBI species name of type species representative already in use: {}: {}, {}'.format(
                                    ncbi_sp, rid, used_sp_names[proposed_gtdb_sp]))
            
            if proposed_gtdb_sp in synonyms:
                self.logger.error('Assigning synonym to type species representative cluster: {}: {}'.format(proposed_gtdb_sp, rid))

            cluster_sp_names[rid] = proposed_gtdb_sp
            used_sp_names[proposed_gtdb_sp] = rid
       
        self.logger.info(f' ... updated name of {num_updates:,} GTDB clusters.')
        self.logger.warning(f' ... GTDB genus did not reflect type species of genus for {num_conflicting_type_sp:,} species clusters.')

    def name_type_strain_clusters(self,
                                    rids_by_naming_priority,
                                    cluster_sp_names,
                                    used_sp_names,
                                    clusters,
                                    prev_genomes, 
                                    cur_genomes, 
                                    gtdb_type_strain_ledger,
                                    synonyms):
        """Assign names to clusters represented by a validally or effectively published type strain genome."""
                                    
        num_updates = 0
        num_conflicting_specific = 0
        num_conflicting_generic = 0
        dup_name_count = 0
        for rid in rids_by_naming_priority:
            if rid in cluster_sp_names:
                continue
                
            if not cur_genomes[rid].is_effective_type_strain():
                continue
                
            ncbi_sp = cur_genomes[rid].ncbi_taxa.species
            if ncbi_sp == 's__':
                # NCBI does define type material that lacks a binomial species name
                continue
                
            if ncbi_sp in synonyms:
                ncbi_sp = synonyms[ncbi_sp]
                
            if rid in gtdb_type_strain_ledger:
                proposed_gtdb_sp = gtdb_type_strain_ledger[rid]
            else:
                # update GTDB generic name only if it is unassigned, and
                # the GTDB specific name unless it looks like only a 
                # small change to the suffix which likely reflect a
                # change in the gender of the genus
                gtdb_genus = cur_genomes[rid].gtdb_taxa.genus
                gtdb_sp = cur_genomes[rid].gtdb_taxa.species
                if gtdb_sp in synonyms:
                    gtdb_sp = synonyms[gtdb_sp]
                    
                genus = gtdb_genus
                if genus == 'g__':
                    genus = 'g__' + generic_name(ncbi_sp)

                sp_epithet = self._determine_sp_epithet(rid, cur_genomes, synonyms)

                proposed_gtdb_sp = '{} {}'.format(genus.replace('g__', 's__'), sp_epithet)
                
                assert genus != 'g__' and sp_epithet != ''

            ncbi_specific = specific_epithet(ncbi_sp)
            proposed_specific = specific_epithet(proposed_gtdb_sp)
            if not self._test_same_epithet(ncbi_specific, proposed_specific):
                num_conflicting_specific += 1
                self.curation_log(rid, 
                                    cur_genomes, 
                                    prev_genomes, 
                                    proposed_gtdb_sp,
                                    True,
                                    "GTDB specific species name does not reflect type strain of species",
                                    rid in gtdb_type_strain_ledger)
            
            if generic_name(proposed_gtdb_sp) != generic_name(ncbi_sp):
                mandatory = False
                issue = "GTDB generic name does not reflect type strain of species"
                if rid in prev_genomes and generic_name(prev_genomes[rid].ncbi_taxa.species) != generic_name(ncbi_sp):
                    issue += '; NCBI generic name has been updated'
                    mandatory = True
                    num_conflicting_generic += 1
                    self.curation_log(rid, 
                                        cur_genomes, 
                                        prev_genomes, 
                                        proposed_gtdb_sp,
                                        mandatory,
                                        issue,
                                        rid in gtdb_type_strain_ledger)
                                        
            if proposed_gtdb_sp != cur_genomes[rid].gtdb_taxa.species:
                if rid in prev_genomes:
                    num_updates += 1
                    self.update_log(rid, cur_genomes, prev_genomes, proposed_gtdb_sp, "Updated to reflect genome assembled from type strain of species.")
                for cid in clusters[rid]:
                    #***cur_genomes[cid].gtdb_taxa.genus = 'g__' + generic_name(proposed_gtdb_sp)
                    cur_genomes[cid].gtdb_taxa.species = proposed_gtdb_sp
            
            if proposed_gtdb_sp in used_sp_names:
                dup_name_count += 1
                self.curation_log(rid, 
                                    cur_genomes, 
                                    prev_genomes, 
                                    proposed_gtdb_sp,
                                    True,
                                    "Assigned species name supported by multiple type strain genomes: {} {}".format(
                                        rid, used_sp_names[proposed_gtdb_sp]),
                                    rid in gtdb_type_strain_ledger)
                self.logger.warning('NCBI species name of type strain representative already in use: {}: {}, {}'.format(
                                    proposed_gtdb_sp, rid, used_sp_names[proposed_gtdb_sp]))
            
            if proposed_gtdb_sp in synonyms:
                self.logger.error('Assigning synonym to type strain representative cluster: {}: {}'.format(proposed_gtdb_sp, rid))
                sys.exit(-1)

            cluster_sp_names[rid] = proposed_gtdb_sp
            used_sp_names[proposed_gtdb_sp] = rid
            
        self.logger.info(f' ... updated name of {num_updates:,} GTDB clusters.')
        self.logger.warning(f' ... GTDB generic name does not reflect type strain of species and NCBI genus assignment changed this release for {num_conflicting_generic:,} species clusters.')
        self.logger.warning(f' ... GTDB specific name does not reflect type strain of species for {num_conflicting_specific:,} species clusters.')
        self.logger.warning(f' ... same name assigned to {dup_name_count:,} species clusters.')
        
    def name_binomial_clusters(self,
                                rids_by_naming_priority,
                                cluster_sp_names,
                                used_sp_names,
                                clusters,
                                prev_genomes, 
                                cur_genomes, 
                                synonyms):
        """Assign names to clusters with suitable binomial Latin names which lack type strain genomes."""
        
        # get species assignment of previous GTDB species clusters
        gtdb_canonical_sp_epithets = defaultdict(set)
        for rid in clusters:
            gtdb_genus = cur_genomes[rid].gtdb_taxa.genus
            gtdb_sp_epithet = cur_genomes[rid].gtdb_taxa.specific_epithet
            gtdb_canonical_sp_epithets[gtdb_genus].add(canonical_taxon(gtdb_sp_epithet))
        
        num_updates = 0
        dup_name_count = 0
        for rid in rids_by_naming_priority:
            if rid in cluster_sp_names:
                continue

            ncbi_sp = cur_genomes[rid].ncbi_taxa.species
            if is_placeholder_taxon(ncbi_sp):
                continue
                
            if ncbi_sp in synonyms:
                ncbi_sp = synonyms[ncbi_sp]
                
            # update GTDB generic name only if unassigned, and
            # GTDB specific name only if no other species has 
            # this assignment in the genus
            gtdb_genus = cur_genomes[rid].gtdb_taxa.genus
            gtdb_sp = cur_genomes[rid].gtdb_taxa.species
            if gtdb_sp in synonyms:
                gtdb_sp = synonyms[gtdb_sp]
                
            genus = gtdb_genus
            if genus == 'g__':
                genus = 'g__' + generic_name(ncbi_sp)

            gtdb_sp_epithet = specific_epithet(gtdb_sp)
            ncbi_sp_epithet = specific_epithet(ncbi_sp)
            
            # *** HACK: this isn't stable since a given epithet could
            # occur twice
            if ncbi_sp_epithet not in gtdb_canonical_sp_epithets[gtdb_genus]:
                # should update to NCBI epithet in order to reflect 
                # recognized species name
                sp_epithet = ncbi_sp_epithet
            else:
                sp_epithet = gtdb_sp_epithet
                if not gtdb_sp_epithet:
                    sp = '{} {}'.format(genus.replace('g__', 's__'), ncbi_sp_epithet)
                    next_suffix = self.sp_name_mngr.taxon_suffix_manager.next_suffix(sp)
                    sp_epithet = '{}_{}'.format(ncbi_sp_epithet, next_suffix)

            proposed_gtdb_sp = '{} {}'.format(genus.replace('g__', 's__'), sp_epithet)
            if proposed_gtdb_sp in used_sp_names:
                # better suffix name
                next_suffix = self.sp_name_mngr.taxon_suffix_manager.next_suffix(proposed_gtdb_sp)
                sp_epithet = '{}_{}'.format(canonical_taxon(specific_epithet(proposed_gtdb_sp)), next_suffix)
                proposed_gtdb_sp = '{} {}'.format(genus.replace('g__', 's__'), sp_epithet)

            assert genus != 'g__' and sp_epithet != ''

            if proposed_gtdb_sp != cur_genomes[rid].gtdb_taxa.species:
                if rid in prev_genomes:
                    num_updates += 1
                    self.update_log(rid, cur_genomes, prev_genomes, proposed_gtdb_sp, "Updated to reflect bionomial NCBI species name.")
                for cid in clusters[rid]:
                    #***cur_genomes[cid].gtdb_taxa.genus = 'g__' + generic_name(proposed_gtdb_sp)
                    cur_genomes[cid].gtdb_taxa.species = proposed_gtdb_sp
            
            if proposed_gtdb_sp in used_sp_names:
                dup_name_count += 1
                self.curation_log(rid, 
                                    cur_genomes, 
                                    prev_genomes, 
                                    proposed_gtdb_sp,
                                    True,
                                    "Assigned species name supported by NCBI binomial species name: {} {}".format(
                                        rid, used_sp_names[proposed_gtdb_sp]),
                                    False)
                self.logger.warning('NCBI binomial species name already in use: {}: {}, {}'.format(
                                    proposed_gtdb_sp, rid, used_sp_names[proposed_gtdb_sp]))
                                    
            if proposed_gtdb_sp in synonyms:
                self.logger.error('Assigning synonym to non-type strain species cluster: {}: {}'.format(ncbi_sp, rid))
                                    
            cluster_sp_names[rid] = proposed_gtdb_sp
            used_sp_names[proposed_gtdb_sp] = rid
            
        self.logger.info(f' ... updated name of {num_updates:,} GTDB clusters.')
        self.logger.warning(f' ... same name assigned to {dup_name_count:,} species clusters.')
        
    def name_placeholder_clusters(self,
                                    rids_by_naming_priority,
                                    cluster_sp_names,
                                    used_sp_names,
                                    clusters,
                                    prev_genomes, 
                                    cur_genomes, 
                                    synonyms):
        """Assign placeholder names to clusters that could not be assigned a binomial Latin name."""
        
        num_updates = 0
        dup_name_count = 0
        for rid in rids_by_naming_priority:
            if rid in cluster_sp_names:
                continue

            gtdb_sp = cur_genomes[rid].gtdb_taxa.species
            if gtdb_sp != 's__' and gtdb_sp not in used_sp_names:
                proposed_gtdb_sp = gtdb_sp
            else:
                genus = 'g__{unresolved}'
                gtdb_genus = cur_genomes[rid].gtdb_taxa.genus
                if gtdb_genus != 'g__':
                    genus = gtdb_genus
                    
                proposed_gtdb_sp = 's__{} sp{}'.format(genus.replace('g__', ''), rid[1:])
            
            if proposed_gtdb_sp != cur_genomes[rid].gtdb_taxa.species:
                if rid in prev_genomes:
                    num_updates += 1
                    self.update_log(rid, cur_genomes, prev_genomes, proposed_gtdb_sp, "Updated to reflect bionomial NCBI species name.")
                    
                for cid in clusters[rid]:
                    #***cur_genomes[cid].gtdb_taxa.genus = 'g__' + generic_name(proposed_gtdb_sp)
                    cur_genomes[cid].gtdb_taxa.species = proposed_gtdb_sp
            
            if proposed_gtdb_sp in used_sp_names:
                dup_name_count += 1
                self.curation_log(rid, 
                                    cur_genomes, 
                                    prev_genomes, 
                                    proposed_gtdb_sp,
                                    True,
                                    "Placeholder name assigned to multiple genomes: {} {}".format(
                                        rid, used_sp_names[proposed_gtdb_sp]),
                                    False)
                self.logger.warning('Placeholder name already in use: {}: {}, {}'.format(
                                    proposed_gtdb_sp, rid, used_sp_names[proposed_gtdb_sp]))
                                    
            if proposed_gtdb_sp in synonyms:
                self.logger.error('Assigning synonym to non-type strain species cluster: {}: {}'.format(ncbi_sp, rid))
                                    
            cluster_sp_names[rid] = proposed_gtdb_sp
            used_sp_names[proposed_gtdb_sp] = rid
            
        self.logger.info(f' ... updated name of {num_updates:,} GTDB clusters.')
        self.logger.warning(f' ... same name assigned to {dup_name_count:,} species clusters.')
        
    def  _sort_by_naming_priority(self, 
                                    clusters,
                                    prev_genomes, 
                                    cur_genomes, 
                                    gtdb_type_strain_ledger,
                                    synonyms):
        """Sort representatives by naming priority."""
        
        # group by naming priority
        type_species = []
        type_strains = []
        binomial = []
        placeholder = []
        for rid in clusters:
            ncbi_sp = cur_genomes[rid].ncbi_taxa.species
            if ncbi_sp in synonyms:
                ncbi_sp = synonyms[ncbi_sp]
                
            gtdb_sp = cur_genomes[rid].gtdb_taxa.species
            if gtdb_sp in synonyms:
                gtdb_sp = synonyms[gtdb_sp]
             
            if cur_genomes[rid].is_gtdb_type_species():
                type_species.append(rid)
            elif cur_genomes[rid].is_effective_type_strain():
                type_strains.append(rid)
            elif not is_placeholder_taxon(ncbi_sp):
                binomial.append(rid)
            else:
                placeholder.append(rid)
                
        assert len(clusters) == len(type_species) + len(type_strains) + len(binomial) + len(placeholder)
                
        # sort groups so previous GTDB representatives
        # are processed first
        rids_by_naming_priority = []
        for d in [type_species, type_strains, binomial, placeholder]:
            prev_reps = []
            new_reps = []
            for rid in d:
                if rid in prev_genomes.sp_clusters:
                    prev_reps.append(rid)
                else:
                    new_reps.append(rid)
                
            rids_by_naming_priority.extend(prev_reps)
            rids_by_naming_priority.extend(new_reps)
            
        assert len(rids_by_naming_priority) == len(clusters)
        
        return rids_by_naming_priority

    def update_species_names(self, 
                                clusters,
                                prev_genomes, 
                                cur_genomes, 
                                gtdb_type_strain_ledger,
                                synonyms):
        """Determine most suitable name for each GTDB species cluster."""
        
        # identify NCBI species with multiple genomes assembled from type strain of species
        self.logger.info('Determining type species genomes in each NCBI species.')
        ncbi_sp_type_strain_genomes = defaultdict(set)
        for gid in cur_genomes:
            if cur_genomes[gid].is_effective_type_strain():
                if gid in gtdb_type_strain_ledger:
                    ncbi_sp = gtdb_type_strain_ledger[gid]
                else:
                    ncbi_sp = cur_genomes[gid].ncbi_taxa.species
                    
                if ncbi_sp != 's__':
                    # yes, NCBI has genomes marked as assembled from type material
                    # that do not actually have a binomial species name
                    ncbi_sp_type_strain_genomes[ncbi_sp].add(gid)
        self.logger.info(' ... identified effective type strain genomes for {:,} NCBI species.'.format(
                            len(ncbi_sp_type_strain_genomes)))
        
        # verify that type genomes for a species are contained in a 
        # single GTDB species cluster
        rid_map = {}
        for rid, gids in clusters.items():
            rid_map[rid] = rid
            for gid in gids:
                rid_map[gid] = rid
                
        for ncbi_sp, type_gids in ncbi_sp_type_strain_genomes.items():
            gtdb_rids = set([rid_map[gid] for gid in type_gids])
            if len(gtdb_rids) > 1:
                self.logger.warning('Type strain genomes from NCBI species {} were assigned to {:,} GTDB species clusters: {}.'.format(
                                    ncbi_sp,
                                    len(gtdb_rids),
                                    [(rid, cur_genomes[rid].gtdb_taxa.species) for rid in gtdb_rids]))
                                    
        # order representatives by naming priority
        rids_by_naming_priority = self._sort_by_naming_priority(clusters,
                                                                prev_genomes, 
                                                                cur_genomes, 
                                                                gtdb_type_strain_ledger,
                                                                synonyms)
        
        # determine names for species clusters represented by a validally
        # published species name with a designated type species of genus
        self.logger.info('Assigning names to cluster represented by type species genomes.')
        cluster_sp_names = {}
        used_sp_names = {}
        self.name_type_species_clusters(rids_by_naming_priority,
                                        cluster_sp_names,
                                        used_sp_names,
                                        clusters,
                                        prev_genomes, 
                                        cur_genomes, 
                                        gtdb_type_strain_ledger,
                                        synonyms)

        num_type_species = len(cluster_sp_names)
        self.logger.info(f' ... assigned {num_type_species:,} name from type species.')
                
        # determine names for species clusters represented by a validally
        # or effectively published species name with a designated type strain
        self.logger.info('Assigning names to cluster represented by type strain genomes.')
        self.name_type_strain_clusters(rids_by_naming_priority,
                                        cluster_sp_names,
                                        used_sp_names,
                                        clusters,
                                        prev_genomes, 
                                        cur_genomes, 
                                        gtdb_type_strain_ledger,
                                        synonyms)
        num_type_strains = len(cluster_sp_names) - num_type_species
        total_sp_names = len(cluster_sp_names)
        self.logger.info(f' ... assigned {num_type_strains:,} names from type strains and {total_sp_names:,} total names.')

        # determine name for species clusters where a NCBI-define 
        # binomial latin name exists, but there is no designated type strain
        # (e.g., Candidatus species)
        self.logger.info('Assigning names to clusters with binomial NCBI species names.')
        self.name_binomial_clusters(rids_by_naming_priority,
                                    cluster_sp_names,
                                    used_sp_names,
                                    clusters,
                                    prev_genomes, 
                                    cur_genomes, 
                                    synonyms)
        num_assigned = len(cluster_sp_names) - total_sp_names
        total_sp_names = len(cluster_sp_names)
        self.logger.info(f' ... assigned {num_assigned:,} additional names and {total_sp_names:,} total names.')
            
        # determine names for remaining clusters which should all have a
        # placeholder name 
        self.logger.info('Assigning names to clusters without NCBI species names.')
        self.name_placeholder_clusters(rids_by_naming_priority,
                                        cluster_sp_names,
                                        used_sp_names,
                                        clusters,
                                        prev_genomes, 
                                        cur_genomes, 
                                        synonyms)
        num_assigned = len(cluster_sp_names) - total_sp_names
        total_sp_names = len(cluster_sp_names)
        self.logger.info(f' ... assigned {num_assigned:,} additional names and {total_sp_names:,} total names.')
        
        # sanity check
        for rid, sp in cluster_sp_names.items():
            assert cur_genomes[rid].gtdb_taxa.species == sp

    def _parse_explicit_taxa_updates(self, gtdb_taxa_updates_ledger):
        """Get explicit updates to previous GTDB taxon names."""
        
        explicit_taxon_updates = {}
        with open(gtdb_taxa_updates_ledger) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                
                prev_taxon = tokens[0].strip()
                assert '__' in prev_taxon
                    
                new_taxon = tokens[1].strip()
                assert '__' in new_taxon
                
                explicit_taxon_updates[prev_taxon] = new_taxon
                
        return explicit_taxon_updates
        
    def write_rep_info(self, 
                        cur_genomes,
                        clusters, 
                        cluster_sp_names, 
                        excluded_from_refseq_note,
                        ani_af,
                        output_file):
        """Write out information about selected representative genomes."""
                                            
        fout = open(output_file, 'w')
        fout.write('Species\tRepresentative\tNCBI assembly level\tNCBI genome category')
        fout.write('\tGenome size (bp)\tQuality score\tCompleteness (%)\tContamination (%)\tNo. scaffolds\tNo. contigs\tN50 contigs\tAmbiguous bases\tSSU count\tSSU length (bp)')
        fout.write('\tNo. genomes in cluster\tMean ANI\tMean AF\tMin ANI\tMin AF\tNCBI exclude from RefSeq\n')
        
        for gid in clusters:
            fout.write('{}\t{}\t{}\t{}'.format(
                        cluster_sp_names[gid], 
                        gid, 
                        cur_genomes[gid].ncbi_assembly_level,
                        cur_genomes[gid].ncbi_genome_category))

            fout.write('\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{}\t{:.1f}\t{}\t{}\t{}'.format(
                            cur_genomes[gid].length,
                            cur_genomes[gid].score_assembly(), 
                            cur_genomes[gid].comp,
                            cur_genomes[gid].cont,
                            cur_genomes[gid].scaffold_count,
                            cur_genomes[gid].contig_count,
                            cur_genomes[gid].contig_n50,
                            cur_genomes[gid].ambiguous_bases,
                            cur_genomes[gid].ssu_count,
                            cur_genomes[gid].ssu_length))
                            
            anis = []
            afs = []
            for cluster_id in clusters[gid]:
                ani, af = symmetric_ani(ani_af, gid, cluster_id)
                anis.append(ani)
                afs.append(af)
            
            if anis:
                fout.write('\t{}\t{:.1f}\t{:.2f}\t{:.1f}\t{:.2f}\t{}\n'.format(len(clusters[gid]),
                                                                    np_mean(anis), np_mean(afs),
                                                                    min(anis), min(afs),
                                                                    excluded_from_refseq_note.get(gid, '')))
            else:
                fout.write('\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(len(clusters[gid]),
                                                            'n/a', 'n/a', 'n/a', 'n/a',
                                                            excluded_from_refseq_note.get(gid, '')))
        fout.close()
        
    def parse_gtdb_type_strain_ledger(self, gtdb_type_strains_ledger, cur_genomes):
        """Read and validate GTDB type strain ledger."""
        
        gtdb_type_strain_ledger = {}
        with open(gtdb_type_strains_ledger) as f:
            header = f.readline().strip().split('\t')
            
            gid_index = header.index('Genome ID')
            sp_name_index = header.index('Proposed species name')
            
            for line in f:
                tokens = line.strip().split('\t')
                
                gid = canonical_gid(tokens[gid_index].strip())
                gtdb_sp_name = tokens[sp_name_index].strip()
                if not gtdb_sp_name.startswith('s__'):
                    gtdb_sp_name = 's__' + gtdb_sp_name
                gtdb_type_strain_ledger[gid] = gtdb_sp_name
                
                # validate assignment
                ncbi_sp = cur_genomes[gid].ncbi_taxa.species
                if ncbi_sp != 's__' and ncbi_sp != gtdb_sp_name:
                    self.logger.warning('GTDB type strain ledger disagrees with NCBI species assignment: {} {} {}'.format(
                                        gid, gtdb_sp_name, ncbi_sp))

        return gtdb_type_strain_ledger

    def run(self, 
                gtdb_clusters_file,
                prev_gtdb_metadata_file,
                prev_genomic_path_file,
                cur_gtdb_metadata_file,
                cur_genomic_path_file,
                uba_genome_paths,
                qc_passed_file,
                gtdbtk_classify_file,
                ncbi_genbank_assembly_file,
                untrustworthy_type_file,
                synonym_file,
                gtdb_type_strains_ledger,
                sp_priority_ledger,
                gtdb_taxa_updates_ledger,
                dsmz_bacnames_file):
        """Perform initial actions required for changed representatives."""

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
        
        # get path to previous and current genomic FASTA files
        self.logger.info('Reading path to previous and current genomic FASTA files.')
        prev_genomes.load_genomic_file_paths(prev_genomic_path_file)
        prev_genomes.load_genomic_file_paths(uba_genome_paths)
        cur_genomes.load_genomic_file_paths(cur_genomic_path_file)
        cur_genomes.load_genomic_file_paths(uba_genome_paths)
        
        # read named GTDB species clusters
        self.logger.info('Reading GTDB species clusters.')
        clusters, rep_radius = read_clusters(gtdb_clusters_file)
        self.logger.info(' ... identified {:,} clusters spanning {:,} genomes.'.format(
                            len(clusters),
                            sum([len(gids) + 1 for gids in clusters.values()])))
                            
        # write out archaeal and bacterial genome files
        self.logger.info('Creating genomic files for archaeal and bacterial reference genomes.')
        for prefix, domain in [('ar', 'd__Archaea'), ('bac', 'd__Bacteria')]:
            fout = open(os.path.join(self.output_dir, '{}_genomes.lst'.format(prefix)), 'w')
            genome_count = 0
            for rid in clusters:
                gtdb_domain = cur_genomes[rid].gtdb_taxa.domain
                if gtdb_domain == domain:
                    fout.write('{}\n'.format(cur_genomes[rid].ncbi_accn))
                    genome_count += 1
                        
            self.logger.info(f' ... identified {genome_count:,} {domain} reference genomes.')
            fout.close()
            
        # read ledger with explicit names for type strain genomes
        self.logger.info('Parsing GTDB type strain ledger.')
        gtdb_type_strain_ledger = self.parse_gtdb_type_strain_ledger(gtdb_type_strains_ledger, cur_genomes)
        self.logger.info(' ... identified {:,} genomes in ledger.'.format(len(gtdb_type_strain_ledger)))

        # get list of synonyms in order to restrict usage of species names
        self.logger.info('Reading GTDB synonyms.')
        synonyms = parse_synonyms(synonym_file)
        self.logger.info(' ... identified {:,} synonyms from {:,} distinct species.'.format(
                            len(synonyms),
                            len(set(synonyms.values()))))
        
        # set current genomes to have same GTDB assignments as in previous
        # GTDB release. This is necessary, since genomes may have different
        # NCBI accession numbers between releases and thus the previous GTDB
        # taxonomy will not be reflected in the latest GTDB database. The 
        # exception is if a genome has changed domains, in which case the
        # previous assignment is invalid.
        self.logger.info('Setting GTDB taxonomy of genomes in current genome set.')
        update_count = 0
        conflicting_domain_count = 0
        for prev_gid in prev_genomes:
            if prev_gid in cur_genomes:
                if prev_genomes[prev_gid].gtdb_taxa != cur_genomes[prev_gid].gtdb_taxa:
                    if prev_genomes[prev_gid].gtdb_taxa.domain == cur_genomes[prev_gid].gtdb_taxa.domain:
                        update_count += 1
                        cur_genomes[prev_gid].gtdb_taxa.update_taxa(prev_genomes[prev_gid].gtdb_taxa)
                    else:
                        conflicting_domain_count += 1
        self.logger.info(f' ... updated {update_count:,} genomes.')
        self.logger.info(f' ... identified {conflicting_domain_count:,} genomes with conflicting domain assignments.')

       # get explicit updates to previous GTDB taxa
        self.logger.info('Reading explicit taxa updates.')
        explicit_taxon_updates = self._parse_explicit_taxa_updates(gtdb_taxa_updates_ledger)
        self.logger.info(f' ... identified {len(explicit_taxon_updates):,} updates.')
        
        self.logger.info('Updating current genomes to reflect explicit taxa updates.')
        update_count = 0
        for cur_taxon, new_taxon in explicit_taxon_updates.items():
            rank_prefix = cur_taxon[0:3]
            rank_index = Taxonomy.rank_prefixes.index(rank_prefix)

            for gid in cur_genomes:
                if cur_genomes[gid].gtdb_taxa.get_taxa(rank_index) == cur_taxon:
                    update_count += 1
                    cur_genomes[gid].gtdb_taxa.set_taxa(rank_index, new_taxon)
                    
                    if rank_prefix == 'g__':
                        # should also update the species name
                        new_sp = cur_genomes[gid].gtdb_taxa.species.replace(cur_taxon[3:], new_taxon[3:])
                        cur_genomes[gid].gtdb_taxa.set_taxa(rank_index+1, new_sp)
                
        self.logger.info(f' ... updated {update_count:,} genomes.')

        # create species name manager
        self.logger.info('Initializing species name manager.')
        self.sp_name_mngr = SpeciesNameManager(prev_genomes, 
                                                    cur_genomes,
                                                    self.fastani)
                                                    
        # initialize species priority manager
        self.sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger,
                                                        dsmz_bacnames_file)
 
        # establish appropriate species names for GTDB clusters with new representatives
        self.update_species_names(clusters,
                                    prev_genomes, 
                                    cur_genomes, 
                                    gtdb_type_strain_ledger,
                                    synonyms)
                                    
        # write out taxonomy files
        bac_taxonomy_out = open(os.path.join(self.output_dir, 'gtdb_bac_taxonomy.tsv'), 'w')
        ar_taxonomy_out = open(os.path.join(self.output_dir, 'gtdb_ar_taxonomy.tsv'), 'w')
        for rid in clusters:
            gtdb_domain = cur_genomes[rid].gtdb_taxa.domain
            fout = bac_taxonomy_out
            if gtdb_domain == 'd__Archaea':
                fout = ar_taxonomy_out
            fout.write('{}\t{}\n'.format(rid, cur_genomes[rid].gtdb_taxa))
            
        bac_taxonomy_out.close()
        ar_taxonomy_out.close()

        # write out cluster information with finalized GTDB cluster names
        if False: #***
            self.logger.info('Writing {:,} species clusters to file.'.format(len(clusters)))
            self.logger.info('Writing {:,} cluster radius information to file.'.format(len(rep_radius)))
            
            write_clusters(clusters, 
                            rep_radius, 
                            cur_genomes,
                            os.path.join(self.output_dir, 'gtdb_clusters_de_novo.tsv'))

            write_rep_radius(rep_radius, 
                                cur_genomes,
                                os.path.join(self.output_dir, 'gtdb_ani_radius_de_novo.tsv'))