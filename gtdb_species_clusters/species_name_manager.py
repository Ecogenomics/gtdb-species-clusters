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

from gtdb_species_clusters.fastani import FastANI
from gtdb_species_clusters.taxon_utils import specific_epithet, generic_name, canonical_taxon
from gtdb_species_clusters.taxon_suffix_manager import TaxonSuffixManager


class SpeciesNameManager(object):
    """Manage assignment and updating of GTDB species names.
    
    This manager handle the following:
     - tracking valid suffixes
     - determining presence of species name in GTDB
     - determining presence of species epithet in a given genus
     - tracking if a species name has been assigned multiple times
    """

    def __init__(self, prev_genomes, cur_genomes, fastani):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')
        self.fastani = fastani
        
        self.taxon_suffix_manager = TaxonSuffixManager()
        
        self.prev_genomes = prev_genomes
        self.cur_genomes = cur_genomes
        
        # get species assignment of previous GTDB species clusters
        self.gtdb_sp_epithets = defaultdict(set)
        self.gtdb_canonical_sp_epithets = defaultdict(set)
        self.sp_epithets_rid = defaultdict(lambda: {})
        for rid, sp in prev_genomes.sp_clusters.species():
            gtdb_genus = prev_genomes[rid].gtdb_taxa.genus
            gtdb_sp_epithet = prev_genomes[rid].gtdb_taxa.specific_epithet
            self.gtdb_sp_epithets[gtdb_genus].add(gtdb_sp_epithet)
            self.gtdb_canonical_sp_epithets[gtdb_genus].add(canonical_taxon(gtdb_sp_epithet))
            self.sp_epithets_rid[gtdb_genus][gtdb_sp_epithet] = rid

    def add_suffix(self, generic, specific):
        """Determine suffix for species."""
        
        sp = 's__{} {}'.format(generic, specific)
        next_suffix = self.taxon_suffix_manager.next_suffix(sp)
        
        canonical_specific = canonical_taxon(specific)
        
        return '{}_{}'.format(canonical_specific, next_suffix)
        
    def _nontype_sp_epithet(self, new_cur_ncbi_sp_epithet, orig_prev_gtdb_genus):
        """Get appropriate and potentially suffixed species epithet for non-type genome."""
        
        if new_cur_ncbi_sp_epithet not in self.gtdb_canonical_sp_epithets[orig_prev_gtdb_genus]:
            # only species in genus with this epithet so it can remain unsuffixed
            epithet = new_cur_ncbi_sp_epithet
        else:
            # other species already have suffixed versions of this epithet so 
            # must also suffix this epithet since the true location of this 
            # species is unclear
            epithet = self.add_suffix(orig_prev_gtdb_genus[3:], new_cur_ncbi_sp_epithet)
        
        return epithet
        
    def update(self,
                orig_rid,
                orig_prev_gtdb_sp,
                orig_prev_ncbi_sp,
                orig_cur_ncbi_sp,
                new_rid,
                new_cur_ncbi_sp):
        """Determine if species name can be updated."""

        new_is_type_strain = (self.cur_genomes[new_rid].is_gtdb_type_strain() 
                                or self.cur_genomes[new_rid].is_ncbi_type_strain())
        
        orig_prev_gtdb_genus = 'g__{}'.format(generic_name(orig_prev_gtdb_sp))
        orig_prev_gtdb_sp_epithet = specific_epithet(orig_prev_gtdb_sp)
        orig_cur_ncbi_sp_epithet = specific_epithet(orig_cur_ncbi_sp)
        new_cur_ncbi_sp_epithet = specific_epithet(new_cur_ncbi_sp)
        
        if new_cur_ncbi_sp_epithet == orig_prev_gtdb_sp_epithet:
            # given that the species epithet has not changed, there is no
            # need to update the GTDB species name
            action = 'UNCHANGED'
            epithet = orig_prev_gtdb_sp_epithet 
        elif new_is_type_strain:
            # new representative is a type strain with a new species epithet so the 
            # GTDB species name needs to be updated
            action = 'UPDATED'
            if new_cur_ncbi_sp_epithet in self.sp_epithets_rid[orig_prev_gtdb_genus]:
                # the new species name already exists in the GTDB
                conflicting_rid = self.sp_epithets_rid[orig_prev_gtdb_genus][new_cur_ncbi_sp_epithet]
                if self.cur_genomes[conflicting_rid].is_gtdb_type_strain():
                    # this is an issue since we have two GTDB clusters claiming
                    # to be represented by a type strain genome
                    action = 'CONFLICT'
                    epithet = new_cur_ncbi_sp_epithet + '-CONFLICT'
                else:
                    # need to update the GTDB species cluster that currently
                    # has this species name since it is not based on the type strain
                    conflicting_epithet = self._nontype_sp_epithet(new_cur_ncbi_sp_epithet, 
                                                                    orig_prev_gtdb_genus)
            else:
                # new species name doesn't exist in GTDB so can update name of GTDB cluster
                epithet = new_cur_ncbi_sp_epithet
        else:
            if new_cur_ncbi_sp_epithet == canonical_taxon(orig_prev_gtdb_sp_epithet):
                # since the new NCBI species epithet matches the canonical name of
                # the current GTDB species the current GTDB species assignment should 
                # be retained
                action = 'UNCHANGED'
                epithet = orig_prev_gtdb_sp_epithet
            elif new_cur_ncbi_sp_epithet == orig_prev_gtdb_sp_epithet:
                # NCBI species epithet associated with GTDB species cluster
                # has not changed, so there should be no need to update
                # the name of the cluster
                pass
            elif new_cur_ncbi_sp_epithet == orig_cur_ncbi_sp_epithet:
                # the previous and new representative have the same NCBI 
                # epithet so the GTDB species cluster should reflect this
                action = 'UPDATED'
                epithet = self._nontype_sp_epithet(new_cur_ncbi_sp_epithet, 
                                                    orig_prev_gtdb_genus)
            
        
        
        
        
        
        
        
        
        
        
        if False:
            conflicting_rid = conflicting_epithet = None
            if new_ncbi_sp == 's__':
                # given a lack of new species information, the current GTDB species
                # name should be retained
                action = 'UNCHANGED'
                epithet = prev_gtdb_sp_epithet
            elif new_ncbi_sp_epithet == prev_gtdb_sp_epithet:
                # given that the species epithet has not changed, there is no
                # need to update the GTDB species name
                action = 'UNCHANGED'
                epithet = prev_gtdb_sp_epithet
            elif (not new_is_type_strain 
                    and new_ncbi_sp_epithet == canonical_taxon(prev_gtdb_sp_epithet)):
                # since the NCBI species epithet matches the canonical name of
                # the current GTDB species epithet and the new genome is not a 
                # type strain the current GTDB species assignment should be retained
                action = 'UNCHANGED'
                epithet = prev_gtdb_sp_epithet
            elif new_is_type_strain:
                # new genomes is a type strain with a new species epithet so the 
                # GTDB species name needs to be updated
                if new_ncbi_sp_epithet in self.sp_epithets_rid[prev_gtdb_genus]:
                    # the new species name already exists in the GTDB
                    conflicting_rid = self.sp_epithets_rid[prev_gtdb_genus][new_ncbi_sp_epithet]
                    if self.cur_genomes[conflicting_rid].is_gtdb_type_strain():
                        # this is an issue since we have two GTDB clusters claiming
                        # to be represented by a type strain genome
                        action = 'CONFLICT'
                        epithet = new_ncbi_sp_epithet + '-CONFLICT'
                    else:
                        # need to update the GTDB species cluster that currently
                        # has this species name since it is not based on the type strain
                        conflicting_epithet = self._nontype_sp_epithet(new_ncbi_sp_epithet, prev_gtdb_genus)
                else:
                    # new species name doesn't exist in GTDB so can update name of GTDB cluster
                    action = 'UPDATED'
                    epithet = new_ncbi_sp_epithet
            elif cur_ncbi_sp_epithet == new_ncbi_sp_epithet:
                # the previous and new representative have the same NCBI 
                # epithet so the GTDB species cluster should reflect this
                action = 'UPDATED'
                epithet = self._nontype_sp_epithet(new_ncbi_sp_epithet, 
                                                            prev_gtdb_genus)
            else:
                # the previous and new representatives do not agree on the
                # species epithet so must decide which epithet to use
            
            
            
                if new_ncbi_sp_epithet not in self.gtdb_sp_epithets[prev_gtdb_genus]:
                    # the proposed NCBI species epithet is not in the GTDB genus so
                    # the GTDB species name can be update accordingly
                    action = 'UPDATED'
                    epithet = self._nontype_sp_epithet(new_ncbi_sp_epithet, 
                                                            prev_gtdb_genus)
                else:
                    # since the NCBI species epithet already exists in the genus
                    # the mostly likely situation is that this new genome is 
                    # misclassified so retain the current GTDB species name
                    action = 'UNCHANGED'
                    epithet = prev_gtdb_sp_epithet

        if False:
            if new_ncbi_sp_epithet not in self.gtdb_sp_epithets[prev_gtdb_genus]:
                # the proposed NCBI species epithet is not in the GTDB genus so
                # the GTDB species name can be update accordingly
                action = 'UPDATED'
                
                if new_is_type_strain:
                    # genome is type strain so species should have unsuffixed epithet
                    epithet = new_ncbi_sp_epithet
                else:
                    # the correct epithet depends on if there are other species with
                    # suffixed versions of this epithet in the genus
                    epithet = self._nontype_sp_epithet(new_ncbi_sp_epithet, 
                                                        prev_gtdb_genus)
            else:
                # check if species cluster should be merged
                ncbi_proposed_rid = self.sp_epithets_rid[prev_gtdb_genus][new_ncbi_sp_epithet]
                ani, af = self.fastani.symmetric_ani_cached(new_rid, ncbi_proposed_rid,
                                                            self.cur_genomes[new_rid].genomic_file, 
                                                            self.cur_genomes[ncbi_proposed_rid].genomic_file)

                if ani >= 95 and af >= 0.65:
                    # since new representative is close to another GTDB representative with the
                    # same species name they should be merged
                    action = 'MERGED'
                    epithet = new_ncbi_sp_epithet
                else:
                    if not new_is_type_strain:
                        # since the NCBI species epithet already exists in the genus
                        # the mostly likely situation is that this new genome is 
                        # misclassified so retain the current GTDB species name
                        action = 'UNCHANGED'
                        epithet = prev_gtdb_sp_epithet
                    else:
                        # this is a direct conflict since the NCBI species name is
                        # from a type strain genome and thus must be reflected in
                        # the species name
                        action = 'CONFLICT'
                        epithet = new_ncbi_sp_epithet + '-CONFLICT'
                    
                    
        return action, 's__{} {}'.format(prev_gtdb_genus[3:], epithet), conflicting_rid, conflicting_epithet
