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
                                                gtdb_merged_genera,
                                                sort_by_naming_priority,
                                                longest_common_prefix,
                                                is_placeholder_taxon,
                                                is_placeholder_sp_epithet)


class UpdateSpeciesFinal(object):
    """Finalize species names based on results of manual curation."""

    def __init__(self, ani_cache_file, cpus, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
        self.fastani = FastANI(ani_cache_file, cpus)
        
        self.sp_name_log = open(os.path.join(self.output_dir, 'sp_name_log.tsv'), 'w')
        self.sp_name_log.write('GTDB domain\tGenome ID\tPrevious GTDB species\tNew GTDB species\tStatus\tReason for update\n')

        self.sp_curation_log = open(os.path.join(self.output_dir, 'sp_curation_log.tsv'), 'w')
        self.sp_curation_log.write('GTDB domain\tGenome ID\tPrevious NCBI species\tCurrent NCBI species')
        self.sp_curation_log.write('\tPrevious GTDB genus\tPrevious GTDB species\tProposed GTDB species')
        self.sp_curation_log.write('\tCase\tIssue\n')

    def update_log(self, 
                    gid, 
                    cur_genomes, 
                    prev_genomes, 
                    sp,
                    status,
                    reason):
        """Add entry to update log."""

        prev_gtdb_sp = 'n/a'
        if gid in prev_genomes:
            prev_gtdb_sp = prev_genomes[gid].gtdb_taxa.species

        self.sp_name_log.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            cur_genomes[gid].gtdb_taxa.domain,
            gid,
            prev_gtdb_sp,
            sp,
            status,
            reason))
            
    def curation_log(self, 
                        gid, 
                        cur_genomes, 
                        prev_genomes, 
                        proposed_gtdb_sp, 
                        case,
                        issue):
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

        self.sp_curation_log.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                    cur_genomes[gid].gtdb_taxa.domain,
                                    gid,
                                    prev_ncbi_sp,
                                    cur_ncbi_sp,
                                    prev_gtdb_genus,
                                    prev_gtdb_sp,
                                    proposed_gtdb_sp,
                                    case,
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
        
    def parse_species_rep_status(self, updated_species_reps_file):
        """Parse file indicating status of previous GTDB species representatives."""
        
        prev_gtdb_rep_status = {}
        new_to_prev_rid = {}
        with open(updated_species_reps_file) as f:
            header = f.readline().strip().split('\t')
            
            prev_rid_index = header.index('Previous representative ID')
            new_rid_index = header.index('New representative ID')
            rep_status_index = header.index('Representative status')
            
            for line in f:
                tokens = line.strip().split('\t')
                
                prev_rid = tokens[prev_rid_index]
                new_rid = tokens[new_rid_index]
                if new_rid.lower() != 'none':
                    if new_rid in prev_gtdb_rep_status:
                        self.logger.error(f'New representative indicated multiple times: {new_rid}')
                        sys.exit(-1)
                        
                    prev_gtdb_rep_status[new_rid] = tokens[rep_status_index]
                    new_to_prev_rid[new_rid] = prev_rid
                else:
                    prev_gtdb_rep_status[prev_rid] = tokens[rep_status_index]
        
        return prev_gtdb_rep_status, new_to_prev_rid
        
    def sanity_check_sp(self,
                            gid,
                            gtdb_taxa_final,
                            prev_genomes,
                            cur_genomes):
        """Sanity check species assignment."""
        
        genus = gtdb_taxa_final[gid][Taxonomy.GENUS_INDEX]
        species = gtdb_taxa_final[gid][Taxonomy.SPECIES_INDEX]
        
        # check for binomial species name
        if len(species.split()) != 2:
            self.curation_log(gid, 
                                cur_genomes, 
                                prev_genomes, 
                                species, 
                                "INVALID_SPECIES_NAME",
                                "GTDB species name is not binomial")
            return

        generic = generic_name(species)
        specific = specific_epithet(species)
        
        ncbi_genus = cur_genomes[gid].ncbi_taxa.genus
        ncbi_sp = cur_genomes[gid].ncbi_taxa.species
        ncbi_specific = specific_epithet(ncbi_sp)
        
        # check generic name is the same as the genus name
        if genus[3:] != generic:
            self.curation_log(gid, 
                                cur_genomes, 
                                prev_genomes, 
                                species, 
                                "GENERIC_NAME_CONFLICT",
                                "Generic name of GTDB species does not reflect GTDB genus assignment")
        
        # genus name must reflect type species of genus
        if cur_genomes[gid].is_gtdb_type_species:
            if ncbi_genus != genus:
                self.curation_log(gid, 
                                    cur_genomes, 
                                    prev_genomes, 
                                    species, 
                                    "TYPE_SPECIES_OF_GENUS",
                                    "GTDB genus assignment does not reflect type species of genus")
                                    
        # specific name must reflect type strain of species
        if cur_genomes[gid].is_gtdb_type_species:
            if ncbi_specific != specific:
                self.curation_log(gid, 
                                    cur_genomes, 
                                    prev_genomes, 
                                    species, 
                                    "TYPE_STRAIN_OF_SPECIES",
                                    "GTDB specific species name does not reflect type strain of species")
                                    
    def key_taxon(self, gid, taxonomy):
        """Get key taxon for determine species assignments."""
        
        genus = taxonomy[gid][Taxonomy.GENUS_INDEX]
        species = taxonomy[gid][Taxonomy.SPECIES_INDEX]
        generic = generic_name(species)
        specific = specific_epithet(species)

        return genus, species, generic, specific

    def set_type_species(self, 
                            gid, 
                            gtdb_taxa_final, 
                            cur_genomes, 
                            prev_merged_gtdb_genera, 
                            case_count):
        """Establish final species name for type species of genus genome."""
        
        note = ''
        
        genus, species, generic, specific = self.key_taxon(gid, gtdb_taxa_final)
        final_sp = 's__{} {}'.format(genus[3:], specific)
        
        ncbi_sp = cur_genomes[gid].ncbi_taxa.species
        
        if final_sp == ncbi_sp:
            return final_sp, note
        else:
            if generic_name(final_sp) != generic_name(ncbi_sp):
                case_count['UNRESOLVED_TYPE_SPECIES_GENERIC'] += 1
            elif specific_epithet(final_sp) != specific_epithet(ncbi_sp):
                case_count['UNRESOLVED_TYPE_SPECIES_SPECIFIC'] += 1
            else:
                self.logger.error('Unexpected error with proposed GTDB type species name: {} {}'.format(
                                    final_sp, ncbi_sp))
                sys.exit(-1)
                
            note = 'unresolved error'

        return final_sp, note
        
    def set_type_strain(self, gid, gtdb_taxa_final, cur_genomes, synonyms, case_count):
        """Establish final species name for type strain of species genome."""
        
        note = ''
        
        genus, species, generic, specific = self.key_taxon(gid, gtdb_taxa_final)
        final_sp = 's__{} {}'.format(genus[3:], specific)
        
        if specific == specific_epithet(cur_genomes[gid].ncbi_taxa.species):
            return final_sp, note
        else:
            case_count['UNRESOLVED_TYPE_STRAIN'] += 1

        return final_sp, note

    def finalize_species_names(self,
                                gtdb_taxa_final, 
                                updated_gtdb_genera, 
                                updated_gtdb_species,
                                clusters,
                                prev_genomes, 
                                cur_genomes, 
                                gtdb_type_strain_ledger,
                                prev_gtdb_rep_status,
                                new_to_prev_rid,
                                synonyms):
        """Establish final species names."""
        
        prev_merged_gtdb_genera = gtdb_merged_genera(prev_genomes, 
                                                        self.sp_priority_mngr, 
                                                        self.output_dir)
        
        # get all previous GTDB species
        prev_gtdb_sp = set()
        for gid in prev_genomes:
            prev_gtdb_sp.add(prev_genomes[gid].gtdb_taxa.species)
            
        # order representatives by naming priority
        rids_by_naming_priority = sort_by_naming_priority(clusters,
                                                            prev_genomes, 
                                                            cur_genomes, 
                                                            gtdb_type_strain_ledger)
        
        # establish final name for each GTDB species cluster
        case_count = defaultdict(int)
        for rid in rids_by_naming_priority:
            if rid in updated_gtdb_species:
                case = 'MANUAL_CURATION'
                note = 'Species name set by manual curation'
                final_sp = gtdb_taxa_final[rid][Taxonomy.SPECIES_INDEX]
            elif cur_genomes[rid].is_gtdb_type_species():
                case = 'TYPE_SPECIES_OF_GENUS'
                final_sp, note = self.set_type_species(rid,
                                                        gtdb_taxa_final, 
                                                        cur_genomes, 
                                                        prev_merged_gtdb_genera,
                                                        case_count)
            elif cur_genomes[rid].is_gtdb_type_strain():
                case = 'TYPE_STRAIN_OF_SPECIES'
                final_sp, note = self.set_type_strain(rid,
                                                        gtdb_taxa_final, 
                                                        cur_genomes, 
                                                        synonyms,
                                                        case_count)

            else:
                pass
                
            # status string should indicate:
            #  status of GTDB species cluster: EXISTING or NEW
            #  status of representative genome: UNCHANGED, REPLACED, or NEW
            #  status of GTDB species name: UNCHANGED, CHANGED, or NEW
            status = 'NEW_SPECIES_CLUSTER:NEW_REPRESENTATIVE:NEW_SPECIES_NAME'
            if rid in prev_gtdb_rep_status:
                status = 'EXISTING_SPECIES_CLUSTER'
                status += ':{}_REPRESENTATIVE'.format(prev_gtdb_rep_status[rid])
                if prev_gtdb_rep_status[rid] not in ['REPLACED', 'UNCHANGED', 'LOST']:
                    self.logger.error(f'Unexpected status for previous GTDB representative: {rid}, {prev_gtdb_rep_status[rid]}')
                    sys.exit(-1)

                prev_sp_rid = new_to_prev_rid.get(rid, rid)
                if final_sp == prev_genomes[prev_sp_rid].gtdb_taxa.species:
                    status += ':UNCHANGED_SPECIES_NAME'
                else:
                    status += ':CHANGED_SPECIES_NAME'
            else:
                # new GTDB species cluster so sanity check species name is new
                if rid in prev_genomes.sp_clusters:
                    self.logger.error(f'New species cluster assigned previous representative genomes: {rid}')
                                    
                if final_sp in prev_gtdb_sp:
                    self.logger.error(f'New species cluster assigned previous GTDB species name: {rid}, {final_sp}')
                
            gtdb_taxa_final[rid][Taxonomy.SPECIES_INDEX] = final_sp
            
            self.sanity_check_sp(rid, 
                                    gtdb_taxa_final,
                                    prev_genomes,
                                    cur_genomes)

        for case, count in case_count.items():
            print('{}\t{}'.format(case, count))
        
    def record_modified_gtdb_sp_names(self,
                                        gtdb_taxa_final, 
                                        updated_gtdb_genera, 
                                        updated_gtdb_species,
                                        clusters,
                                        prev_genomes, 
                                        cur_genomes,
                                        gtdb_type_strain_ledger,
                                        prev_gtdb_rep_status,
                                        new_to_prev_rid,
                                        synonyms):
        """Create record of modified GTDB species names."""

        for rid, status in prev_gtdb_rep_status.items():
            if status == 'LOST' and rid not in new_to_prev_rid:
                # species cluster is effectively retired
                self.update_log(rid, 
                                    cur_genomes, 
                                    prev_genomes, 
                                    'n/a',
                                    'RETIRED',
                                    'Species cluster retired as only representative of species is no longer available.')
            else:
                assert rid in gtdb_taxa_final
                
                prev_sp_rid = new_to_prev_rid[rid]
                prev_gtdb_sp = prev_genomes[prev_sp_rid].gtdb_taxa.species
                sp_final = gtdb_taxa_final[rid][Taxonomy.SPECIES_INDEX]
            
                if prev_gtdb_sp == sp_final:
                    self.update_log(rid, 
                                    cur_genomes, 
                                    prev_genomes, 
                                    sp_final,
                                    'UNCHANGED',
                                    '')
                else:
                    self.update_log(rid, 
                                    cur_genomes, 
                                    prev_genomes, 
                                    sp_final,
                                    'CHANGED',
                                    '')
        
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
        
    def parse_manual_curation(self,
                                bac_curation_file,
                                ar_curation_file,
                                bac_init_taxonomy,
                                ar_init_taxonomy):
        """Parse manually curated taxonomy and identify updated genera and species."""
        
        bac_cur_taxonomy = Taxonomy().read(bac_curation_file)
        ar_cur_taxonomy = Taxonomy().read(ar_curation_file)
        gtdb_taxa_final = bac_cur_taxonomy
        gtdb_taxa_final.update(ar_cur_taxonomy)
        
        bac_init_taxonomy = Taxonomy().read(bac_init_taxonomy)
        ar_init_taxonomy = Taxonomy().read(ar_init_taxonomy)
        init_gtdb_taxonomy = bac_init_taxonomy
        init_gtdb_taxonomy.update(ar_init_taxonomy)
        
        assert set(bac_cur_taxonomy) == set(bac_init_taxonomy)
        assert set(ar_cur_taxonomy) == set(ar_init_taxonomy)
        
        updated_gtdb_genera = set()
        updated_gtdb_species = set()
        for gid in gtdb_taxa_final:
            cur_genus = gtdb_taxa_final[gid][Taxonomy.GENUS_INDEX]
            init_genus = init_gtdb_taxonomy[gid][Taxonomy.GENUS_INDEX]
            if cur_genus != init_genus:
                updated_gtdb_genera.add(gid)
                
            cur_sp = gtdb_taxa_final[gid][Taxonomy.SPECIES_INDEX]
            init_sp = init_gtdb_taxonomy[gid][Taxonomy.SPECIES_INDEX]
            if cur_sp != init_sp:
                updated_gtdb_species.add(gid)
        
        return gtdb_taxa_final, updated_gtdb_genera, updated_gtdb_species

    def run(self,
                bac_curation_file,
                ar_curation_file,
                bac_init_taxonomy,
                ar_init_taxonomy,
                gtdb_clusters_file,
                prev_gtdb_metadata_file,
                cur_gtdb_metadata_file,
                uba_genome_paths,
                qc_passed_file,
                gtdbtk_classify_file,
                updated_species_reps_file,
                ncbi_genbank_assembly_file,
                untrustworthy_type_file,
                synonym_file,
                gtdb_type_strains_ledger,
                sp_priority_ledger,
                gtdb_taxa_updates_ledger,
                dsmz_bacnames_file):
        """Finalize species names based on results of manual curation."""
        
        # identify species and genus names updated during manual curation
        self.logger.info('Parsing manually curated taxonomy.')
        rtn = self.parse_manual_curation(bac_curation_file,
                                                ar_curation_file,
                                                bac_init_taxonomy,
                                                ar_init_taxonomy)
        gtdb_taxa_final, updated_gtdb_genera, updated_gtdb_species = rtn
        self.logger.info(' ... identified taxonomy strings for {:,} genomes.'.format(
                            len(gtdb_taxa_final)))
        self.logger.info(' ... identified {:,} genera and {:,} species names updated by manual curation.'.format(
                            len(updated_gtdb_genera), len(updated_gtdb_species)))

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
                                                untrustworthy_type_ledger=untrustworthy_type_file)
        self.logger.info(f' ... current genome set contains {len(cur_genomes):,} genomes.')
        
        cur_genomes.set_prev_gtdb_classifications(prev_genomes)
        
        # update current genomes with GTDB-Tk classifications
        self.logger.info('Updating current genomes with GTDB-Tk classifications.')
        num_updated, num_ncbi_sp = cur_genomes.set_gtdbtk_classification(gtdbtk_classify_file, 
                                                                            prev_genomes)
        self.logger.info(f' ... set GTDB taxa for {num_updated:,} genomes with {num_ncbi_sp:,} genomes using NCBI genus and species name.')
        
        # read ledger with explicit names for type strain genomes
        self.logger.info('Validating species names given in GTDB type strain ledger.')
        gtdb_type_strain_ledger = self.parse_gtdb_type_strain_ledger(gtdb_type_strains_ledger, cur_genomes)
        self.logger.info(' ... identified {:,} genomes in ledger.'.format(len(gtdb_type_strain_ledger)))
        num_conflict = 0
        for gid in gtdb_type_strain_ledger:
            if gid not in gtdb_taxa_final:
                continue
                
            cur_gtdb_sp = canonical_taxon(gtdb_taxa_final[gid][Taxonomy.SPECIES_INDEX])
            ledger_sp = gtdb_type_strain_ledger[gid]
            if cur_gtdb_sp != ledger_sp:
                num_conflict += 1
                self.logger.error('Assigned species name for {} conflicts with GTDB type strain ledger: {} {}'.format(
                                    gid, cur_gtdb_sp, ledger_sp))
                
        if num_conflict > 0:
            self.exit(-1)
        
        # read named GTDB species clusters
        self.logger.info('Reading GTDB species clusters.')
        clusters, rep_radius = read_clusters(gtdb_clusters_file)
        self.logger.info(' ... identified {:,} clusters spanning {:,} genomes.'.format(
                            len(clusters),
                            sum([len(gids) + 1 for gids in clusters.values()])))
        assert set(clusters) == set(gtdb_taxa_final)

        # get list of synonyms in order to restrict usage of species names
        self.logger.info('Reading GTDB synonyms.')
        synonyms = parse_synonyms(synonym_file)
        self.logger.info(' ... identified {:,} synonyms from {:,} distinct species.'.format(
                            len(synonyms),
                            len(set(synonyms.values()))))

        # create species name manager
        self.logger.info('Initializing species name manager.')
        self.sp_name_mngr = SpeciesNameManager(prev_genomes, 
                                                    cur_genomes,
                                                    self.fastani)
                                                    
        # initialize species priority manager
        self.sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger,
                                                        dsmz_bacnames_file)
                                                        
        # parse file indicating status of previous GTDB species representatives
        self.logger.info('Reading status of previous GTDB species representatives.')
        prev_gtdb_rep_status, new_to_prev_rid = self.parse_species_rep_status(updated_species_reps_file)
        self.logger.info(' ... read status for {:,} representatives.'.format(
                            len(prev_gtdb_rep_status)))
        assert len(prev_gtdb_rep_status) == len(prev_genomes.sp_clusters)
 
        # establish appropriate species names for GTDB clusters with new representatives
        self.finalize_species_names(gtdb_taxa_final, 
                                        updated_gtdb_genera, 
                                        updated_gtdb_species,
                                        clusters,
                                        prev_genomes, 
                                        cur_genomes,
                                        gtdb_type_strain_ledger,
                                        prev_gtdb_rep_status,
                                        new_to_prev_rid,
                                        synonyms)
                                        
        # establish reason for any modified GTDB species names
        self.record_modified_gtdb_sp_names(gtdb_taxa_final, 
                                            updated_gtdb_genera, 
                                            updated_gtdb_species,
                                            clusters,
                                            prev_genomes, 
                                            cur_genomes,
                                            gtdb_type_strain_ledger,
                                            prev_gtdb_rep_status,
                                            new_to_prev_rid,
                                            synonyms)
                                            
        self.sp_name_log.close()
        self.sp_curation_log.close()
                                    
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