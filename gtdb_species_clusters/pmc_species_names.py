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
import copy
import argparse
import logging
import pickle
from collections import defaultdict, Counter

from numpy import (mean as np_mean, std as np_std)

from biolib.taxonomy import Taxonomy

from gtdb_species_clusters.genomes import Genomes

from gtdb_species_clusters.species_name_manager import SpeciesNameManager
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters.specific_epithet_manager import SpecificEpithetManager
from gtdb_species_clusters.ncbi_species_manager import NCBI_SpeciesManager
from gtdb_species_clusters.genome_utils import canonical_gid
from gtdb_species_clusters.type_genome_utils import (read_clusters)
from gtdb_species_clusters.taxon_utils import (generic_name,
                                                specific_epithet,
                                                canonical_taxon,
                                                canonical_species,
                                                taxon_suffix,
                                                sort_by_naming_priority,
                                                is_placeholder_taxon,
                                                is_placeholder_sp_epithet,
                                                is_alphanumeric_taxon,
                                                is_alphanumeric_sp_epithet,
                                                is_suffixed_taxon,
                                                test_same_epithet)


class PMC_SpeciesNames(object):
    """Establish final species names based on manual curation."""

    def __init__(self, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
        self.final_name_log = open(os.path.join(output_dir, 'final_sp_name.log'), 'w')
        self.final_name_log.write('Genome ID\tGTDB species\tCase\tNote\n')
        
    def parse_updated_species_reps(self, updated_species_reps):
        """Get map indicating the updating of GTDB representatives."""
        
        new_to_prev_rid = {}
        with open(updated_species_reps) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                
                prev_rid = tokens[0]
                new_rid = tokens[1]
                
                if new_rid.lower() != 'none':
                    new_to_prev_rid[new_rid] = prev_rid
                    
        return new_to_prev_rid

    def key_taxon(self, gid, taxonomy):
        """Get key taxon for determine species assignments."""
        
        genus = taxonomy[gid][Taxonomy.GENUS_INDEX]
        species = taxonomy[gid][Taxonomy.SPECIES_INDEX]
        generic = generic_name(species)
        specific = specific_epithet(species)

        return genus, species, generic, specific
        
    def create_placeholder_species_name(self, gid, generic, specific):
        """Create most appropriate placeholder species name.
        
        Species names with a Latin specific name are suffixed, other names
        are turned into sp<accn> IDs.
        """
        
        sp = 's__{} {}'.format(generic, specific)
        
        if is_placeholder_sp_epithet(specific):
            specific = self.sp_name_mngr.numeric_placeholder_sp_epithet(gid)
        else:
            specific = self.sp_name_mngr.suffixed_placeholder_sp_epithet(generic, specific)
            
        return 's__{} {}'.format(generic, specific)
        
    def create_numeric_sp_placeholder(self, gid, generic):
        """Explicitly create a numeric sp<accn> placeholder name."""
        
        specific = self.sp_name_mngr.numeric_placeholder_sp_epithet(gid)
        
        return 's__{} {}'.format(generic, specific)

    def set_type_species(self, 
                            rid, 
                            final_taxonomy, 
                            cur_genomes, 
                            unambiguous_ncbi_sp,
                            case_count):
        """Establish final species name for type species of genus genome."""

        gtdb_genus, gtdb_species, gtdb_generic, gtdb_specific = self.key_taxon(rid, 
                                                                                final_taxonomy)
        ncbi_genus = cur_genomes[rid].ncbi_taxa.genus
        ncbi_sp = cur_genomes[rid].ncbi_taxa.species
        
        if ncbi_sp not in unambiguous_ncbi_sp or rid != unambiguous_ncbi_sp[ncbi_sp][0]:
            self.logger.error('NCBI species {} represented by type species genome {} not designated as an unambiguous NCBI species.'.format(
                                ncbi_sp,
                                rid))
            sys.exit(-1)
        
        note = ''
        if gtdb_species != ncbi_sp:
            gtdb_priority = self.sp_priority_mngr.genus_priority(gtdb_genus, ncbi_genus)

            if gtdb_genus != ncbi_genus and (gtdb_priority is None or gtdb_priority == ncbi_genus):
                case_count['INCONGRUENT_TYPE_SPECIES_GENERIC'] += 1
                note = 'incongruent type species'
                print('INCONGRUENT_TYPE_SPECIES_GENERIC', rid, gtdb_genus, ncbi_genus)
            
            if not test_same_epithet(gtdb_specific, specific_epithet(ncbi_sp)):
                case_count['INCONGRUENT_TYPE_SPECIES_SPECIFIC'] += 1
                note = 'incongruent type species'
                print('INCONGRUENT_TYPE_SPECIES_SPECIFIC', rid, gtdb_species, ncbi_sp)

        return gtdb_species, note
        
    def set_type_strain(self, rid, final_taxonomy, cur_genomes, unambiguous_ncbi_sp, case_count):
        """Establish final species name for type strain of species genome."""

        gtdb_genus, gtdb_species, gtdb_generic, gtdb_specific = self.key_taxon(rid, 
                                                                                final_taxonomy)
                                                                                
        ncbi_sp = cur_genomes[rid].ncbi_taxa.species
        ncbi_specific = specific_epithet(ncbi_sp)
        
        if ncbi_sp not in unambiguous_ncbi_sp or rid != unambiguous_ncbi_sp[ncbi_sp][0]:
            self.logger.error('NCBI species {} represented by type strain genome {} not designated as an unambiguous NCBI species.'.format(
                                ncbi_sp,
                                rid))
            sys.exit(-1)

        note = ''
        final_sp = gtdb_species
        if not test_same_epithet(canonical_taxon(gtdb_specific), ncbi_specific):
            # specific epithets disagree which should never occur for type strains,
            # so changing this to the NCBI proposed specific epithet
            case_count['INCONGRUENT_TYPE_STRAIN'] += 1
            note = 'type strain genome had different GTDB and NCBI specific epithets'
            final_sp = 's__{} {}'.format(gtdb_generic, ncbi_specific)
            print('INCONGRUENT_TYPE_STRAIN', rid, gtdb_species, ncbi_sp, final_sp)
        elif is_placeholder_sp_epithet(gtdb_specific):
            # Latin name with a suffix or sp<accn> specific name is by definition an 
            # error since this is a type strain genome
            case_count['ERRONEOUS_SUFFIX_TYPE_STRAIN'] += 1
            note = 'type strain genome had erroneous placeholder suffix'
            final_sp = 's__{} {}'.format(gtdb_generic, ncbi_specific)
            print('ERRONEOUS_SUFFIX_TYPE_STRAIN', rid, gtdb_species, ncbi_sp, final_sp)

        return final_sp, note
        
    def generate_suffixed_final_name(self, gtdb_generic, prev_gtdb_specific, ncbi_specific):
        """Generate Latin name with suffix, using previously defined suffix where possible."""

        if (prev_gtdb_specific 
            and canonical_taxon(prev_gtdb_specific) == ncbi_specific
            and taxon_suffix(prev_gtdb_specific)):
            # must recycle previously used suffix for this GTDB species cluster
            note = 'Using previously assigned suffix'
            prev_suffix = taxon_suffix(prev_gtdb_specific)
            final_species = 's__{} {}_{}'.format(gtdb_generic, ncbi_specific, prev_suffix)
        else:
            # determine next unused suffix associated with this GTDB species
            note = 'Generated new suffix'
            suffixed_specific_name = self.sp_name_mngr.suffixed_placeholder_sp_epithet(gtdb_generic, ncbi_specific)
            final_species = 's__{} {}'.format(gtdb_generic, suffixed_specific_name)
                    
        return final_species, note
        
    def set_nontype_cluster(self,
                            rid,
                            final_taxonomy, 
                            cur_genomes, 
                            prev_genomes,
                            cur_clusters,
                            new_to_prev_rid,
                            misclassified_gids,
                            unambiguous_ncbi_sp,
                            ambiguous_ncbi_sp,
                            ncbi_synonyms,
                            cur_gtdb_sp,
                            case_count):
        """Establish final species name for non-type material genome."""
        
        forbidden_names = set(['cyanobacterium'])

        # get taxon names
        gtdb_genus, gtdb_species, gtdb_generic, gtdb_specific = self.key_taxon(rid, 
                                                                                final_taxonomy)
                                                                                
        prev_gtdb_species = None
        prev_gtdb_specific = None
        if rid in new_to_prev_rid:
            prev_rid = new_to_prev_rid[rid]
            prev_gtdb_species = prev_genomes[prev_rid].gtdb_taxa.species
            prev_gtdb_specific = specific_epithet(prev_gtdb_species)

        # check if GTDB species cluster consists of genomes with one or multiple
        # NCBI species assignments
        ncbi_species_in_cluster = set()
        for cid in cur_clusters[rid]:
            if cid in misclassified_gids:
                # ignore genomes identified as having an erroneous 
                # NCBI species assignment
                continue
                
            cur_ncbi_species = cur_genomes[cid].ncbi_taxa.species
            cur_ncbi_species = ncbi_synonyms.get(cur_ncbi_species, cur_ncbi_species)
            cur_specific = specific_epithet(cur_ncbi_species)
            
            if cur_ncbi_species != 's__' and cur_specific not in forbidden_names:
                ncbi_species_in_cluster.add(cur_ncbi_species)
                
        final_species = gtdb_species
        note = ''
        if len(ncbi_species_in_cluster) == 1:
            # all genomes have the same NCBI species assignment so the name of the
            # GTDB species cluster should reflect this NCBI name
            ncbi_species = ncbi_species_in_cluster.pop()
            ncbi_species = ncbi_synonyms.get(ncbi_species, ncbi_species)
            ncbi_generic = generic_name(ncbi_species)
            ncbi_specific = specific_epithet(ncbi_species)
        
            if ncbi_species in unambiguous_ncbi_sp and unambiguous_ncbi_sp[ncbi_species][0] == rid:
                # all genomes with this NCBI species assignment are contained in this
                # GTDB species cluster so the cluster should be assigned the NCBI 
                # specific epithet without a suffix
                assert unambiguous_ncbi_sp[ncbi_species][1] == 'UNANIMOUS_CONSENSUS'
                note = 'GTDB cluster is unanimous consensus for NCBI species'
                final_species = 's__{} {}'.format(gtdb_generic, ncbi_specific)

                if final_species in cur_gtdb_sp:
                    # name already exists so this collision needs to be resolved
                    prev_sp_rid = cur_gtdb_sp[final_species][0]
                    if (cur_genomes[prev_sp_rid].is_effective_type_strain()
                        and ncbi_generic != gtdb_generic):
                        # previous assignment should unambiguously have this species name,
                        # and this genome must have a conflict since it was transfer into
                        # this genus and happens to have the same specific name
                        # (e.g. Natronohydrobacter thiooxidans transferred into Roseinatronobacter,
                        #  while already contains the species R. thiooxidans)
                        final_species, note = self.generate_suffixed_final_name(gtdb_generic, prev_gtdb_specific, ncbi_specific)
                        note = 'GTDB cluster contains single unambiguous NCBI species; Species transferred to genus with a species containing the same specific name; {}'.format(note)
                    else:
                        self.logger.error('Manual curation require to resolve name of genome {} as proposed name {} was assigned previously.'.format(
                                            rid, final_species))
            elif ncbi_species in unambiguous_ncbi_sp and unambiguous_ncbi_sp[ncbi_species][0] != rid:
                # all genomes in cluster have the name of an NCBI species which can be unambiguously
                # placed due to have a type strain genome for the NCBI species. However, this is a 
                # different cluster so needs to be given a specific epithet with an appropriate suffix
                final_species, note = self.generate_suffixed_final_name(gtdb_generic, prev_gtdb_specific, ncbi_specific)
                note = 'GTDB cluster contains single unambiguous NCBI species, but this is not the type strain cluster; {}'.format(note)
            elif ncbi_species in ambiguous_ncbi_sp:
                # multiple GTDB species clusters contain genomes with this NCBI
                # species assignment so this cluster should be assigned the NCBI
                # specific epithet with an appropriate suffix
                final_species, note = self.generate_suffixed_final_name(gtdb_generic, prev_gtdb_specific, ncbi_specific)
                note = 'GTDB cluster contains single ambiguous NCBI species; {}'.format(note)
            else:
                self.logger.error('NCBI species name is neither unambiguous or ambiguous: {}'.format(ncbi_cluster_name))
                sys.exit(-1)
        else:
            # species cluster must be given a numeric specific epithet, as it contain either
            # no genomes with an NCBI species assignment or genomes spanning 2 or more NCBI
            # species assignments
            if prev_gtdb_specific and is_alphanumeric_sp_epithet(prev_gtdb_specific):
                # must recycle previously used alphanumeric placeholder name for cluster
                note = 'GTDB cluster contains multiple NCBI species; using previous alphanumeric specific name'
                final_species = 's__{} {}'.format(gtdb_generic, prev_gtdb_specific)
            else:
                note = 'GTDB cluster contains multiple NCBI species; generated new alphanumeric specific name'
                final_species = self.create_numeric_sp_placeholder(rid, gtdb_generic)

        return final_species, note
        
    def resolve_specific_epithet_suffixes(self, final_taxonomy, cur_genomes):
        """Resolve cases where specific epithet needs to be modified to account for genus transfer."""
        
        for rid in final_taxonomy:
            gtdb_genus, gtdb_species, gtdb_generic, gtdb_specific = self.key_taxon(rid, final_taxonomy)
            
            canonical_gtdb_species = canonical_taxon(gtdb_species)
            gtdb_species_corr = self.sp_epithet_mngr.translate_species(canonical_gtdb_species)
            gtdb_specific_corr = specific_epithet(gtdb_species_corr)

            if canonical_taxon(gtdb_specific) != canonical_taxon(gtdb_specific_corr):
                suffix = taxon_suffix(gtdb_specific)
                if suffix is None:
                    final_taxonomy[rid][Taxonomy.SPECIES_INDEX] = gtdb_species_corr
                else:
                    final_taxonomy[rid][Taxonomy.SPECIES_INDEX] = '{}_{}'.format(gtdb_species_corr, suffix)

    def resolve_duplicate_type_strain_species(self, mc_taxonomy, cur_genomes):
        """Resolve cases where 2 or more type strain genomes have the same species name."""
        
        sp_gids = defaultdict(list)
        for rid, taxa in mc_taxonomy.items():
            if cur_genomes[rid].is_gtdb_type_strain():
                sp = taxa[Taxonomy.SPECIES_INDEX]
                sp_gids[sp].append(rid)
                
        resolved_gids = set()
        for sp, rids in sp_gids.items():
            if len(rids) > 1:
                # if only 1 genome has a congruent NCBI and GTDB genus
                # genus assignment, the situation can be automatically resolved
                common_genera = []
                for rid in rids:
                    genus = mc_taxonomy[rid][Taxonomy.GENUS_INDEX]
                    if cur_genomes[rid].ncbi_taxa.genus == canonical_taxon(genus):
                        common_genera.append(rid)
                        
                if len(common_genera) == 1:
                    for rid in rids:
                        if rid not in common_genera:
                            resolved_gids.add(rid)
                            new_sp = self.create_placeholder_species_name(rid, generic_name(sp), specific_epithet(sp))
                            mc_taxonomy[rid][Taxonomy.SPECIES_INDEX] = new_sp
                            
                            print('Modified {} from {} to {} as multiple type strain species had the proposed name and this genome was transferred from a different NCBI genera, {}.'.format( 
                                    rid,
                                    sp,
                                    mc_taxonomy[rid][Taxonomy.SPECIES_INDEX],
                                    cur_genomes[rid].ncbi_taxa.species))
                                    
        return resolved_gids

    def finalize_species_name(self, 
                                rid, 
                                final_sp, 
                                case,
                                note,
                                final_taxonomy, 
                                cur_gtdb_sp):
    
        if final_sp in cur_gtdb_sp:
            prev_rid, prev_case = cur_gtdb_sp[final_sp]
            self.logger.error('Finalized GTDB species name {} already exists: rid={} prev_rid={} prev_case={}'.format(
                                final_sp, rid, prev_rid, prev_case))
                                    
        cur_gtdb_sp[final_sp] = (rid, case)
        
        final_taxonomy[rid][Taxonomy.SPECIES_INDEX] = final_sp
        
        self.final_name_log.write('{}\t{}\t{}\t{}\n'.format(
                                    rid,
                                    final_sp,
                                    case,
                                    note))

    def finalize_species_names(self,
                                mc_taxonomy, 
                                mc_species,
                                cur_clusters,
                                new_to_prev_rid,
                                prev_genomes, 
                                cur_genomes, 
                                misclassified_gids,
                                unambiguous_ncbi_sp, 
                                ambiguous_ncbi_sp,
                                ncbi_synonyms):
        """Establish final species names."""
        
        final_taxonomy = copy.deepcopy(mc_taxonomy)

        # order representatives by naming priority
        rtn = sort_by_naming_priority(final_taxonomy,
                                        cur_clusters,
                                        prev_genomes, 
                                        cur_genomes, 
                                        mc_species)
        (manual_curation_rids, type_species_rids, type_strains_rids, binomial_rids, placeholder_rids) = rtn
        sorted_rids = manual_curation_rids + type_species_rids + type_strains_rids + binomial_rids + placeholder_rids

        # establish final name for each GTDB species cluster
        case_count = defaultdict(int)
        cur_gtdb_sp = {}
        
        for rid in manual_curation_rids:
            case = 'MANUAL_CURATION'
            note = 'Species name set by manual curation'
            final_sp = mc_species[rid]
            
            self.finalize_species_name(rid,
                                        final_sp,
                                        case,
                                        note,
                                        final_taxonomy, 
                                        cur_gtdb_sp)

        for rid in type_species_rids:
            if rid == 'G001884735':
                print('!!!HERE!!!', rid, cur_genomes[rid].is_gtdb_type_species())
                
            case = 'TYPE_SPECIES_OF_GENUS'
            final_sp, note = self.set_type_species(rid,
                                                    final_taxonomy, 
                                                    cur_genomes, 
                                                    unambiguous_ncbi_sp,
                                                    case_count)
                                                    
            self.finalize_species_name(rid,
                                        final_sp,
                                        case,
                                        note,
                                        final_taxonomy, 
                                        cur_gtdb_sp)
                                        
        for rid in type_strains_rids:
            case = 'TYPE_STRAIN_OF_SPECIES'
            final_sp, note = self.set_type_strain(rid,
                                                    final_taxonomy, 
                                                    cur_genomes, 
                                                    unambiguous_ncbi_sp,
                                                    case_count)
                                                        
            self.finalize_species_name(rid,
                                        final_sp,
                                        case,
                                        note,
                                        final_taxonomy, 
                                        cur_gtdb_sp)
                                        
        # resolve type strain genomes with the same name
        resolved_type_strain_gids = self.resolve_duplicate_type_strain_species(final_taxonomy, 
                                                                                cur_genomes)

        # establish final names for non-type genomes
        for rid in binomial_rids + placeholder_rids:
            case = 'NONTYPE_SPECIES_CLUSTER'
            final_sp, note = self.set_nontype_cluster(rid,
                                                        final_taxonomy, 
                                                        cur_genomes, 
                                                        prev_genomes,
                                                        cur_clusters,
                                                        new_to_prev_rid,
                                                        misclassified_gids,
                                                        unambiguous_ncbi_sp,
                                                        ambiguous_ncbi_sp,
                                                        ncbi_synonyms,
                                                        cur_gtdb_sp,
                                                        case_count)

            self.finalize_species_name(rid,
                                        final_sp,
                                        case,
                                        note,
                                        final_taxonomy, 
                                        cur_gtdb_sp)

        for case, count in case_count.items():
            print('{}\t{}'.format(case, count))
            
        # resolve specific epithets requiring changes due to genus transfers
        self.resolve_specific_epithet_suffixes(final_taxonomy, cur_genomes)
            
        return final_taxonomy

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

    def resolve_merged_prev_representatives(self, cur_rid, prev_genomes, updated_gtdb_rids, prev_rids_in_cluster):
        """Determine best representative to associated with new species cluster contains multiple previous representative.
        
        This association between previous and new representatives is used in some instances to establish
        the most suitable name for a species cluster.        
        """
        
        if cur_rid in prev_genomes:
            # if current representative was a previous representative, 
            # it is the most natural choice
            return cur_rid
        elif cur_rid in updated_gtdb_rids:
            prev_rid = updated_gtdb_rids[cur_rid]
            if prev_rid in prev_rids_in_cluster:
                # makes sense to use the previous representative that was explicitly
                # updated to the new, current representative
                return prev_rid
            else:
                self.logger.error('Updated representative of cluster {} no longer contains the previous representative of cluster: {}'.format(
                                    cur_rid,
                                    prev_rid))
                sys.exit(-1)
        else:
            self.logger.error('Updated representative of cluster {} contains multiple previous representative and situation could not be resolved: {}'.format(
                                    cur_rid,
                                    prev_rid))
            sys.exit(-1)
        
    def run(self,
                manual_taxonomy,
                manual_sp_names,
                pmc_custom_species,
                gtdb_clusters_file,
                prev_gtdb_metadata_file,
                cur_gtdb_metadata_file,
                uba_genome_paths,
                qc_passed_file,
                ncbi_genbank_assembly_file,
                untrustworthy_type_file,
                synonym_file,
                updated_species_reps,
                gtdb_type_strains_ledger,
                sp_priority_ledger,
                genus_priority_ledger,
                specific_epithet_ledger,
                dsmz_bacnames_file):
        """Finalize species names based on results of manual curation."""

        # read manually-curated taxonomy
        self.logger.info('Parsing manually-curated taxonomy.')
        mc_taxonomy = Taxonomy().read(manual_taxonomy, use_canonical_gid=True)
        self.logger.info(' - identified taxonomy strings for {:,} genomes.'.format(
                            len(mc_taxonomy)))
                            
        # read species names explicitly set via manual curation
        self.logger.info('Parsing manually-curated species.')
        mc_species = {}
        with open(manual_sp_names) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                mc_species[tokens[0]] = tokens[2]
        self.logger.info(' - identified manually-curated species names for {:,} genomes.'.format(
                            len(mc_species)))
                            
        # read post-curation, manually defined species
        self.logger.info('Parsing post-curation, manually-curated species.')
        pmc_species = {}
        with open(pmc_custom_species) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                gid = tokens[0]
                species = tokens[1]
                if gid in mc_species:
                    self.logger.warning('Manually-curated genome {} reassigned from {} to {}.'.format(
                                            gid, mc_species[gid], species))
                pmc_species[gid] = species
                mc_species[gid] = species
        self.logger.info(' - identified post-curation, manually-curated species names for {:,} genomes.'.format(
                            len(pmc_species)))

        # create previous and current GTDB genome sets
        self.logger.info('Creating previous GTDB genome set.')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                uba_genome_file=uba_genome_paths,
                                                ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                                untrustworthy_type_ledger=untrustworthy_type_file)
        self.logger.info(' - previous genome set has {:,} species clusters spanning {:,} genomes.'.format(
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
        self.logger.info(f' - current genome set contains {len(cur_genomes):,} genomes.')
        
        cur_genomes.set_prev_gtdb_classifications(prev_genomes)

        # read ledger with explicit names for type strain genomes
        self.logger.info('Validating species names given in GTDB type strain ledger.')
        gtdb_type_strain_ledger = self.parse_gtdb_type_strain_ledger(gtdb_type_strains_ledger, cur_genomes)
        self.logger.info(' - identified {:,} genomes in ledger.'.format(len(gtdb_type_strain_ledger)))
        for gid in gtdb_type_strain_ledger:
            if gid not in mc_taxonomy:
                continue
                
            cur_gtdb_sp = mc_taxonomy[gid][Taxonomy.SPECIES_INDEX]
            ledger_sp = gtdb_type_strain_ledger[gid]
            if cur_gtdb_sp != ledger_sp:
                self.logger.warning('Assigned species name for {} conflicts with GTDB type strain ledger: {} {}'.format(
                                    gid, cur_gtdb_sp, ledger_sp))
        
        # read named GTDB species clusters
        self.logger.info('Reading GTDB species clusters.')
        cur_clusters, rep_radius = read_clusters(gtdb_clusters_file)
        self.logger.info(' - identified {:,} clusters spanning {:,} genomes.'.format(
                            len(cur_clusters),
                            sum([len(gids) + 1 for gids in cur_clusters.values()])))
        assert len(set(mc_taxonomy) - set(cur_clusters)) == 0

        # get list of synonyms in order to restrict usage of species names
        self.logger.info('Reading GTDB synonyms.')
        ncbi_species_mngr = NCBI_SpeciesManager(cur_genomes, cur_clusters, self.output_dir)
        ncbi_synonyms = ncbi_species_mngr.parse_synonyms_table(synonym_file)
        self.logger.info(' - identified {:,} synonyms from {:,} distinct species.'.format(
                            len(ncbi_synonyms),
                            len(set(ncbi_synonyms.values()))))

        # create specific epithet manager
        self.sp_epithet_mngr = SpecificEpithetManager()
        self.sp_epithet_mngr.parse_specific_epithet_ledger(specific_epithet_ledger)
        self.sp_epithet_mngr.infer_epithet_map(mc_taxonomy, 
                                                mc_species, 
                                                cur_genomes, 
                                                cur_clusters)
        self.sp_epithet_mngr.write_diff_epithet_map(os.path.join(self.output_dir, 'specific_epithet_diff_map.tsv'))
        self.sp_epithet_mngr.write_epithet_map(os.path.join(self.output_dir, 'specific_epithet_map.tsv'))
        
        # create species name manager
        self.logger.info('Initializing species name manager.')
        self.sp_name_mngr = SpeciesNameManager(prev_genomes, 
                                                    cur_genomes,
                                                    None)
                                                    
        # initialize species priority manager
        self.sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger,
                                                        genus_priority_ledger,
                                                        dsmz_bacnames_file)
                                                        
        # establish state of NCBI species
        ncbi_species_mngr = NCBI_SpeciesManager(cur_genomes, cur_clusters, self.output_dir)
        
        # get list of synonyms in order to restrict usage of species names
        self.logger.info('Reading GTDB synonyms.')
        ncbi_synonyms = ncbi_species_mngr.parse_synonyms_table(synonym_file)
        self.logger.info(' - identified {:,} synonyms from {:,} distinct species.'.format(
                            len(ncbi_synonyms),
                            len(set(ncbi_synonyms.values()))))

        # identified genomes with misclassified species assignments at NCBI
        self.logger.info('Identify genomes with misclassified NCBI species assignments.')
        misclassified_gids = ncbi_species_mngr.identify_misclassified_genomes(mc_taxonomy,
                                                                                ncbi_synonyms, 
                                                                                mc_species)
        
        # classify each Latin NCBI species as either unambiguously or ambiguously assigned
        self.logger.info('Classifying NCBI species as unambiguous, ambiguous, or synonym.')
        ncbi_classification_table = os.path.join(self.output_dir, 'ncbi_sp_classification.tsv')
        if os.path.exists(ncbi_classification_table):
            self.logger.warning('Reading classification of NCBI species from existing table: {}'.format(
                                    ncbi_classification_table))
            ncbi_species, unambiguous_ncbi_sp, ambiguous_ncbi_sp = ncbi_species_mngr.parse_ncbi_classification_table(
                                                                    ncbi_classification_table)
        else:
            ncbi_species, unambiguous_ncbi_sp, ambiguous_ncbi_sp = ncbi_species_mngr.classify_ncbi_species(
                                                                        ncbi_synonyms, 
                                                                        misclassified_gids)
        
        unambiguous_classifications = [classification for rid, classification in unambiguous_ncbi_sp.values()]
        unambiguous_type_strain_count = unambiguous_classifications.count('TYPE_STRAIN_GENOME')
        unambiguous_unanimous_consensus_count = unambiguous_classifications.count('UNANIMOUS_CONSENSUS')
        assert unambiguous_type_strain_count + unambiguous_unanimous_consensus_count == len(unambiguous_ncbi_sp)
                
        self.logger.info(' - identified {:,} NCBI species.'.format(len(ncbi_species)))
        self.logger.info(' - identified {:,} ({:.2f}%) as unambiguous assignments.'.format(
                            len(unambiguous_ncbi_sp),
                            len(unambiguous_ncbi_sp) * 100.0 / len(ncbi_species)))
        self.logger.info('   - assigned {:,} ({:.2f}%) by type strain genome.'.format(
                            unambiguous_type_strain_count,
                            unambiguous_type_strain_count * 100.0 / len(ncbi_species)))
        self.logger.info('   - assigned {:,} ({:.2f}%) by unanimous consensus.'.format(
                            unambiguous_unanimous_consensus_count,
                            unambiguous_unanimous_consensus_count * 100.0 / len(ncbi_species)))
        self.logger.info(' - identified {:,} ({:.2f}%) as ambiguous assignments.'.format(
                            len(ambiguous_ncbi_sp),
                            len(ambiguous_ncbi_sp) * 100.0 / len(ncbi_species)))
        self.logger.info(' - identified {:,} ({:.2f}%) as GTDB synonyms.'.format(
                            len(ncbi_synonyms),
                            len(ncbi_synonyms) * 100.0 / len(ncbi_species)))

        # get mapping between new and old representatives. This can't just be taken
        # from the `updated_species_reps` file since the final results of de novo
        # cluster can, in some rare cases, cause small movements in the genomes 
        # associated with each species clusters and the generation or loss of
        # a representative
        self.logger.info('Mapping current GTDB representatives to previous representatives.')
        updated_gtdb_rids = self.parse_updated_species_reps(updated_species_reps)
        
        prev_rids = set(prev_genomes.sp_clusters)
        new_to_prev_rid = {}
        for cur_rid, cur_cids in cur_clusters.items():
            prev_rids_in_cluster = prev_rids.intersection(cur_cids)
            if len(prev_rids_in_cluster) == 1:
                new_to_prev_rid[cur_rid] = prev_rids_in_cluster.pop()
            elif len(prev_rids_in_cluster) > 1:
                resolved_rid = self.resolve_merged_prev_representatives(cur_rid, 
                                                                            prev_genomes, 
                                                                            updated_gtdb_rids,
                                                                            prev_rids_in_cluster)
               
                if resolved_rid:
                    self.logger.info(' - cluster {} contains multiple previous representatives and has been associated with representative {}.'.format(
                                        cur_rid, resolved_rid))
                  
                    new_to_prev_rid[cur_rid] = resolved_rid
                                        
        self.logger.info(' - mapped {:,} current representatives to previous representatives.'.format(len(new_to_prev_rid)))

        # establish appropriate species names for GTDB clusters with new representatives
        final_taxonomy = self.finalize_species_names(mc_taxonomy, 
                                                    mc_species,
                                                    cur_clusters,
                                                    new_to_prev_rid,
                                                    prev_genomes, 
                                                    cur_genomes,
                                                    misclassified_gids,
                                                    unambiguous_ncbi_sp, 
                                                    ambiguous_ncbi_sp,
                                                    ncbi_synonyms)
                                        
        # write out final taxonomy
        fout = open(os.path.join(self.output_dir, 'final_taxonomy.tsv'), 'w')
        for rid, taxa in final_taxonomy.items():
            fout.write('{}\t{}\n'.format(
                        rid,
                        ';'.join(taxa)))
        fout.close()
        
        # write out final taxonomy with additional metadata
        fout = open(os.path.join(self.output_dir, 'final_taxonomy_metadata.tsv'), 'w')
        fout.write('Genome ID\tGTDB taxonomy\tType strain\tEffective type strain\tType species')
        fout.write('\tNCBI species\tNCBI species classification\n')
        for rid, taxa in final_taxonomy.items():
            ncbi_species = cur_genomes[rid].ncbi_taxa.species
            ncbi_sp_classification = 'n/a'
            if ncbi_species in unambiguous_ncbi_sp:
                ncbi_sp_classification = 'UNAMBIGUOUS'
            elif ncbi_species in ambiguous_ncbi_sp:
                ncbi_sp_classification = 'AMBIGUOUS'
            elif ncbi_species in ncbi_synonyms:
                ncbi_sp_classification = 'SYNONYM: {}'.format(ncbi_synonyms[ncbi_species])

            fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        rid,
                        ';'.join(taxa),
                        cur_genomes[rid].is_gtdb_type_strain(),
                        cur_genomes[rid].is_effective_type_strain(),
                        cur_genomes[rid].is_gtdb_type_species(),
                        ncbi_species,
                        ncbi_sp_classification))
        fout.close()
        
        self.final_name_log.close()
                                        