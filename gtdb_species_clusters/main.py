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
import csv
import logging
import random
from collections import defaultdict

from biolib.common import check_file_exists, make_sure_path_exists, is_float
from biolib.external.execute import check_dependencies
from biolib.taxonomy import Taxonomy
from biolib.newick import parse_label

from gtdb_species_clusters.qc_genomes import QcGenomes
from gtdb_species_clusters.mash import Mash
from gtdb_species_clusters.select_type_genomes import SelectTypeGenomes
from gtdb_species_clusters.cluster_named_types import ClusterNamedTypes
from gtdb_species_clusters.cluster_de_novo import ClusterDeNovo
from gtdb_species_clusters.cluster_user import ClusterUser
from gtdb_species_clusters.tree_gids import TreeGIDs
from gtdb_species_clusters.cluster_stats import ClusterStats

from gtdb_species_clusters.update_new_genomes import NewGenomes
from gtdb_species_clusters.update_resolve_types import ResolveTypes
from gtdb_species_clusters.update_gtdbtk import GTDB_Tk
from gtdb_species_clusters.update_rep_changes import RepChanges
from gtdb_species_clusters.update_rep_actions import RepActions
from gtdb_species_clusters.update_select_reps import UpdateSelectRepresentatives
from gtdb_species_clusters.update_cluster_named_reps import UpdateClusterNamedReps
from gtdb_species_clusters.update_synonyms import UpdateSynonyms
from gtdb_species_clusters.update_cluster_de_novo import UpdateClusterDeNovo
from gtdb_species_clusters.update_genus_names import UpdateGenusNames
from gtdb_species_clusters.update_curation_trees import UpdateCurationTrees
from gtdb_species_clusters.update_species_init import UpdateSpeciesInit

from gtdb_species_clusters.pmc_checks import PMC_Checks
from gtdb_species_clusters.pmc_check_type_species import PMC_CheckTypeSpecies
from gtdb_species_clusters.pmc_check_type_strains import PMC_CheckTypeStrains
from gtdb_species_clusters.pmc_species_names import PMC_SpeciesNames
from gtdb_species_clusters.pmc_validation import PMC_Validation
from gtdb_species_clusters.update_summary_stats import UpdateSummaryStats

from gtdb_species_clusters.merge_test import MergeTest
from gtdb_species_clusters.intra_sp_derep import IntraSpeciesDereplication

from gtdb_species_clusters.exceptions import GTDB_Error


class OptionsParser():
    def __init__(self):
        """Initialization"""
        self.logger = logging.getLogger()
        
    def qc_genomes(self, args):
        """Quality check all potential GTDB genomes."""
        
        check_file_exists(args.gtdb_metadata_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.gtdb_domain_report)
        check_file_exists(args.qc_exception_file)
        check_file_exists(args.species_exception_file)
        make_sure_path_exists(args.output_dir)

        try:
            p = QcGenomes(args.ncbi_genbank_assembly_file,
                        args.gtdb_domain_report,
                        args.qc_exception_file,
                        args.species_exception_file,
                        args.min_comp,
                        args.max_cont,
                        args.min_quality,
                        args.sh_exception,
                        args.min_perc_markers,
                        args.max_contigs,
                        args.min_N50,
                        args.max_ambiguous,
                        args.output_dir)
        except GTDB_Error as e:
            print(e.message)
            raise SystemExit

        self.logger.info('Quality checking information written to: %s' % args.output_dir)
        
    def select_type_genomes(self, args):
        """Select representative genomes for named species."""

        check_file_exists(args.qc_file)
        check_file_exists(args.gtdb_metadata_file)
        check_file_exists(args.genome_path_file)
        check_file_exists(args.prev_rep_file)
        check_file_exists(args.ncbi_refseq_assembly_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.gtdb_domain_report)
        check_file_exists(args.species_exception_file)
        check_file_exists(args.gtdb_type_genome_file)
        make_sure_path_exists(args.output_dir)

        try:
            p = SelectTypeGenomes(args.ani_cache_file, args.cpus, args.output_dir)
            p.run(args.qc_file,
                        args.gtdb_metadata_file,
                        args.ltp_blast_file,
                        args.genome_path_file,
                        args.prev_rep_file,
                        args.ncbi_refseq_assembly_file,
                        args.ncbi_genbank_assembly_file,
                        args.gtdb_domain_report,
                        args.species_exception_file,
                        args.gtdb_type_genome_file)
        except GTDB_Error as e:
            print(e.message)
            raise SystemExit

        self.logger.info('GTDB type genomes written to: %s' % args.output_dir)

    def cluster_named_types(self, args):
        """Cluster genomes to selected GTDB type genomes."""

        check_file_exists(args.qc_file)
        check_file_exists(args.gtdb_metadata_file)
        check_file_exists(args.genome_path_file)
        check_file_exists(args.named_type_genome_file)
        check_file_exists(args.type_genome_ani_file)
        check_file_exists(args.species_exception_file)
        make_sure_path_exists(args.output_dir)

        try:
            p = ClusterNamedTypes(args.ani_sp,
                                    args.af_sp,
                                    args.ani_cache_file, 
                                    args.cpus,
                                    args.output_dir)
            p.run(args.qc_file,
                    args.gtdb_metadata_file,
                    args.genome_path_file,
                    args.named_type_genome_file,
                    args.type_genome_ani_file,
                    args.mash_sketch_file,
                    args.species_exception_file)
        except GTDB_Error as e:
            print(e.message)
            raise SystemExit

        self.logger.info('Clustering results written to: %s' % args.output_dir)
        
    def cluster_de_novo(self, args):
        """Infer de novo species clusters and type genomes for remaining genomes."""

        check_file_exists(args.qc_file)
        check_file_exists(args.gtdb_metadata_file)
        check_file_exists(args.gtdb_user_genomes_file)
        check_file_exists(args.genome_path_file)
        check_file_exists(args.type_genome_cluster_file)
        check_file_exists(args.type_genome_synonym_file)
        check_file_exists(args.ncbi_refseq_assembly_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.ani_af_nontype_vs_type)
        check_file_exists(args.species_exception_file)
        make_sure_path_exists(args.output_dir)

        try:
            p = ClusterDeNovo(args.ani_sp,
                                    args.af_sp,
                                    args.ani_cache_file, 
                                    args.cpus,
                                    args.output_dir)
            p.run(args.qc_file,
                        args.gtdb_metadata_file,
                        args.gtdb_user_genomes_file,
                        args.genome_path_file,
                        args.type_genome_cluster_file,
                        args.type_genome_synonym_file,
                        args.ncbi_refseq_assembly_file,
                        args.ncbi_genbank_assembly_file,
                        args.ani_af_nontype_vs_type,
                        args.species_exception_file,
                        args.rnd_type_genome)
        except GTDB_Error as e:
            print(e.message)
            raise SystemExit

        self.logger.info('Clustering results written to: %s' % args.output_dir)
        
    def cluster_user(self, args):
        """Cluster User genomes to GTDB species clusters."""

        check_file_exists(args.gtdb_metadata_file)
        check_file_exists(args.genome_path_file)
        check_file_exists(args.final_cluster_file)
        make_sure_path_exists(args.output_dir)

        try:
            p = ClusterUser(args.ani_cache_file, 
                                args.cpus,
                                args.output_dir)
            p.run(args.gtdb_metadata_file,
                        args.genome_path_file,
                        args.final_cluster_file)
        except GTDB_Error as e:
            print(e.message)
            raise SystemExit

        self.logger.info('Clustering results written to: %s' % args.output_dir)
        
    def tree_gids(self, args):
        """Determine genome IDs for test/validation tree."""

        check_file_exists(args.qc_file)
        check_file_exists(args.gtdb_metadata_file)
        check_file_exists(args.gtdb_final_clusters)
        check_file_exists(args.species_exception_file)

        try:
            p = TreeGIDs()
            p.run(args.qc_file,
                    args.gtdb_metadata_file,
                    args.gtdb_final_clusters,
                    args.species_exception_file,
                    args.output_dir)
        except GTDB_Error as e:
            print(e.message)
            raise SystemExit

        self.logger.info('Results written to: %s' % args.output_dir)
            
    def u_new_genomes(self, args):
        """Identify new and updated genomes."""

        check_file_exists(args.prev_gtdb_metadata_file)
        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.cur_genome_paths)
        check_file_exists(args.ncbi_assembly_summary_genbank)

        make_sure_path_exists(args.output_dir)
        
        try:
            p = NewGenomes(args.output_dir)
            p.run(args.prev_gtdb_metadata_file,
                    args.cur_gtdb_metadata_file,
                    args.cur_genome_paths,
                    args.ncbi_assembly_summary_genbank)
        except GTDB_Error as e:
            print(e.message)
            raise SystemExit
            
        self.logger.info('Done.')
        
    def u_qc_genomes(self, args):
        """Quality check new and updated genomes."""
        
        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.cur_uba_gid_file)
        check_file_exists(args.cur_genbank_assembly_file)
        check_file_exists(args.cur_gtdb_domain_report)
        check_file_exists(args.qc_exception_file)
        make_sure_path_exists(args.output_dir)

        try:
            p = QcGenomes()
            p.run(args.cur_gtdb_metadata_file,
                        args.cur_uba_gid_file,
                        args.cur_genbank_assembly_file,
                        args.cur_gtdb_domain_report,
                        args.qc_exception_file,
                        args.min_comp,
                        args.max_cont,
                        args.min_quality,
                        args.sh_exception,
                        args.min_perc_markers,
                        args.max_contigs,
                        args.min_N50,
                        args.max_ambiguous,
                        args.output_dir)
        except GTDB_Error as e:
            print(e.message)
            raise SystemExit

        self.logger.info('Quality checking information written to: %s' % args.output_dir)
        
    def u_resolve_types(self, args):
        """Resolve cases where a species has multiple genomes assembled from the type strain."""

        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.cur_genomic_path_file)
        check_file_exists(args.qc_passed_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.ltp_taxonomy_file)
        check_file_exists(args.gtdb_type_strains_ledger)
        check_file_exists(args.untrustworthy_type_ledger)
        make_sure_path_exists(args.output_dir)
        
        p = ResolveTypes(args.ani_cache_file, args.cpus, args.output_dir)
        p.run(args.cur_gtdb_metadata_file,
                args.cur_genomic_path_file,
                args.qc_passed_file,
                args.ncbi_genbank_assembly_file,
                args.ltp_taxonomy_file,
                args.gtdb_type_strains_ledger,
                args.untrustworthy_type_ledger)
        
        self.logger.info('Done.')
        
    def u_gtdbtk(self, args):
        """Perform initial classification of new and updated genomes using GTDB-Tk."""
        
        check_file_exists(args.genomes_new_updated_file)
        check_file_exists(args.qc_passed_file)
        make_sure_path_exists(args.output_dir)
        
        p = GTDB_Tk(args.cpus, args.output_dir)
        p.run(args.genomes_new_updated_file,
                args.qc_passed_file,
                args.batch_size)

        self.logger.info('Done.')
            
    def u_rep_changes(self, args):
        """Identify species representatives that have changed from previous release."""

        check_file_exists(args.prev_gtdb_metadata_file)
        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.cur_uba_gid_file)
        check_file_exists(args.genomes_new_updated_file)
        check_file_exists(args.qc_passed_file)
        check_file_exists(args.gtdbtk_classify_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.untrustworthy_type_file)
        check_file_exists(args.gtdb_type_strains_ledger)
        make_sure_path_exists(args.output_dir)
        
        p = RepChanges(args.output_dir)
        p.run(args.prev_gtdb_metadata_file,
                args.cur_gtdb_metadata_file,
                args.cur_uba_gid_file,
                args.genomes_new_updated_file,
                args.qc_passed_file,
                args.gtdbtk_classify_file,
                args.ncbi_genbank_assembly_file,
                args.untrustworthy_type_file,
                args.gtdb_type_strains_ledger)
        
        self.logger.info('Done.')
        
    def u_rep_actions(self, args):
        """Perform initial actions required for changed representatives."""
        
        check_file_exists(args.rep_change_summary_file)
        check_file_exists(args.prev_genomic_path_file)
        check_file_exists(args.cur_genomic_path_file)
        check_file_exists(args.uba_genome_paths)
        check_file_exists(args.genomes_new_updated_file)
        check_file_exists(args.qc_passed_file)
        check_file_exists(args.gtdbtk_classify_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.untrustworthy_type_file)
        check_file_exists(args.gtdb_type_strains_ledger)
        check_file_exists(args.sp_priority_ledger)
        make_sure_path_exists(args.output_dir)
        
        p = RepActions(args.ani_cache_file, args.cpus, args.output_dir)
        p.run(args.rep_change_summary_file,
                args.prev_gtdb_metadata_file,
                args.prev_genomic_path_file,
                args.cur_gtdb_metadata_file,
                args.cur_genomic_path_file,
                args.uba_genome_paths,
                args.genomes_new_updated_file,
                args.qc_passed_file,
                args.gtdbtk_classify_file,
                args.ncbi_genbank_assembly_file,
                args.untrustworthy_type_file,
                args.gtdb_type_strains_ledger,
                args.sp_priority_ledger)
        
        self.logger.info('Done.')
        
    def u_sel_reps(self, args):
        """Select representatives for all named species at NCBI."""
        
        check_file_exists(args.updated_sp_cluster_file)
        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.cur_genomic_path_file)
        check_file_exists(args.uba_genome_paths)
        check_file_exists(args.qc_passed_file)
        check_file_exists(args.gtdbtk_classify_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.untrustworthy_type_file)
        check_file_exists(args.gtdb_type_strains_ledger)
        check_file_exists(args.sp_priority_ledger)
        make_sure_path_exists(args.output_dir)
        
        p = UpdateSelectRepresentatives(args.ani_cache_file, 
                                    args.cpus, 
                                    args.output_dir)
        p.run(args.updated_sp_cluster_file,
                args.cur_gtdb_metadata_file,
                args.cur_genomic_path_file,
                args.uba_genome_paths,
                args.qc_passed_file,
                args.gtdbtk_classify_file,
                args.ncbi_genbank_assembly_file,
                args.untrustworthy_type_file,
                args.gtdb_type_strains_ledger,
                args.sp_priority_ledger)
        
        self.logger.info('Done.')
        
    def u_cluster_named_reps(self, args):
        """Cluster genomes to selected GTDB representatives."""
        
        check_file_exists(args.named_rep_file)
        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.cur_genomic_path_file)
        check_file_exists(args.uba_genome_paths)
        check_file_exists(args.qc_passed_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.untrustworthy_type_file)
        check_file_exists(args.rep_mash_sketch_file)
        check_file_exists(args.rep_ani_file)
        check_file_exists(args.gtdb_type_strains_ledger)
        make_sure_path_exists(args.output_dir)
        
        p = UpdateClusterNamedReps(args.ani_sp,
                                    args.af_sp,
                                    args.ani_cache_file, 
                                    args.cpus, 
                                    args.output_dir)
        p.run(args.named_rep_file,
                args.cur_gtdb_metadata_file,
                args.cur_genomic_path_file,
                args.uba_genome_paths,
                args.qc_passed_file,
                args.ncbi_genbank_assembly_file,
                args.untrustworthy_type_file,
                args.rep_mash_sketch_file,
                args.rep_ani_file,
                args.gtdb_type_strains_ledger)
        
        self.logger.info('Done.')
        
    def u_synonyms(self, args):
        """Determine synonyms for validly or effectively published species."""
        
        check_file_exists(args.gtdb_clusters_file)
        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.uba_genome_paths)
        check_file_exists(args.qc_passed_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.untrustworthy_type_file)
        check_file_exists(args.ani_af_rep_vs_nonrep)
        check_file_exists(args.gtdb_type_strains_ledger)
        check_file_exists(args.sp_priority_ledger)
        check_file_exists(args.genus_priority_ledger)
        check_file_exists(args.dsmz_bacnames_file)
        make_sure_path_exists(args.output_dir)
        
        p = UpdateSynonyms(args.output_dir)
        p.run(args.gtdb_clusters_file,
                args.cur_gtdb_metadata_file,
                args.uba_genome_paths,
                args.qc_passed_file,
                args.ncbi_genbank_assembly_file,
                args.untrustworthy_type_file,
                args.ani_af_rep_vs_nonrep,
                args.gtdb_type_strains_ledger,
                args.sp_priority_ledger,
                args.genus_priority_ledger,
                args.dsmz_bacnames_file)
        
        self.logger.info('Done.')
        
    def u_cluster_de_novo(self, args):
        """Infer de novo species clusters and representatives for remaining genomes."""
        
        check_file_exists(args.named_cluster_file)
        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.cur_genomic_path_file)
        check_file_exists(args.uba_genome_paths)
        check_file_exists(args.qc_passed_file)
        check_file_exists(args.gtdbtk_classify_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.untrustworthy_type_file)
        check_file_exists(args.ani_af_rep_vs_nonrep)
        check_file_exists(args.gtdb_type_strains_ledger)
        make_sure_path_exists(args.output_dir)
        
        p = UpdateClusterDeNovo(args.ani_sp,
                                    args.af_sp,
                                    args.ani_cache_file, 
                                    args.cpus, 
                                    args.output_dir)
        p.run(args.named_cluster_file,
                args.cur_gtdb_metadata_file,
                args.cur_genomic_path_file,
                args.uba_genome_paths,
                args.qc_passed_file,
                args.gtdbtk_classify_file,
                args.ncbi_genbank_assembly_file,
                args.untrustworthy_type_file,
                args.ani_af_rep_vs_nonrep,
                args.gtdb_type_strains_ledger)
        
        self.logger.info('Done.')
        
    def u_genus_names(self, args):
        """Update genus names as a precursor for establish binomial species names."""
        
        check_file_exists(args.gtdb_clusters_file)
        check_file_exists(args.prev_gtdb_metadata_file)
        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.uba_genome_paths)
        check_file_exists(args.qc_passed_file)
        check_file_exists(args.gtdbtk_classify_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.untrustworthy_type_file)
        check_file_exists(args.synonym_file)
        check_file_exists(args.gtdb_type_strains_ledger)
        check_file_exists(args.sp_priority_ledger)
        check_file_exists(args.gtdb_taxa_updates_ledger)
        check_file_exists(args.dsmz_bacnames_file)
        make_sure_path_exists(args.output_dir)

        p = UpdateGenusNames(args.output_dir)
        p.run(args.gtdb_clusters_file,
                args.prev_gtdb_metadata_file,
                args.cur_gtdb_metadata_file,
                args.uba_genome_paths,
                args.qc_passed_file,
                args.gtdbtk_classify_file,
                args.ncbi_genbank_assembly_file,
                args.untrustworthy_type_file,
                args.synonym_file,
                args.gtdb_type_strains_ledger,
                args.sp_priority_ledger,
                args.gtdb_taxa_updates_ledger,
                args.dsmz_bacnames_file)
        
        self.logger.info('Done.')
        
    def u_curation_trees(self, args):
        """Produce curation trees highlighting new NCBI taxa."""
        
        check_file_exists(args.gtdb_clusters_file)
        check_file_exists(args.prev_gtdb_metadata_file)
        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.uba_genome_paths)
        check_file_exists(args.qc_passed_file)

        p = UpdateCurationTrees(args.output_dir, args.output_prefix)
        p.run(args.gtdb_clusters_file,
                args.prev_gtdb_metadata_file,
                args.cur_gtdb_metadata_file,
                args.uba_genome_paths,
                args.qc_passed_file)
        
        self.logger.info('Done.')
        
    def u_species_init(self, args):
        """Produce initial best guess at GTDB species clusters."""
        
        check_file_exists(args.gtdb_clusters_file)
        check_file_exists(args.prev_gtdb_metadata_file)
        check_file_exists(args.prev_genomic_path_file)
        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.cur_genomic_path_file)
        check_file_exists(args.uba_genome_paths)
        check_file_exists(args.qc_passed_file)
        check_file_exists(args.gtdbtk_classify_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.untrustworthy_type_file)
        check_file_exists(args.synonym_file)
        check_file_exists(args.gtdb_type_strains_ledger)
        check_file_exists(args.sp_priority_ledger)
        check_file_exists(args.genus_priority_ledger)
        check_file_exists(args.gtdb_taxa_updates_ledger)
        check_file_exists(args.dsmz_bacnames_file)
        make_sure_path_exists(args.output_dir)

        p = UpdateSpeciesInit(args.ani_cache_file, 
                                args.cpus, 
                                args.output_dir)
        p.run(args.gtdb_clusters_file,
                args.prev_gtdb_metadata_file,
                args.prev_genomic_path_file,
                args.cur_gtdb_metadata_file,
                args.cur_genomic_path_file,
                args.uba_genome_paths,
                args.qc_passed_file,
                args.gtdbtk_classify_file,
                args.ncbi_genbank_assembly_file,
                args.untrustworthy_type_file,
                args.synonym_file,
                args.gtdb_type_strains_ledger,
                args.sp_priority_ledger,
                args.genus_priority_ledger,
                args.gtdb_taxa_updates_ledger,
                args.dsmz_bacnames_file)
        
        self.logger.info('Done.')
        
    def pmc_manual_species(self, args):
        """Identify species names manually set by curators."""
        
        check_file_exists(args.init_taxonomy)
        check_file_exists(args.manually_curated_tree)
        make_sure_path_exists(args.output_dir)

        p = PMC_Checks(args.output_dir)
        p.manual_species(args.init_taxonomy, 
                            args.manually_curated_tree)
        
        self.logger.info('Done.')
        
    def pmc_replace_generic(self, args):
        """Replace generic names with genus assignment."""
        
        check_file_exists(args.manual_species_names)
        check_file_exists(args.manual_taxonomy)
        make_sure_path_exists(args.output_dir)

        p = PMC_Checks(args.output_dir)
        p.replace_generic(args.manual_species_names, 
                            args.manual_taxonomy)
        
        self.logger.info('Done.')
        
    def pmc_check_type_species(self, args):
        """Check for agreement between GTDB genera and genomes assembled from type species of genus."""

        check_file_exists(args.manual_taxonomy)
        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.uba_genome_paths)
        check_file_exists(args.qc_passed_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.untrustworthy_type_file)
        check_file_exists(args.synonym_file)
        check_file_exists(args.gtdb_type_strains_ledger)
        check_file_exists(args.sp_priority_ledger)
        check_file_exists(args.genus_priority_ledger)
        check_file_exists(args.dsmz_bacnames_file)
        make_sure_path_exists(args.output_dir)

        p = PMC_CheckTypeSpecies(args.output_dir)
        p.run(args.manual_taxonomy,
                args.cur_gtdb_metadata_file,
                args.uba_genome_paths,
                args.qc_passed_file,
                args.ncbi_genbank_assembly_file,
                args.untrustworthy_type_file,
                args.synonym_file,
                args.gtdb_type_strains_ledger,
                args.sp_priority_ledger,
                args.genus_priority_ledger,
                args.dsmz_bacnames_file)
        
        self.logger.info('Done.')
        
    def pmc_check_type_strains(self, args):
        """Check for agreement between GTDB species and genomes assembled from type strain of species."""

        check_file_exists(args.manual_taxonomy)
        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.uba_genome_paths)
        check_file_exists(args.qc_passed_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.untrustworthy_type_file)
        check_file_exists(args.synonym_file)
        check_file_exists(args.gtdb_type_strains_ledger)
        check_file_exists(args.sp_priority_ledger)
        check_file_exists(args.genus_priority_ledger)
        check_file_exists(args.dsmz_bacnames_file)
        make_sure_path_exists(args.output_dir)

        p = PMC_CheckTypeStrains(args.output_dir)
        p.run(args.manual_taxonomy,
                args.cur_gtdb_metadata_file,
                args.uba_genome_paths,
                args.qc_passed_file,
                args.ncbi_genbank_assembly_file,
                args.untrustworthy_type_file,
                args.synonym_file,
                args.gtdb_type_strains_ledger,
                args.sp_priority_ledger,
                args.genus_priority_ledger,
                args.dsmz_bacnames_file)
        
        self.logger.info('Done.')
        
    def pmc_species_names(self, args):
        """Establish final species names based on manual curation."""
        
        check_file_exists(args.manual_taxonomy)
        check_file_exists(args.manual_sp_names)
        check_file_exists(args.pmc_custom_species)
        check_file_exists(args.gtdb_clusters_file)
        check_file_exists(args.prev_gtdb_metadata_file)
        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.uba_genome_paths)
        check_file_exists(args.qc_passed_file)
        check_file_exists(args.updated_species_reps_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.untrustworthy_type_file)
        check_file_exists(args.synonym_file)
        check_file_exists(args.gtdb_type_strains_ledger)
        check_file_exists(args.sp_priority_ledger)
        check_file_exists(args.genus_priority_ledger)
        check_file_exists(args.specific_epithet_ledger)
        check_file_exists(args.dsmz_bacnames_file)
        make_sure_path_exists(args.output_dir)

        p = PMC_SpeciesNames(args.output_dir)
        p.run(args.manual_taxonomy,
                args.manual_sp_names,
                args.pmc_custom_species,
                args.gtdb_clusters_file,
                args.prev_gtdb_metadata_file,
                args.cur_gtdb_metadata_file,
                args.uba_genome_paths,
                args.qc_passed_file,
                args.updated_species_reps_file,
                args.ncbi_genbank_assembly_file,
                args.untrustworthy_type_file,
                args.synonym_file,
                args.gtdb_type_strains_ledger,
                args.sp_priority_ledger,
                args.genus_priority_ledger,
                args.specific_epithet_ledger,
                args.dsmz_bacnames_file)
        
        self.logger.info('Done.')
        
    def pmc_validate(self, args):
        """Validate final species names."""
        
        check_file_exists(args.final_taxonomy)
        check_file_exists(args.final_scaled_tree)
        check_file_exists(args.manual_sp_names)
        check_file_exists(args.pmc_custom_species)
        check_file_exists(args.gtdb_clusters_file)
        check_file_exists(args.prev_gtdb_metadata_file)
        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.uba_genome_paths)
        check_file_exists(args.qc_passed_file)
        check_file_exists(args.updated_species_reps_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.untrustworthy_type_file)
        check_file_exists(args.synonym_file)
        check_file_exists(args.gtdb_type_strains_ledger)
        check_file_exists(args.sp_priority_ledger)
        check_file_exists(args.genus_priority_ledger)
        check_file_exists(args.specific_epithet_ledger)
        check_file_exists(args.dsmz_bacnames_file)
        check_file_exists(args.ground_truth_test_cases)
        make_sure_path_exists(args.output_dir)

        p = PMC_Validation(args.output_dir)
        p.run(args.final_taxonomy,
                args.final_scaled_tree,
                args.manual_sp_names,
                args.pmc_custom_species,
                args.gtdb_clusters_file,
                args.prev_gtdb_metadata_file,
                args.cur_gtdb_metadata_file,
                args.uba_genome_paths,
                args.qc_passed_file,
                args.updated_species_reps_file,
                args.ncbi_genbank_assembly_file,
                args.untrustworthy_type_file,
                args.synonym_file,
                args.gtdb_type_strains_ledger,
                args.sp_priority_ledger,
                args.genus_priority_ledger,
                args.specific_epithet_ledger,
                args.dsmz_bacnames_file,
                args.ground_truth_test_cases)
        
        self.logger.info('Done.')
        
    def u_summary_stats(self, args):
        """Summary statistics indicating changes to GTDB species clusters."""
        
        check_file_exists(args.updated_sp_rep_file)
        check_file_exists(args.gtdb_clusters_file)
        check_file_exists(args.prev_gtdb_metadata_file)
        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.uba_genome_paths)
        check_file_exists(args.qc_passed_file)
        check_file_exists(args.gtdbtk_classify_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.untrustworthy_type_file)
        check_file_exists(args.synonym_file)
        check_file_exists(args.gtdb_type_strains_ledger)
        make_sure_path_exists(args.output_dir)

        p = UpdateSummaryStats(args.output_dir)
        p.run(args.updated_sp_rep_file,
                args.gtdb_clusters_file,
                args.prev_gtdb_metadata_file,
                args.cur_gtdb_metadata_file,
                args.uba_genome_paths,
                args.qc_passed_file,
                args.gtdbtk_classify_file,
                args.ncbi_genbank_assembly_file,
                args.untrustworthy_type_file,
                args.synonym_file,
                args.gtdb_type_strains_ledger)
        
        self.logger.info('Done.')
        
    def merge_test(self, args):
        """Produce information relevant to merging two sister species."""
        
        check_file_exists(args.gtdb_metadata_file)
        check_file_exists(args.genome_path_file)
        
        make_sure_path_exists(args.output_dir)
        
        p = MergeTest(args.ani_cache_file, args.cpus, args.output_dir)
        p.run(args.gtdb_metadata_file,
                args.genome_path_file,
                args.species1,
                args.species2)
        
        self.logger.info('Done.')
        
    def intra_sp_derep(self, args):
        """Dereplicate GTDB species clusters using ANI/AF criteria."""
        
        check_file_exists(args.gtdb_clusters_file)
        check_file_exists(args.gtdb_metadata_file)
        check_file_exists(args.genomic_path_file)
        check_file_exists(args.uba_gid_table)
        
        make_sure_path_exists(args.output_dir)
        
        p = IntraSpeciesDereplication(args.derep_ani,
                                        args.derep_af,
                                        args.max_genomes_per_sp,
                                        args.ani_cache_file, 
                                        args.cpus, 
                                        args.output_dir)
        p.run(args.gtdb_clusters_file,
                args.gtdb_metadata_file,
                args.genomic_path_file,
                args.uba_gid_table)
        
        self.logger.info('Done.')

    def rep_compare(self, args):
        """Compare current and previous representatives."""

        check_file_exists(args.cur_metadata_file)
        check_file_exists(args.prev_metadata_file)
        
        # get representatives in current taxonomy
        cur_gids = set()
        cur_species = set()
        cur_genera = set()
        cur_reps_taxa = {}
        cur_rep_species = set()
        cur_rep_genera = set()
        header = True
        for row in csv.reader(open(args.cur_metadata_file)):
            if header:
                header = False
                gtdb_rep_index = row.index('gtdb_representative')
                gtdb_taxonomy_index = row.index('gtdb_taxonomy')
            else:
                gid = row[0]
                cur_gids.add(gid)
                
                gtdb_taxonomy = row[gtdb_taxonomy_index]
                if gtdb_taxonomy:
                    gtdb_taxa = [t.strip() for t in row[gtdb_taxonomy_index].split(';')]
                    if gtdb_taxa[6] != 's__':
                        cur_species.add(gtdb_taxa[6])
                    if gtdb_taxa[5] != 'g__':
                        cur_genera.add(gtdb_taxa[5])

                if row[gtdb_rep_index] == 't':
                    cur_reps_taxa[gid] = gtdb_taxa
                    
                    if gtdb_taxa[6] != 's__':
                        cur_rep_species.add(gtdb_taxa[6])
                        
                    if gtdb_taxa[5] != 'g__':
                        cur_rep_genera.add(gtdb_taxa[5])
                    
        # get representatives in previous taxonomy
        prev_reps_taxa = {}
        prev_rep_species = set()
        prev_rep_genera = set()
        header = True
        for row in csv.reader(open(args.prev_metadata_file)):
            if header:
                header = False
                gtdb_rep_index = row.index('gtdb_representative')
                gtdb_taxonomy_index = row.index('gtdb_taxonomy')
            else:
                if row[gtdb_rep_index] == 't':
                    gid = row[0]
                    gtdb_taxonomy = row[gtdb_taxonomy_index]
                    if gtdb_taxonomy:
                        gtdb_taxa = [t.strip() for t in row[gtdb_taxonomy_index].split(';')]
                            
                        prev_reps_taxa[gid] = gtdb_taxa
                        
                        if gtdb_taxa[6] != 's__':
                            prev_rep_species.add(gtdb_taxa[6])
                            
                        if gtdb_taxa[5] != 'g__':
                            prev_rep_genera.add(gtdb_taxa[5])
                    
        # summarize differences
        print('No. current representatives: %d' % len(cur_reps_taxa))
        print('No. previous representatives: %d' % len(prev_reps_taxa))
        
        print('')
        print('No. current species with representatives: %d' % len(cur_rep_species))
        print('No. previous species with representatives: %d' % len(prev_rep_species))
        
        print('')
        print('No. new representatives: %d' % len(set(cur_reps_taxa) - set(prev_reps_taxa)))
        print('No. retired representatives: %d' % len(set(prev_reps_taxa) - set(cur_reps_taxa)))
        
        print('')
        print('No. new species with representative: %d' % len(cur_rep_species - prev_rep_species))
        print('No. new genera with representative: %d' % len(cur_rep_genera - prev_rep_genera))
        
        print('')
        missing_sp_reps = prev_rep_species.intersection(cur_species) - cur_rep_species
        print('No. species that no longer have a representative: %d' % len(missing_sp_reps))
        for sp in missing_sp_reps:
            print('  ' + sp)
        
        print('')
        missing_genera_reps = prev_rep_genera.intersection(cur_genera) - cur_rep_genera
        print('No. genera that no longer have a representative: %d' % len(missing_genera_reps))
        for g in missing_genera_reps:
            print('  ' + g)
        
        print('')
        deprecated_reps = set(prev_reps_taxa).intersection(cur_gids) - set(cur_reps_taxa)
        print('No. deprecated previous representatives: %d' % len(deprecated_reps))
        
    def cluster_stats(self, args):
        """Calculate statistics for species cluster."""

        check_file_exists(args.cluster_file)
        check_file_exists(args.genome_path_file)
        check_file_exists(args.gtdb_metadata_file)
        
        p = ClusterStats(args.af_sp,
                            args.max_genomes,
                            args.ani_cache_file,
                            args.cpus, 
                            args.output_dir)
        p.run(args.cluster_file, 
                args.genome_path_file,
                args.gtdb_metadata_file)

    def run(self, args):
        """Parse user arguments and call the correct pipeline(s)"""

        logging.basicConfig(format='', level=logging.INFO)

        if args.subparser_name == 'qc_genomes':
            self.qc_genomes(args)
        elif args.subparser_name == 'mash_dist':
            self.mash_dist(args)
        elif args.subparser_name == 'select_type_genomes':
            self.select_type_genomes(args)
        elif args.subparser_name == 'cluster_named_types':
            self.cluster_named_types(args)
        elif args.subparser_name == 'cluster_de_novo':
            self.cluster_de_novo(args)
        elif args.subparser_name == 'cluster_user':
            self.cluster_user(args)
        elif args.subparser_name == 'tree_gids':
            self.tree_gids(args)
        elif args.subparser_name == 'u_new_genomes':
            self.u_new_genomes(args)
        elif args.subparser_name == 'u_qc_genomes':
            self.u_qc_genomes(args)
        elif args.subparser_name == 'u_resolve_types':
            self.u_resolve_types(args)
        elif args.subparser_name == 'u_gtdbtk':
            self.u_gtdbtk(args)
        elif args.subparser_name == 'u_rep_changes':
            self.u_rep_changes(args)
        elif args.subparser_name == 'u_rep_actions':
            self.u_rep_actions(args)
        elif args.subparser_name == 'u_sel_reps':
            self.u_sel_reps(args)
        elif args.subparser_name == 'u_cluster_named_reps':
            self.u_cluster_named_reps(args)
        elif args.subparser_name == 'u_synonyms':
            self.u_synonyms(args)
        elif args.subparser_name == 'u_cluster_de_novo':
            self.u_cluster_de_novo(args)
        elif args.subparser_name == 'u_genus_names':
            self.u_genus_names(args)
        elif args.subparser_name == 'u_curation_trees':
            self.u_curation_trees(args)
        elif args.subparser_name == 'u_species_init':
            self.u_species_init(args)
        elif args.subparser_name == 'pmc_manual_species':
            self.pmc_manual_species(args)
        elif args.subparser_name == 'pmc_replace_generic':
            self.pmc_replace_generic(args)
        elif args.subparser_name == 'pmc_check_type_species':
            self.pmc_check_type_species(args)
        elif args.subparser_name == 'pmc_check_type_strains':
            self.pmc_check_type_strains(args)
        elif args.subparser_name == 'pmc_species_names':
            self.pmc_species_names(args)
        elif args.subparser_name == 'pmc_validate':
            self.pmc_validate(args)
        elif args.subparser_name == 'u_summary_stats':
            self.u_summary_stats(args)
        elif args.subparser_name == 'merge_test':
            self.merge_test(args)
        elif args.subparser_name == 'intra_sp_derep':
            self.intra_sp_derep(args)
        elif args.subparser_name == 'rep_compare':
            self.rep_compare(args)
        elif args.subparser_name == 'cluster_stats':
            self.cluster_stats(args)
        else:
            self.logger.error('Unknown gtdb_species_clusters command: ' + args.subparser_name + '\n')
            sys.exit()

        return 0
