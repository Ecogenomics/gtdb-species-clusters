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

from gtdblib.util.shell.gtdbshutil import check_file_exists, make_sure_path_exists

from gtdb_species_clusters.update_new_genomes import NewGenomes
from gtdb_species_clusters.update_qc_genomes import QcGenomes, QcCriteria
from gtdb_species_clusters.update_lpsn_ssu_types import LPSN_SSU_Types
from gtdb_species_clusters.update_resolve_types import ResolveTypes
from gtdb_species_clusters.update_gtdbtk import GTDB_Tk
from gtdb_species_clusters.update_rep_changes import RepChanges
from gtdb_species_clusters.update_rep_actions import RepActions
from gtdb_species_clusters.update_select_reps import UpdateSelectRepresentatives
from gtdb_species_clusters.update_cluster_named_reps import UpdateClusterNamedReps
from gtdb_species_clusters.update_erroneous_ncbi import UpdateErroneousNCBI
from gtdb_species_clusters.update_synonyms import UpdateSynonyms
from gtdb_species_clusters.update_cluster_de_novo import UpdateClusterDeNovo
from gtdb_species_clusters.update_cluster_stats import UpdateClusterStats
from gtdb_species_clusters.update_curation_trees import UpdateCurationTrees
from gtdb_species_clusters.update_species_init import UpdateSpeciesInit

from gtdb_species_clusters.pmc_checks import PMC_Checks
from gtdb_species_clusters.pmc_check_type_species import PMC_CheckTypeSpecies
from gtdb_species_clusters.pmc_check_type_strains import PMC_CheckTypeStrains
from gtdb_species_clusters.pmc_species_names import PMC_SpeciesNames
from gtdb_species_clusters.pmc_validation import PMC_Validation
from gtdb_species_clusters.pmc_cluster_stats import PMC_ClusterStats
from gtdb_species_clusters.pmc_cleanup import PMC_Cleanup

from gtdb_species_clusters.merge_test import MergeTest
from gtdb_species_clusters.ani_sp_pair import ANI_SpeciesPair
from gtdb_species_clusters.intra_sp_derep import IntraSpeciesDereplication
from gtdb_species_clusters.intra_genus_ani import IntraGenusANI
from gtdb_species_clusters.rep_genomic_similarity import RepGenomicSimilarity

from gtdb_species_clusters.inspect_genomes import InspectGenomes

from gtdb_species_clusters.sandbox import Sandbox


if sys.version_info >= (3, 11):
    import tomllib

    def parse_toml_file(toml_file: str):
        with open(toml_file, 'rb') as f:
            return tomllib.load(f)
else:
    import toml as tomllib

    def parse_toml_file(toml_file: str):
        with open(toml_file, 'r') as f:
            return tomllib.load(f)


class OptionsParser():
    """Validate and execute command-line interface."""

    def __init__(self):
        """Initialization"""
        self.log = logging.getLogger()

    def u_new_genomes(self, args):
        """Identify new and updated genomes."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        check_file_exists(params['data_files']['prev_gtdb_metadata_file'])
        check_file_exists(params['data_files']['cur_gtdb_metadata_file'])
        check_file_exists(params['data_files']['cur_genome_path_file'])
        check_file_exists(params['data_files']['ncbi_assembly_summary_genbank_file'])

        make_sure_path_exists(params['output_dirs']['u_new_genomes'])
        self.log.info(f"Output directory: {params['output_dirs']['u_new_genomes']}")

        p = NewGenomes(params['output_dirs']['u_new_genomes'])
        p.run(params['data_files']['prev_gtdb_metadata_file'],
              params['data_files']['cur_gtdb_metadata_file'],
              params['data_files']['cur_genome_path_file'],
              params['data_files']['ncbi_assembly_summary_genbank_file'])

    def u_qc_genomes(self, args):
        """Quality check new and updated genomes."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        check_file_exists(params['data_files']['prev_gtdb_metadata_file'])
        check_file_exists(params['data_files']['cur_gtdb_metadata_file'])
        check_file_exists(params['data_files']['ncbi_assembly_summary_genbank_file'])
        check_file_exists(params['data_files']['cur_gtdb_domain_report'])
        check_file_exists(params['ledgers']['gtdb_type_strains_ledger'])
        check_file_exists(params['ledgers']['qc_exception_ledger'])
        check_file_exists(params['ledgers']['ncbi_env_bioproject_ledger'])

        make_sure_path_exists(params['output_dirs']['u_qc_genomes'])
        self.log.info(f"Output directory: {params['output_dirs']['u_qc_genomes']}")

        qc_criteria = QcCriteria(
            args.min_comp,
            args.max_cont,
            args.min_quality,
            args.sh_exception,
            args.min_perc_markers,
            args.max_contigs,
            args.min_N50,
            args.max_ambiguous,
            args.max_cm_xor_contigs)

        p = QcGenomes(params['output_dirs']['u_qc_genomes'])
        p.run(params['data_files']['prev_gtdb_metadata_file'],
              params['data_files']['cur_gtdb_metadata_file'],
              params['data_files']['ncbi_assembly_summary_genbank_file'],
              params['data_files']['cur_gtdb_domain_report'],
              params['ledgers']['gtdb_type_strains_ledger'],
              params['ledgers']['qc_exception_ledger'],
              params['ledgers']['ncbi_env_bioproject_ledger'],
              qc_criteria)

    def u_gtdbtk(self, args):
        """Perform initial classification of new and updated genomes using GTDB-Tk."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        genomes_new_updated = os.path.join(params['output_dirs']['u_new_genomes'], 'genomes_new_updated.tsv')
        qc_passed = os.path.join(params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')

        check_file_exists(genomes_new_updated)
        check_file_exists(qc_passed)

        make_sure_path_exists(params['output_dirs']['u_gtdbtk'])
        self.log.info(f"Output directory: {params['output_dirs']['u_gtdbtk']}")

        p = GTDB_Tk(args.cpus, params['output_dirs']['u_gtdbtk'])
        p.run(genomes_new_updated,
              qc_passed,
              params['data_files']['mash_db_file'],
              args.batch_size)

    def u_lpsn_rna_types(self, args):
        """Identify type genomes based on type 16S rRNA sequences indicated at LPSN."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        qc_passed = os.path.join(params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')

        check_file_exists(params['data_files']['lpsn_metadata_file'])
        check_file_exists(params['data_files']['cur_gtdb_metadata_file'])
        check_file_exists(params['data_files']['cur_genome_path_file'])
        check_file_exists(qc_passed)
        check_file_exists(params['data_files']['ncbi_assembly_summary_genbank_file'])
        check_file_exists(params['ledgers']['gtdb_type_strains_ledger'])
        check_file_exists(params['ledgers']['untrustworthy_type_ledger'])

        make_sure_path_exists(params['output_dirs']['u_lpsn_rna_types'])
        self.log.info(f"Output directory: {params['output_dirs']['u_lpsn_rna_types']}")

        p = LPSN_SSU_Types(args.cpus, params['output_dirs']['u_lpsn_rna_types'])
        p.run(params['data_files']['lpsn_metadata_file'],
              params['data_files']['cur_gtdb_metadata_file'],
              params['data_files']['cur_genome_path_file'],
              qc_passed,
              params['data_files']['ncbi_assembly_summary_genbank_file'],
              params['ledgers']['gtdb_type_strains_ledger'],
              params['ledgers']['untrustworthy_type_ledger'])

    def u_resolve_types(self, args):
        """Resolve cases where a species has multiple genomes assembled from the type strain."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        qc_passed = os.path.join(params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')

        check_file_exists(qc_passed)
        check_file_exists(params['data_files']['cur_gtdb_metadata_file'])
        check_file_exists(params['data_files']['cur_genome_path_file'])
        check_file_exists(params['data_files']['ncbi_assembly_summary_genbank_file'])
        check_file_exists(params['ltp']['ltp_taxonomy_file'])
        check_file_exists(params['ledgers']['gtdb_type_strains_ledger'])
        check_file_exists(params['ledgers']['untrustworthy_type_ledger'])
        check_file_exists(params['ledgers']['ncbi_env_bioproject_ledger'])

        make_sure_path_exists(params['output_dirs']['u_resolve_types'])
        self.log.info(f"Output directory: {params['output_dirs']['u_resolve_types']}")

        p = ResolveTypes(args.cpus, params['ltp']['ltp_dir'], params['output_dirs']['u_resolve_types'])
        p.run(params['data_files']['cur_gtdb_metadata_file'],
              params['data_files']['cur_genome_path_file'],
              qc_passed,
              params['data_files']['ncbi_assembly_summary_genbank_file'],
              params['ltp']['ltp_taxonomy_file'],
              params['ledgers']['gtdb_type_strains_ledger'],
              params['ledgers']['untrustworthy_type_ledger'],
              params['ledgers']['ncbi_env_bioproject_ledger'])

    def u_rep_changes(self, args):
        """Identify species representatives that have changed from previous release."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        genomes_new_updated = os.path.join(params['output_dirs']['u_new_genomes'], 'genomes_new_updated.tsv')
        qc_passed = os.path.join(params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')
        gtdbtk_classify_file = os.path.join(params['output_dirs']['u_gtdbtk'], 'gtdbtk_classify.tsv')
        untrustworthy_type_file = os.path.join(params['output_dirs']['u_resolve_types'], 'untrustworthy_type_material.tsv')

        check_file_exists(params['data_files']['prev_gtdb_metadata_file'])
        check_file_exists(params['data_files']['cur_gtdb_metadata_file'])
        check_file_exists(genomes_new_updated)
        check_file_exists(qc_passed)
        check_file_exists(gtdbtk_classify_file)
        check_file_exists(params['data_files']['ncbi_assembly_summary_genbank_file'])
        check_file_exists(untrustworthy_type_file)
        check_file_exists(params['ledgers']['disband_cluster_ledger'])
        check_file_exists(params['ledgers']['gtdb_type_strains_ledger'])
        check_file_exists(params['ledgers']['ncbi_env_bioproject_ledger'])

        make_sure_path_exists(params['output_dirs']['u_rep_changes'])
        self.log.info(f"Output directory: {params['output_dirs']['u_rep_changes']}")

        p = RepChanges(params['output_dirs']['u_rep_changes'])
        p.run(params['data_files']['prev_gtdb_metadata_file'],
              params['data_files']['cur_gtdb_metadata_file'],
              genomes_new_updated,
              qc_passed,
              gtdbtk_classify_file,
              params['data_files']['ncbi_assembly_summary_genbank_file'],
              untrustworthy_type_file,
              params['ledgers']['disband_cluster_ledger'],
              params['ledgers']['gtdb_type_strains_ledger'],
              params['ledgers']['ncbi_env_bioproject_ledger'])

    def u_rep_actions(self, args):
        """Perform initial actions required for changed representatives."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        genomes_new_updated = os.path.join(params['output_dirs']['u_new_genomes'], 'genomes_new_updated.tsv')
        qc_passed = os.path.join(params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')
        gtdbtk_classify_file = os.path.join(params['output_dirs']['u_gtdbtk'], 'gtdbtk_classify.tsv')
        untrustworthy_type_file = os.path.join(params['output_dirs']['u_resolve_types'], 'untrustworthy_type_material.tsv')
        rep_change_summary_file = os.path.join(params['output_dirs']['u_rep_changes'], 'rep_change_summary.tsv')

        check_file_exists(rep_change_summary_file)
        check_file_exists(params['data_files']['prev_genome_path_file'])
        check_file_exists(params['data_files']['cur_genome_path_file'])
        check_file_exists(genomes_new_updated)
        check_file_exists(qc_passed)
        check_file_exists(gtdbtk_classify_file)
        check_file_exists(params['data_files']['ncbi_assembly_summary_genbank_file'])
        check_file_exists(untrustworthy_type_file)
        check_file_exists(params['ledgers']['gtdb_type_strains_ledger'])
        check_file_exists(params['ledgers']['sp_priority_ledger'])
        check_file_exists(params['ledgers']['genus_priority_ledger'])
        check_file_exists(params['ledgers']['ncbi_env_bioproject_ledger'])
        check_file_exists(params['data_files']['lpsn_gss_file'])

        make_sure_path_exists(params['output_dirs']['u_rep_actions'])
        self.log.info(f"Output directory: {params['output_dirs']['u_rep_actions']}")

        p = RepActions(args.cpus, params['output_dirs']['u_rep_actions'])
        p.run(rep_change_summary_file,
              params['data_files']['prev_gtdb_metadata_file'],
              params['data_files']['prev_genome_path_file'],
              params['data_files']['cur_gtdb_metadata_file'],
              params['data_files']['cur_genome_path_file'],
              genomes_new_updated,
              qc_passed,
              gtdbtk_classify_file,
              params['data_files']['ncbi_assembly_summary_genbank_file'],
              untrustworthy_type_file,
              params['ledgers']['gtdb_type_strains_ledger'],
              params['ledgers']['sp_priority_ledger'],
              params['ledgers']['genus_priority_ledger'],
              params['ledgers']['ncbi_env_bioproject_ledger'],
              params['data_files']['lpsn_gss_file'])

    def u_sel_reps(self, args):
        """Select representatives for all named species at NCBI."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        qc_passed = os.path.join(params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')
        untrustworthy_type_file = os.path.join(params['output_dirs']['u_resolve_types'], 'untrustworthy_type_material.tsv')
        updated_sp_cluster_file = os.path.join(params['output_dirs']['u_rep_actions'], 'updated_sp_clusters.tsv')

        check_file_exists(updated_sp_cluster_file)
        check_file_exists(params['data_files']['cur_gtdb_metadata_file'])
        check_file_exists(params['data_files']['cur_genome_path_file'])
        check_file_exists(qc_passed)
        check_file_exists(params['data_files']['ncbi_assembly_summary_genbank_file'])
        check_file_exists(untrustworthy_type_file)
        check_file_exists(params['ledgers']['gtdb_type_strains_ledger'])
        check_file_exists(params['ledgers']['ncbi_untrustworthy_sp_ledger'])
        check_file_exists(params['ledgers']['sp_priority_ledger'])
        check_file_exists(params['ledgers']['genus_priority_ledger'])
        check_file_exists(params['ledgers']['ncbi_env_bioproject_ledger'])
        check_file_exists(params['data_files']['lpsn_gss_file'])

        make_sure_path_exists(params['output_dirs']['u_sel_reps'])
        self.log.info(f"Output directory: {params['output_dirs']['u_sel_reps']}")

        p = UpdateSelectRepresentatives(args.cpus,
                                        params['output_dirs']['u_sel_reps'])
        p.run(updated_sp_cluster_file,
              params['data_files']['cur_gtdb_metadata_file'],
              params['data_files']['cur_genome_path_file'],
              qc_passed,
              params['data_files']['ncbi_assembly_summary_genbank_file'],
              untrustworthy_type_file,
              params['ledgers']['gtdb_type_strains_ledger'],
              params['ledgers']['ncbi_untrustworthy_sp_ledger'],
              params['ledgers']['sp_priority_ledger'],
              params['ledgers']['genus_priority_ledger'],
              params['ledgers']['ncbi_env_bioproject_ledger'],
              params['data_files']['lpsn_gss_file'])

    def u_cluster_named_reps(self, args):
        """Cluster genomes to selected GTDB representatives."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        qc_passed = os.path.join(params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')
        untrustworthy_type_file = os.path.join(params['output_dirs']['u_resolve_types'], 'untrustworthy_type_material.tsv')
        named_rep_file = os.path.join(params['output_dirs']['u_sel_reps'], 'gtdb_named_reps_final.tsv')
        rep_ani_file = os.path.join(params['output_dirs']['u_sel_reps'], 'gtdb_rep_pairwise_ani.tsv')

        check_file_exists(named_rep_file)
        check_file_exists(params['data_files']['cur_gtdb_metadata_file'])
        check_file_exists(params['data_files']['cur_genome_path_file'])
        check_file_exists(qc_passed)
        check_file_exists(params['data_files']['ncbi_assembly_summary_genbank_file'])
        check_file_exists(untrustworthy_type_file)
        check_file_exists(rep_ani_file)
        check_file_exists(params['ledgers']['gtdb_type_strains_ledger'])
        check_file_exists(params['ledgers']['ncbi_env_bioproject_ledger'])

        make_sure_path_exists(params['output_dirs']['u_cluster_named_reps'])
        self.log.info(f"Output directory: {params['output_dirs']['u_cluster_named_reps']}")

        p = UpdateClusterNamedReps(args.ani_sp,
                                   args.af_sp,
                                   args.cpus,
                                   params['output_dirs']['u_cluster_named_reps'])
        p.run(named_rep_file,
              params['data_files']['cur_gtdb_metadata_file'],
              params['data_files']['cur_genome_path_file'],
              qc_passed,
              params['data_files']['ncbi_assembly_summary_genbank_file'],
              untrustworthy_type_file,
              rep_ani_file,
              params['ledgers']['gtdb_type_strains_ledger'],
              params['ledgers']['ncbi_env_bioproject_ledger'])

    def u_cluster_de_novo(self, args):
        """Infer de novo species clusters and representatives for remaining genomes."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        qc_passed = os.path.join(params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')
        untrustworthy_type_file = os.path.join(params['output_dirs']['u_resolve_types'], 'untrustworthy_type_material.tsv')
        named_cluster_file = os.path.join(params['output_dirs']['u_cluster_named_reps'], 'gtdb_named_rep_clusters.tsv')
        ani_af_nonrep_vs_rep = os.path.join(params['output_dirs']['u_cluster_named_reps'], 'ani_af_nonrep_vs_rep.pkl')

        check_file_exists(named_cluster_file)
        check_file_exists(params['data_files']['cur_gtdb_metadata_file'])
        check_file_exists(params['data_files']['cur_genome_path_file'])
        check_file_exists(qc_passed)
        check_file_exists(params['data_files']['ncbi_assembly_summary_genbank_file'])
        check_file_exists(untrustworthy_type_file)
        check_file_exists(ani_af_nonrep_vs_rep)
        check_file_exists(params['ledgers']['gtdb_type_strains_ledger'])
        check_file_exists(params['ledgers']['ncbi_env_bioproject_ledger'])

        make_sure_path_exists(params['output_dirs']['u_cluster_de_novo'])
        self.log.info(f"Output directory: {params['output_dirs']['u_cluster_de_novo']}")

        p = UpdateClusterDeNovo(args.ani_sp,
                                args.af_sp,
                                args.cpus,
                                params['output_dirs']['u_cluster_de_novo'])
        p.run(named_cluster_file,
              params['data_files']['cur_gtdb_metadata_file'],
              params['data_files']['cur_genome_path_file'],
              qc_passed,
              params['data_files']['ncbi_assembly_summary_genbank_file'],
              untrustworthy_type_file,
              ani_af_nonrep_vs_rep,
              params['ledgers']['gtdb_type_strains_ledger'],
              params['ledgers']['ncbi_env_bioproject_ledger'])

    def u_cluster_stats(self, args):
        """Summary statistics indicating changes to GTDB species cluster membership."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        qc_passed = os.path.join(params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')
        untrustworthy_type_file = os.path.join(params['output_dirs']['u_resolve_types'], 'untrustworthy_type_material.tsv')
        gtdb_clusters_file = os.path.join(params['output_dirs']['u_cluster_de_novo'], 'gtdb_clusters_de_novo.tsv')

        check_file_exists(gtdb_clusters_file)
        check_file_exists(params['data_files']['prev_gtdb_metadata_file'])
        check_file_exists(params['data_files']['cur_gtdb_metadata_file'])
        check_file_exists(qc_passed)
        check_file_exists(params['data_files']['ncbi_assembly_summary_genbank_file'])
        check_file_exists(untrustworthy_type_file)
        check_file_exists(params['ledgers']['gtdb_type_strains_ledger'])
        check_file_exists(params['ledgers']['ncbi_env_bioproject_ledger'])

        make_sure_path_exists(params['output_dirs']['u_cluster_stats'])
        self.log.info(f"Output directory: {params['output_dirs']['u_cluster_stats']}")

        p = UpdateClusterStats(params['output_dirs']['u_cluster_stats'])
        p.run(gtdb_clusters_file,
              params['data_files']['prev_gtdb_metadata_file'],
              params['data_files']['cur_gtdb_metadata_file'],
              qc_passed,
              params['data_files']['ncbi_assembly_summary_genbank_file'],
              untrustworthy_type_file,
              params['ledgers']['gtdb_type_strains_ledger'],
              params['ledgers']['ncbi_env_bioproject_ledger'])

    def u_ncbi_erroneous(self, args):
        """Identify genomes with erroneous NCBI species assignments."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        qc_passed = os.path.join(params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')
        untrustworthy_type_file = os.path.join(params['output_dirs']['u_resolve_types'], 'untrustworthy_type_material.tsv')
        gtdb_clusters_file = os.path.join(params['output_dirs']['u_cluster_de_novo'], 'gtdb_clusters_de_novo.tsv')

        check_file_exists(gtdb_clusters_file)
        check_file_exists(params['data_files']['cur_gtdb_metadata_file'])
        check_file_exists(params['data_files']['cur_genome_path_file'])
        check_file_exists(qc_passed)
        check_file_exists(params['data_files']['ncbi_assembly_summary_genbank_file'])
        check_file_exists(untrustworthy_type_file)
        check_file_exists(params['ledgers']['gtdb_type_strains_ledger'])
        check_file_exists(params['ledgers']['ncbi_untrustworthy_sp_ledger'])
        check_file_exists(params['ledgers']['ncbi_env_bioproject_ledger'])

        make_sure_path_exists(params['output_dirs']['u_ncbi_erroneous'])
        self.log.info(f"Output directory: {params['output_dirs']['u_ncbi_erroneous']}")

        p = UpdateErroneousNCBI(params['output_dirs']['u_ncbi_erroneous'])
        p.run(gtdb_clusters_file,
              params['data_files']['cur_gtdb_metadata_file'],
              params['data_files']['cur_genome_path_file'],
              qc_passed,
              params['data_files']['ncbi_assembly_summary_genbank_file'],
              untrustworthy_type_file,
              params['ledgers']['gtdb_type_strains_ledger'],
              params['ledgers']['ncbi_untrustworthy_sp_ledger'],
              params['ledgers']['ncbi_env_bioproject_ledger'])

    def u_synonyms(self, args):
        """Determine synonyms for validly or effectively published species."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        qc_passed = os.path.join(params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')
        untrustworthy_type_file = os.path.join(params['output_dirs']['u_resolve_types'], 'untrustworthy_type_material.tsv')
        gtdb_clusters_file = os.path.join(params['output_dirs']['u_cluster_de_novo'], 'gtdb_clusters_de_novo.tsv')
        ncbi_misclassified_file = os.path.join(
            params['output_dirs']['u_ncbi_erroneous'], 'ncbi_misclassified_sp.gtdb_clustering.tsv')
        ani_af_nonrep_vs_rep = os.path.join(params['output_dirs']['u_cluster_named_reps'], 'ani_af_nonrep_vs_rep.pkl')

        check_file_exists(gtdb_clusters_file)
        check_file_exists(params['data_files']['cur_gtdb_metadata_file'])
        check_file_exists(qc_passed)
        check_file_exists(ncbi_misclassified_file)
        check_file_exists(params['data_files']['ncbi_assembly_summary_genbank_file'])
        check_file_exists(untrustworthy_type_file)
        check_file_exists(ani_af_nonrep_vs_rep)
        check_file_exists(params['ledgers']['gtdb_type_strains_ledger'])
        check_file_exists(params['ledgers']['sp_priority_ledger'])
        check_file_exists(params['ledgers']['genus_priority_ledger'])
        check_file_exists(params['ledgers']['ncbi_untrustworthy_sp_ledger'])
        check_file_exists(params['ledgers']['ncbi_env_bioproject_ledger'])
        check_file_exists(params['data_files']['lpsn_gss_file'])

        make_sure_path_exists(params['output_dirs']['u_synonyms'])
        self.log.info(f"Output directory: {params['output_dirs']['u_synonyms']}")

        p = UpdateSynonyms(params['output_dirs']['u_synonyms'])
        p.run(gtdb_clusters_file,
              params['data_files']['cur_gtdb_metadata_file'],
              qc_passed,
              ncbi_misclassified_file,
              params['data_files']['ncbi_assembly_summary_genbank_file'],
              untrustworthy_type_file,
              ani_af_nonrep_vs_rep,
              params['ledgers']['gtdb_type_strains_ledger'],
              params['ledgers']['sp_priority_ledger'],
              params['ledgers']['genus_priority_ledger'],
              params['ledgers']['ncbi_untrustworthy_sp_ledger'],
              params['ledgers']['ncbi_env_bioproject_ledger'],
              params['data_files']['lpsn_gss_file'])

    def u_curation_trees(self, args):
        """Produce curation trees highlighting new NCBI taxa."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        qc_passed = os.path.join(params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')
        untrustworthy_type_file = os.path.join(params['output_dirs']['u_resolve_types'], 'untrustworthy_type_material.tsv')
        gtdb_clusters_file = os.path.join(params['output_dirs']['u_cluster_de_novo'], 'gtdb_clusters_de_novo.tsv')

        check_file_exists(gtdb_clusters_file)
        check_file_exists(params['data_files']['prev_gtdb_metadata_file'])
        check_file_exists(params['data_files']['cur_gtdb_metadata_file'])
        check_file_exists(qc_passed)
        check_file_exists(params['data_files']['ncbi_assembly_summary_genbank_file'])
        check_file_exists(untrustworthy_type_file)
        check_file_exists(params['ledgers']['gtdb_type_strains_ledger'])
        check_file_exists(params['ledgers']['ncbi_untrustworthy_sp_ledger'])
        check_file_exists(params['ledgers']['ncbi_env_bioproject_ledger'])

        make_sure_path_exists(params['output_dirs']['u_curation_trees'])
        self.log.info(f"Output directory: {params['output_dirs']['u_curation_trees']}")

        p = UpdateCurationTrees(params['output_dirs']['u_curation_trees'], f"gtdb_{params['release']}")
        p.run(gtdb_clusters_file,
              params['data_files']['prev_gtdb_metadata_file'],
              params['data_files']['cur_gtdb_metadata_file'],
              qc_passed,
              params['data_files']['ncbi_assembly_summary_genbank_file'],
              untrustworthy_type_file,
              params['ledgers']['gtdb_type_strains_ledger'],
              params['ledgers']['ncbi_untrustworthy_sp_ledger'],
              params['ledgers']['ncbi_env_bioproject_ledger'])

    def u_species_init(self, args):
        """Produce initial best guess at GTDB species clusters."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        qc_passed = os.path.join(params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')
        untrustworthy_type_file = os.path.join(params['output_dirs']['u_resolve_types'], 'untrustworthy_type_material.tsv')
        gtdb_clusters_file = os.path.join(params['output_dirs']['u_cluster_de_novo'], 'gtdb_clusters_de_novo.tsv')
        gtdbtk_classify_file = os.path.join(params['output_dirs']['u_gtdbtk'], 'gtdbtk_classify.tsv')
        synonym_file = os.path.join(params['output_dirs']['u_synonyms'], 'synonyms.tsv')

        check_file_exists(gtdb_clusters_file)
        check_file_exists(params['data_files']['prev_gtdb_metadata_file'])
        check_file_exists(params['data_files']['cur_gtdb_metadata_file'])
        check_file_exists(qc_passed)
        check_file_exists(gtdbtk_classify_file)
        check_file_exists(params['data_files']['ncbi_assembly_summary_genbank_file'])
        check_file_exists(untrustworthy_type_file)
        check_file_exists(synonym_file)
        check_file_exists(params['ledgers']['gtdb_type_strains_ledger'])
        check_file_exists(params['ledgers']['sp_priority_ledger'])
        check_file_exists(params['ledgers']['genus_priority_ledger'])
        check_file_exists(params['ledgers']['ncbi_untrustworthy_sp_ledger'])
        check_file_exists(params['ledgers']['ncbi_env_bioproject_ledger'])
        check_file_exists(params['data_files']['lpsn_gss_file'])

        make_sure_path_exists(params['output_dirs']['u_species_init'])
        self.log.info(f"Output directory: {params['output_dirs']['u_species_init']}")

        p = UpdateSpeciesInit(params['output_dirs']['u_species_init'])
        p.run(gtdb_clusters_file,
              params['data_files']['prev_gtdb_metadata_file'],
              params['data_files']['cur_gtdb_metadata_file'],
              qc_passed,
              gtdbtk_classify_file,
              params['data_files']['ncbi_assembly_summary_genbank_file'],
              untrustworthy_type_file,
              synonym_file,
              params['ledgers']['gtdb_type_strains_ledger'],
              params['ledgers']['sp_priority_ledger'],
              params['ledgers']['genus_priority_ledger'],
              params['ledgers']['ncbi_untrustworthy_sp_ledger'],
              params['ledgers']['ncbi_env_bioproject_ledger'],
              params['data_files']['lpsn_gss_file'])

    def u_pmc_species_names(self, args):
        """Refine species names using post-manual curation rules."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        qc_passed = os.path.join(params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')
        untrustworthy_type_file = os.path.join(params['output_dirs']['u_resolve_types'], 'untrustworthy_type_material.tsv')
        gtdb_clusters_file = os.path.join(params['output_dirs']['u_cluster_de_novo'], 'gtdb_clusters_de_novo.tsv')
        synonym_file = os.path.join(params['output_dirs']['u_synonyms'], 'synonyms.tsv')
        ncbi_misclassified_file = os.path.join(
            params['output_dirs']['u_ncbi_erroneous'], 'ncbi_misclassified_sp.gtdb_clustering.tsv')
        updated_species_reps_file = os.path.join(params['output_dirs']['u_rep_actions'], 'updated_species_reps.tsv')
        taxonomy_init_ar = os.path.join(params['output_dirs']['u_species_init'], 'gtdb_ar_taxonomy.tsv')
        taxonomy_init_bac = os.path.join(params['output_dirs']['u_species_init'], 'gtdb_bac_taxonomy.tsv')

        check_file_exists(taxonomy_init_ar)
        check_file_exists(taxonomy_init_bac)
        check_file_exists(params['data_files']['seqcode_file'])
        check_file_exists(gtdb_clusters_file)
        check_file_exists(params['data_files']['prev_gtdb_metadata_file'])
        check_file_exists(params['data_files']['cur_gtdb_metadata_file'])
        check_file_exists(qc_passed)
        check_file_exists(ncbi_misclassified_file)
        check_file_exists(params['data_files']['ncbi_assembly_summary_genbank_file'])
        check_file_exists(untrustworthy_type_file)
        check_file_exists(synonym_file)
        check_file_exists(updated_species_reps_file)
        check_file_exists(params['ledgers']['gtdb_type_strains_ledger'])
        check_file_exists(params['ledgers']['species_classification_ledger'])
        check_file_exists(params['ledgers']['sp_priority_ledger'])
        check_file_exists(params['ledgers']['genus_priority_ledger'])
        check_file_exists(params['ledgers']['specific_epithet_ledger_ar'])
        check_file_exists(params['ledgers']['specific_epithet_ledger_bac'])
        check_file_exists(params['ledgers']['ncbi_untrustworthy_sp_ledger'])
        check_file_exists(params['ledgers']['ncbi_env_bioproject_ledger'])
        check_file_exists(params['data_files']['lpsn_gss_file'])
        check_file_exists(args.manual_sp_file)

        # create archaeal species names
        ar_out_dir = os.path.join(params['output_dirs']['u_pmc_species_names'], 'ar')
        make_sure_path_exists(ar_out_dir)
        self.log.info(f"Output directory (ar): {ar_out_dir}")

        p = PMC_SpeciesNames(ar_out_dir)
        p.run(None,
              taxonomy_init_ar,
              params['data_files']['seqcode_file'],
              None,
              args.manual_sp_file,
              gtdb_clusters_file,
              params['data_files']['prev_gtdb_metadata_file'],
              params['data_files']['cur_gtdb_metadata_file'],
              qc_passed,
              ncbi_misclassified_file,
              params['data_files']['ncbi_assembly_summary_genbank_file'],
              untrustworthy_type_file,
              synonym_file,
              updated_species_reps_file,
              params['ledgers']['gtdb_type_strains_ledger'],
              params['ledgers']['species_classification_ledger'],
              params['ledgers']['sp_priority_ledger'],
              params['ledgers']['genus_priority_ledger'],
              params['ledgers']['specific_epithet_ledger_ar'],
              params['ledgers']['ncbi_untrustworthy_sp_ledger'],
              params['ledgers']['ncbi_env_bioproject_ledger'],
              params['data_files']['lpsn_gss_file'])

        # create bacterial species names
        bac_out_dir = os.path.join(params['output_dirs']['u_pmc_species_names'], 'bac')
        make_sure_path_exists(bac_out_dir)
        self.log.info(f"Output directory (bac): {bac_out_dir}")

        p = PMC_SpeciesNames(bac_out_dir)
        p.run(None,
              taxonomy_init_bac,
              params['data_files']['seqcode_file'],
              None,
              args.manual_sp_file,
              gtdb_clusters_file,
              params['data_files']['prev_gtdb_metadata_file'],
              params['data_files']['cur_gtdb_metadata_file'],
              qc_passed,
              ncbi_misclassified_file,
              params['data_files']['ncbi_assembly_summary_genbank_file'],
              untrustworthy_type_file,
              synonym_file,
              updated_species_reps_file,
              params['ledgers']['gtdb_type_strains_ledger'],
              params['ledgers']['species_classification_ledger'],
              params['ledgers']['sp_priority_ledger'],
              params['ledgers']['genus_priority_ledger'],
              params['ledgers']['specific_epithet_ledger_bac'],
              params['ledgers']['ncbi_untrustworthy_sp_ledger'],
              params['ledgers']['ncbi_env_bioproject_ledger'],
              params['data_files']['lpsn_gss_file'])

    def pmc_manual_species(self, args):
        """Identify species names manually set by curators."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        if args.domain == 'Archaea':
            init_taxonomy = os.path.join(params['root_dir'],
                                         params['output_dirs']['u_pmc_species_names'],
                                         'ar',
                                         'final_taxonomy.tsv')
        else:
            init_taxonomy = os.path.join(params['root_dir'],
                                         params['output_dirs']['u_pmc_species_names'],
                                         'bac',
                                         'final_taxonomy.tsv')

        check_file_exists(init_taxonomy)
        check_file_exists(args.manually_curated_tree)
        make_sure_path_exists(args.output_dir)

        p = PMC_Checks(args.output_dir)
        p.manual_species(init_taxonomy,
                         args.manually_curated_tree)

    def pmc_replace_generic(self, args):
        """Replace generic names with genus assignment."""

        check_file_exists(args.manual_species_names)
        check_file_exists(args.manual_taxonomy)
        make_sure_path_exists(args.output_dir)

        p = PMC_Checks(args.output_dir)
        p.replace_generic(args.manual_species_names,
                          args.manual_taxonomy)

    def pmc_check_type_species(self, args):
        """Check for agreement between GTDB genera and genomes assembled from type species of genus."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        cur_gtdb_metadata_file = os.path.join(params['root_dir'], params['data_files']['cur_gtdb_metadata_file'])
        qc_passed_file = os.path.join(params['root_dir'], params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')
        ncbi_assembly_summary_genbank_file = os.path.join(
            params['root_dir'], params['data_files']['ncbi_assembly_summary_genbank_file'])
        untrustworthy_type_file = os.path.join(params['root_dir'], params['output_dirs']['u_resolve_types'], 'untrustworthy_type_material.tsv')
        gtdb_type_strains_ledger = os.path.join(params['root_dir'], params['ledgers']['gtdb_type_strains_ledger'])
        sp_priority_ledger = os.path.join(params['root_dir'], params['ledgers']['sp_priority_ledger'])
        genus_priority_ledger = os.path.join(params['root_dir'], params['ledgers']['genus_priority_ledger'])
        ncbi_env_bioproject_ledger = os.path.join(params['root_dir'], params['ledgers']['ncbi_env_bioproject_ledger'])
        lpsn_gss_file = os.path.join(params['root_dir'], params['data_files']['lpsn_gss_file'])

        check_file_exists(args.manual_taxonomy)
        check_file_exists(args.manual_sp_file)
        check_file_exists(args.pmc_custom_species)
        check_file_exists(cur_gtdb_metadata_file)
        check_file_exists(qc_passed_file)
        check_file_exists(ncbi_assembly_summary_genbank_file)
        check_file_exists(untrustworthy_type_file)
        check_file_exists(gtdb_type_strains_ledger)
        check_file_exists(sp_priority_ledger)
        check_file_exists(genus_priority_ledger)
        check_file_exists(ncbi_env_bioproject_ledger)
        check_file_exists(lpsn_gss_file)
        make_sure_path_exists(args.output_dir)

        p = PMC_CheckTypeSpecies(args.output_dir)
        p.run(args.manual_taxonomy,
                args.manual_sp_file,
                args.pmc_custom_species,
                cur_gtdb_metadata_file,
                qc_passed_file,
                ncbi_assembly_summary_genbank_file,
                untrustworthy_type_file,
                gtdb_type_strains_ledger,
                sp_priority_ledger,
                genus_priority_ledger,
                ncbi_env_bioproject_ledger,
                lpsn_gss_file)

    def pmc_check_type_strains(self, args):
        """Check for agreement between GTDB species and genomes assembled from type strain of species."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        cur_gtdb_metadata_file = os.path.join(params['root_dir'], params['data_files']['cur_gtdb_metadata_file'])
        qc_passed_file = os.path.join(params['root_dir'], params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')
        ncbi_assembly_summary_genbank_file = os.path.join(
            params['root_dir'], params['data_files']['ncbi_assembly_summary_genbank_file'])
        untrustworthy_type_file = os.path.join(params['root_dir'], params['output_dirs']['u_resolve_types'], 'untrustworthy_type_material.tsv')
        gtdb_type_strains_ledger = os.path.join(params['root_dir'], params['ledgers']['gtdb_type_strains_ledger'])
        seqcode_file = os.path.join(params['root_dir'], params['data_files']['seqcode_file'])
        ncbi_env_bioproject_ledger = os.path.join(params['root_dir'], params['ledgers']['ncbi_env_bioproject_ledger'])
        
        check_file_exists(args.manual_taxonomy)
        check_file_exists(args.manual_sp_file)
        check_file_exists(args.pmc_custom_species)
        check_file_exists(cur_gtdb_metadata_file)
        check_file_exists(qc_passed_file)
        check_file_exists(ncbi_assembly_summary_genbank_file)
        check_file_exists(untrustworthy_type_file)
        check_file_exists(gtdb_type_strains_ledger)
        check_file_exists(seqcode_file)
        check_file_exists(ncbi_env_bioproject_ledger)
        make_sure_path_exists(args.output_dir)

        p = PMC_CheckTypeStrains(args.output_dir)
        p.run(args.manual_taxonomy,
                args.manual_sp_file,
                args.pmc_custom_species,
                cur_gtdb_metadata_file,
                qc_passed_file,
                ncbi_assembly_summary_genbank_file,
                untrustworthy_type_file,
                gtdb_type_strains_ledger,
                seqcode_file,
                ncbi_env_bioproject_ledger)

    def pmc_species_names(self, args):
        """Establish final species names based on manual curation."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        if args.domain == 'Archaea':
            specific_epithet_ledger = os.path.join(params['root_dir'], params['ledgers']['specific_epithet_ledger_ar'])
        else:
            specific_epithet_ledger = os.path.join(params['root_dir'], params['ledgers']['specific_epithet_ledger_bac'])

        seqcode_file = os.path.join(params['root_dir'], params['data_files']['seqcode_file'])
        gtdb_clusters_file = os.path.join(params['root_dir'], params['output_dirs']['u_cluster_de_novo'], 'gtdb_clusters_de_novo.tsv')
        prev_gtdb_metadata_file = os.path.join(params['root_dir'], params['data_files']['prev_gtdb_metadata_file'])
        cur_gtdb_metadata_file = os.path.join(params['root_dir'], params['data_files']['cur_gtdb_metadata_file'])
        qc_passed_file = os.path.join(params['root_dir'], params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')
        ncbi_misclassified_file = os.path.join(
            params['root_dir'], params['output_dirs']['u_ncbi_erroneous'], 'ncbi_misclassified_sp.gtdb_clustering.tsv')
        ncbi_genbank_assembly_file = os.path.join(
            params['root_dir'], params['data_files']['ncbi_assembly_summary_genbank_file'])
        untrustworthy_type_file = os.path.join(params['root_dir'], params['output_dirs']['u_resolve_types'], 'untrustworthy_type_material.tsv')
        synonym_file = os.path.join(params['root_dir'], params['output_dirs']['u_synonyms'], 'synonyms.tsv')
        updated_species_reps = os.path.join(params['root_dir'], params['output_dirs']
                                            ['u_rep_actions'], 'updated_species_reps.tsv')
        gtdb_type_strains_ledger = os.path.join(params['root_dir'], params['ledgers']['gtdb_type_strains_ledger'])
        species_classification_ledger = os.path.join(params['root_dir'], params['ledgers']['species_classification_ledger'])
        sp_priority_ledger = os.path.join(params['root_dir'], params['ledgers']['sp_priority_ledger'])
        genus_priority_ledger = os.path.join(params['root_dir'], params['ledgers']['genus_priority_ledger'])
        ncbi_untrustworthy_sp_ledger = os.path.join(params['root_dir'], params['ledgers']['ncbi_untrustworthy_sp_ledger'])
        ncbi_env_bioproject_ledger = os.path.join(params['root_dir'], params['ledgers']['ncbi_env_bioproject_ledger'])
        lpsn_gss_file = os.path.join(params['root_dir'], params['data_files']['lpsn_gss_file'])

        check_file_exists(args.curation_tree)
        check_file_exists(args.manual_taxonomy)
        check_file_exists(seqcode_file)
        check_file_exists(args.manual_sp_names)
        check_file_exists(args.pmc_custom_species)
        check_file_exists(gtdb_clusters_file)
        check_file_exists(prev_gtdb_metadata_file)
        check_file_exists(cur_gtdb_metadata_file)
        check_file_exists(qc_passed_file)
        check_file_exists(ncbi_misclassified_file)
        check_file_exists(ncbi_genbank_assembly_file)
        check_file_exists(untrustworthy_type_file)
        check_file_exists(synonym_file)
        check_file_exists(updated_species_reps)
        check_file_exists(gtdb_type_strains_ledger)
        check_file_exists(species_classification_ledger)
        check_file_exists(sp_priority_ledger)
        check_file_exists(genus_priority_ledger)
        check_file_exists(specific_epithet_ledger)
        check_file_exists(ncbi_untrustworthy_sp_ledger)
        check_file_exists(ncbi_env_bioproject_ledger)
        check_file_exists(lpsn_gss_file)
        make_sure_path_exists(args.output_dir)

        p = PMC_SpeciesNames(args.output_dir)
        p.run(args.curation_tree,
              args.manual_taxonomy,
              seqcode_file,
              args.manual_sp_names,
              args.pmc_custom_species,
              gtdb_clusters_file,
              prev_gtdb_metadata_file,
              cur_gtdb_metadata_file,
              qc_passed_file,
              ncbi_misclassified_file,
              ncbi_genbank_assembly_file,
              untrustworthy_type_file,
              synonym_file,
              updated_species_reps,
              gtdb_type_strains_ledger,
              species_classification_ledger,
              sp_priority_ledger,
              genus_priority_ledger,
              specific_epithet_ledger,
              ncbi_untrustworthy_sp_ledger,
              ncbi_env_bioproject_ledger,
              lpsn_gss_file)

    def pmc_validate(self, args):
        """Validate final species names."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        if args.domain == 'Archaea':
            specific_epithet_ledger = os.path.join(params['root_dir'], params['ledgers']['specific_epithet_ledger_ar'])
        else:
            specific_epithet_ledger = os.path.join(params['root_dir'], params['ledgers']['specific_epithet_ledger_bac'])

        gtdb_clusters_file = os.path.join(params['root_dir'], params['output_dirs']['u_cluster_de_novo'], 'gtdb_clusters_de_novo.tsv')
        prev_gtdb_metadata_file = os.path.join(params['root_dir'], params['data_files']['prev_gtdb_metadata_file'])
        cur_gtdb_metadata_file = os.path.join(params['root_dir'], params['data_files']['cur_gtdb_metadata_file'])
        qc_passed_file = os.path.join(params['root_dir'], params['output_dirs']['u_qc_genomes'], 'qc_passed.tsv')
        ncbi_misclassified_file = os.path.join(params['root_dir'], params['output_dirs']['u_ncbi_erroneous'], 'ncbi_misclassified_sp.gtdb_clustering.tsv')
        ncbi_genbank_assembly_file = os.path.join(params['root_dir'], params['data_files']['ncbi_assembly_summary_genbank_file'])
        untrustworthy_type_file = os.path.join(params['root_dir'], params['output_dirs']['u_resolve_types'], 'untrustworthy_type_material.tsv')
        synonym_file = os.path.join(params['root_dir'], params['output_dirs']['u_synonyms'], 'synonyms.tsv')
        updated_species_reps = os.path.join(params['root_dir'], params['output_dirs']['u_rep_actions'], 'updated_species_reps.tsv')
        gtdb_type_strains_ledger = os.path.join(params['root_dir'], params['ledgers']['gtdb_type_strains_ledger'])
        species_classification_ledger = os.path.join(params['root_dir'], params['ledgers']['species_classification_ledger'])
        sp_priority_ledger = os.path.join(params['root_dir'], params['ledgers']['sp_priority_ledger'])
        genus_priority_ledger = os.path.join(params['root_dir'], params['ledgers']['genus_priority_ledger'])
        ncbi_env_bioproject_ledger = os.path.join(params['root_dir'], params['ledgers']['ncbi_env_bioproject_ledger'])
        lpsn_gss_metadata_file = os.path.join(params['root_dir'], params['data_files']['lpsn_gss_file'])
        lpsn_gss_file = os.path.join(params['root_dir'], params['data_files']['lpsn_metadata_file'])
        seqcode_file = os.path.join(params['root_dir'], params['data_files']['seqcode_file'])

        check_file_exists(args.final_taxonomy)

        if args.final_scaled_tree.lower() != 'none':
            check_file_exists(args.final_scaled_tree)

        if args.manual_sp_names.lower() != 'none':
            check_file_exists(args.manual_sp_names)

        if args.pmc_custom_species.lower() != 'none':
            check_file_exists(args.pmc_custom_species)

        check_file_exists(gtdb_clusters_file)
        check_file_exists(prev_gtdb_metadata_file)
        check_file_exists(cur_gtdb_metadata_file)
        check_file_exists(qc_passed_file)
        check_file_exists(ncbi_misclassified_file)
        check_file_exists(ncbi_genbank_assembly_file)
        check_file_exists(untrustworthy_type_file)
        check_file_exists(synonym_file)
        check_file_exists(updated_species_reps)
        check_file_exists(gtdb_type_strains_ledger)
        check_file_exists(species_classification_ledger)
        check_file_exists(sp_priority_ledger)
        check_file_exists(genus_priority_ledger)
        check_file_exists(specific_epithet_ledger)
        check_file_exists(ncbi_env_bioproject_ledger)
        check_file_exists(lpsn_gss_metadata_file)
        check_file_exists(lpsn_gss_file)
        check_file_exists(seqcode_file)

        if args.ground_truth_test_cases.lower() != 'none':
            check_file_exists(args.ground_truth_test_cases)

        make_sure_path_exists(args.output_dir)

        p = PMC_Validation(args.output_dir)
        p.run(args.domain,
              args.final_taxonomy,
              args.final_scaled_tree,
              args.manual_sp_names,
              args.pmc_custom_species,
              gtdb_clusters_file,
              prev_gtdb_metadata_file,
              cur_gtdb_metadata_file,
              qc_passed_file,
              ncbi_misclassified_file,
              ncbi_genbank_assembly_file,
              untrustworthy_type_file,
              synonym_file,
              updated_species_reps,
              gtdb_type_strains_ledger,
              species_classification_ledger,
              sp_priority_ledger,
              genus_priority_ledger,
              specific_epithet_ledger,
              ncbi_env_bioproject_ledger,
              lpsn_gss_metadata_file,
              lpsn_gss_file,
              seqcode_file,
              args.ground_truth_test_cases,
              args.skip_full_taxonomy_checks,
              args.skip_genus_checks)

    def pmc_cluster_stats(self, args):
        """Calculate final statistics for species cluster."""

        check_file_exists(args.sp_cluster_file)
        check_file_exists(args.genome_path_file)
        check_file_exists(args.gtdb_metadata_file)

        p = PMC_ClusterStats(args.af_sp,
                             args.max_genomes,
                             args.ani_cache_file,
                             args.cpus,
                             args.output_dir)
        p.run(args.sp_cluster_file,
              args.genome_path_file,
              args.gtdb_metadata_file)

    def pmc_cleanup(self, args):
        """Remove temporary files and compress large files."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        p = PMC_Cleanup()
        p.run(params['root_dir'], params['output_dirs'])

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

    def ani_sp_pair(self, args):
        """Calculate all pairwise ANI/AF values between genomes in two species."""

        check_file_exists(args.gtdb_metadata_file)
        check_file_exists(args.genome_path_file)

        make_sure_path_exists(args.output_dir)

        p = ANI_SpeciesPair(args.ani_cache_file, args.cpus, args.output_dir)
        p.run(args.gtdb_metadata_file,
              args.genome_path_file,
              args.species1,
              args.species2)

    def intra_sp_derep(self, args):
        """Dereplicate GTDB species clusters using ANI/AF criteria."""

        check_file_exists(args.gtdb_metadata_file)
        check_file_exists(args.genomic_path_file)

        make_sure_path_exists(args.output_dir)

        p = IntraSpeciesDereplication(args.derep_ani,
                                      args.derep_af,
                                      args.max_genomes_per_sp,
                                      args.ani_cache_file,
                                      args.cpus,
                                      args.output_dir)
        p.run(args.gtdb_metadata_file,
              args.genomic_path_file)

    def intra_genus_ani(self, args):
        """Calculate intra-genus ANI/AF values between GTDB representative genomes."""

        check_file_exists(args.gtdb_metadata_file)
        check_file_exists(args.genomic_path_file)

        make_sure_path_exists(args.output_dir)

        p = IntraGenusANI(args.ani_cache_file,
                          args.cpus,
                          args.output_dir)

        p.run(args.target_genus,
              args.gtdb_metadata_file,
              args.genomic_path_file)

    def ani_af_reps(self, args):
        """Calculate ANI/AF between GTDB representative genomes with the same genus."""

        check_file_exists(args.gtdb_metadata_file)
        check_file_exists(args.genomic_path_file)

        make_sure_path_exists(args.output_dir)

        p = RepGenomicSimilarity(args.ani_cache_file,
                                 args.cpus,
                                 args.output_dir)

        p.run(args.gtdb_metadata_file,
              args.genomic_path_file)

    def type_status(self, args):
        """Report information related to a genome being type material."""

        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.qc_passed_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.untrustworthy_type_file)
        check_file_exists(args.gtdb_type_strains_ledger)
        check_file_exists(args.ncbi_env_bioproject_ledger)

        p = InspectGenomes()
        p.type_status(args.cur_gtdb_metadata_file,
                      args.qc_passed_file,
                      args.ncbi_genbank_assembly_file,
                      args.untrustworthy_type_file,
                      args.gtdb_type_strains_ledger,
                      args.ncbi_env_bioproject_ledger,
                      args.genome_ids)

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

        reader = csv.reader(open(args.cur_metadata_file))
        header = next(reader, None)
        if header:
            gtdb_rep_index = header.index('gtdb_representative')
            gtdb_taxonomy_index = header.index('gtdb_taxonomy')
        else:
            self.log.error(f"No header found in {args.cur_metadata_file}")
            sys.exit(1)

        for row in reader:
            gid = row[0]
            cur_gids.add(gid)

            gtdb_taxa = ['']*7
            gtdb_taxonomy = row[gtdb_taxonomy_index]
            if gtdb_taxonomy:
                gtdb_taxa = [t.strip()
                             for t in row[gtdb_taxonomy_index].split(';')]
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

        reader = csv.reader(open(args.prev_metadata_file))
        header = next(reader, None)
        if header:
            gtdb_rep_index = header.index('gtdb_representative')
            gtdb_taxonomy_index = header.index('gtdb_taxonomy')
        else:
            self.log.error(f"No header found in {args.prev_metadata_file}")
            sys.exit(1)

        for row in reader:
            if row[gtdb_rep_index] == 't':
                gid = row[0]
                gtdb_taxonomy = row[gtdb_taxonomy_index]
                if gtdb_taxonomy:
                    gtdb_taxa = [t.strip()
                                 for t in row[gtdb_taxonomy_index].split(';')]

                    prev_reps_taxa[gid] = gtdb_taxa

                    if gtdb_taxa[6] != 's__':
                        prev_rep_species.add(gtdb_taxa[6])

                    if gtdb_taxa[5] != 'g__':
                        prev_rep_genera.add(gtdb_taxa[5])

        # summarize differences
        print('No. current representatives: %d' % len(cur_reps_taxa))
        print('No. previous representatives: %d' % len(prev_reps_taxa))

        print('')
        print('No. current species with representatives: %d' %
              len(cur_rep_species))
        print('No. previous species with representatives: %d' %
              len(prev_rep_species))

        print('')
        print('No. new representatives: %d' %
              len(set(cur_reps_taxa) - set(prev_reps_taxa)))
        print('No. retired representatives: %d' %
              len(set(prev_reps_taxa) - set(cur_reps_taxa)))

        print('')
        print('No. new species with representative: %d' %
              len(cur_rep_species - prev_rep_species))
        print('No. new genera with representative: %d' %
              len(cur_rep_genera - prev_rep_genera))

        print('')
        missing_sp_reps = prev_rep_species.intersection(
            cur_species) - cur_rep_species
        print('No. species that no longer have a representative: %d' %
              len(missing_sp_reps))
        for sp in missing_sp_reps:
            print('  ' + sp)

        print('')
        missing_genera_reps = prev_rep_genera.intersection(
            cur_genera) - cur_rep_genera
        print('No. genera that no longer have a representative: %d' %
              len(missing_genera_reps))
        for g in missing_genera_reps:
            print('  ' + g)

        print('')
        deprecated_reps = set(prev_reps_taxa).intersection(
            cur_gids) - set(cur_reps_taxa)
        print('No. deprecated previous representatives: %d' %
              len(deprecated_reps))

    def sandbox(self, args):
        """Play to explore new ideas or calculate one off statistics relted to species clusters."""

        check_file_exists(args.prev_gtdb_metadata_file)
        check_file_exists(args.cur_gtdb_metadata_file)
        check_file_exists(args.qc_passed_file)
        check_file_exists(args.ncbi_genbank_assembly_file)
        check_file_exists(args.untrustworthy_type_file)
        check_file_exists(args.gtdb_type_strains_ledger)
        check_file_exists(args.sp_priority_ledger)
        check_file_exists(args.genus_priority_ledger)
        check_file_exists(args.ncbi_env_bioproject_ledger)
        check_file_exists(args.lpsn_gss_file)
        make_sure_path_exists(args.output_dir)

        p = Sandbox(args.output_dir)
        p.run(args.prev_gtdb_metadata_file,
              args.cur_gtdb_metadata_file,
              args.qc_passed_file,
              args.ncbi_genbank_assembly_file,
              args.untrustworthy_type_file,
              args.gtdb_type_strains_ledger,
              args.sp_priority_ledger,
              args.genus_priority_ledger,
              args.ncbi_env_bioproject_ledger,
              args.lpsn_gss_file)

    def run(self, args):
        """Parse user arguments and call the correct pipeline(s)"""

        logging.basicConfig(format='', level=logging.INFO)

        if args.subparser_name == 'u_new_genomes':
            self.u_new_genomes(args)
        elif args.subparser_name == 'u_qc_genomes':
            self.u_qc_genomes(args)
        elif args.subparser_name == 'u_lpsn_rna_types':
            self.u_lpsn_rna_types(args)
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
        elif args.subparser_name == 'u_ncbi_erroneous':
            self.u_ncbi_erroneous(args)
        elif args.subparser_name == 'u_synonyms':
            self.u_synonyms(args)
        elif args.subparser_name == 'u_cluster_de_novo':
            self.u_cluster_de_novo(args)
        elif args.subparser_name == 'u_cluster_stats':
            self.u_cluster_stats(args)
        elif args.subparser_name == 'u_curation_trees':
            self.u_curation_trees(args)
        elif args.subparser_name == 'u_species_init':
            self.u_species_init(args)
        elif args.subparser_name == 'u_pmc_species_names':
            self.u_pmc_species_names(args)
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
        elif args.subparser_name == 'pmc_cluster_stats':
            self.pmc_cluster_stats(args)
        elif args.subparser_name == 'pmc_cleanup':
            self.pmc_cleanup(args)
        elif args.subparser_name == 'merge_test':
            self.merge_test(args)
        elif args.subparser_name == 'ani_sp_pair':
            self.ani_sp_pair(args)
        elif args.subparser_name == 'intra_sp_derep':
            self.intra_sp_derep(args)
        elif args.subparser_name == 'intra_genus_ani':
            self.intra_genus_ani(args)
        elif args.subparser_name == 'ani_af_reps':
            self.ani_af_reps(args)
        elif args.subparser_name == 'type_status':
            self.type_status(args)
        elif args.subparser_name == 'rep_compare':
            self.rep_compare(args)
        elif args.subparser_name == 'sandbox':
            self.sandbox(args)
        else:
            self.log.error(
                f'Unknown gtdb_species_clusters command: {args.subparser_name}\n')
            sys.exit()

        return 0
