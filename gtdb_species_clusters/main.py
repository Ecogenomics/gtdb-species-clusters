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
from gtdb_species_clusters.update_gtdbtk import GTDB_Tk
from gtdb_species_clusters.update_rep_changes import RepChanges
from gtdb_species_clusters.update_rep_actions import RepActions

from gtdb_species_clusters.merge_test import MergeTest

from gtdb_species_clusters.exceptions import GTDB_Error


class OptionsParser():
    def __init__(self):
        """Initialization"""
        self.logger = logging.getLogger()
        
    def qc_genomes(self, options):
        """Quality check all potential GTDB genomes."""
        
        check_file_exists(options.gtdb_metadata_file)
        check_file_exists(options.ncbi_genbank_assembly_file)
        check_file_exists(options.gtdb_domain_report)
        check_file_exists(options.qc_exception_file)
        check_file_exists(options.species_exception_file)
        make_sure_path_exists(options.output_dir)

        try:
            p = QcGenomes(options.ncbi_genbank_assembly_file,
                        options.gtdb_domain_report,
                        options.qc_exception_file,
                        options.species_exception_file,
                        options.min_comp,
                        options.max_cont,
                        options.min_quality,
                        options.sh_exception,
                        options.min_perc_markers,
                        options.max_contigs,
                        options.min_N50,
                        options.max_ambiguous,
                        options.output_dir)
        except GTDB_Error as e:
            print(e.message)
            raise SystemExit

        self.logger.info('Quality checking information written to: %s' % options.output_dir)
        
    def select_type_genomes(self, options):
        """Select representative genomes for named species."""

        check_file_exists(options.qc_file)
        check_file_exists(options.gtdb_metadata_file)
        check_file_exists(options.genome_path_file)
        check_file_exists(options.prev_rep_file)
        check_file_exists(options.ncbi_refseq_assembly_file)
        check_file_exists(options.ncbi_genbank_assembly_file)
        check_file_exists(options.gtdb_domain_report)
        check_file_exists(options.species_exception_file)
        check_file_exists(options.gtdb_type_genome_file)
        make_sure_path_exists(options.output_dir)

        try:
            p = SelectTypeGenomes(options.ani_cache_file, options.cpus, options.output_dir)
            p.run(options.qc_file,
                        options.gtdb_metadata_file,
                        options.ltp_blast_file,
                        options.genome_path_file,
                        options.prev_rep_file,
                        options.ncbi_refseq_assembly_file,
                        options.ncbi_genbank_assembly_file,
                        options.gtdb_domain_report,
                        options.species_exception_file,
                        options.gtdb_type_genome_file)
        except GTDB_Error as e:
            print(e.message)
            raise SystemExit

        self.logger.info('GTDB type genomes written to: %s' % options.output_dir)

    def cluster_named_types(self, options):
        """Cluster genomes to selected GTDB type genomes."""

        check_file_exists(options.qc_file)
        check_file_exists(options.gtdb_metadata_file)
        check_file_exists(options.genome_path_file)
        check_file_exists(options.named_type_genome_file)
        check_file_exists(options.type_genome_ani_file)
        check_file_exists(options.species_exception_file)
        make_sure_path_exists(options.output_dir)

        try:
            p = ClusterNamedTypes(options.ani_sp,
                                    options.af_sp,
                                    options.ani_cache_file, 
                                    options.cpus,
                                    options.output_dir)
            p.run(options.qc_file,
                    options.gtdb_metadata_file,
                    options.genome_path_file,
                    options.named_type_genome_file,
                    options.type_genome_ani_file,
                    options.mash_sketch_file,
                    options.species_exception_file)
        except GTDB_Error as e:
            print(e.message)
            raise SystemExit

        self.logger.info('Clustering results written to: %s' % options.output_dir)
        
    def cluster_de_novo(self, options):
        """Infer de novo species clusters and type genomes for remaining genomes."""

        check_file_exists(options.qc_file)
        check_file_exists(options.gtdb_metadata_file)
        check_file_exists(options.gtdb_user_genomes_file)
        check_file_exists(options.genome_path_file)
        check_file_exists(options.type_genome_cluster_file)
        check_file_exists(options.type_genome_synonym_file)
        check_file_exists(options.ncbi_refseq_assembly_file)
        check_file_exists(options.ncbi_genbank_assembly_file)
        check_file_exists(options.ani_af_nontype_vs_type)
        check_file_exists(options.species_exception_file)
        make_sure_path_exists(options.output_dir)

        try:
            p = ClusterDeNovo(options.ani_sp,
                                    options.af_sp,
                                    options.ani_cache_file, 
                                    options.cpus,
                                    options.output_dir)
            p.run(options.qc_file,
                        options.gtdb_metadata_file,
                        options.gtdb_user_genomes_file,
                        options.genome_path_file,
                        options.type_genome_cluster_file,
                        options.type_genome_synonym_file,
                        options.ncbi_refseq_assembly_file,
                        options.ncbi_genbank_assembly_file,
                        options.ani_af_nontype_vs_type,
                        options.species_exception_file,
                        options.rnd_type_genome)
        except GTDB_Error as e:
            print(e.message)
            raise SystemExit

        self.logger.info('Clustering results written to: %s' % options.output_dir)
        
    def cluster_user(self, options):
        """Cluster User genomes to GTDB species clusters."""

        check_file_exists(options.gtdb_metadata_file)
        check_file_exists(options.genome_path_file)
        check_file_exists(options.final_cluster_file)
        make_sure_path_exists(options.output_dir)

        try:
            p = ClusterUser(options.ani_cache_file, 
                                options.cpus,
                                options.output_dir)
            p.run(options.gtdb_metadata_file,
                        options.genome_path_file,
                        options.final_cluster_file)
        except GTDB_Error as e:
            print(e.message)
            raise SystemExit

        self.logger.info('Clustering results written to: %s' % options.output_dir)
        
    def tree_gids(self, options):
        """Determine genome IDs for test/validation tree."""

        check_file_exists(options.qc_file)
        check_file_exists(options.gtdb_metadata_file)
        check_file_exists(options.gtdb_final_clusters)
        check_file_exists(options.species_exception_file)

        try:
            p = TreeGIDs()
            p.run(options.qc_file,
                    options.gtdb_metadata_file,
                    options.gtdb_final_clusters,
                    options.species_exception_file,
                    options.output_dir)
        except GTDB_Error as e:
            print(e.message)
            raise SystemExit

        self.logger.info('Results written to: %s' % options.output_dir)
            
    def u_new_genomes(self, options):
        """Identify new and updated genomes."""

        check_file_exists(options.prev_gtdb_metadata_file)
        check_file_exists(options.cur_gtdb_metadata_file)
        check_file_exists(options.cur_genome_paths)
        check_file_exists(options.ncbi_assembly_summary_genbank)

        make_sure_path_exists(options.output_dir)
        
        try:
            p = NewGenomes(options.output_dir)
            p.run(options.prev_gtdb_metadata_file,
                    options.cur_gtdb_metadata_file,
                    options.cur_genome_paths,
                    options.ncbi_assembly_summary_genbank)
        except GTDB_Error as e:
            print(e.message)
            raise SystemExit
            
        self.logger.info('Done.')
        
    def u_qc_genomes(self, options):
        """Quality check new and updated genomes."""
        
        check_file_exists(options.cur_gtdb_metadata_file)
        check_file_exists(options.cur_uba_gid_file)
        check_file_exists(options.cur_genbank_assembly_file)
        check_file_exists(options.cur_gtdb_domain_report)
        check_file_exists(options.qc_exception_file)
        check_file_exists(options.species_exception_file)
        check_file_exists(options.genus_exception_file)
        make_sure_path_exists(options.output_dir)

        try:
            p = QcGenomes()
            p.run(options.cur_gtdb_metadata_file,
                        options.cur_uba_gid_file,
                        options.cur_genbank_assembly_file,
                        options.cur_gtdb_domain_report,
                        options.qc_exception_file,
                        options.species_exception_file,
                        options.genus_exception_file,
                        options.min_comp,
                        options.max_cont,
                        options.min_quality,
                        options.sh_exception,
                        options.min_perc_markers,
                        options.max_contigs,
                        options.min_N50,
                        options.max_ambiguous,
                        options.output_dir)
        except GTDB_Error as e:
            print(e.message)
            raise SystemExit

        self.logger.info('Quality checking information written to: %s' % options.output_dir)
        
    def u_gtdbtk(self, options):
        """Perform initial classification of new and updated genomes using GTDB-Tk."""
        
        check_file_exists(options.genomes_new_updated_file)
        check_file_exists(options.qc_passed_file)
        make_sure_path_exists(options.output_dir)
        
        p = GTDB_Tk(options.cpus, options.output_dir)
        p.run(options.genomes_new_updated_file,
                options.qc_passed_file,
                options.batch_size)

        self.logger.info('Done.')
            
    def u_rep_changes(self, options):
        """Identify species representatives that have changed from previous release."""

        check_file_exists(options.prev_gtdb_metadata_file)
        check_file_exists(options.cur_gtdb_metadata_file)
        check_file_exists(options.cur_uba_gid_file)
        check_file_exists(options.genomes_new_updated_file)
        check_file_exists(options.qc_passed_file)
        check_file_exists(options.gtdbtk_classify_file)
        check_file_exists(options.species_exception_file)
        check_file_exists(options.genus_exception_file)
        make_sure_path_exists(options.output_dir)
        
        p = RepChanges(options.output_dir)
        p.run(options.prev_gtdb_metadata_file,
                options.cur_gtdb_metadata_file,
                options.cur_uba_gid_file,
                options.genomes_new_updated_file,
                options.qc_passed_file,
                options.gtdbtk_classify_file,
                options.species_exception_file,
                options.genus_exception_file)
        
        self.logger.info('Done.')
        
    def u_rep_actions(self, options):
        """Perform initial actions required for changed representatives."""
        
        check_file_exists(options.rep_change_summary_file)
        check_file_exists(options.prev_genomic_path_file)
        check_file_exists(options.cur_genomic_path_file)
        check_file_exists(options.uba_genome_paths)
        check_file_exists(options.genomes_new_updated_file)
        check_file_exists(options.qc_passed_file)
        check_file_exists(options.gtdbtk_classify_file)
        check_file_exists(options.species_exception_file)
        check_file_exists(options.genus_exception_file)
        make_sure_path_exists(options.output_dir)
        
        p = RepActions(options.ani_cache_file, options.output_dir)
        p.run(options.rep_change_summary_file,
                options.prev_gtdb_metadata_file,
                options.prev_genomic_path_file,
                options.cur_gtdb_metadata_file,
                options.cur_genomic_path_file,
                options.uba_genome_paths,
                options.genomes_new_updated_file,
                options.qc_passed_file,
                options.gtdbtk_classify_file,
                options.species_exception_file,
                options.genus_exception_file)
        
        self.logger.info('Done.')
        
    def merge_test(self, options):
        """Produce information relevant to merging two sister species."""
        
        check_file_exists(options.gtdb_metadata_file)
        check_file_exists(options.genome_path_file)
        
        make_sure_path_exists(options.output_dir)
        
        p = MergeTest(options.ani_cache_file, options.cpus, options.output_dir)
        p.run(options.gtdb_metadata_file,
                options.genome_path_file,
                options.species1,
                options.species2)
        
        self.logger.info('Done.')

    def rep_compare(self, options):
        """Compare current and previous representatives."""

        check_file_exists(options.cur_metadata_file)
        check_file_exists(options.prev_metadata_file)
        
        # get representatives in current taxonomy
        cur_gids = set()
        cur_species = set()
        cur_genera = set()
        cur_reps_taxa = {}
        cur_rep_species = set()
        cur_rep_genera = set()
        header = True
        for row in csv.reader(open(options.cur_metadata_file)):
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
        for row in csv.reader(open(options.prev_metadata_file)):
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
        
    def cluster_stats(self, options):
        """Calculate statistics for species cluster."""

        check_file_exists(options.cluster_file)
        check_file_exists(options.genome_path_file)
        check_file_exists(options.gtdb_metadata_file)
        
        p = ClusterStats(options.af_sp,
                            options.max_genomes,
                            options.ani_cache_file,
                            options.cpus, 
                            options.output_dir)
        p.run(options.cluster_file, 
                options.genome_path_file,
                options.gtdb_metadata_file)

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        logging.basicConfig(format='', level=logging.INFO)

        if options.subparser_name == 'qc_genomes':
            self.qc_genomes(options)
        elif options.subparser_name == 'mash_dist':
            self.mash_dist(options)
        elif options.subparser_name == 'select_type_genomes':
            self.select_type_genomes(options)
        elif options.subparser_name == 'cluster_named_types':
            self.cluster_named_types(options)
        elif options.subparser_name == 'cluster_de_novo':
            self.cluster_de_novo(options)
        elif options.subparser_name == 'cluster_user':
            self.cluster_user(options)
        elif options.subparser_name == 'tree_gids':
            self.tree_gids(options)
        elif options.subparser_name == 'u_new_genomes':
            self.u_new_genomes(options)
        elif options.subparser_name == 'u_qc_genomes':
            self.u_qc_genomes(options)
        elif options.subparser_name == 'u_gtdbtk':
            self.u_gtdbtk(options)
        elif options.subparser_name == 'u_rep_changes':
            self.u_rep_changes(options)
        elif options.subparser_name == 'u_rep_actions':
            self.u_rep_actions(options)
        elif options.subparser_name == 'merge_test':
            self.merge_test(options)
        elif options.subparser_name == 'rep_compare':
            self.rep_compare(options)
        elif options.subparser_name == 'cluster_stats':
            self.cluster_stats(options)
        else:
            self.logger.error('Unknown gtdb_species_clusters command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
