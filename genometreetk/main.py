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

import dendropy

from biolib.common import check_file_exists, make_sure_path_exists
from biolib.external.execute import check_dependencies
from biolib.taxonomy import Taxonomy

from genometreetk.exceptions import GenomeTreeTkError
from genometreetk.trusted_genome_workflow import TrustedGenomeWorkflow
from genometreetk.dereplication_workflow import DereplicationWorkflow
from genometreetk.marker_workflow import MarkerWorkflow
from genometreetk.infer_workflow import InferWorkflow
from genometreetk.ssu_workflow import SSU_Workflow
from genometreetk.bootstrap import Bootstrap
from genometreetk.jackknife_markers import JackknifeMarkers
from genometreetk.jackknife_taxa import JackknifeTaxa
from genometreetk.combine_support import CombineSupport
from genometreetk.reroot_tree import RerootTree
from genometreetk.representatives import Representatives
from genometreetk.cluster import Cluster
from genometreetk.common import read_gtdb_metadata


class OptionsParser():
    def __init__(self):
        """Initialization"""
        self.logger = logging.getLogger()

    def _read_config_file(self):
        """Read configuration info."""

        cfg_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'genometreetk.cfg')

        d = {}
        for line in open(cfg_file):
            key, value = line.split('=')
            d[key] = value.strip()

        return d

    def trusted(self, options):
        """Determine trusted genomes."""

        check_file_exists(options.metadata_file)

        trusted_genome_workflow = TrustedGenomeWorkflow()

        trusted_genome_workflow.run(options.metadata_file,
                                        options.trusted_comp,
                                        options.trusted_cont,
                                        options.max_contigs,
                                        options.min_N50,
                                        options.trusted_genomes_file)

        self.logger.info('Trusted genome list written to: %s' % options.trusted_genomes_file)

    def dereplicate(self, options):
        """Dereplicate genomes based on taxonomy."""

        if options.trusted_genomes_file:
            check_file_exists(options.trusted_genomes_file)

        try:
            dereplication_workflow = DereplicationWorkflow()

            dereplication_workflow.run(options.max_species,
                                       options.trusted_genomes_file,
                                       None,
                                       None,
                                       options.derep_genome_file)
        except GenomeTreeTkError as e:
            print e.message
            raise SystemExit

        self.logger.info('Dereplicated genome list written to: %s' % options.derep_genome_file)

    def markers(self, options):
        """Determine marker genes."""

        make_sure_path_exists(options.output_dir)

        config_data = self._read_config_file()

        self.logger.error('NEED TO FIX ISSUE WITH GENOME DIR FILES!')
        sys.exit(-1)

        marker_workflow = MarkerWorkflow(config_data['genome_dir_file'],
                                            config_data['pfam_model_file'],
                                            config_data['tigrfams_model_dir'],
                                            options.cpus)

        phylo_hmm_out = marker_workflow.run(options.ingroup_file,
                                                options.ubiquity,
                                                options.single_copy,
                                                options.redundancy,
                                                options.min_support,
                                                options.min_per_taxa,
                                                options.perc_markers,
                                                options.restict_marker_list,
                                                options.output_dir)

        self.logger.info('Marker genes written to: %s' % phylo_hmm_out)

    def infer(self, options):
        """Infer genome tree."""

        check_file_exists(options.genome_id_file)
        check_file_exists(options.marker_id_file)
        make_sure_path_exists(options.output_dir)

        self.logger.error('NEED TO FIX ISSUE WITH GENOME DIR FILES!')
        sys.exit(-1)

        config_data = self._read_config_file()
        infer_workflow = InferWorkflow(config_data['genome_dir_file'],
                                            config_data['pfam_model_file'],
                                            config_data['tigrfams_model_dir'],
                                            options.cpus)

        infer_workflow.run(options.genome_id_file,
                                options.marker_id_file,
                                options.model,
                                options.output_dir)

        self.logger.info('Results written to: %s' % options.output_dir)

    def ssu_tree(self, options):
        """Infer 16S tree spanning GTDB genomes."""

        check_dependencies(['ssu-align', 'ssu-mask', 'FastTreeMP'])

        check_file_exists(options.gtdb_metadata_file)
        check_file_exists(options.gtdb_dir_file)
        make_sure_path_exists(options.output_dir)

        ssu_workflow = SSU_Workflow(options.gtdb_metadata_file, options.gtdb_dir_file)
        ssu_workflow.run(options.min_ssu_length,
                         options.min_ssu_contig,
                         options.min_quality,
                         options.max_contigs,
                         options.min_N50,
                         options.ncbi_rep_only,
                         options.user_genomes,
                         options.output_dir)

        self.logger.info('Results written to: %s' % options.output_dir)

    def bootstrap(self, options):
        """Bootstrap multiple sequence alignment."""

        check_file_exists(options.input_tree)
        check_file_exists(options.msa_file)
        make_sure_path_exists(options.output_dir)

        bootstrap = Bootstrap(options.cpus)
        output_tree = bootstrap.run(options.input_tree,
                                    options.msa_file,
                                    options.num_replicates,
                                    options.model,
                                    options.base_type,
                                    options.fraction,
                                    options.output_dir)

        self.logger.info('Bootstrapped tree written to: %s' % output_tree)

    def jk_markers(self, options):
        """Jackknife marker genes."""

        check_file_exists(options.input_tree)
        check_file_exists(options.msa_file)
        make_sure_path_exists(options.output_dir)

        jackknife_markers = JackknifeMarkers(options.cpus)
        output_tree = jackknife_markers.run(options.input_tree,
                                                options.msa_file,
                                                options.gene_length_file,
                                                options.perc_markers,
                                                options.num_replicates,
                                                options.model,
                                                options.gamma,
                                                options.output_dir)

        self.logger.info('Jackknifed marker tree written to: %s' % output_tree)

    def jk_taxa(self, options):
        """Jackknife taxa."""

        check_file_exists(options.input_tree)
        check_file_exists(options.msa_file)
        make_sure_path_exists(options.output_dir)

        jackknife_taxa = JackknifeTaxa(options.cpus)
        output_tree = jackknife_taxa.run(options.input_tree,
                                            options.msa_file,
                                            options.outgroup_ids,
                                            options.perc_taxa,
                                            options.num_replicates,
                                            options.model,
                                            options.gamma,
                                            options.output_dir)

        self.logger.info('Jackknifed taxa tree written to: %s' % output_tree)

    def combine(self, options):
        """Combine support values into a single tree."""

        combineSupport = CombineSupport()
        combineSupport.run(options.support_type,
                            options.bootstrap_tree,
                            options.jk_marker_tree,
                            options.jk_taxa_tree,
                            options.output_tree)

    def support_wf(self, options):
        """"Perform entire tree support workflow."""

        self.bootstrap(options)
        self.jk_markers(options)
        self.jk_taxa(options)
        self.combine(options)

    def midpoint(self, options):
        """"Midpoint root tree."""

        reroot = RerootTree()
        reroot.midpoint(options.input_tree, options.output_tree)

    def outgroup(self, options):
        """Reroot tree with outgroup."""

        check_file_exists(options.taxonomy_file)

        self.logger.info('Identifying genomes from the specified outgroup.')
        outgroup = set()
        for genome_id, taxa in Taxonomy().read(options.taxonomy_file).iteritems():
            if options.outgroup_taxon in taxa:
                outgroup.add(genome_id)
        self.logger.info('Identifying %d genomes in the outgroup.' % len(outgroup))

        reroot = RerootTree()
        reroot.root_with_outgroup(options.input_tree, options.output_tree, outgroup)

    def refseq_representatives(self, options):
        """Determine representative genomes in RefSeq."""

        check_file_exists(options.metadata_file)

        try:
            rep_workflow = DereplicationWorkflow()

            rep_workflow.run(options.max_species,
                                       None,
                                       options.metadata_file,
                                       options.min_rep_comp,
                                       options.max_rep_cont,
                                       options.rep_genome_file)
        except GenomeTreeTkError as e:
            print e.message
            raise SystemExit

        self.logger.info('RefSeq representative genomes written to: %s' % options.rep_genome_file)

    def representatives(self, options):
        """Determine additional representatives genomes."""

        check_file_exists(options.ar_msa_file)
        check_file_exists(options.bac_msa_file)
        check_file_exists(options.refseq_representatives)
        check_file_exists(options.metadata_file)

        try:
            representatives = Representatives()
            representatives.run(options.refseq_representatives,
                                    options.ar_msa_file,
                                    options.bac_msa_file,
                                    options.threshold,
                                    options.min_rep_comp,
                                    options.max_rep_cont,
                                    options.metadata_file,
                                    options.rep_genome_file)

            self.logger.info('Representative genomes written to: %s' % options.rep_genome_file)

        except GenomeTreeTkError as e:
            print e.message
            raise SystemExit

    def aai_cluster(self, options):
        """Cluster genomes based on AAI."""

        check_file_exists(options.ar_msa_file)
        check_file_exists(options.bac_msa_file)
        check_file_exists(options.representatives)
        check_file_exists(options.metadata_file)

        try:
            cluster = Cluster(options.cpus)
            cluster.run(options.representatives,
                        options.ar_msa_file,
                        options.bac_msa_file,
                        options.threshold,
                        options.metadata_file,
                        options.cluster_file)

            self.logger.info('Clustering information written to: %s' % options.cluster_file)

        except GenomeTreeTkError as e:
            print e.message
            raise SystemExit

    def validate(self, options):
        """Check taxonomy file is formatted as expected."""

        check_file_exists(options.input_taxonomy)

        taxonomy = Taxonomy()
        t = taxonomy.read(options.input_taxonomy)

        taxonomy.validate(t,
                          check_prefixes=True,
                          check_ranks=True,
                          check_hierarchy=True,
                          check_species=True,
                          report_errors=True)

        self.logger.info('Finished performing validation tests.')

    def fill_ranks(self, options):
        """Ensure taxonomy strings contain all 7 canonical ranks."""

        check_file_exists(options.input_taxonomy)

        fout = open(options.output_taxonomy, 'w')
        taxonomy = Taxonomy()
        t = taxonomy.read(options.input_taxonomy)

        for genome_id, taxon_list in t.iteritems():
            full_taxon_list = taxonomy.fill_missing_ranks(taxon_list)

            taxonomy_str = ';'.join(full_taxon_list)
            if not taxonomy.check_full(taxonomy_str):
                sys.exit(-1)

            fout.write('%s\t%s\n' % (genome_id, taxonomy_str))

        fout.close()

        self.logger.info('Revised taxonomy written to: %s' % options.output_taxonomy)

    def binomial(self, options):
        """Ensure species are designated using binomial nomenclature."""

        check_file_exists(options.input_taxonomy)

        fout = open(options.output_taxonomy, 'w')
        taxonomy = Taxonomy()
        t = taxonomy.read(options.input_taxonomy)

        for genome_id, taxon_list in t.iteritems():
            taxonomy_str = ';'.join(taxon_list)
            if not taxonomy.check_full(taxonomy_str):
                sys.exit(-1)

            genus = taxon_list[5][3:]
            species = taxon_list[6][3:]
            if species and genus not in species:
                taxon_list[6] = 's__' + genus + ' ' + species
                taxonomy_str = ';'.join(taxon_list)

            fout.write('%s\t%s\n' % (genome_id, taxonomy_str))

        fout.close()


        self.logger.info('Revised taxonomy written to: %s' % options.output_taxonomy)

    def propagate(self, options):
        """Propagate labels to all genomes in a cluster."""

        check_file_exists(options.input_taxonomy)
        check_file_exists(options.metadata_file)

        # get representative genome information
        rep_metadata = read_gtdb_metadata(options.metadata_file, ['gtdb_representative',
                                                                  'gtdb_clustered_genomes'])

        taxonomy = Taxonomy()
        t = taxonomy.read(options.input_taxonomy)
        expanded_taxonomy = {}
        for genome_id, taxon_list in t.iteritems():
            if genome_id in expanded_taxonomy:
                # this should never happen
                self.logger.error('Genome has both a representative and explicit taxonomy string: %s' % cluster_genome_id)
                sys.exit(-1)

            taxonomy_str = ';'.join(taxon_list)
            expanded_taxonomy[genome_id] = taxonomy_str

            # propagate taxonomy strings from representatives. Note that a genome may
            # not have metadata as it is possible a User has removed a genome that
            # is in the provided taxonomy file.
            _rep_genome, clustered_genomes = rep_metadata.get(genome_id, (None, None))
            if clustered_genomes:
                for cluster_genome_id in clustered_genomes.split(';'):
                    if cluster_genome_id == genome_id:
                        continue

                    if cluster_genome_id in expanded_taxonomy:
                        # this should never happen
                        self.logger.error('Genome has both a representative and explicit taxonomy string: %s' % cluster_genome_id)
                        sys.exit(-1)

                    expanded_taxonomy[cluster_genome_id] = taxonomy_str

        fout = open(options.output_taxonomy, 'w')
        for genome_id, taxonomy_str in expanded_taxonomy.iteritems():
            fout.write('%s\t%s\n' % (genome_id, taxonomy_str))
        fout.close()

        self.logger.info('Taxonomy written to: %s' % options.output_taxonomy)

    def strip(self, options):
        """Remove taxonomic labels from tree."""

        check_file_exists(options.input_tree)

        outgroup_in_tree = set()
        tree = dendropy.Tree.get_from_path(options.input_tree,
                                            schema='newick',
                                            rooting='force-rooted',
                                            preserve_underscores=True)

        for node in tree.internal_nodes():
            if node.label:
                if ':' in node.label:
                    support, _taxa = node.label.split(':')
                    node.label = support

        tree.write_to_path(options.output_tree,
                            schema='newick',
                            suppress_rooting=True,
                            unquoted_underscores=True)

        self.logger.info('Stripped tree written to: %s' % options.output_tree)

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        logging.basicConfig(format='', level=logging.INFO)

        check_dependencies(('FastTree', 'hmmsearch'))

        if options.subparser_name == 'trusted':
            self.trusted(options)
        elif options.subparser_name == 'dereplicate':
            self.dereplicate(options)
        elif options.subparser_name == 'markers':
            self.markers(options)
        elif options.subparser_name == 'infer':
            self.infer(options)
        elif options.subparser_name == 'ssu_tree':
            self.ssu_tree(options)
        elif options.subparser_name == 'bootstrap':
            self.bootstrap(options)
        elif options.subparser_name == 'jk_markers':
            self.jk_markers(options)
        elif options.subparser_name == 'jk_taxa':
            self.jk_taxa(options)
        elif options.subparser_name == 'combine':
            self.combine(options)
        elif options.subparser_name == 'midpoint':
            self.midpoint(options)
        elif options.subparser_name == 'outgroup':
            self.outgroup(options)
        elif options.subparser_name == 'refseq_reps':
            self.refseq_representatives(options)
        elif options.subparser_name == 'reps':
            self.representatives(options)
        elif options.subparser_name == 'aai_cluster':
            self.aai_cluster(options)
        elif options.subparser_name == 'validate':
            self.validate(options)
        elif options.subparser_name == 'binomial':
            self.binomial(options)
        elif options.subparser_name == 'binomial':
            self.binomial(options)
        elif options.subparser_name == 'propagate':
            self.propagate(options)
        elif options.subparser_name == 'fill_ranks':
            self.fill_ranks(options)
        elif options.subparser_name == 'strip':
            self.strip(options)
        else:
            self.logger.error('  [Error] Unknown GenomeTreeTk command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
