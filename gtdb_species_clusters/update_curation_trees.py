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

from biolib.taxonomy import Taxonomy

from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.type_genome_utils import read_clusters


class UpdateCurationTrees(object):
    """Produce curation trees highlighting new NCBI taxa."""

    def __init__(self, output_dir, output_prefix):
        """Initialization."""
        
        self.output_dir = output_dir
        self.output_prefix = output_prefix
        
        self.logger = logging.getLogger('timestamp')
         
    def new_ncbi_taxa(self, prev_genomes, cur_genomes, cur_clusters):
        """Determine new NCBI taxa as these should be considered by curators."""
        
        for rank_index in range(1, 7):
            rank_label = Taxonomy.rank_labels[rank_index]
            self.logger.info(f'Determining new NCBI taxa at rank of {rank_label}.')

            # get NCBI taxa in previous GTDB release
            prev_ncbi_taxa = set()
            for gid in prev_genomes:
                ncbi_taxa = prev_genomes[gid].ncbi_taxa.get_taxa(rank_index)
                prev_ncbi_taxa.add(ncbi_taxa)
            self.logger.info(' - identified {:,} NCBI taxa in previous GTDB release.'.format(
                                len(prev_ncbi_taxa)))
            
            # get NCBI taxa in current GTDB release
            cur_ncbi_taxa = set()
            taxon_domain = {}
            for rid in cur_clusters:
                ncbi_taxa = cur_genomes[rid].ncbi_taxa.get_taxa(rank_index)
                cur_ncbi_taxa.add(ncbi_taxa)
                
                taxon_domain[ncbi_taxa] = cur_genomes[rid].ncbi_taxa.domain
            self.logger.info(' - identified {:,} NCBI taxa in current GTDB release.'.format(
                                len(cur_ncbi_taxa)))
            
            # determine new NCBI taxa
            new_ncbi_taxa = cur_ncbi_taxa - prev_ncbi_taxa
            self.logger.info(' - identified {:,} NCBI taxa that are new to the current GTDB release.'.format(
                                len(new_ncbi_taxa)))
                                
            # determine genomes from new NCBI taxa
            new_taxa_rids = defaultdict(list)
            for rid in cur_clusters:
                ncbi_taxa = cur_genomes[rid].ncbi_taxa.get_taxa(rank_index)
                if ncbi_taxa in new_ncbi_taxa:
                    new_taxa_rids[ncbi_taxa].append(rid)
                    
            # report each new NCBI taxa
            out_table = os.path.join(self.output_dir, f'{self.output_prefix}_new_ncbi_{rank_label}.tsv')
            fout = open(out_table, 'w')
            fout.write('NCBI domain\tNCBI taxa\tNo. representatives\tRepresentative ID(s)\n')
            new_rids = set()
            for new_taxa in new_ncbi_taxa:
                new_rids.update(new_taxa_rids[new_taxa])
                fout.write('{}\t{}\t{}\t{}\n'.format(
                            taxon_domain[new_taxa],
                            new_taxa,
                            len(new_taxa_rids[new_taxa]),
                            ', '.join(new_taxa_rids[new_taxa])))
            fout.close()
            
            out_tree = os.path.join(self.output_dir, f'{self.output_prefix}_new_ncbi_{rank_label}.tree')
            fout_tree = open(out_tree, 'w')
            fout_tree.write('({});\n'.format(','.join(new_rids)))
            fout_tree.close()

    def run(self, 
                gtdb_clusters_file,
                prev_gtdb_metadata_file,
                cur_gtdb_metadata_file,
                qc_passed_file,
                ncbi_genbank_assembly_file,
                untrustworthy_type_file,
                gtdb_type_strains_ledger,
                ncbi_env_bioproject_ledger):
        """Perform initial actions required for changed representatives."""

        # create previous and current GTDB genome sets
        self.logger.info('Creating previous GTDB genome set.')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                                untrustworthy_type_ledger=untrustworthy_type_file,
                                                ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)
        self.logger.info(' - previous genome set has {:,} species clusters spanning {:,} genomes.'.format(
                            len(prev_genomes.sp_clusters),
                            prev_genomes.sp_clusters.total_num_genomes()))

        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                create_sp_clusters=False,
                                                qc_passed_file=qc_passed_file,
                                                ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                                untrustworthy_type_ledger=untrustworthy_type_file,
                                                ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)

        # read named GTDB species clusters
        self.logger.info('Reading GTDB species clusters.')
        cur_clusters, _ = read_clusters(gtdb_clusters_file)
        self.logger.info(' - identified {:,} clusters spanning {:,} genomes.'.format(
                            len(cur_clusters),
                            sum([len(gids) + 1 for gids in cur_clusters.values()])))

        # create curate tree and table indicating new NCBI taxa as these
        # should be considered by GTDB curators
        self.new_ncbi_taxa(prev_genomes, 
                                cur_genomes, 
                                cur_clusters)
 