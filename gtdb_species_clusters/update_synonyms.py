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

from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.ncbi_species_manager import NCBI_SpeciesManager
                                    
from gtdb_species_clusters.type_genome_utils import (read_clusters)


class UpdateSynonyms(object):
    """Determine synonyms for validly or effectively published species."""

    def __init__(self, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir

        self.logger = logging.getLogger('timestamp')

    def run(self, gtdb_clusters_file,
                    cur_gtdb_metadata_file,
                    uba_genome_paths,
                    qc_passed_file,
                    ncbi_misclassified_file,
                    ncbi_genbank_assembly_file,
                    untrustworthy_type_file,
                    ani_af_rep_vs_nonrep,
                    gtdb_type_strains_ledger,
                    sp_priority_ledger,
                    genus_priority_ledger,
                    dsmz_bacnames_file):
        """Cluster genomes to selected GTDB representatives."""
        
        # create current GTDB genome sets
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
        
        # read named GTDB species clusters
        self.logger.info('Reading named and previous placeholder GTDB species clusters.')
        cur_clusters, rep_radius = read_clusters(gtdb_clusters_file)
        self.logger.info(' - identified {:,} clusters spanning {:,} genomes.'.format(
                            len(cur_clusters),
                            sum([len(gids) + 1 for gids in cur_clusters.values()])))

        # identified genomes with misclassified species assignments at NCBI
        self.logger.info('Identify genomes with misclassified NCBI species assignments.')
        ncbi_species_mngr = NCBI_SpeciesManager(cur_genomes, 
                                                cur_clusters, 
                                                {}, # no manually-curated GTDB species assignments at this point
                                                self.output_dir)
        ncbi_misclassified_gids = ncbi_species_mngr.parse_ncbi_misclassified_table(ncbi_misclassified_file)
        self.logger.info(' - identified {:,} genomes with erroneous NCBI species assignments'.format(
                            len(ncbi_misclassified_gids)))
                            
        # identify NCBI species considered to be synonyms under the GTDB
        type_strain_synonyms = ncbi_species_mngr.identify_type_strain_synonyms(ncbi_misclassified_gids)
        consensus_synonyms = ncbi_species_mngr.identify_consensus_synonyms(ncbi_misclassified_gids)

        # read ANI and AF between representatives and non-representative genomes
        self.logger.info('Reading ANI and AF between representative and non-representative genomes.')
        ani_af = pickle.load(open(ani_af_rep_vs_nonrep, 'rb'))
        
        # write out synonyms
        ncbi_species_mngr.write_synonym_table(type_strain_synonyms,
                                                consensus_synonyms,
                                                ani_af,
                                                sp_priority_ledger,
                                                genus_priority_ledger,
                                                dsmz_bacnames_file)
