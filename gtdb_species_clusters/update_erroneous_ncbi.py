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
from gtdb_species_clusters.genome_utils import canonical_gid
from gtdb_species_clusters.type_genome_utils import (read_clusters, symmetric_ani)
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


class UpdateErroneousNCBI(object):
    """Identify genomes with erroneous NCBI species assignments."""

    def __init__(self, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
    def identify_misclassified_genomes_cluster(self, cur_genomes, cur_clusters):
        """Identify genomes with erroneous NCBI species assignments, based on GTDB clustering of type strain genomes."""
        
        forbidden_names = set(['cyanobacterium'])
        
        # get mapping from genomes to their representatives
        gid_to_rid = {}
        for rid, cids in cur_clusters.items():
            for cid in cids:
                gid_to_rid[cid] = rid
                
        # get genomes with NCBI species assignment
        ncbi_sp_gids = defaultdict(list)
        for gid in cur_genomes:
            ncbi_species = cur_genomes[gid].ncbi_taxa.species
            ncbi_specific = specific_epithet(ncbi_species)
            
            if ncbi_species != 's__' and ncbi_specific not in forbidden_names:
                ncbi_sp_gids[ncbi_species].append(gid)
                
        # get NCBI species anchored by a type strain genome
        ncbi_type_anchored_species = {}
        for rid, cids in cur_clusters.items():
            for cid in cids:
                if cur_genomes[cid].is_effective_type_strain():
                    ncbi_type_species = cur_genomes[cid].ncbi_taxa.species
                    ncbi_specific = specific_epithet(ncbi_species)
                    if ncbi_type_species != 's__' and ncbi_specific not in forbidden_names:
                        if (ncbi_type_species in ncbi_type_anchored_species
                            and rid != ncbi_type_anchored_species[ncbi_type_species]):
                            self.logger.error('NCBI species {} has multiple effective type strain genomes in different clusters.'.format(
                                                ncbi_type_species))
                            sys.exit(-1)
                            
                        ncbi_type_anchored_species[ncbi_type_species] = rid
        self.logger.info(' - identified {:,} NCBI species anchored by a type strain genome.'.format(
                            len(ncbi_type_anchored_species)))
        
        # identify genomes with erroneous NCBI species assignments
        fout = open(os.path.join(self.output_dir, 'ncbi_misclassified_sp.gtdb_clustering.tsv'), 'w')
        fout.write('Genome ID\tNCBI species\tGenome cluster\tType species cluster\n')
        
        misclassified_gids = set()
        for idx, (ncbi_species, species_gids) in enumerate(ncbi_sp_gids.items()):
            if ncbi_species not in ncbi_type_anchored_species:
                continue
                
            # find genomes with NCBI species assignments that are in a
            # different cluster than the type strain genome
            type_rid = ncbi_type_anchored_species[ncbi_species]
            for gid in species_gids:
                cur_rid = gid_to_rid[gid]
                if type_rid != cur_rid:
                    misclassified_gids.add(gid)
                    fout.write('{}\t{}\t{}\t{}\t\n'.format(
                                gid,
                                ncbi_species,
                                cur_rid,
                                type_rid))
        
        sys.stdout.write('\n')
        fout.close()

        misclassified_species = set([cur_genomes[gid].ncbi_taxa.species for gid in misclassified_gids])
        self.logger.info(' - identified {:,} genomes from {:,} species as having misclassified NCBI species assignments.'.format(
                            len(misclassified_gids),
                            len(misclassified_species)))
                            
        return misclassified_gids
        
    def run(self, gtdb_clusters_file,
                    cur_gtdb_metadata_file,
                    cur_genomic_path_file,
                    uba_genome_paths,
                    qc_passed_file,
                    ncbi_genbank_assembly_file,
                    untrustworthy_type_file,
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
        self.logger.info(f' ... current genome set contains {len(cur_genomes):,} genomes.')
        
        # get path to previous and current genomic FASTA files
        self.logger.info('Reading path to current genomic FASTA files.')
        cur_genomes.load_genomic_file_paths(cur_genomic_path_file)
        cur_genomes.load_genomic_file_paths(uba_genome_paths)
        
        # read named GTDB species clusters
        self.logger.info('Reading named and previous placeholder GTDB species clusters.')
        cur_clusters, rep_radius = read_clusters(gtdb_clusters_file)
        self.logger.info(' ... identified {:,} clusters spanning {:,} genomes.'.format(
                            len(cur_clusters),
                            sum([len(gids) + 1 for gids in cur_clusters.values()])))
        
        # identify genomes with erroneous NCBI species assignments
        self.logger.info('Identifying genomes with erroneous NCBI species assignments as established by GTDB cluster of type strain genomes.')
        self.identify_misclassified_genomes_cluster(cur_genomes, cur_clusters)
