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
import operator
import shutil
import tempfile
import ntpath
import pickle
from itertools import combinations
from collections import defaultdict, namedtuple

from gtdb_species_clusters.genome import Genome
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_clusters import SpeciesClusters
                                    
from gtdb_species_clusters.type_genome_utils import (GenomeRadius,
                                                        symmetric_ani,
                                                        read_clusters)
                                                        


class UpdateSynonyms(object):
    """Determine synonyms for validly or effectively published species."""

    def __init__(self, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir

        self.logger = logging.getLogger('timestamp')
        
    def write_synonym_table(self,
                            synonyms,
                            cur_genomes,
                            ani_af):
        """Create table indicating species names that should be considered synonyms."""

        out_file = os.path.join(self.output_dir, 'synonyms.tsv')
        fout = open(out_file, 'w')
        fout.write('NCBI species\tGTDB species\tRepresentative\tStrain IDs\tRepresentative type sources\tPriority year\tGTDB type species\tGTDB type strain\tNCBI assembly type')
        fout.write('\tNCBI synonym\tGTDB synonym\tSynonym genome\tSynonym strain IDs\tSynonym type sources\tPriority year\tGTDB type species\tGTDB type strain\tSynonym NCBI assembly type')
        fout.write('\tANI\tAF\tWarnings\n')
        
        ambiguous_priority = 0
        incorrect_priority = 0
        incorrect_type_species_of_genus = 0
        for rid, synonym_ids in synonyms.items():
            for gid in synonym_ids:
                ani, af = symmetric_ani(ani_af, rid, gid)

                fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                            cur_genomes[rid].ncbi_taxa.species,
                            cur_genomes[rid].gtdb_taxa.species,
                            rid,
                            ','.join(sorted(cur_genomes[rid].strain_ids())),
                            ','.join(sorted(cur_genomes[rid].gtdb_type_sources())).upper().replace('STRAININFO', 'StrainInfo'),
                            cur_genomes[rid].year_of_priority(),
                            cur_genomes[rid].is_gtdb_type_species(),
                            cur_genomes[rid].is_gtdb_type_strain(),
                            cur_genomes[rid].ncbi_type_material))
                fout.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                            cur_genomes[gid].ncbi_taxa.species,
                            cur_genomes[gid].gtdb_taxa.species,
                            gid,
                            ','.join(sorted(cur_genomes[gid].strain_ids())),
                            ','.join(sorted(cur_genomes[gid].gtdb_type_sources())).upper().replace('STRAININFO', 'StrainInfo'),
                            cur_genomes[gid].year_of_priority(),
                            cur_genomes[gid].is_gtdb_type_species(),
                            cur_genomes[gid].is_gtdb_type_strain(),
                            cur_genomes[gid].ncbi_type_material))
                fout.write('\t{:.3f}\t{:.4f}'.format(ani, af))
                
                if (cur_genomes[rid].year_of_priority() == cur_genomes[gid].year_of_priority()
                    and cur_genomes[rid].year_of_priority() != Genome.NO_PRIORITY_YEAR
                    and not cur_genomes[rid].is_gtdb_type_species()):
                    ambiguous_priority += 1
                    fout.write('\tAmbiguous priority')
                    
                if (cur_genomes[rid].year_of_priority() > cur_genomes[gid].year_of_priority()
                    and not cur_genomes[rid].is_gtdb_type_species()):
                    incorrect_priority += 1
                    fout.write('\tIncorrect priority of publication dates')
                    
                if cur_genomes[gid].is_gtdb_type_species():
                    incorrect_type_species_of_genus += 1
                    fout.write('\tIncorrect priority of type species of genus')
                    
                fout.write('\n')
        
        self.logger.info(f'Identified {ambiguous_priority:,} synonyms with ambiguous priority based on year of publication.')
        
        if incorrect_priority > 0:
            self.logger.info(f'Identified {incorrect_priority:,} synonyms with incorrect priority based on year of publication.')
            
        if incorrect_type_species_of_genus > 0:
            self.logger.info(f'Identified {incorrect_type_species_of_genus:,} synonyms with incorrect priority based on type species of genus.')

    def run(self, named_cluster_file,
                    cur_gtdb_metadata_file,
                    uba_genome_paths,
                    qc_passed_file,
                    gtdbtk_classify_file,
                    ncbi_genbank_assembly_file,
                    untrustworthy_type_file,
                    ani_af_rep_vs_nonrep,
                    gtdb_type_strains_ledger):
        """Cluster genomes to selected GTDB representatives."""
        
        # create current GTDB genome sets
        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                create_sp_clusters=False,
                                                uba_genome_file=uba_genome_paths,
                                                qc_passed_file=qc_passed_file,
                                                gtdbtk_classify_file=gtdbtk_classify_file,
                                                ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                                untrustworthy_type_ledger=untrustworthy_type_file)
        self.logger.info(f' ... current genome set contains {len(cur_genomes):,} genomes.')
        
        # read named GTDB species clusters
        self.logger.info('Reading named and previous placeholder GTDB species clusters.')
        clusters, rep_radius = read_clusters(named_cluster_file)
        self.logger.info(' ... identified {:,} clusters spanning {:,} genomes.'.format(
                            len(clusters),
                            sum([len(gids) + 1 for gids in clusters.values()])))

        # identify NCBI species with multiple genomes assembled from type strain of species
        self.logger.info('Determining number of type strain genomes in each NCBI species.')
        sp_type_strain_genomes = defaultdict(set)
        for gid in cur_genomes:
            if cur_genomes[gid].is_effective_type_strain():
                ncbi_sp = cur_genomes[gid].ncbi_taxa.species
                if ncbi_sp != 's__':
                    # yes, NCBI has genomes marked as assembled from type material
                    # that do not actually have a binomial species name
                    sp_type_strain_genomes[ncbi_sp].add(gid)
        
        # verify that type genomes for a species are contained in a 
        # single GTDB species cluster
        rid_map = {}
        for rid, gids in clusters.items():
            rid_map[rid] = rid
            for gid in gids:
                rid_map[gid] = rid
                
        for ncbi_sp, type_gids in sp_type_strain_genomes.items():
            gtdb_rids = set([rid_map[gid] for gid in type_gids if gid in rid_map])
            if len(gtdb_rids) > 1:
                self.logger.warning('Type strain genomes from NCBI species {} were assigned to {:,} GTDB species clusters: {}.'.format(
                                    ncbi_sp,
                                    len(gtdb_rids),
                                    [(rid, cur_genomes[rid].gtdb_taxa.species) for rid in gtdb_rids]))

        # identify synonyms
        self.logger.info('Identifying synonyms.')
        synonyms = defaultdict(list)
        for rid, gids in clusters.items():
            if not cur_genomes[rid].is_effective_type_strain():
                continue
                
            rep_ncbi_sp = cur_genomes[rid].ncbi_taxa.species
            
            # find species that are a synonym to the current representative,
            # using the best quality genome for each species to establish
            # synonym statistics such as ANI and AF
            type_gids = [gid for gid in gids if cur_genomes[gid].is_effective_type_strain()]
            q = {gid:cur_genomes[gid].score_type_strain() for gid in type_gids}
            q_sorted = sorted(q.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
            processed_sp = set()
            for gid, _quality in q_sorted:
                cur_ncbi_sp = cur_genomes[gid].ncbi_taxa.species
                
                if cur_ncbi_sp in processed_sp:
                    continue
                
                if cur_ncbi_sp != rep_ncbi_sp:
                    synonyms[rid].append(gid)
                    processed_sp.add(cur_ncbi_sp)
                    
        self.logger.info(' ... identified {:,} GTDB representatives resulting in {:,} synonyms.'.format(
                            len(synonyms),
                            sum([len(gids) for gids in synonyms.values()])))

        # read ANI and AF between representatives and non-representative genomes
        self.logger.info('Reading ANI and AF between representative and non-representative genomes.')
        ani_af = pickle.load(open(ani_af_rep_vs_nonrep, 'rb'))
        
        # write out synonyms
        self.write_synonym_table(synonyms,
                                    cur_genomes,
                                    ani_af)
        
