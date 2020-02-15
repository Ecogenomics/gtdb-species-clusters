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
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
                                    
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
                            ani_af,
                            sp_priority_ledger):
        """Create table indicating species names that should be considered synonyms."""
        
        sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger)

        out_file = os.path.join(self.output_dir, 'synonyms.tsv')
        fout = open(out_file, 'w')
        fout.write('NCBI species\tGTDB species\tRepresentative\tStrain IDs\tRepresentative type sources\tPriority year\tGTDB type species\tGTDB type strain\tNCBI assembly type')
        fout.write('\tNCBI synonym\tGTDB synonym\tSynonym genome\tSynonym strain IDs\tSynonym type sources\tPriority year\tGTDB type species\tGTDB type strain\tSynonym NCBI assembly type')
        fout.write('\tANI\tAF\tWarnings\n')
        
        incorrect_priority = 0
        failed_type_strain_priority = 0
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
                
                if cur_genomes[rid].is_gtdb_type_strain() and cur_genomes[gid].is_gtdb_type_strain():
                        priority_gid, note = sp_priority_mngr.priority(cur_genomes, rid, gid)
                        if priority_gid != rid:
                            incorrect_priority += 1
                            fout.write('\tIncorrect priority: {}'.format(note))
                elif not cur_genomes[rid].is_gtdb_type_strain() and cur_genomes[gid].is_gtdb_type_strain():
                        failed_type_strain_priority += 1
                        fout.write('\tFailed to prioritize type strain of species')

                    
                fout.write('\n')
        
        if incorrect_priority:
            self.logger.warning(f'Identified {incorrect_priority:,} synonyms with incorrect priority.')
            
        if failed_type_strain_priority:
            self.logger.warning(f'Identified {failed_type_strain_priority:,} synonyms that failed to priotize the type strain of the species.')

    def run(self, named_cluster_file,
                    cur_gtdb_metadata_file,
                    uba_genome_paths,
                    qc_passed_file,
                    gtdbtk_classify_file,
                    ncbi_genbank_assembly_file,
                    untrustworthy_type_file,
                    ani_af_rep_vs_nonrep,
                    gtdb_type_strains_ledger,
                    sp_priority_ledger):
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
        self.logger.info('Determining effective type strain genomes in each NCBI species.')
        ncbi_sp_type_strain_genomes = cur_genomes.ncbi_sp_effective_type_genomes()
        self.logger.info(' ... identified effective type strain genomes for {:,} NCBI species.'.format(
                            len(ncbi_sp_type_strain_genomes)))
        
        # verify that type genomes for a species are contained in a 
        # single GTDB species cluster
        rid_map = {}
        for rid, gids in clusters.items():
            rid_map[rid] = rid
            for gid in gids:
                rid_map[gid] = rid
                
        for ncbi_sp, type_gids in ncbi_sp_type_strain_genomes.items():
            gtdb_rids = set([rid_map[gid] for gid in type_gids if gid in rid_map])
            if len(gtdb_rids) > 1:
                self.logger.warning('Type strain genomes from NCBI species {} were assigned to {:,} GTDB species clusters: {}.'.format(
                                    ncbi_sp,
                                    len(gtdb_rids),
                                    [(rid, cur_genomes[rid].gtdb_taxa.species) for rid in gtdb_rids]))

        # identify synonyms
        self.logger.info('Identifying synonyms.')
        synonyms = defaultdict(list)
        failed_type_strain_priority = 0
        for rid, gids in clusters.items():
            rep_ncbi_sp = cur_genomes[rid].ncbi_taxa.species
            
            # find species that are a synonym to the current representative,
            # using the best quality genome for each species to establish
            # synonym statistics such as ANI and AF
            type_gids = [gid for gid in gids if cur_genomes[gid].is_effective_type_strain()]
            if not cur_genomes[rid].is_effective_type_strain() and len(type_gids) > 0:
                failed_type_strain_priority += 1
                continue

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
        
        if failed_type_strain_priority:
            self.logger.warning(f'Identified {failed_type_strain_priority:,} non-type strain representatives that failed to priotize an effective type strain genome.')

        # read ANI and AF between representatives and non-representative genomes
        self.logger.info('Reading ANI and AF between representative and non-representative genomes.')
        ani_af = pickle.load(open(ani_af_rep_vs_nonrep, 'rb'))
        
        # write out synonyms
        self.write_synonym_table(synonyms,
                                    cur_genomes,
                                    ani_af,
                                    sp_priority_ledger)
        
