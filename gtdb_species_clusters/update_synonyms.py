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
                                                        

from gtdb_species_clusters.taxon_utils import (specific_epithet)


class UpdateSynonyms(object):
    """Determine synonyms for validly or effectively published species."""

    def __init__(self, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir

        self.logger = logging.getLogger('timestamp')
        
    def write_synonym_table(self,
                            type_strain_synonyms,
                            majority_vote_synonyms,
                            cur_genomes,
                            ani_af,
                            sp_priority_ledger,
                            genus_priority_ledger,
                            dsmz_bacnames_file):
        """Create table indicating species names that should be considered synonyms."""
        
        sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger,
                                                    genus_priority_ledger,
                                                    dsmz_bacnames_file)

        out_file = os.path.join(self.output_dir, 'synonyms.tsv')
        fout = open(out_file, 'w')
        fout.write('Synonym type\tNCBI species\tGTDB representative\tStrain IDs\tType sources\tPriority year')
        fout.write('\tGTDB type species\tGTDB type strain\tNCBI assembly type')
        fout.write('\tNCBI synonym\tHighest-quality synonym genome\tSynonym strain IDs\tSynonym type sources\tSynonym priority year')
        fout.write('\tSynonym GTDB type species\tSynonym GTDB type strain\tSynonym NCBI assembly type')
        fout.write('\tANI\tAF\tWarnings\n')
        
        incorrect_priority = 0
        failed_type_strain_priority = 0
        for synonyms, synonym_type in [(type_strain_synonyms, 'TYPE_STRAIN_SYNONYM'), 
                                        (majority_vote_synonyms, 'MAJORITY_VOTE_SYNONYM')]:
            for rid, synonym_ids in synonyms.items():
                for gid in synonym_ids:
                    ani, af = symmetric_ani(ani_af, rid, gid)

                    fout.write(synonym_type)
                    fout.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                                cur_genomes[rid].ncbi_taxa.species,
                                rid,
                                ','.join(sorted(cur_genomes[rid].strain_ids())),
                                ','.join(sorted(cur_genomes[rid].gtdb_type_sources())).upper().replace('STRAININFO', 'StrainInfo'),
                                cur_genomes[rid].year_of_priority(),
                                cur_genomes[rid].is_gtdb_type_species(),
                                cur_genomes[rid].is_gtdb_type_strain(),
                                cur_genomes[rid].ncbi_type_material))
                    
                    synonym_priority_year = cur_genomes[gid].year_of_priority()
                    if synonym_priority_year == Genome.NO_PRIORITY_YEAR:
                        synonym_priority_year = 'n/a'
                    
                    fout.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                                cur_genomes[gid].ncbi_taxa.species,
                                gid,
                                ','.join(sorted(cur_genomes[gid].strain_ids())),
                                ','.join(sorted(cur_genomes[gid].gtdb_type_sources())).upper().replace('STRAININFO', 'StrainInfo'),
                                cur_genomes[gid].year_of_priority(),
                                cur_genomes[gid].is_gtdb_type_species(),
                                cur_genomes[gid].is_gtdb_type_strain(),
                                cur_genomes[gid].ncbi_type_material))
                    fout.write('\t{:.3f}\t{:.4f}'.format(ani, af))
                    
                    if cur_genomes[rid].is_gtdb_type_strain() and cur_genomes[gid].is_gtdb_type_strain():
                            priority_gid, note = sp_priority_mngr.species_priority(cur_genomes, rid, gid)
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

    def identify_type_strain_synonyms(self, cur_genomes, clusters, ncbi_sp_type_strain_genomes):
        """Identify synonyms arising from multiple type strain genomes residing in the same GTDB cluster."""

        # verify that type genomes for a species are contained in a single GTDB species cluster
        self.logger.info('Verifying that all type strain genomes for a NCBI species occur in a single GTDB cluster.')
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

        # identify type strain synonyms
        self.logger.info('Identifying type strain synonyms.')
        type_strain_synonyms = defaultdict(list)
        failed_type_strain_priority = 0
        for rid, gids in clusters.items():
            if not cur_genomes[rid].is_effective_type_strain():
                continue 
                
            rep_ncbi_sp = cur_genomes[rid].ncbi_taxa.species
            
            # find species that are a type strain synonym to the current representative,
            # using the best quality genome for each species to establish
            # synonym statistics such as ANI and AF
            type_gids = [gid for gid in gids if cur_genomes[gid].is_effective_type_strain()]
            q = {gid:cur_genomes[gid].score_type_strain() for gid in type_gids}
            q_sorted = sorted(q.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
            processed_sp = set([rep_ncbi_sp])
            for gid, _quality in q_sorted:
                cur_ncbi_sp = cur_genomes[gid].ncbi_taxa.species
                
                if cur_ncbi_sp in processed_sp:
                    continue
                
                type_strain_synonyms[rid].append(gid)
                processed_sp.add(cur_ncbi_sp)
                    
        self.logger.info(' ... identified {:,} GTDB representatives resulting in {:,} type strain synonyms.'.format(
                            len(type_strain_synonyms),
                            sum([len(gids) for gids in type_strain_synonyms.values()])))
        
        if failed_type_strain_priority:
            self.logger.warning(f'Identified {failed_type_strain_priority:,} non-type strain representatives that failed to priotize an effective type strain genome.')
            
        return type_strain_synonyms
        
    def identify_majority_vote_synonyms(self, cur_genomes, clusters, ncbi_sp_type_strain_genomes):
        """Identify synonyms arising from >50% of genomes with an NCBI species classification 
            which lack a type strain genome being contained in a GTDB cluster defined by 
            a type strain genome."""
           
        forbidden_specific_names = set(['cyanobacterium'])
        
        # get genomes with NCBI species assignment
        ncbi_species_gids = defaultdict(list)
        ncbi_all_sp_gids = set()
        for gid in cur_genomes:
            ncbi_species = cur_genomes[gid].ncbi_taxa.species
            ncbi_specific = specific_epithet(ncbi_species)
            if ncbi_species != 's__' and ncbi_specific not in forbidden_specific_names:
                ncbi_species_gids[ncbi_species].append(gid)
                ncbi_all_sp_gids.add(gid)
        
        # identify majority vote synonyms
        majority_vote_synonyms = defaultdict(list)
        for ncbi_species, ncbi_sp_gids in ncbi_species_gids.items():
            if ncbi_species in ncbi_sp_type_strain_genomes:
                # NCBI species define by type strain genome so should 
                # be either the representative of a GTDB species cluster
                # or a type strain synonym
                continue

            # determine if >50% of genomes in NCBI species are contained
            # in a single GTDB species represented by a type strain genome
            for gtdb_rid, cids in clusters.items():
                assert gtdb_rid in cids
                
                if not cur_genomes[gtdb_rid].is_effective_type_strain():
                    continue
                
                ncbi_cur_sp_in_cluster = cids.intersection(ncbi_sp_gids)
                ncbi_cur_sp_in_cluster_perc = len(ncbi_cur_sp_in_cluster)*100.0/len(ncbi_sp_gids)
                if ncbi_cur_sp_in_cluster_perc == 100:
                    # using the best quality genome in NCBI species to establish
                    # synonym statistics such as ANI and AF
                    q = {gid:cur_genomes[gid].score_type_strain() for gid in ncbi_sp_gids}
                    q_sorted = sorted(q.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
                    majority_vote_synonyms[gtdb_rid].append(q_sorted[0][0])
                    break
                        
        self.logger.info(' ... identified {:,} GTDB representatives resulting in {:,} majority vote synonyms.'.format(
                            len(majority_vote_synonyms),
                            sum([len(gids) for gids in majority_vote_synonyms.values()])))
                        
        return majority_vote_synonyms
            
    def run(self, gtdb_clusters_file,
                    cur_gtdb_metadata_file,
                    uba_genome_paths,
                    qc_passed_file,
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
        self.logger.info(f' ... current genome set contains {len(cur_genomes):,} genomes.')
        
        # read named GTDB species clusters
        self.logger.info('Reading named and previous placeholder GTDB species clusters.')
        clusters, rep_radius = read_clusters(gtdb_clusters_file)
        self.logger.info(' ... identified {:,} clusters spanning {:,} genomes.'.format(
                            len(clusters),
                            sum([len(gids) + 1 for gids in clusters.values()])))
                            
        # identify NCBI species with genomes assembled from type strain of species
        self.logger.info('Determining effective type strain genomes in each NCBI species.')
        ncbi_sp_type_strain_genomes = cur_genomes.ncbi_sp_effective_type_genomes()
        self.logger.info(' ... identified effective type strain genomes for {:,} NCBI species.'.format(
                            len(ncbi_sp_type_strain_genomes)))

        # identify type strain synonyms
        type_strain_synonyms = self.identify_type_strain_synonyms(cur_genomes, 
                                                                    clusters,
                                                                    ncbi_sp_type_strain_genomes)
        
        # identify majority vote synonyms
        majority_vote_synonyms = self.identify_majority_vote_synonyms(cur_genomes, 
                                                                        clusters,
                                                                        ncbi_sp_type_strain_genomes)

        # read ANI and AF between representatives and non-representative genomes
        self.logger.info('Reading ANI and AF between representative and non-representative genomes.')
        ani_af = pickle.load(open(ani_af_rep_vs_nonrep, 'rb'))
        
        # write out synonyms
        self.write_synonym_table(type_strain_synonyms,
                                    majority_vote_synonyms,
                                    cur_genomes,
                                    ani_af,
                                    sp_priority_ledger,
                                    genus_priority_ledger,
                                    dsmz_bacnames_file)
        
