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
import argparse
import logging
import pickle
from collections import defaultdict, Counter

from numpy import (mean as np_mean, std as np_std)

from biolib.taxonomy import Taxonomy

from gtdb_species_clusters.fastani import FastANI
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_clusters import SpeciesClusters
from gtdb_species_clusters.species_name_manager import SpeciesNameManager
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters.genome_utils import canonical_gid
from gtdb_species_clusters.type_genome_utils import (symmetric_ani, 
                                                        read_clusters, 
                                                        write_clusters, 
                                                        write_rep_radius)
from gtdb_species_clusters.taxon_utils import (generic_name,
                                                specific_epithet,
                                                parse_synonyms,
                                                longest_common_prefix,
                                                is_placeholder_taxon,
                                                is_placeholder_sp_epithet)


class UpdateGenusNames(object):
    """Update genus names as a precursor for establish binomial species names."""

    def __init__(self, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
    def _parse_explicit_taxa_updates(self, gtdb_taxa_updates_ledger):
        """Get explicit updates to previous GTDB taxon names."""
        
        explicit_taxon_updates = {}
        with open(gtdb_taxa_updates_ledger) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                
                prev_taxon = tokens[0].strip()
                assert '__' in prev_taxon
                    
                new_taxon = tokens[1].strip()
                assert '__' in new_taxon
                
                explicit_taxon_updates[prev_taxon] = new_taxon
                
        return explicit_taxon_updates
        
    def new_ncbi_genera(self, prev_genomes, cur_genomes, cur_clusters, gtdbtk_classify_file):
        """Determine new NCBI genera that likely need to be considered by curators."""
        
        self.logger.info('Determining new NCBI genera for consideration by curators.')
        
        # read GTDB-Tk classification information
        gtdbtk = {}
        with open(gtdbtk_classify_file) as f:
            header = f.readline().strip().split('\t')
            
            classification_index = header.index('classification')
            red_index = header.index('red_value')
            note_index = header.index('note')
            
            for line in f:
                tokens = line.strip().split('\t')
                gid = tokens[0]
                gid = canonical_gid(gid)
                
                classification = [t.strip() for t in tokens[classification_index].split(';')]
                red = tokens[red_index]
                if red != 'N/A':
                    red = float(tokens[red_index])
                
                gtdbtk[gid] = (classification, red, tokens[note_index])
                
        # get NCBI genera in previous GTDB release
        prev_ncbi_genera = set()
        for gid in prev_genomes:
            ncbi_genera = prev_genomes[gid].ncbi_taxa.genus
            prev_ncbi_genera.add(ncbi_genera)
        self.logger.info(' ... identified {:,} NCBI genera in previous GTDB release.'.format(
                            len(prev_ncbi_genera)))
        
        # get NCBI genera in current GTDB release
        cur_ncbi_genera = set()
        for rid in cur_clusters:
            ncbi_genera = cur_genomes[rid].ncbi_taxa.genus
            cur_ncbi_genera.add(ncbi_genera)
        self.logger.info(' ... identified {:,} NCBI genera in current GTDB release.'.format(
                            len(cur_ncbi_genera)))
        
        # determine new NCBI genera
        new_ncbi_genera = cur_ncbi_genera - prev_ncbi_genera
        self.logger.info(' ... identified {:,} NCBI genera that are new to the current GTDB release.'.format(
                            len(new_ncbi_genera)))
                            
        # determine genomes from new NCBI genera
        new_genera_rids = defaultdict(list)
        for rid in cur_clusters:
            ncbi_genera = cur_genomes[rid].ncbi_taxa.genus
            if ncbi_genera in new_ncbi_genera:
                new_genera_rids[ncbi_genera].append(rid)
                
        # report each new NCBI genera
        fout = open(os.path.join(self.output_dir, 'new_ncbi_genera.tsv'), 'w')
        fout.write('NCBI genus\tRepresentative ID(s)\tNew genomes\tGTDB-Tk classification\tNCBI genus conflict\tRED\tGTDB-Tk note\n')
        for new_genus in new_ncbi_genera:
            new_genomes = []
            red_values = []
            gtdbtk_classification = []
            gtdbtk_notes = []
            ncbi_genus_conflict = False
            for rid in new_genera_rids[new_genus]:
                new_genomes.append(str(rid not in prev_genomes))
                if rid in gtdbtk:
                    classification, red, note = gtdbtk[rid]
                    gtdbtk_notes.append(note)
                    
                    if red != 'N/A':
                        red_values.append('{:.2f}'.format(red))
                    else:
                        red_values.append('n/a')
                        
                    if classification[6] != 's__':
                        gtdbtk_classification.append(classification[6])
                    else:
                        gtdbtk_classification.append(classification[5])
                    
                    gtdbtk_genus = classification[5]
                    if gtdbtk_genus != 'g__' and gtdbtk_genus != new_genus:
                        ncbi_genus_conflict = True
                else:
                    red_values.append('n/a')
                    gtdbtk_classification.append('n/a')
                    gtdbtk_notes.append('n/a')

            fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        new_genus,
                        ', '.join(new_genera_rids[new_genus]),
                        ', '.join(new_genomes),
                        ', '.join(gtdbtk_classification),
                        str(ncbi_genus_conflict),
                        ', '.join(red_values),
                        ', '.join(gtdbtk_notes)))
        
        fout.close()
        
    def new_ncbi_families(self, prev_genomes, cur_genomes, cur_clusters, gtdbtk_classify_file):
        """Determine new NCBI families that likely need to be considered by curators."""
        
        self.logger.info('Determining new NCBI families for consideration by curators.')
        
        # read GTDB-Tk classification information
        gtdbtk = {}
        with open(gtdbtk_classify_file) as f:
            header = f.readline().strip().split('\t')
            
            classification_index = header.index('classification')
            red_index = header.index('red_value')
            note_index = header.index('note')
            
            for line in f:
                tokens = line.strip().split('\t')
                gid = tokens[0]
                gid = canonical_gid(gid)
                
                classification = [t.strip() for t in tokens[classification_index].split(';')]
                red = tokens[red_index]
                if red != 'N/A':
                    red = float(tokens[red_index])
                
                gtdbtk[gid] = (classification, red, tokens[note_index])
                
        # get NCBI families in previous GTDB release
        prev_ncbi_families = set()
        for gid in prev_genomes:
            ncbi_family = prev_genomes[gid].ncbi_taxa.family
            prev_ncbi_families.add(ncbi_family)
        self.logger.info(' ... identified {:,} NCBI families in previous GTDB release.'.format(
                            len(prev_ncbi_families)))
        
        # get NCBI families in current GTDB release
        cur_ncbi_families = set()
        for rid in cur_clusters:
            ncbi_family = cur_genomes[rid].ncbi_taxa.family
            cur_ncbi_families.add(ncbi_family)
        self.logger.info(' ... identified {:,} NCBI families in current GTDB release.'.format(
                            len(cur_ncbi_families)))
        
        # determine new NCBI families
        new_ncbi_families = cur_ncbi_families - prev_ncbi_families
        self.logger.info(' ... identified {:,} NCBI families that are new to the current GTDB release.'.format(
                            len(new_ncbi_families)))
                            
        # determine genomes from new NCBI families
        new_family_rids = defaultdict(list)
        for rid in cur_clusters:
            ncbi_family = cur_genomes[rid].ncbi_taxa.family
            if ncbi_family in new_ncbi_families:
                new_family_rids[ncbi_family].append(rid)
                
        # report each new NCBI families
        fout = open(os.path.join(self.output_dir, 'new_ncbi_families.tsv'), 'w')
        fout.write('NCBI family\tRepresentative ID(s)\tNew genomes\tGTDB-Tk classification\tNCBI family conflict\tRED\tGTDB-Tk note\n')
        for new_family in new_ncbi_families:
            new_genomes = []
            red_values = []
            gtdbtk_classification = []
            gtdbtk_notes = []
            ncbi_family_conflict = False
            for rid in new_family_rids[new_family]:
                new_genomes.append(str(rid not in prev_genomes))
                if rid in gtdbtk:
                    classification, red, note = gtdbtk[rid]
                    gtdbtk_notes.append(note)
                    
                    if red != 'N/A':
                        red_values.append('{:.2f}'.format(red))
                    else:
                        red_values.append('n/a')
                        
                    if classification[6] != 's__':
                        gtdbtk_classification.append(classification[6])
                    elif classification[5] != 'g__':
                        gtdbtk_classification.append(classification[5])
                    else:
                        gtdbtk_classification.append(classification[4])
                    
                    gtdbtk_family = classification[4]
                    if gtdbtk_family != 'f__' and gtdbtk_family != new_family:
                        ncbi_family_conflict = True
                else:
                    red_values.append('n/a')
                    gtdbtk_classification.append('n/a')
                    gtdbtk_notes.append('n/a')

            fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        new_family,
                        ', '.join(new_family_rids[new_family]),
                        ', '.join(new_genomes),
                        ', '.join(gtdbtk_classification),
                        str(ncbi_family_conflict),
                        ', '.join(red_values),
                        ', '.join(gtdbtk_notes)))
        
        fout.close()

    def run(self, 
                gtdb_clusters_file,
                prev_gtdb_metadata_file,
                cur_gtdb_metadata_file,
                uba_genome_paths,
                qc_passed_file,
                gtdbtk_classify_file,
                ncbi_genbank_assembly_file,
                untrustworthy_type_file,
                synonym_file,
                gtdb_type_strains_ledger,
                sp_priority_ledger,
                gtdb_taxa_updates_ledger,
                dsmz_bacnames_file):
        """Perform initial actions required for changed representatives."""

        # create previous and current GTDB genome sets
        self.logger.info('Creating previous GTDB genome set.')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                uba_genome_file=uba_genome_paths,
                                                ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                                untrustworthy_type_ledger=untrustworthy_type_file)
        self.logger.info(' ... previous genome set has {:,} species clusters spanning {:,} genomes.'.format(
                            len(prev_genomes.sp_clusters),
                            prev_genomes.sp_clusters.total_num_genomes()))

        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                create_sp_clusters=False,
                                                uba_genome_file=uba_genome_paths,
                                                qc_passed_file=qc_passed_file,
                                                ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                                untrustworthy_type_ledger=untrustworthy_type_file,
                                                gtdbtk_classify_file=gtdbtk_classify_file)
        self.logger.info(f' ... current genome set contains {len(cur_genomes):,} genomes.')
        
        # read named GTDB species clusters
        self.logger.info('Reading GTDB species clusters.')
        cur_clusters, _ = read_clusters(gtdb_clusters_file)
        self.logger.info(' ... identified {:,} clusters spanning {:,} genomes.'.format(
                            len(cur_clusters),
                            sum([len(gids) + 1 for gids in cur_clusters.values()])))

        # get list of synonyms in order to restrict usage of species names
        self.logger.info('Reading GTDB synonyms.')
        synonyms = parse_synonyms(synonym_file)
        self.logger.info(' ... identified {:,} synonyms from {:,} distinct species.'.format(
                            len(synonyms),
                            len(set(synonyms.values()))))
        
        # set current genomes to have same GTDB assignments as in previous
        # GTDB release. This is necessary, since genomes may have different
        # NCBI accession numbers between releases and thus the previous GTDB
        # taxonomy will not be reflected in the latest GTDB database. The 
        # exception is if a genome has changed domains, in which case the
        # previous assignment is invalid.
        self.logger.info('Setting GTDB taxonomy of genomes in current genome set.')
        update_count = 0
        conflicting_domain_count = 0
        for prev_gid in prev_genomes:
            if prev_gid in cur_genomes:
                if prev_genomes[prev_gid].gtdb_taxa != cur_genomes[prev_gid].gtdb_taxa:
                    if prev_genomes[prev_gid].gtdb_taxa.domain == cur_genomes[prev_gid].gtdb_taxa.domain:
                        update_count += 1
                        cur_genomes[prev_gid].gtdb_taxa.update_taxa(prev_genomes[prev_gid].gtdb_taxa)
                    else:
                        conflicting_domain_count += 1
        self.logger.info(f' ... updated {update_count:,} genomes.')
        self.logger.info(f' ... identified {conflicting_domain_count:,} genomes with conflicting domain assignments.')

        # get explicit updates to previous GTDB taxa
        self.logger.info('Reading explicit taxa updates.')
        explicit_taxon_updates = self._parse_explicit_taxa_updates(gtdb_taxa_updates_ledger)
        self.logger.info(f' ... identified {len(explicit_taxon_updates):,} updates.')
        
        self.logger.info('Updating current genomes to reflect explicit taxa updates.')
        update_count = 0
        for cur_taxon, new_taxon in explicit_taxon_updates.items():
            rank_prefix = cur_taxon[0:3]
            rank_index = Taxonomy.rank_prefixes.index(rank_prefix)

            for gid in cur_genomes:
                if cur_genomes[gid].gtdb_taxa.get_taxa(rank_index) == cur_taxon:
                    update_count += 1
                    cur_genomes[gid].gtdb_taxa.set_taxa(rank_index, new_taxon)
                    
                    if rank_prefix == 'g__':
                        # should also update the species name
                        new_sp = cur_genomes[gid].gtdb_taxa.species.replace(cur_taxon[3:], new_taxon[3:])
                        cur_genomes[gid].gtdb_taxa.set_taxa(rank_index+1, new_sp)
                
        self.logger.info(f' ... updated {update_count:,} genomes.')

        # initialize species priority manager
        self.sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger,
                                                        dsmz_bacnames_file)
                                                        
        # create table with new NCBI genera that likely need to be incorporated into
        # this release of the GTDB
        self.new_ncbi_genera(prev_genomes, 
                                cur_genomes, 
                                cur_clusters, 
                                gtdbtk_classify_file)
                                
        self.new_ncbi_families(prev_genomes, 
                                cur_genomes, 
                                cur_clusters, 
                                gtdbtk_classify_file)
 