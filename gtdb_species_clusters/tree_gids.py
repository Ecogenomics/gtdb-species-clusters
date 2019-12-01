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
import random
import operator
from collections import defaultdict

from gtdb_species_clusters.common import read_gtdb_metadata
                                    
from gtdb_species_clusters.genome_utils import read_qc_file
                                    
from gtdb_species_clusters.taxon_utils import (read_gtdb_ncbi_taxonomy,
                                            read_gtdb_taxonomy)
                                    
from gtdb_species_clusters.type_genome_utils import (read_clusters, 
                                                        read_quality_metadata, 
                                                        quality_score)


class TreeGIDs(object):
    """Quality check all potential GTDB genomes."""

    def __init__(self):
        """Initialization."""
        
        self.logger = logging.getLogger('timestamp')
        
        
    def _incongruent_specific_names(self, species, ncbi_taxonomy, prev_gtdb_taxonomy, type_metadata, output_dir):
        """Identify species with incongruent specific names between GTDB releases."""
        
        fout = open(os.path.join(output_dir, 'incongruent_specific_name.tsv'), 'w')
        fout.write('Genome ID\tProposed species\tPrevious species')
        fout.write('\tgtdb_type_designation\tgtdb_type_designation_sources\tgtdb_type_species_of_genus')
        fout.write('\tGTDB taxonomy\tNCBI taxonomy\n')
        num_diff_specific = 0
        for rid, sp in species.items():
            prev_gtdb_sp = prev_gtdb_taxonomy[rid][6]
            if prev_gtdb_sp == 's__':
                continue
                
            prev_genus, prev_specific = prev_gtdb_sp.split()
            cur_sp = sp.replace('Candidatus ', '')
            cur_genus, cur_specific = cur_sp.split()

            if cur_specific != prev_specific:
                num_diff_specific += 1
                
                fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                                    rid, 
                                    sp, 
                                    prev_gtdb_sp,
                                    type_metadata[rid].gtdb_type_designation,
                                    type_metadata[rid].gtdb_type_designation_sources,
                                    type_metadata[rid].gtdb_type_species_of_genus,
                                    '; '.join(prev_gtdb_taxonomy[rid]),
                                    '; '.join(ncbi_taxonomy[rid])))
                
        fout.close()
                        
        self.logger.info('Identified %d genomes with incongruent specific species names.' % num_diff_specific)
        
    def _incongruent_genus_names(self, species, ncbi_taxonomy, prev_gtdb_taxonomy, type_metadata, output_dir):
        """Identify species with incongruent genus names between GTDB releases."""
        
        fout = open(os.path.join(output_dir, 'incongruent_genus_name.tsv'), 'w')
        fout.write('Genome ID\tProposed species\tPrevious species')
        fout.write('\tgtdb_type_designation\tgtdb_type_designation_sources\tgtdb_type_species_of_genus')
        fout.write('\tGTDB taxonomy\tNCBI taxonomy\n')
        num_diff_genus = 0
        for rid, sp in species.items():
            prev_gtdb_sp = prev_gtdb_taxonomy[rid][6]
            if prev_gtdb_sp == 's__':
                continue
                
            prev_genus, prev_specific = prev_gtdb_sp.split()
            cur_sp = sp.replace('Candidatus ', '')
            cur_genus, cur_specific = cur_sp.split()

            if prev_genus != cur_genus:
                num_diff_genus += 1
                
                fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                                    rid, 
                                    sp, 
                                    prev_gtdb_sp,
                                    type_metadata[rid].gtdb_type_designation,
                                    type_metadata[rid].gtdb_type_designation_sources,
                                    type_metadata[rid].gtdb_type_species_of_genus,
                                    '; '.join(prev_gtdb_taxonomy[rid]),
                                    '; '.join(ncbi_taxonomy[rid])))
                
        fout.close()
                        
        self.logger.info('Identified %d genomes with incongruent genus names.' % num_diff_genus)

    def run(self, 
                qc_file,
                gtdb_metadata_file,
                gtdb_final_clusters,
                species_exception_file,
                output_dir):
        """Quality check all potential GTDB genomes."""
        
        # identify genomes failing quality criteria
        self.logger.info('Reading QC file.')
        passed_qc = read_qc_file(qc_file)
        self.logger.info('Identified %d genomes passing QC.' % len(passed_qc))
        
        # get GTDB and NCBI taxonomy strings for each genome
        self.logger.info('Reading NCBI and GTDB taxonomy from GTDB metadata file.')
        ncbi_taxonomy, ncbi_update_count = read_gtdb_ncbi_taxonomy(gtdb_metadata_file, species_exception_file)
        prev_gtdb_taxonomy = read_gtdb_taxonomy(gtdb_metadata_file)
        self.logger.info('Read NCBI taxonomy for %d genomes with %d manually defined updates.' % (len(ncbi_taxonomy), ncbi_update_count))
        self.logger.info('Read GTDB taxonomy for %d genomes.' % len(prev_gtdb_taxonomy))
        
        # get GTDB metadata
        type_metadata = read_gtdb_metadata(gtdb_metadata_file, ['gtdb_type_designation',
                                                                    'gtdb_type_designation_sources',
                                                                    'gtdb_type_species_of_genus'])
                                                                    
        quality_metadata = read_quality_metadata(gtdb_metadata_file)

        # read species clusters
        sp_clusters, species, _rep_radius = read_clusters(gtdb_final_clusters)
        self.logger.info('Read %d species clusters.' % len(sp_clusters))
        
        # sanity check species clusters all defined by genomes passing QC
        for gid in sp_clusters:
            if gid not in passed_qc:
                self.logger.error('Genome %s defines a species cluster, but fails QC.' % gid)
                sys.exit(-1)
                
        # modify GTDB taxonomy to reflect new species clustering and report incongruencies
        self.logger.info('Identifying species with incongruent specific names.')
        self._incongruent_specific_names(species, 
                                            ncbi_taxonomy,
                                            prev_gtdb_taxonomy, 
                                            type_metadata, 
                                            output_dir)
        
        self._incongruent_genus_names(species, 
                                            ncbi_taxonomy,
                                            prev_gtdb_taxonomy, 
                                            type_metadata, 
                                            output_dir)
                                            

        # get GIDs for canonical and validation trees
        fout_bac_can_gtdb = open(os.path.join(output_dir, 'bac_can_taxonomy.tsv'), 'w')
        fout_bac_val_gtdb = open(os.path.join(output_dir, 'bac_val_taxonomy.tsv'), 'w')
        fout_ar_can_gtdb = open(os.path.join(output_dir, 'ar_can_taxonomy.tsv'), 'w')
        fout_ar_val_gtdb = open(os.path.join(output_dir, 'ar_val_taxonomy.tsv'), 'w')
            
        fout_bac_val = open(os.path.join(output_dir, 'gids_bac_validation.lst'), 'w')
        fout_ar_val = open(os.path.join(output_dir, 'gids_ar_validation.lst'), 'w')
        fout_bac_can = open(os.path.join(output_dir, 'gids_bac_canonical.lst'), 'w')
        fout_ar_can = open(os.path.join(output_dir, 'gids_ar_canonical.lst'), 'w')
        fout_bac_val.write('#Accession\tSpecies\tNote\n')
        fout_ar_val.write('#Accession\tSpecies\tNote\n')
        fout_bac_can.write('#Accession\tSpecies\tNote\n')
        fout_ar_can.write('#Accession\tSpecies\tNote\n')
   
        for rid in sp_clusters:
            domain = prev_gtdb_taxonomy[rid][0]
            if domain == 'd__Bacteria':
                fout_val = fout_bac_val
                fout_can = fout_bac_can
                
                fout_can_gtdb = fout_bac_can_gtdb
                fout_val_gtdb = fout_bac_val_gtdb
            elif domain == 'd__Archaea':
                fout_val = fout_ar_val
                fout_can = fout_ar_can
                fout_can_gtdb = fout_ar_can_gtdb
                fout_val_gtdb = fout_ar_val_gtdb
            else:
                self.logger.error('Genome %s has no GTDB domain assignment.' % rid)
                sys.exit(-1)
            
            # substitute proposed species name into GTDB taxonomy
            taxa = prev_gtdb_taxonomy[rid][0:6] + [species[rid]]
            new_gtdb_str = '; '.join(taxa)
            fout_can_gtdb.write('%s\t%s\n' % (rid, new_gtdb_str))
            fout_val_gtdb.write('%s\t%s\n' % (rid, new_gtdb_str))
            
            fout_val.write('%s\t%s\t%s\n' % (rid, species[rid], 'GTDB type or representative genome'))
            fout_can.write('%s\t%s\t%s\n' % (rid, species[rid], 'GTDB type or representative genome'))
            
            cluster_gids = set(sp_clusters[rid])
            for gid in cluster_gids:
                if gid not in passed_qc:
                    self.logger.error('Genome %s is in a species cluster, but fails QC.' % gid)
                    sys.exit(-1)
                    
            if len(cluster_gids) > 0:
                # select highest-quality genome
                q = quality_score(cluster_gids, quality_metadata)
                gid = max(q.items(), key=operator.itemgetter(1))[0]
                
                taxa = prev_gtdb_taxonomy[gid][0:6] + [species[rid]]
                new_gtdb_str = '; '.join(taxa)
    
                fout_val.write('%s\t%s\t%s\n' % (gid, species[rid], 'selected highest-quality genome (Q=%.2f)' % q[gid]))
                fout_val_gtdb.write('%s\t%s\n' % (gid, new_gtdb_str))
                    
        fout_bac_val.close()
        fout_ar_val.close()
        fout_bac_can.close()
        fout_ar_can.close()
        
        fout_bac_can_gtdb.close()
        fout_bac_val_gtdb.close()
        fout_ar_can_gtdb.close()
        fout_ar_val_gtdb.close()
        