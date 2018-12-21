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

from genometreetk.common import (binomial_species,
                                    read_gtdb_metadata,
                                    read_gtdb_ncbi_taxonomy)
                                    
from genometreetk.type_genome_utils import (exclude_from_refseq, 
                                            ncbi_type_strain_of_species,
                                            gtdb_type_strain_of_species,
                                            pass_qc)


class QcTypeGenomes(object):
    """Quality check potential type genomes."""

    def __init__(self):
        """Initialization."""
        
        self.logger = logging.getLogger('timestamp')
        
        self.true_str = ['t', 'T', 'true', 'True']

    def run(self, metadata_file,
                ncbi_refseq_assembly_file,
                ncbi_genbank_assembly_file,
                min_comp,
                max_cont,
                min_quality,
                sh_exception,
                max_contigs,
                min_N50,
                max_ambiguous,
                output_dir):
        """Quality check potential type genomes."""

        # get GTDB and NCBI taxonomy strings for each genome
        self.logger.info('Reading NCBI taxonomy from GTDB metadata file.')
        ncbi_taxonomy = read_gtdb_ncbi_taxonomy(metadata_file)
        ncbi_species = binomial_species(ncbi_taxonomy)
        self.logger.info('Reading NCBI taxonomy from GTDB metadata file.')
        
        # calculate quality score for genomes
        self.logger.info('Parsing QC statistics for each genome.')
        quality_metadata = read_gtdb_metadata(metadata_file, ['checkm_completeness',
                                                                'checkm_contamination',
                                                                'checkm_strain_heterogeneity',
                                                                'contig_count',
                                                                'n50_contigs',
                                                                'ambiguous_bases',
                                                                'genome_size'])
                                                                
        # parse NCBI assembly files
        self.logger.info('Parsing NCBI assembly files.')
        excluded_from_refseq_note = exclude_from_refseq(ncbi_refseq_assembly_file, ncbi_genbank_assembly_file)

        # get type material designations for each genome
        self.logger.info('Reading type material designations for genomes from GTDB metadata file.')
        type_metadata = read_gtdb_metadata(metadata_file, ['ncbi_type_material_designation',
                                                                'gtdb_type_designation',
                                                                'gtdb_type_designation_sources'])
                                                                
        ncbi_tsp = ncbi_type_strain_of_species(type_metadata)
        gtdb_tsp = gtdb_type_strain_of_species(type_metadata)
                                                                
        # QC genomes in each named species
        self.logger.info('Performing QC of type genome for each of the %d species.' % len(ncbi_species))
        
        fout_type_fail = open(os.path.join(output_dir, 'type_genomes_fail_qc.tsv'), 'w')
        fout_type_fail.write('Species\tAccession\tType sources\tNCBI assembly type\tGenome size (bp)')
        fout_type_fail.write('\tCompleteness (%)\tContamination (%)\tQuality\tStrain heterogeneity')
        fout_type_fail.write('\tNo. contigs\tN50 contigs\tAmbiguous bases\tNCBI exclude from RefSeq\tLost species\n')
        
        fout_fail = open(os.path.join(output_dir, 'species_fail_qc.tsv'), 'w')
        fout_fail.write('Species\tAccession\tAssembled from type material\tGenome size (bp)')
        fout_fail.write('\tCompleteness (%)\tContamination (%)\tQuality\tStrain heterogeneity')
        fout_fail.write('\tNo. contigs\tN50 contigs\tAmbiguous bases\tNCBI exclude from RefSeq\n')
        
        fout_sp_lost = open(os.path.join(output_dir, 'species_lost.tsv'), 'w')
        fout_sp_lost.write('Species\tNo. genomes\tNo. type genomes')
        fout_sp_lost.write('\tFail completeness\tFail contamination\tFail quality\tFail no. contigs\tFail N50 contigs\tFail ambiguous bases\n')
        
        lost_type = 0
        lost_sp = 0
        for sp, gids in ncbi_species.items():
            type_pass = set()
            type_fail = set()
            other_pass = set()
            other_fail = set()
            
            failed_tests = defaultdict(int)
            for gid in gids:
                if gid in gtdb_tsp or gid in ncbi_tsp:
                    if pass_qc(quality_metadata[gid], 
                                min_comp,
                                max_cont,
                                min_quality,
                                sh_exception,
                                max_contigs,
                                min_N50,
                                max_ambiguous,
                                failed_tests):
                        type_pass.add(gid)
                    else:
                        type_fail.add(gid)
                else:
                    if pass_qc(quality_metadata[gid], 
                                min_comp,
                                max_cont,
                                min_quality,
                                sh_exception,
                                max_contigs,
                                min_N50,
                                max_ambiguous,
                                failed_tests):
                        other_pass.add(gid)
                    else:
                        other_fail.add(gid)

            if len(type_pass) >= 1:
                # great: one or more type genomes pass QC and will be selected as the type genome
                continue 
            
            if len(type_fail):
                # all potential type genomes for species failed QC so report these for manual inspection
                lost_type += 1
                for gid in type_fail:
                    fout_type_fail.write('%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d\t%s\t%s\n' % (
                                            sp,
                                            gid,
                                            type_metadata[gid].gtdb_type_designation_sources,
                                            type_metadata[gid].ncbi_type_material_designation,
                                            float(quality_metadata[gid].genome_size)/1e6,
                                            quality_metadata[gid].checkm_completeness,
                                            quality_metadata[gid].checkm_contamination,
                                            quality_metadata[gid].checkm_completeness-5*quality_metadata[gid].checkm_contamination,
                                            quality_metadata[gid].checkm_strain_heterogeneity,
                                            quality_metadata[gid].contig_count,
                                            quality_metadata[gid].n50_contigs,
                                            quality_metadata[gid].ambiguous_bases,
                                            excluded_from_refseq_note[gid],
                                            len(other_pass) == 0))
                
            if len(other_pass) == 0:
                # no genomes for species pass QC so report loss of species
                lost_sp += 1
                fout_sp_lost.write('%s\t%d\t%d' % (sp, len(gids), len(type_fail)))
                fout_sp_lost.write('\t%d\t%d\t%d\t%d\t%d\t%d\n' % (
                                    failed_tests['comp'],
                                    failed_tests['cont'],
                                    failed_tests['qual'],
                                    failed_tests['contig_count'],
                                    failed_tests['N50'],
                                    failed_tests['ambig']))
                                    
                for gid in type_fail.union(other_fail):
                    fout_fail.write('%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d\t%s\n' % (sp,
                                                                                            gid,
                                                                                            gid in type_fail,
                                                                                            float(quality_metadata[gid].genome_size)/1e6,
                                                                                            quality_metadata[gid].checkm_completeness,
                                                                                            quality_metadata[gid].checkm_contamination,
                                                                                            quality_metadata[gid].checkm_completeness-5*quality_metadata[gid].checkm_contamination,
                                                                                            quality_metadata[gid].checkm_strain_heterogeneity,
                                                                                            quality_metadata[gid].contig_count,
                                                                                            quality_metadata[gid].n50_contigs,
                                                                                            quality_metadata[gid].ambiguous_bases,
                                                                                            excluded_from_refseq_note[gid]))

        fout_type_fail.close()
        fout_fail.close()
        fout_sp_lost.close()

        self.logger.info('Identified %d species with type genomes failing QC and %d total species failing QC.' % (lost_type, lost_sp))
