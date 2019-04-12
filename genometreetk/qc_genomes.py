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
                                    read_gtdb_ncbi_taxonomy,
                                    read_gtdb_taxonomy)
                                    
from genometreetk.type_genome_utils import (exclude_from_refseq, 
                                            ncbi_type_strain_of_species,
                                            gtdb_type_strain_of_species,
                                            parse_marker_percentages,
                                            pass_qc)


class QcGenomes(object):
    """Quality check all potential GTDB genomes."""

    def __init__(self):
        """Initialization."""
        
        self.logger = logging.getLogger('timestamp')
        
        self.true_str = ['t', 'T', 'true', 'True']
        
    def _gtdb_user_genomes(self, gtdb_user_genomes_file, metadata_file):
        """Get map between GTDB User genomes and GenBank accessions."""
        
        uba_to_genbank = {}
        for line in open(gtdb_user_genomes_file):
            line_split = line.strip().split('\t')
            gb_acc = 'GB_' + line_split[0]
            uba_id = line_split[4]
            uba_to_genbank[uba_id] = gb_acc
        
        user_to_genbank = {}
        m = read_gtdb_metadata(metadata_file, ['organism_name'])
        for gid, metadata in m.items():
            if '(UBA' in str(metadata.organism_name):
                uba_id = metadata.organism_name[metadata.organism_name.find('(')+1:-1]
                if uba_id in uba_to_genbank:
                    if uba_to_genbank[uba_id] in m:
                        # use only NCBI accessioned version of genome
                        continue
                    else:
                        user_to_genbank[gid] = uba_to_genbank[uba_id]

        return user_to_genbank

    def run(self, metadata_file,
                gtdb_user_genomes_file,
                gtdb_user_reps,
                ncbi_refseq_assembly_file,
                ncbi_genbank_assembly_file,
                gtdb_domain_report,
                qc_exception_file,
                species_exception_file,
                min_comp,
                max_cont,
                min_quality,
                sh_exception,
                min_perc_markers,
                max_contigs,
                min_N50,
                max_ambiguous,
                output_dir):
        """Quality check all potential GTDB genomes."""

        # get GTDB and NCBI taxonomy strings for each genome
        self.logger.info('Reading NCBI taxonomy from GTDB metadata file.')
        ncbi_taxonomy, ncbi_update_count = read_gtdb_ncbi_taxonomy(metadata_file, species_exception_file)
        ncbi_species = binomial_species(ncbi_taxonomy)
        gtdb_taxonomy = read_gtdb_taxonomy(metadata_file)
        self.logger.info('Read NCBI taxonomy for %d genomes with %d manually defined updates.' % (len(ncbi_taxonomy), ncbi_update_count))
        self.logger.info('Read GTDB taxonomy for %d genomes.' % len(gtdb_taxonomy))
        
        # determine User genomes to retain for consideration
        gtdb_user_to_genbank = self._gtdb_user_genomes(gtdb_user_genomes_file, metadata_file)
        self.logger.info('Identified %d GTDB User genomes with GenBank accessions to retain for potential inclusion in GTDB.' % len(gtdb_user_to_genbank))
        
        user_genomes = 0
        for line in open(gtdb_user_reps):
            line_split = line.strip().split('\t')
            gid, taxonomy = line_split
            if gid not in gtdb_user_to_genbank:
                if 'd__Bacteria' in taxonomy:
                    self.logger.warning('Bacterial genome %s has no NCBI accession and is being skipped.' % gid)
                else:
                    gtdb_user_to_genbank[gid] = gid
                    user_genomes += 1
        self.logger.info('Identified %d archaeal GTDB User genome WITHOUT GenBank accessions to retain for potential inclusion in GTDB.' % user_genomes)

        # parse genomes flagged as exceptions from QC
        qc_exceptions = set()
        for line in open(qc_exception_file):
            qc_exceptions.add(line.split('\t')[0].strip())
        self.logger.info('Identified %d genomes flagged as exceptions from QC.' % len(qc_exceptions))
        
        # calculate quality score for genomes
        self.logger.info('Parsing QC statistics for each genome.')
        quality_metadata = read_gtdb_metadata(metadata_file, ['checkm_completeness',
                                                                'checkm_contamination',
                                                                'checkm_strain_heterogeneity_100',
                                                                'contig_count',
                                                                'n50_contigs',
                                                                'ambiguous_bases',
                                                                'genome_size'])
                                                                
        marker_perc = parse_marker_percentages(gtdb_domain_report)
                                                                
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
        
        # QC all genomes
        self.logger.info('Validating genomes.')
        fout_retained = open(os.path.join(output_dir, 'qc_passed.tsv'), 'w')
        fout_failed = open(os.path.join(output_dir, 'qc_failed.tsv'), 'w')
        
        header = 'Accession\tNCBI species'
        header += '\tCompleteness (%)\tContamination (%)\tQuality\tStrain heterogeneity at 100%'
        header += '\tMarkers (%)\tNo. contigs\tN50 contigs\tAmbiguous bases'
        
        fout_retained.write(header + '\tNote\n')
        fout_failed.write(header)
        fout_failed.write('\tFailed completeness\tFailed contamination\tFailed quality')
        fout_failed.write('\tFailed marker percentage\tFailed no. contigs\tFailed N50 contigs\tFailed ambiguous bases\n')

        num_retained = 0
        num_filtered = 0
        for gid in quality_metadata:
            if gid.startswith('U_') and gid not in gtdb_user_to_genbank:
                # skip user genomes not marked for retention
                continue

            failed_tests = defaultdict(int)
            passed_qc = pass_qc(quality_metadata[gid], 
                                    marker_perc[gid],
                                    min_comp,
                                    max_cont,
                                    min_quality,
                                    sh_exception,
                                    min_perc_markers,
                                    max_contigs,
                                    min_N50,
                                    max_ambiguous,
                                    failed_tests)

            if passed_qc or gid in qc_exceptions:
                num_retained += 1
                fout_retained.write('%s\t%s' % (gid, ncbi_taxonomy[gid][6]))
                fout_retained.write('\t%.2f\t%.2f\t%.2f\t%s\t%.2f\t%d\t%d\t%d\t%s\n' % (
                                        quality_metadata[gid].checkm_completeness,
                                        quality_metadata[gid].checkm_contamination,
                                        quality_metadata[gid].checkm_completeness-5*quality_metadata[gid].checkm_contamination,
                                        ('%.2f' % quality_metadata[gid].checkm_strain_heterogeneity_100) if quality_metadata[gid].checkm_strain_heterogeneity_100 else '-',
                                        marker_perc[gid],
                                        quality_metadata[gid].contig_count,
                                        quality_metadata[gid].n50_contigs,
                                        quality_metadata[gid].ambiguous_bases,
                                        'Passed QC' if passed_qc else 'Flagged as exception'))
            else:
                num_filtered += 1 
                fout_failed.write('%s\t%s' % (gid, ncbi_taxonomy[gid][6]))
                fout_failed.write('\t%.2f\t%.2f\t%.2f\t%s\t%.2f\t%d\t%d\t%d' % (
                                        quality_metadata[gid].checkm_completeness,
                                        quality_metadata[gid].checkm_contamination,
                                        quality_metadata[gid].checkm_completeness-5*quality_metadata[gid].checkm_contamination,
                                        ('%.2f' % quality_metadata[gid].checkm_strain_heterogeneity_100) if quality_metadata[gid].checkm_strain_heterogeneity_100 else '-',
                                        marker_perc[gid],
                                        quality_metadata[gid].contig_count,
                                        quality_metadata[gid].n50_contigs,
                                        quality_metadata[gid].ambiguous_bases))
                fout_failed.write('\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n' % (
                                    failed_tests['comp'],
                                    failed_tests['cont'],
                                    failed_tests['qual'],
                                    failed_tests['marker_perc'],
                                    failed_tests['contig_count'],
                                    failed_tests['N50'],
                                    failed_tests['ambig']))
        fout_retained.close()
        fout_failed.close()
        
        self.logger.info('Retained %d genomes and filtered %d genomes.' % (num_retained, num_filtered))
                                                                
        # QC genomes in each named species
        self.logger.info('Performing QC of type genome for each of the %d NCBI species.' % len(ncbi_species))
        
        fout_type_fail = open(os.path.join(output_dir, 'type_genomes_fail_qc.tsv'), 'w')
        fout_type_fail.write('Species\tAccession\tGTDB taxonomy\tNCBI taxonomy\tType sources\tNCBI assembly type\tGenome size (bp)')
        fout_type_fail.write('\tCompleteness (%)\tContamination (%)\tQuality\tStrain heterogeneity at 100%')
        fout_type_fail.write('\tMarkers (%)\tNo. contigs\tN50 contigs\tAmbiguous bases\tNCBI exclude from RefSeq\tLost species\n')
        
        fout_fail_sp = open(os.path.join(output_dir, 'species_fail_qc.tsv'), 'w')
        fout_fail_sp.write('Species\tAccession\tGTDB taxonomy\tNCBI taxonomy\tAssembled from type material\tGenome size (bp)')
        fout_fail_sp.write('\tCompleteness (%)\tContamination (%)\tQuality\tStrain heterogeneity at 100%')
        fout_fail_sp.write('\tMarkers (%)\tNo. contigs\tN50 contigs\tAmbiguous bases')
        fout_fail_sp.write('\tFailed completeness\tFailed contamination\tFailed quality')
        fout_fail_sp.write('\tFailed marker percentage\tFailed no. contigs\tFailed N50 contigs\tFailed ambiguous bases')
        fout_fail_sp.write('\tNCBI exclude from RefSeq\n')
        
        fout_sp_lost = open(os.path.join(output_dir, 'species_lost.tsv'), 'w')
        fout_sp_lost.write('Species\tNo. genomes\tNo. type genomes')
        fout_sp_lost.write('\tFail completeness\tFail contamination\tFail quality\tFailed percent markers')
        fout_sp_lost.write('\tFail no. contigs\tFail N50 contigs\tFail ambiguous bases\n')
        
        lost_type = 0
        lost_sp = 0
        filtered_genomes = 0
        failed_tests_cumulative = defaultdict(int)
        for sp, gids in ncbi_species.items():
            type_pass = set()
            type_fail = set()
            other_pass = set()
            other_fail = set()
            
            failed_tests_gids = {}
            for gid in gids:
                failed_tests = defaultdict(int)
                passed_qc = pass_qc(quality_metadata[gid], 
                                    marker_perc[gid],
                                    min_comp,
                                    max_cont,
                                    min_quality,
                                    sh_exception,
                                    min_perc_markers,
                                    max_contigs,
                                    min_N50,
                                    max_ambiguous,
                                    failed_tests)
                                    
                failed_tests_gids[gid] = failed_tests

                if gid in gtdb_tsp or gid in ncbi_tsp:
                    if passed_qc:
                        type_pass.add(gid)
                    else:
                        type_fail.add(gid)
                        filtered_genomes += 1
                else:
                    if passed_qc:
                        other_pass.add(gid)
                    else:
                        other_fail.add(gid)
                        filtered_genomes += 1
                        
            # tally failed species
            for test, count in failed_tests.items():
                failed_tests_cumulative[test] += count

            if len(type_pass) >= 1:
                # great: one or more type genomes pass QC and will be selected as the type genome
                continue 
            
            if len(type_fail):
                # all potential type genomes for species failed QC so report these for manual inspection
                lost_type += 1
                for gid in type_fail:
                    fout_type_fail.write('%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d\t%s\t%s\n' % (
                                            sp,
                                            gid,
                                            '; '.join(gtdb_taxonomy[gid]),
                                            '; '.join(ncbi_taxonomy[gid]),
                                            type_metadata[gid].gtdb_type_designation_sources,
                                            type_metadata[gid].ncbi_type_material_designation,
                                            float(quality_metadata[gid].genome_size)/1e6,
                                            quality_metadata[gid].checkm_completeness,
                                            quality_metadata[gid].checkm_contamination,
                                            quality_metadata[gid].checkm_completeness-5*quality_metadata[gid].checkm_contamination,
                                            quality_metadata[gid].checkm_strain_heterogeneity_100,
                                            marker_perc[gid],
                                            quality_metadata[gid].contig_count,
                                            quality_metadata[gid].n50_contigs,
                                            quality_metadata[gid].ambiguous_bases,
                                            excluded_from_refseq_note[gid],
                                            len(other_pass) == 0))
                
            if len(other_pass) == 0:
                # no genomes for species pass QC so report loss of species
                lost_sp += 1
                fout_sp_lost.write('%s\t%d\t%d' % (sp, len(gids), len(type_fail)))
                fout_sp_lost.write('\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n' % (
                                    sum([failed_tests_gids[gid]['comp'] for gid in gids]),
                                    sum([failed_tests_gids[gid]['cont'] for gid in gids]),
                                    sum([failed_tests_gids[gid]['qual'] for gid in gids]),
                                    sum([failed_tests_gids[gid]['marker_perc'] for gid in gids]),
                                    sum([failed_tests_gids[gid]['contig_count'] for gid in gids]),
                                    sum([failed_tests_gids[gid]['N50'] for gid in gids]),
                                    sum([failed_tests_gids[gid]['ambig'] for gid in gids])))
                                    
                for gid in type_fail.union(other_fail):
                    fout_fail_sp.write('%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d' % (
                                            sp,
                                            gid,
                                            '; '.join(gtdb_taxonomy[gid]),
                                            '; '.join(ncbi_taxonomy[gid]),
                                            gid in type_fail,
                                            float(quality_metadata[gid].genome_size)/1e6,
                                            quality_metadata[gid].checkm_completeness,
                                            quality_metadata[gid].checkm_contamination,
                                            quality_metadata[gid].checkm_completeness-5*quality_metadata[gid].checkm_contamination,
                                            quality_metadata[gid].checkm_strain_heterogeneity_100,
                                            marker_perc[gid],
                                            quality_metadata[gid].contig_count,
                                            quality_metadata[gid].n50_contigs,
                                            quality_metadata[gid].ambiguous_bases))
                    fout_fail_sp.write('\t%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                                        failed_tests_gids[gid]['comp'],
                                        failed_tests_gids[gid]['cont'],
                                        failed_tests_gids[gid]['qual'],
                                        failed_tests_gids[gid]['marker_perc'],
                                        failed_tests_gids[gid]['contig_count'],
                                        failed_tests_gids[gid]['N50'],
                                        failed_tests_gids[gid]['ambig']))
                    fout_fail_sp.write('\t%s\n' % excluded_from_refseq_note[gid])

        fout_type_fail.close()
        fout_fail_sp.close()
        fout_sp_lost.close()
        
        self.logger.info('Genomes filtered for each criterion:')
        for test in sorted(failed_tests_cumulative):
            self.logger.info('%s: %d' % (test, failed_tests_cumulative[test]))

        self.logger.info('Filtered %d genomes assigned to NCBI species.' % filtered_genomes)
        self.logger.info('Identified %d species with type genomes failing QC and %d total species failing QC.' % (lost_type, lost_sp))
