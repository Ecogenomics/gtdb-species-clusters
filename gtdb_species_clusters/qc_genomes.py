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

from gtdb_species_clusters.genome_utils import (canonical_gid,
                                                exclude_from_refseq)

from gtdb_species_clusters.genomes import Genomes


class QcGenomes(object):
    """Quality check all potential GTDB genomes."""

    def __init__(self):
        """Initialization."""
        
        self.logger = logging.getLogger('timestamp')
        
    def read_marker_percentages(self, gtdb_domain_report, cur_genomes):
        """Parse percentage of marker genes for each genome."""
        
        marker_perc = {}
        with open(gtdb_domain_report, encoding='utf-8') as f:
            header = f.readline().rstrip().split('\t')
            
            bac_marker_perc_index = header.index('Bacterial Marker Percentage')
            ar_marker_perc_index = header.index('Archaeal Marker Percentage')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gid = canonical_gid(line_split[0])
                gid = cur_genomes.user_uba_id_map.get(gid, gid)
                bac_perc = float(line_split[bac_marker_perc_index])
                ar_perc = float(line_split[ar_marker_perc_index])

                marker_perc[gid] = max(bac_perc, ar_perc)
              
        return marker_perc
        
    def check_domain_assignments(self, gtdb_domain_report, cur_genomes, pass_qc_gids):
        """Check GTDB domain assignment."""
        
        with open(gtdb_domain_report, encoding='utf-8') as f:
            header = f.readline().rstrip().split('\t')
            
            domain_index = header.index('Predicted domain')
            bac_marker_perc_index = header.index('Bacterial Marker Percentage')
            ar_marker_perc_index = header.index('Archaeal Marker Percentage')
            ncbi_taxonomy_index = header.index('NCBI taxonomy')
            gtdb_taxonomy_index = header.index('GTDB taxonomy')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gid = canonical_gid(line_split[0])
                gid = cur_genomes.user_uba_id_map.get(gid, gid)
                
                if gid not in pass_qc_gids:
                    continue
                
                domain = line_split[domain_index]
                bac_perc = float(line_split[bac_marker_perc_index])
                ar_perc = float(line_split[ar_marker_perc_index])
                ncbi_domain = [t.strip() for t in line_split[ncbi_taxonomy_index].split(';')][0]
                gtdb_domain = [t.strip() for t in line_split[gtdb_taxonomy_index].split(';')][0]

                if not gid.startswith('U'):
                    if ncbi_domain != gtdb_domain and ncbi_domain != 'None':
                        print(f'[WARNING] NCBI ({ncbi_domain}) and GTDB ({gtdb_domain}) domains disagree in domain report (Bac = {bac_perc:.1f}%; Ar = {ar_perc:.1f}%): {gid}')

                    if domain != gtdb_domain and domain != 'None':
                        print(f'[WARNING] GTDB and predicted domain (Bac = {bac_perc:.1f}%; Ar = {ar_perc:.1f}%) disagree in domain report: {gid}')

    def run(self, 
                metadata_file,
                cur_uba_gid_file,
                ncbi_genbank_assembly_file,
                gtdb_domain_report,
                qc_exception_file,
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

        # create current GTDB genome sets
        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(metadata_file,
                                                create_sp_clusters=False,
                                                uba_genome_file=cur_uba_gid_file)
        self.logger.info(f' ...current genome set contains {len(cur_genomes):,} genomes.')

        # parse genomes flagged as exceptions from QC
        qc_exceptions = set()
        with open(qc_exception_file, encoding='utf-8') as f:
            f.readline()
            for line in f:
                gid = canonical_gid(line.split('\t')[0].strip())
                qc_exceptions.add(gid)
        self.logger.info(f'Identified {len(qc_exceptions):,} genomes flagged as exceptions from QC.')
        
        # get percentage of bac120 or ar122 marker genes
        marker_perc = self.read_marker_percentages(gtdb_domain_report, 
                                                    cur_genomes)

        # parse NCBI assembly files
        self.logger.info('Parsing NCBI assembly files.')
        excluded_from_refseq_note = exclude_from_refseq(ncbi_genbank_assembly_file)

        # QC all genomes
        self.logger.info('Validating genomes.')
        fout_retained = open(os.path.join(output_dir, 'qc_passed.tsv'), 'w')
        fout_failed = open(os.path.join(output_dir, 'qc_failed.tsv'), 'w')
        
        header = 'Accession\tNCBI species\tGTDB taxonomy'
        header += '\tCompleteness (%)\tContamination (%)\tQuality\tStrain heterogeneity at 100%'
        header += '\tMarkers (%)\tNo. contigs\tN50 contigs\tAmbiguous bases'
        
        fout_retained.write(header + '\tNote\n')
        fout_failed.write(header)
        fout_failed.write('\tFailed completeness\tFailed contamination\tFailed quality')
        fout_failed.write('\tFailed marker percentage\tFailed no. contigs\tFailed N50 contigs\tFailed ambiguous bases\n')

        pass_qc_gids = set()
        failed_qc_gids = set()
        for gid in cur_genomes:
            failed_tests = defaultdict(int)
            passed_qc = cur_genomes[gid].pass_qc(marker_perc[gid],
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
                pass_qc_gids.add(gid)
                fout_retained.write('%s\t%s\t%s' % (gid, cur_genomes[gid].ncbi_taxa.species, cur_genomes[gid].gtdb_taxa))
                fout_retained.write('\t%.2f\t%.2f\t%.2f\t%s\t%.2f\t%d\t%d\t%d\t%s\n' % (
                                        cur_genomes[gid].comp,
                                        cur_genomes[gid].cont,
                                        cur_genomes[gid].comp-5*cur_genomes[gid].cont,
                                        ('%.2f' % cur_genomes[gid].strain_heterogeneity_100) if cur_genomes[gid].strain_heterogeneity_100 else '-',
                                        marker_perc[gid],
                                        cur_genomes[gid].contig_count,
                                        cur_genomes[gid].contig_n50,
                                        cur_genomes[gid].ambiguous_bases,
                                        'Passed QC' if passed_qc else 'Flagged as exception'))
            else:
                failed_qc_gids.add(gid) 
                fout_failed.write('%s\t%s\t%s' % (gid, cur_genomes[gid].ncbi_taxa.species, cur_genomes[gid].gtdb_taxa))
                fout_failed.write('\t%.2f\t%.2f\t%.2f\t%s\t%.2f\t%d\t%d\t%d' % (
                                        cur_genomes[gid].comp,
                                        cur_genomes[gid].cont,
                                        cur_genomes[gid].comp-5*cur_genomes[gid].cont,
                                        ('%.2f' % cur_genomes[gid].strain_heterogeneity_100) if cur_genomes[gid].strain_heterogeneity_100 else '-',
                                        marker_perc[gid],
                                        cur_genomes[gid].contig_count,
                                        cur_genomes[gid].contig_n50,
                                        cur_genomes[gid].ambiguous_bases))
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
        
        self.logger.info('Retained {:,} ({:.2f}%) genomes and filtered {:,} ({:.2f}%) genomes.'.format(
                            len(pass_qc_gids),
                            len(pass_qc_gids)*100.0/len(cur_genomes),
                            len(failed_qc_gids),
                            len(failed_qc_gids)*100.0/len(cur_genomes)))
        
        # check domain assignment of genomes passing QC
        # report potential issues
        self.check_domain_assignments(gtdb_domain_report, 
                                        cur_genomes,
                                        pass_qc_gids)
                                                                
        # QC genomes in each named species
        named_ncbi_species = cur_genomes.named_ncbi_species()
        self.logger.info(f'Performing QC of type genome for each of the {len(named_ncbi_species):,} NCBI species.')
        
        fout_type_fail = open(os.path.join(output_dir, 'type_genomes_fail_qc.tsv'), 'w')
        fout_type_fail.write('NCBI species\tAccession\tGTDB taxonomy\tNCBI taxonomy\tType sources\tNCBI assembly type\tGenome size (bp)')
        fout_type_fail.write('\tCompleteness (%)\tContamination (%)\tQuality\tStrain heterogeneity at 100%')
        fout_type_fail.write('\tMarkers (%)\tNo. contigs\tN50 contigs\tAmbiguous bases\tNCBI exclude from RefSeq\tLost species\n')
        
        fout_fail_sp = open(os.path.join(output_dir, 'species_fail_qc.tsv'), 'w')
        fout_fail_sp.write('NCBI species\tAccession\tGTDB taxonomy\tNCBI taxonomy\tAssembled from type material\tGenome size (bp)')
        fout_fail_sp.write('\tCompleteness (%)\tContamination (%)\tQuality\tStrain heterogeneity at 100%')
        fout_fail_sp.write('\tMarkers (%)\tNo. contigs\tN50 contigs\tAmbiguous bases')
        fout_fail_sp.write('\tFailed completeness\tFailed contamination\tFailed quality')
        fout_fail_sp.write('\tFailed marker percentage\tFailed no. contigs\tFailed N50 contigs\tFailed ambiguous bases')
        fout_fail_sp.write('\tNCBI exclude from RefSeq\n')
        
        fout_sp_lost = open(os.path.join(output_dir, 'species_lost.tsv'), 'w')
        fout_sp_lost.write('NCBI species\tNo. genomes\tNo. type genomes')
        fout_sp_lost.write('\tFail completeness\tFail contamination\tFail quality\tFailed percent markers')
        fout_sp_lost.write('\tFail no. contigs\tFail N50 contigs\tFail ambiguous bases\n')
        
        lost_type = 0
        lost_sp = 0
        filtered_genomes = 0
        failed_tests_cumulative = defaultdict(int)
        for sp, gids in named_ncbi_species.items():
            type_pass = set()
            type_fail = set()
            other_pass = set()
            other_fail = set()
            
            failed_tests_gids = {}
            for gid in gids:
                failed_tests = defaultdict(int)
                passed_qc = cur_genomes[gid].pass_qc(marker_perc[gid],
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

                if cur_genomes[gid].is_gtdb_type_strain() or cur_genomes[gid].is_ncbi_type_strain():
                    if passed_qc or gid in qc_exceptions:
                        type_pass.add(gid)
                    else:
                        type_fail.add(gid)
                        filtered_genomes += 1
                else:
                    if passed_qc or gid in qc_exceptions:
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
                                            cur_genomes[gid].gtdb_taxa,
                                            cur_genomes[gid].ncbi_taxa,
                                            cur_genomes[gid].gtdb_type_designation_sources,
                                            cur_genomes[gid].ncbi_type_material,
                                            float(cur_genomes[gid].length)/1e6,
                                            cur_genomes[gid].comp,
                                            cur_genomes[gid].cont,
                                            cur_genomes[gid].comp-5*cur_genomes[gid].cont,
                                            cur_genomes[gid].strain_heterogeneity_100,
                                            marker_perc[gid],
                                            cur_genomes[gid].contig_count,
                                            cur_genomes[gid].contig_n50,
                                            cur_genomes[gid].ambiguous_bases,
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
                                            cur_genomes[gid].gtdb_taxa,
                                            cur_genomes[gid].ncbi_taxa,
                                            gid in type_fail,
                                            float(cur_genomes[gid].length)/1e6,
                                            cur_genomes[gid].comp,
                                            cur_genomes[gid].cont,
                                            cur_genomes[gid].comp-5*cur_genomes[gid].cont,
                                            cur_genomes[gid].strain_heterogeneity_100,
                                            marker_perc[gid],
                                            cur_genomes[gid].contig_count,
                                            cur_genomes[gid].contig_n50,
                                            cur_genomes[gid].ambiguous_bases))
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
        
        self.logger.info(f'Filtered {filtered_genomes:,} genomes assigned to NCBI species.')
        self.logger.info(f'Identified {lost_type:,} species with type genomes failing QC and {lost_sp:,} total species failing QC.')
        self.logger.info('Genomes from NCBI species filtered by each criterion:')
        for test in sorted(failed_tests_cumulative):
            self.logger.info(f'{test}: {failed_tests_cumulative[test]:,}')
