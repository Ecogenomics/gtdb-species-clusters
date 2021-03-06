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
import logging
from collections import defaultdict

from gtdb_species_clusters.genome_utils import (canonical_gid,
                                                exclude_from_refseq)

from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.genome_utils import same_assembly_version


class QcGenomes():
    """Quality check all potential GTDB genomes."""

    def __init__(self, output_dir):
        """Initialization."""

        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')

    def parse_marker_percentages(self, gtdb_domain_report):
        """Parse percentage of marker genes for each genome."""

        marker_perc = {}
        with open(gtdb_domain_report, encoding='utf-8') as f:
            header = f.readline().rstrip().split('\t')

            bac_marker_perc_index = header.index('Bacterial Marker Percentage')
            ar_marker_perc_index = header.index('Archaeal Marker Percentage')

            for line in f:
                line_split = line.strip().split('\t')

                gid = canonical_gid(line_split[0])
                bac_perc = float(line_split[bac_marker_perc_index])
                ar_perc = float(line_split[ar_marker_perc_index])

                marker_perc[gid] = max(bac_perc, ar_perc)

        return marker_perc

    def check_domain_assignments(self, gtdb_domain_report, passed_qc_gids):
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
                if gid not in passed_qc_gids:
                    continue

                domain = line_split[domain_index]
                bac_perc = float(line_split[bac_marker_perc_index])
                ar_perc = float(line_split[ar_marker_perc_index])
                ncbi_domain = [t.strip()
                               for t in line_split[ncbi_taxonomy_index].split(';')][0]
                gtdb_domain = [t.strip()
                               for t in line_split[gtdb_taxonomy_index].split(';')][0]

                if not gid.startswith('U'):
                    if ncbi_domain != gtdb_domain and ncbi_domain != 'None':
                        self.logger.warning(
                            f'NCBI ({ncbi_domain}) and GTDB ({gtdb_domain}) domains disagree in domain report (Bac = {bac_perc:.1f}%; Ar = {ar_perc:.1f}%): {gid}')

                    if domain != gtdb_domain and domain != 'None':
                        self.logger.warning(
                            f'GTDB and predicted domain (Bac = {bac_perc:.1f}%; Ar = {ar_perc:.1f}%) disagree in domain report: {gid}')

    def parse_qc_exception_file(self, qc_exception_file):
        """Parse file indicating genomes flagged as exceptions from QC."""

        qc_exceptions = set()
        with open(qc_exception_file, encoding='utf-8') as f:
            f.readline()
            for line in f:
                gid = canonical_gid(line.split('\t')[0].strip())
                qc_exceptions.add(gid)

        return qc_exceptions

    def qc_genomes(self,
                   cur_genomes,
                   marker_perc,
                   qc_exceptions,
                   excluded_from_refseq_note,
                   min_comp,
                   max_cont,
                   min_quality,
                   sh_exception,
                   min_perc_markers,
                   max_contigs,
                   min_N50,
                   max_ambiguous):
        """Apply quality-control to all genomes."""

        fout_passed = open(os.path.join(self.output_dir, 'qc_passed.tsv'), 'w')
        fout_failed = open(os.path.join(self.output_dir, 'qc_failed.tsv'), 'w')

        header = 'Accession\tNCBI species\tGTDB taxonomy'
        header += '\tCompleteness (%)\tContamination (%)\tQuality\tStrain heterogeneity at 100%'
        header += '\tMarkers (%)\tNo. contigs\tN50 contigs\tAmbiguous bases'

        fout_passed.write(
            header + '\tScore\tNote\tNCBI exclude from RefSeq\n')
        fout_failed.write(header)
        fout_failed.write(
            '\tFailed completeness\tFailed contamination\tFailed quality')
        fout_failed.write(
            '\tFailed marker percentage\tFailed no. contigs\tFailed N50 contigs')
        fout_failed.write(
            '\tFailed ambiguous bases\tNCBI exclude from RefSeq\n')

        passed_qc_gids = set()
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
                passed_qc_gids.add(gid)
                fout_passed.write('{}\t{}\t{}'.format(
                    gid, cur_genomes[gid].ncbi_taxa.species, cur_genomes[gid].gtdb_taxa))
                fout_passed.write('\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{:.2f}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\n'.format(
                    cur_genomes[gid].comp,
                    cur_genomes[gid].cont,
                    cur_genomes[gid].comp - 5*cur_genomes[gid].cont,
                    f'{cur_genomes[gid].strain_heterogeneity_100:.2f}' if cur_genomes[gid].strain_heterogeneity_100 else '-',
                    marker_perc[gid],
                    cur_genomes[gid].contig_count,
                    cur_genomes[gid].contig_n50,
                    cur_genomes[gid].ambiguous_bases,
                    cur_genomes[gid].score_type_strain(),
                    'Passed QC' if passed_qc else 'Flagged as exception',
                    excluded_from_refseq_note.get(gid, '')))
            else:
                failed_qc_gids.add(gid)
                fout_failed.write('{}\t{}\t{}'.format(
                    gid, cur_genomes[gid].ncbi_taxa.species, cur_genomes[gid].gtdb_taxa))
                fout_failed.write('\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{:.2f}\t{}\t{}\t{}'.format(
                    cur_genomes[gid].comp,
                    cur_genomes[gid].cont,
                    cur_genomes[gid].comp-5*cur_genomes[gid].cont,
                    f'{cur_genomes[gid].strain_heterogeneity_100:.2f}' if cur_genomes[gid].strain_heterogeneity_100 else '-',
                    marker_perc[gid],
                    cur_genomes[gid].contig_count,
                    cur_genomes[gid].contig_n50,
                    cur_genomes[gid].ambiguous_bases))
                fout_failed.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    failed_tests['comp'],
                    failed_tests['cont'],
                    failed_tests['qual'],
                    failed_tests['marker_perc'],
                    failed_tests['contig_count'],
                    failed_tests['N50'],
                    failed_tests['ambig'],
                    excluded_from_refseq_note.get(gid, '')))
        fout_passed.close()
        fout_failed.close()

        self.logger.info('Retained {:,} ({:.2f}%) genomes and filtered {:,} ({:.2f}%) genomes.'.format(
            len(passed_qc_gids),
            len(passed_qc_gids)*100.0/len(cur_genomes),
            len(failed_qc_gids),
            len(failed_qc_gids)*100.0/len(cur_genomes)))

        return passed_qc_gids, failed_qc_gids

    def check_qc_of_ncbi_species(self,
                                 cur_genomes,
                                 marker_perc,
                                 qc_exceptions,
                                 excluded_from_refseq_note,
                                 min_comp,
                                 max_cont,
                                 min_quality,
                                 sh_exception,
                                 min_perc_markers,
                                 max_contigs,
                                 min_N50,
                                 max_ambiguous):
        """Report impact of QC filtering on NCBI species."""

        # These files can be inspected to determine species where
        # all avaliable genomes assemblies have been filtered, and
        # why they were filtered. These are good candidated for
        # having a genome selected that is an exception to QC, but
        # this must be balanced with ensuring assemblies are of
        # sufficient quality as to not negatively impact tree
        # inference.

        named_ncbi_species = cur_genomes.named_ncbi_species()
        self.logger.info(
            f'Performing QC of type genome for each of the {len(named_ncbi_species):,} NCBI species.')

        fout_type_fail = open(os.path.join(
            self.output_dir, 'type_genomes_fail_qc.tsv'), 'w')
        fout_type_fail.write(
            'NCBI species\tAccession\tGTDB taxonomy\tNCBI taxonomy\tType sources\tNCBI assembly type\tGenome size (bp)')
        fout_type_fail.write(
            '\tCompleteness (%)\tContamination (%)\tQuality\tStrain heterogeneity at 100%')
        fout_type_fail.write(
            '\tMarkers (%)\tNo. contigs\tN50 contigs\tAmbiguous bases\tNCBI exclude from RefSeq\tLost species\n')

        fout_fail_sp = open(os.path.join(
            self.output_dir, 'species_fail_qc.tsv'), 'w')
        fout_fail_sp.write(
            'NCBI species\tAccession\tGTDB taxonomy\tNCBI taxonomy\tAssembled from type material\tGenome size (bp)')
        fout_fail_sp.write(
            '\tCompleteness (%)\tContamination (%)\tQuality\tStrain heterogeneity at 100%')
        fout_fail_sp.write(
            '\tMarkers (%)\tNo. contigs\tN50 contigs\tAmbiguous bases')
        fout_fail_sp.write(
            '\tFailed completeness\tFailed contamination\tFailed quality')
        fout_fail_sp.write(
            '\tFailed marker percentage\tFailed no. contigs\tFailed N50 contigs\tFailed ambiguous bases')
        fout_fail_sp.write('\tNCBI exclude from RefSeq\n')

        fout_sp_lost = open(os.path.join(
            self.output_dir, 'species_lost.tsv'), 'w')
        fout_sp_lost.write('NCBI species\tNo. genomes\tNo. type genomes')
        fout_sp_lost.write(
            '\tFail completeness\tFail contamination\tFail quality\tFailed percent markers')
        fout_sp_lost.write(
            '\tFail no. contigs\tFail N50 contigs\tFail ambiguous bases\n')

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

            if type_fail:
                # all potential type genomes for species failed QC so report these for manual inspection
                lost_type += 1
                for gid in type_fail:
                    fout_type_fail.write('{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{}\t{}\t{}\t{}\n'.format(
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
                fout_sp_lost.write('{}\t{}\t{}'.format(
                    sp, len(gids), len(type_fail)))
                fout_sp_lost.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    sum([failed_tests_gids[gid]['comp'] for gid in gids]),
                    sum([failed_tests_gids[gid]['cont'] for gid in gids]),
                    sum([failed_tests_gids[gid]['qual'] for gid in gids]),
                    sum([failed_tests_gids[gid]['marker_perc']
                         for gid in gids]),
                    sum([failed_tests_gids[gid]['contig_count']
                         for gid in gids]),
                    sum([failed_tests_gids[gid]['N50'] for gid in gids]),
                    sum([failed_tests_gids[gid]['ambig'] for gid in gids])))

                for gid in type_fail.union(other_fail):
                    fout_fail_sp.write('{}\t{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{}\t{}'.format(
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
                    fout_fail_sp.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                        failed_tests_gids[gid]['comp'],
                        failed_tests_gids[gid]['cont'],
                        failed_tests_gids[gid]['qual'],
                        failed_tests_gids[gid]['marker_perc'],
                        failed_tests_gids[gid]['contig_count'],
                        failed_tests_gids[gid]['N50'],
                        failed_tests_gids[gid]['ambig']))
                    fout_fail_sp.write('\t{}\n'.format(
                        excluded_from_refseq_note[gid]))

        fout_type_fail.close()
        fout_fail_sp.close()
        fout_sp_lost.close()

        self.logger.info(
            f'Filtered {filtered_genomes:,} genomes assigned to NCBI species.')
        self.logger.info(
            f'Identified {lost_type:,} species with type genomes failing QC and {lost_sp:,} total species failing QC.')
        self.logger.info(
            'Genomes from NCBI species filtered by each criterion:')
        for test in sorted(failed_tests_cumulative):
            self.logger.info(f'{test}: {failed_tests_cumulative[test]:,}')

    def run(self,
            prev_gtdb_metadata_file,
            cur_gtdb_metadata_file,
            ncbi_genbank_assembly_file,
            gtdb_domain_report,
            gtdb_type_strains_ledger,
            qc_exception_file,
            ncbi_env_bioproject_ledger,
            min_comp,
            max_cont,
            min_quality,
            sh_exception,
            min_perc_markers,
            max_contigs,
            min_N50,
            max_ambiguous):
        """Quality check all potential GTDB genomes."""

        # create previous and current GTDB genome sets
        self.logger.info('Creating previous GTDB genome set.')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                             gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                             ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                             ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)
        self.logger.info(
            f' - previous genome set contains {len(prev_genomes):,} genomes.')
        self.logger.info(' - previous genome set has {:,} species clusters spanning {:,} genomes.'.format(
            len(prev_genomes.sp_clusters),
            prev_genomes.sp_clusters.total_num_genomes()))

        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                            gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                            create_sp_clusters=False,
                                            ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                            ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)

        # parse genomes flagged as exceptions from QC
        qc_exceptions = self.parse_qc_exception_file(qc_exception_file)
        self.logger.info(
            f'Identified {len(qc_exceptions):,} genomes flagged as exceptions from QC.')

        # get percentage of bac120 or ar122 marker genes
        marker_perc = self.parse_marker_percentages(gtdb_domain_report)

        # parse NCBI assembly files
        self.logger.info('Parsing NCBI assembly files.')
        excluded_from_refseq_note = exclude_from_refseq(
            ncbi_genbank_assembly_file)

        # QC all genomes
        self.logger.info('Validating genomes.')
        passed_qc_gids, failed_qc_gids = self.qc_genomes(
            cur_genomes,
            marker_perc,
            qc_exceptions,
            excluded_from_refseq_note,
            min_comp,
            max_cont,
            min_quality,
            sh_exception,
            min_perc_markers,
            max_contigs,
            min_N50,
            max_ambiguous)

        # check domain assignment of genomes passing QC
        # and report potential issues
        self.check_domain_assignments(gtdb_domain_report,
                                      passed_qc_gids)

        # report results of QC on genomes from each NCBI species
        self.check_qc_of_ncbi_species(
            cur_genomes,
            marker_perc,
            qc_exceptions,
            excluded_from_refseq_note,
            min_comp,
            max_cont,
            min_quality,
            sh_exception,
            min_perc_markers,
            max_contigs,
            min_N50,
            max_ambiguous)

        # sanity check QC results by identifying any genomes that passed QC last release, but
        # have now been flagged as failing QC. This should rarely, if ever, happen unless the
        # genomic data of the assembly has been updated.
        unexpected_qc_fail = []
        for gid in prev_genomes:
            if gid in cur_genomes:
                if not same_assembly_version(prev_genomes[gid].ncbi_accn, cur_genomes[gid].ncbi_accn):
                    # genome assembly has changed so QC status is not expected to be the same
                    continue

                if gid in failed_qc_gids:
                    unexpected_qc_fail.append(gid)

        if len(unexpected_qc_fail) > 0:
            self.logger.warning('Identified {:,} genomes that passed QC in previous GTDB release, that failed QC in this release.'.format(
                len(unexpected_qc_fail)))
            self.logger.warning(
                ' - examples: {}'.format(','.join(unexpected_qc_fail[0:10])))
