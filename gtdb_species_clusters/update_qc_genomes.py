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
from dataclasses import dataclass
from typing import Dict, Set, Tuple

from gtdblib.util.bio.accession import canonical_gid, is_same_accn_version

from gtdb_species_clusters.genome_utils import exclude_from_refseq
from gtdb_species_clusters.genomes import Genomes


@dataclass
class QcCriteria:
    """Criteria for passing QC.

    :param min_comp: Minimum estimated genome completeness to pass QC.
    :param max_cont: Maximum estimated genome contamination to pass QC.
    :param min_quality: Minimum estimated genome quality to pass QC, defined as completeness - 5*contamination.
    :param sh_exception: Minimum strain heterogenity to retain genomes with upto 20% contamination.
    :param min_perc_markers: Minimum per4centage of marker genes to pass QC.
    :param max_contigs: Maximum number of contigs to pass QC.
    :param min_N50: Minimum N50 to pass QC.
    :param max_ambiguous: Maximum number of ambiguous bases to pass QC.
    """
    min_comp: float
    max_cont: float
    min_quality: float
    sh_exception: float
    min_perc_markers: float
    max_contigs: int
    min_N50: int
    max_ambiguous: int


class QcGenomes():
    """Quality check all potential GTDB genomes."""

    def __init__(self, output_dir: str):
        """Initialization.

        :param output_dir: Output directory.
        """

        self.output_dir = output_dir
        self.log = logging.getLogger('rich')

    def parse_marker_percentages(self, gtdb_domain_report: str) -> Dict[str, float]:
        """Parse percentage of marker genes for each genome.

        :param gtdb_domain_report: File indicating percent of bac120 and ar53 gene sets identified in each genome.

        :return: Percent of bac120 or ar53 markers, whichever is larger, identified in each genome.
        """

        marker_perc = {}
        with open(gtdb_domain_report, encoding='utf-8') as f:
            header = f.readline().rstrip().split('\t')

            bac_marker_perc_index = header.index('Bac120 Markers (%)')
            ar_marker_perc_index = header.index('Ar53 Markers (%)')

            for line in f:
                tokens = line.strip().split('\t')

                gid = canonical_gid(tokens[0])
                bac_perc = float(tokens[bac_marker_perc_index])
                ar_perc = float(tokens[ar_marker_perc_index])

                marker_perc[gid] = max(bac_perc, ar_perc)

        return marker_perc

    def check_domain_assignments(self, gtdb_domain_report: str, passed_qc_gids: Set[str], cur_genomes: Genomes) -> None:
        """Check GTDB domain assignments and report incongruencies and errors.

        This check is performed as we have previously observed genomes at NCBI that
        are classified to the incorrect domain (e.g. a bacterium classified as Archaea).
        Genomes should have a GTDB domain assignment that is congruent with the domain marker
        set with the largest percentage of identified genes. This GTDB domain assignment is
        handled by the GTDB database codebase, but was previously observed to be in error so we
        now do an explicit check for this case and report instances where the NCBI and GTDB
        domain assignments are in conflict so this can be manually inspected.

        :param gtdb_domain_report: File indicating percent of bac120 and ar53 gene sets identified in each genome.
        :param passed_qc_gids: Identifier of genomes passing QC criteria.
        :param cur_genomes: Metadata for genomes in the current GTDB release.
        """

        with open(gtdb_domain_report, encoding='utf-8') as f:
            header = f.readline().rstrip().split('\t')

            domain_index = header.index('Predicted domain')
            bac_marker_perc_index = header.index('Bac120 Markers (%)')
            ar_marker_perc_index = header.index('Ar53 Markers (%)')
            ncbi_taxonomy_index = header.index('NCBI taxonomy')
            gtdb_taxonomy_index = header.index('GTDB taxonomy')

            for line in f:
                tokens = line.strip().split('\t')

                gid = canonical_gid(tokens[0])
                if gid not in passed_qc_gids:
                    continue

                domain = tokens[domain_index]
                bac_perc = float(tokens[bac_marker_perc_index])
                ar_perc = float(tokens[ar_marker_perc_index])
                ncbi_domain = [t.strip()
                               for t in tokens[ncbi_taxonomy_index].split(';')][0]
                gtdb_domain = [t.strip()
                               for t in tokens[gtdb_taxonomy_index].split(';')][0]

                if ncbi_domain != gtdb_domain and ncbi_domain != 'None':
                    self.log.info(
                        f'NCBI ({ncbi_domain}) and GTDB ({gtdb_domain}) domains disagree in domain report (Bac = {bac_perc:.1f}%; Ar = {ar_perc:.1f}%): {gid}')
                    if cur_genomes[gid].gtdb_taxa.domain != gtdb_domain:
                        self.log.error(
                            f'GTDB ({gtdb_domain}) domains in domain report differs from domain in GTDB taxonomy string ({cur_genomes[gid].gtdb_taxa.domain}): {gid} [THIS MUST BE FIXED BEFORE PROCEEDING]')

                if domain != gtdb_domain and domain != 'None':
                    self.log.error(
                        f'GTDB and predicted domain (Bac = {bac_perc:.1f}%; Ar = {ar_perc:.1f}%) disagree in domain report: {gid} [THIS MUST BE FIXED BEFORE PROCEEDING]')

    def parse_qc_exception_file(self, qc_exception_file: str) -> Tuple[Set[str], Set[str]]:
        """Parse file indicating genomes flagged as exceptions from QC.

        A small number of genomes with nomenclatural importance that would otherwise
        fail QC criteria are maintained. We have also identified a small set of genomes
        that pass our QC criteria, but have been identified as problematic using other
        means so are explicitly flagged as needing to be filtered.

        :param qc_exception_file: File indicating genomes that should pass or fail QC irrespective of QC criteria.
        :return: Sets indicating genomes that must be retained or filtered.
        """

        retain_exceptions = set()
        filter_exceptions = set()
        with open(qc_exception_file, encoding='utf-8') as f:
            header = f.readline().strip().split('\t')
            include_idx = header.index('Include')

            for line in f:
                tokens = line.strip().split('\t')

                gid = canonical_gid(tokens[0].strip())
                retain = tokens[include_idx].lower().startswith('t')

                if retain:
                    retain_exceptions.add(gid)
                else:
                    filter_exceptions.add(gid)

        return retain_exceptions, filter_exceptions

    def qc_genomes(self,
                   cur_genomes: Genomes,
                   marker_perc: float,
                   retain_exceptions: Set[str],
                   filter_exceptions: Set[str],
                   excluded_from_refseq_note: Dict[str, str],
                   qc_criteria: QcCriteria) -> Tuple[Set[str], Set[str]]:
        """Apply QC criteria to all genomes and report genomes passing or failing QC.

        The QC criteria is straight forward except for the use of the CheckM
        strain heterogeneity estimate which is used to relax the contamination
        requirement to 20%. Use of this exception should be revisited and perhaps
        the small number of genomes impacted by this criteria simply explicitly
        marked as exceptions to QC.

        :param cur_genomes: Metadata for genomes in the current GTDB release.
        :param marker_perc: Percent of bac120 or ar53 markers, whichever is larger, identified in each genome.
        :param retain_exceptions: Genomes to retain regardless of QC criteria.
        :param filter_expections: Genomes to filter regardless of QC criteria.
        :param excluded_from_refseq_note: Notes about genomes provided in NCBI RefSeq database.
        :param qc_criteria: Criteria for passing QC.
        :return: Sets indicating genomes that pass or fail QC.
        """

        fout_passed = open(os.path.join(self.output_dir, 'qc_passed.tsv'), 'w')
        fout_failed = open(os.path.join(self.output_dir, 'qc_failed.tsv'), 'w')

        header = 'Accession\tNCBI species\tGTDB taxonomy'
        header += '\tCompleteness (%)\tContamination (%)\tQuality\tStrain heterogeneity at 100%'
        header += '\tMarkers (%)\tNo. contigs\tN50 contigs\tAmbiguous bases'

        fout_passed.write(
            header + '\tScore\tNote\tNCBI exclude from RefSeq\tMarked in GTDB ledger\n')
        fout_failed.write(header)
        fout_failed.write(
            '\tFailed completeness\tFailed contamination\tFailed quality')
        fout_failed.write(
            '\tFailed marker percentage\tFailed no. contigs\tFailed N50 contigs')
        fout_failed.write(
            '\tFailed ambiguous bases\tNCBI exclude from RefSeq\tMarked in GTDB ledger\n')

        passed_qc_gids = set()
        failed_qc_gids = set()
        for gid in cur_genomes:
            failed_tests = defaultdict(int)
            passed_qc = cur_genomes[gid].pass_qc(marker_perc[gid],
                                                 qc_criteria,
                                                 failed_tests)

            if (passed_qc or gid in retain_exceptions) and gid not in filter_exceptions:
                passed_qc_gids.add(gid)
                fout_passed.write('{}\t{}\t{}'.format(
                    gid, cur_genomes[gid].ncbi_taxa.species, cur_genomes[gid].gtdb_taxa))
                fout_passed.write('\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{:.2f}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\n'.format(
                    cur_genomes[gid].comp,
                    cur_genomes[gid].cont,
                    cur_genomes[gid].comp - 5*cur_genomes[gid].cont,
                    f'{cur_genomes[gid].strain_heterogeneity_100:.2f}' if cur_genomes[gid].strain_heterogeneity_100 else '-',
                    marker_perc[gid],
                    cur_genomes[gid].contig_count,
                    cur_genomes[gid].contig_n50,
                    cur_genomes[gid].ambiguous_bases,
                    cur_genomes[gid].score_type_strain(),
                    'Passed QC' if passed_qc else 'Retained as exception',
                    excluded_from_refseq_note.get(gid, ''),
                    gid in retain_exceptions))
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
                fout_failed.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    failed_tests['comp'],
                    failed_tests['cont'],
                    failed_tests['qual'],
                    failed_tests['marker_perc'],
                    failed_tests['contig_count'],
                    failed_tests['N50'],
                    failed_tests['ambig'],
                    excluded_from_refseq_note.get(gid, ''),
                    gid in filter_exceptions))
        fout_passed.close()
        fout_failed.close()

        self.log.info('Retained {:,} ({:.2f}%) genomes and filtered {:,} ({:.2f}%) genomes.'.format(
            len(passed_qc_gids),
            len(passed_qc_gids)*100.0/len(cur_genomes),
            len(failed_qc_gids),
            len(failed_qc_gids)*100.0/len(cur_genomes)))

        return passed_qc_gids, failed_qc_gids

    def check_qc_of_ncbi_species(self,
                                 cur_genomes: Genomes,
                                 marker_perc: float,
                                 retain_exceptions: Set[str],
                                 filter_exceptions: Set[str],
                                 excluded_from_refseq_note: Dict[str, str],
                                 qc_criteria: QcCriteria) -> None:
        """Report impact of QC filtering on NCBI species.

        Produces files that can be inspected to determine species where
        all avaliable genomes assemblies have been filtered and why they
        were filtered. These are good candidated for having a genome selected
        that is an exception to QC, but this must be balanced with ensuring assemblies
        are of sufficient quality as to not negatively impact tree inference.

        :param cur_genomes: Metadata for genomes in the current GTDB release.
        :param marker_perc: Percent of bac120 or ar53 markers, whichever is larger, identified in each genome.
        :param retain_exceptions: Genomes to retain regardless of QC criteria.
        :param filter_expections: Genomes to filter regardless of QC criteria.
        :param excluded_from_refseq_note: Notes about genomes provided in NCBI RefSeq database.
        :param qc_criteria: Criteria for passing QC.
        """

        named_ncbi_species = cur_genomes.named_ncbi_species()
        self.log.info(
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
                                                     qc_criteria,
                                                     failed_tests)

                failed_tests_gids[gid] = failed_tests

                if cur_genomes[gid].is_gtdb_type_strain() or cur_genomes[gid].is_ncbi_type_strain():
                    if (passed_qc or gid in retain_exceptions) and gid not in filter_exceptions:
                        type_pass.add(gid)
                    else:
                        type_fail.add(gid)
                        filtered_genomes += 1
                else:
                    if (passed_qc or gid in retain_exceptions) and gid not in filter_exceptions:
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

        self.log.info(
            f'Filtered {filtered_genomes:,} genomes assigned to NCBI species.')
        self.log.info(
            f'Identified {lost_type:,} species with type genomes failing QC and {lost_sp:,} total species failing QC.')
        self.log.info(
            'Genomes from NCBI species filtered by each criterion:')
        for test in sorted(failed_tests_cumulative):
            self.log.info(f'{test}: {failed_tests_cumulative[test]:,}')

    def run(self,
            prev_gtdb_metadata_file: str,
            cur_gtdb_metadata_file: str,
            ncbi_genbank_assembly_file: str,
            gtdb_domain_report: str,
            gtdb_type_strains_ledger: str,
            qc_exception_file: str,
            ncbi_env_bioproject_ledger: str,
            qc_criteria: QcCriteria) -> None:
        """Quality check genomes being considered for GTDB.

        :param prev_gtdb_metadata_file: File with metadata for genomes in previous GTDB release.
        :param cur_gtdb_metadata_file: File with metadata for genomes in current GTDB release.
        :param ncbi_genbank_assembly_file: File from NCBI with metadata for GenBank assemblies.
        :param gtdb_domain_report: File indicating percent of bac120 and ar53 gene sets identified in each genome.
        :param gtdb_type_strains_ledger: File indicating genomes that should be considered the type strain of a species.
        :param qc_exception_file: File indicating genomes that should fail or pass QC irrespective of the QC criteria.
        :param ncbi_env_bioproject_ledger: File indicating genome that should be considered MAGs.
        :param qc_criteria: Criteria for passing QC.
        """

        # create previous and current GTDB genome sets
        self.log.info('Creating previous GTDB genome set:')
        prev_genomes = Genomes()
        prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                             gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                             ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                             ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)
        self.log.info(
            f' - previous genome set contains {len(prev_genomes):,} genomes')
        self.log.info(' - previous genome set has {:,} species clusters spanning {:,} genomes'.format(
            len(prev_genomes.sp_clusters),
            prev_genomes.sp_clusters.total_num_genomes()))

        self.log.info('Creating current GTDB genome set:')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                            gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                            create_sp_clusters=False,
                                            ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                            ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)

        # parse genomes flagged as exceptions from QC
        self.log.info('Parsing QC ledger:')
        retain_exceptions, filter_exceptions = self.parse_qc_exception_file(
            qc_exception_file)
        self.log.info(
            f' - identified {len(retain_exceptions):,} genomes flagged for retention')
        self.log.info(
            f' - identified {len(filter_exceptions):,} genomes flagged for removal')

        # get percentage of bac120 or ar122 marker genes
        marker_perc = self.parse_marker_percentages(gtdb_domain_report)

        # parse NCBI assembly files
        self.log.info('Parsing NCBI assembly files.')
        excluded_from_refseq_note = exclude_from_refseq(
            ncbi_genbank_assembly_file)

        # QC all genomes
        self.log.info('Validating genomes.')
        passed_qc_gids, failed_qc_gids = self.qc_genomes(
            cur_genomes,
            marker_perc,
            retain_exceptions,
            filter_exceptions,
            excluded_from_refseq_note,
            qc_criteria)

        # check domain assignment of genomes passing QC
        # and report potential issues
        self.check_domain_assignments(gtdb_domain_report,
                                      passed_qc_gids,
                                      cur_genomes)

        # report results of QC on genomes from each NCBI species
        self.check_qc_of_ncbi_species(
            cur_genomes,
            marker_perc,
            retain_exceptions,
            filter_exceptions,
            excluded_from_refseq_note,
            qc_criteria)

        # sanity check QC results by identifying any genomes that passed QC last release, but
        # have now been flagged as failing QC. This should rarely, if ever, happen unless the
        # genomic data of the assembly has been updated.
        fout = open(os.path.join(self.output_dir,
                                 'previously_passed_qc.tsv'), 'w')
        fout.write(
            'Accession\tGTDB taxonomy\tGTDB species representative\tMarker gene (%)\n')
        unexpected_qc_fail = []
        for gid in prev_genomes:
            if gid in cur_genomes:
                if not is_same_accn_version(prev_genomes[gid].ncbi_accn, cur_genomes[gid].ncbi_accn):
                    # genome assembly has been updates so QC status is not expected to be the same
                    continue

                if gid in failed_qc_gids:
                    unexpected_qc_fail.append(gid)
                    fout.write(
                        f"{gid}\t{cur_genomes[gid].gtdb_taxa}\t{cur_genomes[gid].is_gtdb_sp_rep()}\t{marker_perc[gid]:.1f}\n")

        fout.close()

        if len(unexpected_qc_fail) > 0:
            self.log.warning('Identified {:,} genomes that passed QC in previous GTDB release, that failed QC in this release.'.format(
                len(unexpected_qc_fail)))
            self.log.warning(
                ' - examples: {}'.format(','.join(unexpected_qc_fail[0:10])))
