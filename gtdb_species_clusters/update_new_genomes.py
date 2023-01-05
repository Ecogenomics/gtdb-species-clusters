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
from typing import Dict

from gtdblib.util.bio.accession import canonical_gid


class NewGenomes():
    """Identify new or modified genomes."""

    def __init__(self, output_dir: str):
        """Initialization.

        :param output_dir: Output directory.
        """

        self.output_dir = output_dir
        self.log = logging.getLogger("rich")

    def same_genome_accn(self, accn1: str, accn2: str, identical_accns: Dict[str, str]):
        """Check if NCBI genome accessions are the same.

        :param accn1: First accession to compare.
        :param accn2: Second accession to compare.
        :param identical_accns: Map indicating GenBank and RefSeq accessions considered to be identical.
        """

        if accn1 == accn2:
            return True

        if canonical_gid(accn1) != canonical_gid(accn2):
            self.log.error('Genomes have different canonical genome IDs: {}, {}, {}, {}'.format(
                accn1,
                canonical_gid(accn1),
                accn2,
                canonical_gid(accn2)))
            sys.exit(-1)

        accn1 = accn1.replace('RS_', '').replace('GB_', '')
        accn2 = accn2.replace('RS_', '').replace('GB_', '')

        if identical_accns.get(accn1, None) == accn2:
            return True

        return False

    def run(self,
            prev_gtdb_metadata_file: str,
            cur_gtdb_metadata_file: str,
            cur_genome_paths: str,
            ncbi_assembly_summary_genbank: str):
        """Identify new or modified genomes.

        :param prev_gtdb_metadata_file: File with genome metadata for previous GTDB release.
        :param cur_gtdb_metadata_file: File with genome metadata for current GTDB release.
        :param cur_genome_paths: File with path to data files for genomes in current GTDB release.
        :param ncbi_assembly_summary_genbank: File with NCBI assembly summary metadata.
        """

        # get data for genomes in previous GTDB release
        self.log.info('Reading previous GTDB genomes:')
        prev_accns = {}
        gtdb_taxonomy = {}
        gtdb_rep = {}
        ncbi_genome_category = {}
        with open(prev_gtdb_metadata_file, encoding='utf-8') as f:
            header = f.readline().strip().split('\t')

            gtdb_index = header.index('gtdb_taxonomy')
            gtdb_rep_index = header.index('gtdb_representative')
            ncbi_genome_cat_index = header.index('ncbi_genome_category')

            for line in f:
                tokens = line.strip().split('\t')

                gid = tokens[0]
                cid = canonical_gid(gid)
                prev_accns[cid] = gid
                gtdb_taxonomy[cid] = tokens[gtdb_index]
                gtdb_rep[cid] = tokens[gtdb_rep_index]
                ncbi_genome_category[cid] = tokens[ncbi_genome_cat_index]

        self.log.info(f' - identified {len(prev_accns):,} genomes')

        # get genomes in current release
        self.log.info('Reading current GTDB genomes:')
        cur_accns = {}
        with open(cur_gtdb_metadata_file, encoding='utf-8') as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                gid = tokens[0]
                if gid.startswith('U'):  # only concerned with genomes at NCBI
                    continue

                cur_accns[canonical_gid(gid)] = gid
        self.log.info(f' - identified {len(cur_accns):,} genomes')

        # get equivalent GenBank and RefSeq genome assemblies
        self.log.info(
            'Determining identical GenBank and RefSeq accessions.')
        identical_accns = {}
        with open(ncbi_assembly_summary_genbank, encoding='utf-8') as f:
            for line in f:
                if line.startswith('# assembly_accession'):
                    header = line.strip().split('\t')
                    gb_accn_index = header.index('# assembly_accession')
                    rs_accn_index = header.index('gbrs_paired_asm')
                    paired_asm_index = header.index('paired_asm_comp')
                elif line.startswith('#'):
                    continue
                else:
                    tokens = line.strip().split('\t')

                    paired_asm = tokens[paired_asm_index]
                    if paired_asm.lower() == 'identical':
                        gb_accn = tokens[gb_accn_index]
                        rs_accn = tokens[rs_accn_index]
                        identical_accns[gb_accn] = rs_accn
                        identical_accns[rs_accn] = gb_accn

        # identify new and modified genome IDs
        self.log.info('Identifying new or modified genome IDs:')
        new_gids = set()
        updated_gids = set()
        for cur_gid in cur_accns:
            if cur_gid in prev_accns:
                if not self.same_genome_accn(
                        cur_accns[cur_gid],
                        prev_accns[cur_gid],
                        identical_accns):
                    updated_gids.add(cur_gid)
            else:
                # genome not present in previous GTDB release
                new_gids.add(cur_gid)

        lost_gids = set(prev_accns) - set(cur_accns)
        num_lost_gtdb_reps = sum(
            [1 for gid in lost_gids if gtdb_rep[gid] == 't'])
        self.log.info(
            f' - identified {len(new_gids):,} new, {len(updated_gids):,} updated, and {len(lost_gids):,} lost genomes')
        self.log.info(
            f' - {num_lost_gtdb_reps:,} lost genomes were GTDB representatives')

        # get path to current GTDB genome directories
        self.log.info(
            'Identifying path to genomic files for current GTDB genomes:')
        cur_genome_files = {}
        skipped_genomes = 0
        fout = open(os.path.join(self.output_dir, 'skipped_genomes.tsv'), 'w')
        with open(cur_genome_paths) as f:
            for line in f:
                tokens = line.strip().split('\t')
                accn = tokens[0]
                genome_path = tokens[1]

                if len(tokens) == 3:
                    gid = tokens[2]
                    assert canonical_gid(accn) == gid
                else:
                    gid = canonical_gid(accn)

                if gid not in cur_accns:
                    # a genome may not be part of the GTDB release
                    # (e.g., genome has no NCBI taxonomy information
                    # or Prodigal failed to call genes such as for
                    # GCA_000716285.1)
                    skipped_genomes += 1
                    fout.write(gid + '\n')
                    continue

                assembly_id = os.path.basename(os.path.normpath(genome_path))
                genomic_file = os.path.join(
                    genome_path, assembly_id + '_genomic.fna')
                cur_genome_files[gid] = genomic_file

        fout.close()

        assert len(cur_genome_files) == len(cur_accns)

        self.log.info(
            f' - identified genomic file for {len(cur_genome_files):,} genomes')
        self.log.info(
            f' - skipped {skipped_genomes:,} genomes without GTDB metadata')

        # write out new or modified genome IDs
        self.log.info(
            'Writing out path to new and updated genomic FASTA files.')
        output_file = os.path.join(self.output_dir, 'genomes_new_updated.tsv')
        fout = open(output_file, 'w')
        fout.write('Genome ID\tNCBI accession\tStatus\tGenomic file\n')
        for type_str, gids in [('NEW', new_gids), ('UPDATED', updated_gids)]:
            for gid in gids:
                genomic_file = cur_genome_files[gid]
                fout.write('{}\t{}\t{}\t{}\n'.format(
                    gid,
                    cur_accns[gid],
                    type_str,
                    genomic_file))
        fout.close()

        # write out lost genomes
        output_file = os.path.join(self.output_dir, 'genomes_lost.tsv')
        fout = open(output_file, 'w')
        fout.write(
            'Genome ID\tGTDB taxonomy\tGTDB representative\tNCBI genome category\n')
        for gid in lost_gids:
            fout.write('{}\t{}\t{}\t{}\n'.format(
                gid,
                gtdb_taxonomy[gid],
                gtdb_rep[gid],
                ncbi_genome_category[gid]))
        fout.close()
