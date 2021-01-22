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
import ntpath

from biolib.external.execute import check_dependencies, run

from gtdb_species_clusters.genome_utils import read_qc_file


class GTDB_Tk():
    """Perform initial classification of new and updated genomes using GTDB-Tk."""

    def __init__(self, cpus, output_dir):
        """Initialization."""

        check_dependencies(['gtdbtk'])

        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus
        if self.cpus > 64:
            self.logger.error(
                'Testing indicates pplacer will stale if used with more than 64 CPUs.')
            sys.exit(-1)

    def run(self,
            genomes_new_updated_file,
            qc_passed_file,
            batch_size):
        """Perform initial classification of new and updated genomes using GTDB-Tk."""

        # get list of genomes passing QC
        self.logger.info('Reading genomes passing QC.')
        gids_pass_qc = read_qc_file(qc_passed_file)
        self.logger.info(f' - identified {len(gids_pass_qc):,} genomes.')

        # get path to genomes passing QC
        self.logger.info(
            'Reading path to genomic file for new/updated genomes passing QC.')
        genomic_files = []
        new_updated_gids = set()
        total_count = 0
        with open(genomes_new_updated_file, encoding='utf-8') as f:
            header = f.readline().strip().split('\t')

            genomic_file_index = header.index('Genomic file')

            for line in f:
                tokens = line.strip().split('\t')

                gid = tokens[0]
                total_count += 1
                if gid in gids_pass_qc:
                    gf = tokens[genomic_file_index]
                    genomic_files.append((gid, gf))
                    new_updated_gids.add(gid)
        self.logger.info(
            f' - identified {len(genomic_files):,} of {total_count:,} genomes as passing QC.')

        # create batch files
        genome_batch_files = []
        batch_dir = os.path.join(self.output_dir, 'genome_batch_files')
        if os.path.exists(batch_dir):
            self.logger.warning(
                f'Using existing genome batch files in {batch_dir}.')
            for f in os.listdir(batch_dir):
                genome_batch_files.append(os.path.join(batch_dir, f))

            # check if there are genomes not already in a batch file. Ideally,
            # this would never happen, but sometimes we process past this step
            # and then identify genomes missing in the database. These need to
            # be put into a batch file for processing.
            missing_gids = set(new_updated_gids)
            last_batch_idx = 0
            for batch_file in os.listdir(batch_dir):
                idx = int(batch_file.split('_')[1].replace('.lst', ''))
                if idx > last_batch_idx:
                    last_batch_idx = idx

                with open(os.path.join(batch_dir, batch_file)) as f:
                    for line in f:
                        tokens = line.strip().split('\t')
                        missing_gids.discard(tokens[1])

            if len(missing_gids) > 0:
                genome_batch_file = os.path.join(
                    batch_dir, f'genomes_{last_batch_idx+1}.lst')
                genome_batch_files.append(genome_batch_file)
                self.logger.info('Added the batch file {} with {:,} genomes.'.format(
                    genome_batch_file,
                    len(missing_gids)))

                fout = open(genome_batch_file, 'w')
                for gid, gf in genomic_files:
                    if gid in missing_gids:
                        fout.write('{}\t{}\n'.format(gf, gid))
                fout.close()
        else:
            os.makedirs(batch_dir)
            for batch_idx, start in enumerate(range(0, len(genomic_files), batch_size)):
                genome_batch_file = os.path.join(
                    batch_dir, f'genomes_{batch_idx}.lst')
                genome_batch_files.append(genome_batch_file)

                fout = open(genome_batch_file, 'w')
                for i in range(start, min(start+batch_size, len(genomic_files))):
                    gid, gf = genomic_files[i]
                    fout.write('{}\t{}\n'.format(gf, gid))
                fout.close()

        # process genomes with GTDB-Tk in batches
        for genome_batch_file in genome_batch_files:
            batch_idx = ntpath.basename(genome_batch_file).split('_')[
                1].replace('.lst', '')
            out_dir = os.path.join(self.output_dir, f'gtdbtk_batch{batch_idx}')
            if os.path.exists(out_dir):
                self.logger.warning(
                    f'Skipping genome batch {batch_idx} as output directory already exists.')
                continue

            os.makedirs(out_dir)
            cmd = 'gtdbtk classify_wf --cpus {} --force --batchfile {} --out_dir {}'.format(
                self.cpus,
                genome_batch_file,
                out_dir)
            print(cmd)
            run(cmd)

        # combine summary files
        fout = open(os.path.join(self.output_dir, 'gtdbtk_classify.tsv'), 'w')
        bHeader = True
        gtdbtk_processed = set()
        for batch_dir in os.listdir(self.output_dir):
            if not batch_dir.startswith('gtdbtk_batch'):
                continue

            batch_dir = os.path.join(self.output_dir, batch_dir)
            ar_summary = os.path.join(batch_dir, 'gtdbtk.ar122.summary.tsv')
            bac_summary = os.path.join(batch_dir, 'gtdbtk.bac120.summary.tsv')

            for summary_file in [ar_summary, bac_summary]:
                with open(summary_file, encoding='utf-8') as f:
                    header = f.readline()

                    if bHeader:
                        fout.write(header)
                        bHeader = False

                    for line in f:
                        tokens = line.strip().split('\t')
                        gid = tokens[0]
                        if gid in new_updated_gids:
                            # Ideally, this shouldn't be necessary, but
                            # sometimes we process past this step and then
                            # identify genomes missing in the database. This
                            # can result in GTDB-Tk having been applied to
                            # genomes that looked like they were "new", but
                            # really were just erroneously missing from the
                            # database.
                            fout.write(line)
                            gtdbtk_processed.add(gid)

        fout.close()

        self.logger.info(
            'Identified {:,} genomes as being processed by GTDB-Tk.'.format(len(gtdbtk_processed)))
        skipped_gids = new_updated_gids - gtdbtk_processed
        if len(skipped_gids) > 0:
            self.logger.warning('Identified {:,} genomes as being skipped by GTDB-Tk.'.format(
                len(skipped_gids)))
