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
import time
import logging
import datetime
import multiprocessing as mp
from collections import defaultdict

from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.genome_utils import canonical_gid
from gtdb_species_clusters.taxon_utils import generic_name, specific_epithet

from biolib.taxonomy import Taxonomy


class LPSN_SSU_Types():
    """Identify type genomes based on type 16S rRNA sequences indicated at LPSN.

    It should be appreciated that a genome can be assembled from the type strain
    of the species and have a different 16S rRNA gene annotation. The 16S rRNA
    gene indicated at LPSN is simply the gene originally deposited as an examplar
    of the type strain and in many cases is not associated with a genome assembly.
    In fact, for many important species the type is indicated by a partial 16S rRNA
    sequence which, as expected, is not associated with any complete genome assembly.

    That said, if the 16S rRNA sequence indicated at LPSN can be associated with a
    genome assembly, then this assembly should be considered the best choice for the
    GTDB representative unless there is evidence to the contrary. This is especially
    helpful for Candidatus species where there is no notion of a type strain. It is
    also useful for resolving cases where two or more genomes are all assembled from
    the type strain.
    """

    def __init__(self, cpus, output_dir):
        """Initialization."""

        self.cpus = cpus
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')

    def parse_ncbi_assem_report(self, assem_report):
        """Parse sequence information from NCBI Assembly Report."""

        seq_ids = set()
        with open(assem_report) as f:
            for line in f:
                if line.startswith('# Sequence-Name'):
                    header = line.strip().split('\t')
                    gb_accn_idx = header.index('GenBank-Accn')
                    rs_accn_idx = header.index('RefSeq-Accn')
                if line[0] == '#':
                    continue

                tokens = line.strip().split('\t')
                gb_accn = tokens[gb_accn_idx]
                rs_accn = tokens[rs_accn_idx]

                # remove version number as LPSN does not record this
                gb_accn = gb_accn[:gb_accn.rfind('.')]
                rs_accn = rs_accn[:rs_accn.rfind('.')]

                seq_ids.add(gb_accn)
                seq_ids.add(rs_accn)

        return seq_ids

    def parse_lpsn_ssu_metadata(self, lpsn_metadata_file):
        """Parse type 16S rRNA information from LPSN metadata."""

        lpsn_sp_type_ssu = {}
        with open(lpsn_metadata_file) as f:
            header = f.readline().strip().split('\t')

            rank_idx = header.index('Rank')
            name_idx = header.index('Name')
            rRNA_idx = header.index('16S rRNA gene')
            type_strain_idx = header.index('Type strain')

            for line in f:
                tokens = line.strip().split('\t')

                rank = tokens[rank_idx]
                if rank == 'species':
                    species = 's__' + tokens[name_idx]
                    rRNA = tokens[rRNA_idx]
                    strain_ids = tokens[type_strain_idx]

                    if rRNA.lower() != 'n/a':
                        lpsn_sp_type_ssu[species] = rRNA
                        assert '.' not in rRNA

        return lpsn_sp_type_ssu

    def _worker(self, cur_genomes, ncbi_sp_gids, ncbi_candidatus, ncbi_assem_report, queue_in, queue_out):
        """Process each data item in parallel."""

        while True:
            data = queue_in.get(block=True, timeout=None)
            if data is None:
                break

            lpsn_sp, rRNA = data

            # Sort genome ids by type material status as they are the most
            # likely to contain the 16S rRNA sequence. Only type strain
            # genomes are considered for species with >250. Species
            # with large number of genome assemblies are of sufficient interest
            # that it is highly likely the type strain has been sequenced. This
            # was manually verified to be the case on Feb. 16, 2021.
            sp_gids = ncbi_sp_gids.get(lpsn_sp, [])
            type_gids = [gid
                         for gid in sp_gids
                         if cur_genomes[gid].is_gtdb_type_strain() or cur_genomes[gid].is_ncbi_type_strain()]

            if len(sp_gids) < 250:
                sorted_gids = type_gids + list(set(sp_gids) - set(type_gids))
            else:
                sorted_gids = type_gids

            found_ssu = False
            for gid in sorted_gids:
                seq_ids = self.parse_ncbi_assem_report(ncbi_assem_report[gid])

                if rRNA in seq_ids:
                    found_ssu = True
                    queue_out.put((found_ssu, gid, gid in ncbi_candidatus, lpsn_sp, rRNA))

                    # safe to break since a given rRNA sequence
                    # can appear in at most one genome
                    break

            if not found_ssu:
                queue_out.put((found_ssu, None, None, lpsn_sp, rRNA))

    def _writer(self, cur_genomes, gtdb_type_strains, num_sp, queue_writer):
        """Store or write results of worker threads in a single thread."""

        fout = open(os.path.join(self.output_dir, 'lpsn_ssu_type_genomes.tsv'), 'w')
        fout.write('Genome ID\tCandidatus\tSpecies\trRNA\tIs GTDB type genome\tIs NCBI type genome\tGTDB type strain genomes\n')

        processed = 0
        while True:
            results = queue_writer.get(block=True, timeout=None)
            if results is None:
                break

            (found_ssu, gid, is_candidatus, lpsn_sp, rRNA) = results

            if found_ssu:
                fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    gid,
                    is_candidatus,
                    lpsn_sp,
                    rRNA,
                    cur_genomes[gid].is_gtdb_type_strain(),
                    cur_genomes[gid].is_ncbi_type_strain(),
                    ','.join(gtdb_type_strains[lpsn_sp])))

            processed += 1
            statusStr = '-> Processing {:,} of {:,} ({:.2f}%) species.'.format(
                processed,
                num_sp,
                float(processed*100)/num_sp).ljust(86)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

        fout.close()

    def run(self,
            lpsn_metadata_file,
            cur_gtdb_metadata_file,
            cur_genomic_path_file,
            qc_passed_file,
            ncbi_genbank_assembly_file,
            gtdb_type_strains_ledger,
            untrustworthy_type_ledger):
        """Identify type genomes based on type 16S rRNA sequences indicated at LPSN."""

        # create current GTDB genome sets
        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                            gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                            create_sp_clusters=False,
                                            qc_passed_file=qc_passed_file,
                                            ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                            untrustworthy_type_ledger=untrustworthy_type_ledger)
        cur_genomes.load_genomic_file_paths(cur_genomic_path_file)

        # get LPSN species names with specified sequence type material
        self.logger.info('Parsing LPSN type 16S rRNA data.')
        lpsn_sp_type_ssu = self.parse_lpsn_ssu_metadata(lpsn_metadata_file)
        self.logger.info(f' - identified {len(lpsn_sp_type_ssu):,} species with type 16S rRNA sequence.')

        # get NCBI species assignments for genomes and genomes marked as being
        # type strain genomes
        ncbi_candidatus = set()
        ncbi_sp_gids = defaultdict(set)
        ncbi_assem_report = {}
        gtdb_type_strains = defaultdict(set)
        for gid in cur_genomes:
            ncbi_sp = cur_genomes[gid].ncbi_taxa.species
            ncbi_sp_gids[ncbi_sp].add(gid)

            if 'Candidatus' in cur_genomes[gid].ncbi_unfiltered_taxa.species:
                ncbi_candidatus.add(gid)

            if cur_genomes[gid].is_gtdb_type_strain():
                gtdb_type_strains[ncbi_sp].add(gid)

            ncbi_assem_report[gid] = cur_genomes.genomic_files[gid].replace('_genomic.fna', '_assembly_report.txt')

        # match LPSN species with type rRNA sequences to genomes
        # with the same NCBI species classification
        self.logger.info('Identifying type genomes through LPSN type 16S rRNA sequences.')

        worker_queue = mp.Queue()
        writer_queue = mp.Queue()

        for lpsn_sp, rRNA in lpsn_sp_type_ssu.items():
            worker_queue.put((lpsn_sp, rRNA))

        for _ in range(self.cpus):
            worker_queue.put(None)

        try:
            worker_proc = [mp.Process(target=self._worker, args=(cur_genomes,
                                                                 ncbi_sp_gids,
                                                                 ncbi_candidatus,
                                                                 ncbi_assem_report,
                                                                 worker_queue,
                                                                 writer_queue)) for _ in range(self.cpus)]
            write_proc = mp.Process(target=self._writer, args=(cur_genomes,
                                                               gtdb_type_strains,
                                                               len(lpsn_sp_type_ssu),
                                                               writer_queue))

            write_proc.start()

            for p in worker_proc:
                p.start()

            for p in worker_proc:
                p.join()

            writer_queue.put(None)
            write_proc.join()
        except:
            for p in worker_proc:
                p.terminate()
            write_proc.terminate()
