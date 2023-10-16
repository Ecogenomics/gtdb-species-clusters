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
import multiprocessing as mp
from collections import defaultdict
from typing import Set, Dict

from gtdb_species_clusters.genomes import Genomes


class LPSN_SSU_Types():
    """Identify type genomes based on type 16S rRNA sequences indicated at LPSN.

    Example: LPSN indicates AB042061 as "type 16S rRNA sequence" for Bacillus subtilis. Presumably,
    a genome containing AB042061 is assembled from the type strain of B. subtilis.

    It should be appreciated that a genome can be assembled from the type strain
    of the species and have a different 16S rRNA gene annotation. The 16S rRNA
    gene indicated at LPSN is simply the gene originally deposited as an examplar
    of the type strain and in many cases is not associated with a genome assembly.

    That said, if the 16S rRNA sequence indicated at LPSN can be associated with a
    genome assembly, then this assembly should be considered the best choice for the
    GTDB representative unless there is evidence to the contrary. This is especially
    helpful for Candidatus species where there is no notion of a type strain. It is
    also useful for resolving cases where two or more genomes are all assembled from
    the type strain.
    """

    def __init__(self, cpus: int, output_dir: str):
        """Initialization.

        :param cpus: Number of CPUs to be used by GTDB-Tk.
        :param output_dir: Output directory.
        """

        self.cpus = cpus
        self.output_dir = output_dir
        self.log = logging.getLogger('rich')

    def parse_ncbi_assem_report(self, assem_report: str) -> Set[str]:
        """Parse sequence information from NCBI Assembly Report.

        :param assem_report: NCBI assembly report file with assembly information.
        :return: Sequence accessions specified in the assembly report.
        """

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

    def parse_lpsn_ssu_metadata(self, lpsn_metadata_file: str) -> Dict[str, str]:
        """Parse 16S rRNA accession for type strain of species as indicated at LPSN.

        :param lpsn_metadata_file: LPSN metadata file indicating 16S accessions for type strains.
        :return: Mapping indicating 16S rRNA accession for type strain of species.
        """

        lpsn_sp_type_ssu = {}
        with open(lpsn_metadata_file) as f:
            header = f.readline().strip().split('\t')

            rank_idx = header.index('Rank')
            name_idx = header.index('Name')
            rRNA_idx = header.index('16S rRNA gene')

            for line in f:
                tokens = line.strip().split('\t')

                rank = tokens[rank_idx]
                if rank == 'species':
                    species = 's__' + tokens[name_idx]
                    rRNA = tokens[rRNA_idx]

                    if rRNA.lower() != 'n/a':
                        lpsn_sp_type_ssu[species] = rRNA
                        assert '.' not in rRNA

        return lpsn_sp_type_ssu

    def _worker(self,
                cur_genomes: Genomes,
                ncbi_sp_gids: Dict[str, Set[str]],
                ncbi_candidatus: Set[str],
                ncbi_assem_report: Dict[str, str],
                queue_in: mp.Queue,
                queue_out: mp.Queue) -> None:
        """Match LPSN species with type rRNA sequences to genomes with the same NCBI species classification.

        :param cur_genomes: Metadata for genomes in the current GTDB release.
        :param ncbi_sp_gids: Mapping indicating genomes with a NCBI species classification.
        :param ncbi_candidatus: Genomes annotated as `Candidatus` at NCBI.
        :param ncbi_assem_report: Mapping indicating NCBI assembly report file for a genome.
        :param queue_in: Input queue of LPSN species to process in parallel.
        :param queue_out: Output queue for reporting results.
        """

        while True:
            data = queue_in.get(block=True, timeout=None)
            if data is None:
                break

            lpsn_sp, rRNA = data

            # Sort genome ids by type material status as they are the most
            # likely to contain the 16S rRNA sequence. Only type strain
            # genomes are considered for species with >250 to reduce computational
            # requirements. Species with large number of genome assemblies are of sufficient
            # interest that it is highly likely the type strain has been sequenced. This
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
            gid_with_ssu = None
            for gid in sorted_gids:
                seq_ids = self.parse_ncbi_assem_report(ncbi_assem_report[gid])

                if rRNA in seq_ids:
                    found_ssu = True
                    gid_with_ssu = gid

                    # safe to break since a given rRNA sequence
                    # can appear in at most one genome
                    break

            queue_out.put(
                (found_ssu, gid_with_ssu, gid_with_ssu in ncbi_candidatus, lpsn_sp, rRNA))

    def _writer(self,
                cur_genomes: Genomes,
                gtdb_type_strains: Dict[str, Set[str]],
                num_sp: int,
                queue_writer: mp.Queue) -> None:
        """Write results of identifying types strain genomes using 16S rRNA information at LPSN.

        :param cur_genomes: Metadata for genomes in the current GTDB release.
        :param gtdb_type_strains: Mapping indicating genomes identified by GTDB as being type strain genomes of a species.
        :param num_sp: Number of LPSN species being processed.
        :param queue_writer: Queue on which results are placed.
        """

        fout = open(os.path.join(self.output_dir,
                    'lpsn_ssu_type_genomes.tsv'), 'w')
        fout.write(
            'Genome ID\tCandidatus\tSpecies\trRNA\tIs GTDB type genome\tIs NCBI type genome\tGTDB type strain genomes\n')

        processed = 0
        missing_type_designation = 0
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

                if not cur_genomes[gid].is_gtdb_type_strain():
                    missing_type_designation += 1

            processed += 1
            statusStr = '-> Processing {:,} of {:,} ({:.2f}%) species.'.format(
                processed,
                num_sp,
                float(processed*100)/num_sp).ljust(86)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

        fout.close()

        if missing_type_designation > 0:
            self.log.warning(
                f"[IMPORTANT]: add the {missing_type_designation:,} genomes where `Is GTDB type genome` is FALSE to the `gtdb_type_strains` ledger.")

    def run(self,
            lpsn_metadata_file: str,
            cur_gtdb_metadata_file: str,
            cur_genomic_path_file: str,
            qc_passed_file: str,
            ncbi_genbank_assembly_file: str,
            gtdb_type_strains_ledger: str,
            untrustworthy_type_ledger: str) -> None:
        """Identify type genomes based on type 16S rRNA sequences indicated at LPSN.

        :param lpsn_metadata_file: LPSN metadata file indicating 16S accessions for type strains.
        :param cur_gtdb_metadata_file: File with metadata for genomes in current GTDB release.
        :param cur_genomic_path_file: File with path to current genomic FASTA files for genomes.
        :param qc_passed_file: File indicating genomes that have passed QC.
        :param ncbi_genbank_assembly_file: File from NCBI with metadata for GenBank assemblies.
        :param gtdb_type_strains_ledger: File indicating genomes that should be considered assembled from the type strain.
        :param untrustworthy_type_ledger: File listing genomes that should be considered untrustworthy as type material.
        """

        # create current GTDB genome sets
        self.log.info('Creating current GTDB genome set:')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                            gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                            create_sp_clusters=False,
                                            qc_passed_file=qc_passed_file,
                                            ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                            untrustworthy_type_ledger=untrustworthy_type_ledger)
        cur_genomes.load_genomic_file_paths(cur_genomic_path_file)

        # get LPSN species names with specified sequence type material
        self.log.info('Parsing LPSN type 16S rRNA data:')
        lpsn_sp_type_ssu = self.parse_lpsn_ssu_metadata(lpsn_metadata_file)
        self.log.info(
            f' - identified {len(lpsn_sp_type_ssu):,} species with type 16S rRNA sequence.')

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

            gf = cur_genomes.genomic_files[gid]
            assem_report = gf.replace('_genomic.fna.gz', '_assembly_report.txt')
            ncbi_assem_report[gid] = assem_report

        # match LPSN species with type rRNA sequences to genomes
        # with the same NCBI species classification
        self.log.info(
            'Identifying type genomes through LPSN type 16S rRNA sequences.')

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
