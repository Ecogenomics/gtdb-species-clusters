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

import logging

from genometreetk.common import read_gtdb_taxonomy, read_gtdb_ncbi_taxonomy
import genometreetk.ncbi as ncbi

import csv


class TrustedGenomeWorkflow(object):
    """Determine trusted genomes based on genome statistics."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger()

    def _trusted_genomes(self, metadata_file,
                                trusted_comp,
                                trusted_cont,
                                max_contigs,
                                min_N50):
        """Identify trusted genomes.

        Parameters
        ----------
        metadata_file : str
            Metadata, including CheckM estimates, for all genomes.
        trusted_comp : float [0, 100]
            Minimum completeness of trusted genomes.
        trusted_cont : float [0, 100]
            Maximum contamination of trusted genomes.
        max_contigs : int
            Maximum number of contigs within trusted genomes.
        min_N50 : int
            Minimum N50 of trusted genomes.

        Returns
        -------
        dict : d[assembly_accession] -> [comp, cont, # contigs, N50]
            Genome statistics for trusted genomes.
        """

        trusted_genome_stats = {}

        csv_reader = csv.reader(open(metadata_file, 'rt'))
        bHeader = True
        for row in csv_reader:
            if bHeader:
                headers = row
                genome_index = headers.index('genome')
                comp_index = headers.index('checkm_completeness')
                cont_index = headers.index('checkm_contamination')
                num_contigs_index = headers.index('contig_count')
                n50_index = headers.index('n50_contigs')
                bHeader = False
            else:
                genome_id = row[genome_index]
                comp = float(row[comp_index])
                cont = float(row[cont_index])
                num_contigs = int(row[num_contigs_index])
                N50 = float(row[n50_index])

                if (comp >= trusted_comp and
                    cont <= trusted_cont and
                    num_contigs <= max_contigs and
                    N50 >= min_N50):
                        trusted_genome_stats[genome_id] = [comp, cont, num_contigs, N50]

        return trusted_genome_stats

    def run(self, metadata_file,
                    trusted_comp,
                    trusted_cont,
                    max_contigs,
                    min_N50,
                    refseq_rep,
                    output_file):
        """Determine trusted genomes based on genome statistics.

        Parameters
        ----------
        metadata_file : str
            Metadata, including CheckM estimates, for all genomes.
        trusted_comp : float [0, 100]
            Minimum completeness to trust genome for marker set inference.
        trusted_cont : float [0, 100]
            Maximum contamination to trust genome for marker set inference.
        max_contigs : int
            Maximum number of contigs within trusted genomes.
        min_N50 : int
            Minimum N50 of trusted genomes.
        refseq_rep : boolean
            If true, consider only RefSeq representative and reference genomes.
        output_file : str
            Output file to contain list of trusted genomes.
        """

        representative_genomes = None
        if refseq_rep:
            _accession_to_taxid, complete_genomes, representative_genomes = ncbi.read_refseq_metadata(metadata_file, keep_db_prefix=True)

        gtdb_taxonomy = read_gtdb_taxonomy(metadata_file)
        ncbi_taxonomy = read_gtdb_ncbi_taxonomy(metadata_file)

        trusted_genomes_stats = self._trusted_genomes(metadata_file,
                                                      trusted_comp, trusted_cont,
                                                      max_contigs, min_N50)
        if representative_genomes:
            self.logger.info('Limiting genomes to RefSeq representative.')
            for genome_id in trusted_genomes_stats.keys():
                if genome_id not in representative_genomes:
                    del trusted_genomes_stats[genome_id]

        self.logger.info('Identified %d trusted genomes.' % len(trusted_genomes_stats))

        fout = open(output_file, 'w')
        fout.write('# Selection criteria:\n')
        fout.write('# Trusted completeness: %f\n' % trusted_comp)
        fout.write('# Trusted contamination: %f\n' % trusted_cont)
        fout.write('# Maximum contigs: %d\n' % max_contigs)
        fout.write('# Minimum N50: %d\n' % min_N50)
        fout.write('#\n')
        fout.write('# Genome Id\tCompleteness,Contamination,Contig count,N50\tGTDB Taxonomy\tNCBI Taxonomy\n')

        for assembly_accession, stats in trusted_genomes_stats.iteritems():
            fout.write('%s\t%s\t%s\t%s\n' % (assembly_accession,
                                                 ','.join(map(str, stats)),
                                                 ';'.join(gtdb_taxonomy.get(assembly_accession, ['none'])),
                                                 ';'.join(ncbi_taxonomy.get(assembly_accession, ['none']))))

        fout.close()
