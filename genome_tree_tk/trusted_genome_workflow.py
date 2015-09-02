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

from biolib.taxonomy import Taxonomy


class TrustedGenomeWorkflow(object):
    """Determine trusted genomes based on genome statistics."""

    def __init__(self, assembly_metadata_file, checkm_stats_file, taxonomy_file):
        """Initialization.

        Parameters
        ----------
        assembly_metadata_file : str
            File specifying assembly statistics of genomes.
        checkm_stats_file : str
            File specifying CheckM statistics of genomes.
        taxonomy_file : str
            File specifying 7 rank taxonomy of genomes.
        """

        self.logger = logging.getLogger()

        self.assembly_metadata_file = assembly_metadata_file
        self.checkm_stats_file = checkm_stats_file
        self.taxonomy_file = taxonomy_file

    def _read_metadata(self, assembly_metadata_file):
        """Parse assembly metadata file.

        Parameters
        ----------
        assembly_metadata_file : str
            File specifying assembly statistics of genomes.

        Returns
        -------
        dict : d[assembly_accession] -> [# contigs, N50]
            Number of contigs and N50 stats for assemblies.
        """

        assembly_stats = {}
        with open(assembly_metadata_file) as f:
            header = f.readline().split('\t')
            num_contigs_index = header.index('contig-count')
            contig_N50_index = header.index('contig-N50')
            total_length_index = header.index('total-length')

            for line in f:
                line_split = line.split('\t')

                assembly_accession = line_split[0]

                # the number of contigs and N50 is not specified,
                # if the genome is complete
                num_contigs = 1
                if line_split[num_contigs_index]:
                    num_contigs = int(line_split[num_contigs_index])

                N50 = int(line_split[total_length_index])
                if line_split[contig_N50_index]:
                    N50 = int(line_split[contig_N50_index])

                assembly_stats[assembly_accession] = [num_contigs, N50]

        return assembly_stats

    def _trusted_genomes(self, checkm_stats_file,
                                assembly_metadata,
                                trusted_comp,
                                trusted_cont,
                                max_contigs,
                                min_N50,
                                taxonomy,
                                allow_partial_taxonomy):
        """Identify trusted genomes.

        Parameters
        ----------
        checkm_stats_file : str
            File specifying CheckM statistics of genomes.
        assembly_metadata : d[assembly_accession] -> [# contigs, N50]
            Number of contigs and N50 stats for assemblies.
        trusted_comp : float
            Minimum completeness of trusted genomes.
        trusted_cont : float
            Maximum contamination of trusted genomes.
        max_contigs : int
            Maximum number of contigs within trusted genomes.
        min_N50 : int
            Minimum N50 of trusted genomes.
        taxonomy : d[assembly_accession] -> [d__<taxon>, ..., s__<taxon>]
            Taxa indexed by unique ids.
        allow_partial_taxonomy : bool
            Flag indicating if genomes with incomplete taxonomic information should be retained.

        Returns
        -------
        dict : d[assembly_accession] -> [comp, cont, # contigs, N50]
            Genome statistics for trusted genomes.
        """

        trusted_genome_stats = {}
        with open(checkm_stats_file) as f:
            header = f.readline().split('\t')
            comp_index = header.index('Completeness')
            cont_index = header.index('Contamination')

            for line in f:
                line_split = line.split('\t')

                genome_id = line_split[0]
                assembly_accession = genome_id[0:genome_id.find('_', 4)]

                if not allow_partial_taxonomy and assembly_accession not in taxonomy:
                    continue

                comp = float(line_split[comp_index])
                cont = float(line_split[cont_index])
                num_contigs, N50 = assembly_metadata[assembly_accession]

                if (comp >= (trusted_comp * 100) and
                    cont <= (trusted_cont * 100) and
                    num_contigs <= max_contigs and
                    N50 >= min_N50):
                        trusted_genome_stats[assembly_accession] = [comp, cont, num_contigs, N50]

        return trusted_genome_stats

    def run(self, trusted_comp,
                    trusted_cont,
                    max_contigs,
                    min_N50,
                    allow_partial_taxonomy,
                    output_file):
        """Determine trusted genomes based on genome statistics.

        Parameters
        ----------
        trusted_comp : float
            Minimum completeness to trust genome for marker set inference [0, 1].
        trusted_cont : float
            Maximum contamination to trust genome for marker set inference  [0, 1].
        max_contigs : int
            Maximum number of contigs within trusted genomes.
        min_N50 : int
            Minimum N50 of trusted genomes.
        allow_partial_taxonomy : bool
            Flag indicating if genomes with incomplete taxonomic information should be retained.
        output_file : str
            Output file to contain list of trusted genomes.
        """

        taxonomy = Taxonomy().read(self.taxonomy_file)

        assembly_metadata = self._read_metadata(self.assembly_metadata_file)

        trusted_genomes_stats = self._trusted_genomes(self.checkm_stats_file,
                                                      assembly_metadata,
                                                      trusted_comp, trusted_cont,
                                                      max_contigs, min_N50,
                                                      taxonomy, allow_partial_taxonomy)

        self.logger.info('')
        self.logger.info('  Identified %d trusted genomes.' % len(trusted_genomes_stats))

        fout = open(output_file, 'w')
        fout.write('# Selection criteria:\n')
        fout.write('# Trusted completeness: %f\n' % trusted_comp)
        fout.write('# Trusted contamination: %f\n' % trusted_cont)
        fout.write('# Maximum contigs: %d\n' % max_contigs)
        fout.write('# Minimum N50: %d\n' % min_N50)
        fout.write('#\n')
        fout.write('# Genome Id\tCompleteness,Contamination,# contigs,N50\tTaxonomy\n')

        for assembly_accession, stats in trusted_genomes_stats.iteritems():
            fout.write('%s\t%s\t%s\n' % (assembly_accession, ','.join(map(str, stats)), ';'.join(taxonomy.get(assembly_accession, ['none']))))

        fout.close()
