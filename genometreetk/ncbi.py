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

import csv

"""Functions for working with NCBI genomes and metadata."""


def read_genome_dir(genome_dir_file):
    """Parse genome dir file.

    """

    genome_dirs = {}
    for line in open(genome_dir_file):
        line_split = line.rstrip().split('\t')
        genome_dirs[line_split[0]] = line_split[1]

    return genome_dirs


def read_refseq_metadata(metadata_file):
    """Parse metadata for RefSeq genomes from GTDB metadata file.

    Parameters
    ----------
    metadata_file : str
        File specifying metadata for all genomes.

    Returns
    -------
    dict : d[assembly_accession] -> taxonomy id
        Taxonomy id of assemblies.
    set
        Set of complete genomes based on NCBI assembly level.
    set
        Set of reference or representative genomes based on NCBI RefSeq category.
    """

    refseq_genomes = set()
    complete_genomes = set()
    representative_genomes = set()

    csv_reader = csv.reader(open(metadata_file, 'rt'))
    bHeader = True
    for row in csv_reader:
        if bHeader:
            genome_index = row.index('accession')
            assembly_level_index = row.index('ncbi_assembly_level')
            refseq_category_index = row.index('ncbi_refseq_category')
            bHeader = False
        else:
            genome_id = row[genome_index]
            if genome_id.startswith('RS_'):
                refseq_genomes.add(genome_id)

                assembly_level = row[assembly_level_index].lower()
                if assembly_level == 'complete genome':
                    complete_genomes.add(genome_id)

                refseq_category = row[refseq_category_index].lower()
                if 'reference' in refseq_category or 'representative' in refseq_category:
                    representative_genomes.add(genome_id)

    return refseq_genomes, complete_genomes, representative_genomes

