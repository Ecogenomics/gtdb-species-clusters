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


"""Functions for working with NCBI genomes and metadata."""


def read_metadata(assembly_metadata_file):
    """Parse assembly metadata file.

    Parameters
    ----------
    assembly_metadata_file : str
        File specifying assembly statistics of genomes.

    Returns
    -------
    dict : d[assembly_accession] -> taxonomy id
        Taxonomy id of assemblies.
    set
        Set of complete genomes based on NCBI assembly level.
    """

    accession_to_taxid = {}
    complete_genomes = set()
    with open(assembly_metadata_file) as f:
        header = f.readline().split('\t')
        taxid_index = header.index('ncbi_taxid')
        assembly_level_index = header.index('ncbi_assembly_level')

        for line in f:
            line_split = line.split('\t')

            assembly_accession = line_split[0]
            accession_to_taxid[assembly_accession] = line_split[taxid_index]

            assembly_level = line_split[assembly_level_index]
            if assembly_level == 'Complete Genome':
                complete_genomes.add(assembly_accession)

    return accession_to_taxid, complete_genomes


def read_type_strains(type_strain_file):
    """Parse type strain file.

    Parameters
    ----------
    type_strain_file : str
        File specifying NCBI taxonomy identifiers representing type strains.

    Returns
    -------
    set
        NCBI taxonomy identifiers of type strains.
    """

    type_strain_taxids = set()
    with open(type_strain_file) as f:
        for line in f:
            line_split = line.split('\t')
            type_strain_taxids.add(line_split[0])

    return type_strain_taxids


def get_type_strains(genome_ids, accession_to_taxid, type_strain_taxids):
    """Get genomes representing a type strain.

    Parameters
    ----------
    genome_ids : set
        Genomes to check.
    accession_to_taxid : d[assembly_accession] -> taxonomy id
        NCBI taxonomy identifier for each genomes.
    type_strain_taxids : set
        NCBI taxonomy identifiers of type strains.

    Returns
    -------
    set
        Genomes representing a type strain.
    """

    type_strains = set()
    for genome_id in genome_ids:
        taxid = accession_to_taxid[genome_id]

        if taxid in type_strain_taxids:
            type_strains.add(genome_id)

    return type_strains
