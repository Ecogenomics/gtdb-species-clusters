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

from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.genome_utils import canonical_gid

from gtdb_species_clusters.prettytable import PrettyTable


class InspectGenomes(object):
    """Report information about a set of genomes."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

    def type_status(self,
                    cur_gtdb_metadata_file,
                    qc_passed_file,
                    ncbi_genbank_assembly_file,
                    untrustworthy_type_file,
                    gtdb_type_strains_ledger,
                    ncbi_env_bioproject_ledger,
                    genome_ids):
        """Report information related to a genome being type material."""

        # create current GTDB genome sets
        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                            gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                            create_sp_clusters=False,
                                            qc_passed_file=qc_passed_file,
                                            ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                            untrustworthy_type_ledger=untrustworthy_type_file,
                                            ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)
        self.logger.info(
            f' - current genome set contains {len(cur_genomes):,} genomes.')

        # report information
        pt = PrettyTable()
        pt.field_names = ['Genome ID', 'GTDB representative', 'GTDB type strain', 'GTDB untrustworthy as type',
                          'NCBI type strain', 'NCBI untrustworthy as type', 'GTDB species', 'NCBI species', 'NCBI strain IDs\n']
        for gid in genome_ids:
            gid = canonical_gid(gid)
            if gid not in cur_genomes:
                self.logger.warning(f'Genome {gid} not in current genome set.')
                continue

            pt.add_row([gid,
                        cur_genomes[gid].is_gtdb_sp_rep(),
                        cur_genomes[gid].is_gtdb_type_strain(),
                        cur_genomes[gid].is_gtdb_untrustworthy_as_type(),
                        cur_genomes[gid].is_ncbi_type_strain(),
                        cur_genomes[gid].is_ncbi_untrustworthy_as_type(),
                        cur_genomes[gid].gtdb_taxa.species,
                        cur_genomes[gid].ncbi_taxa.species,
                        cur_genomes[gid].ncbi_strain_identifiers])

        print(pt)
