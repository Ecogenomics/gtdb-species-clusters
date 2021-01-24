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

from biolib.taxonomy import Taxonomy

from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager


class PMC_CheckTypeSpecies(object):
    """Check for agreement between GTDB genera and genomes assembled from type species of genus."""

    def __init__(self, output_dir):
        """Initialization."""

        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')

    def run(self,
            manual_taxonomy,
            cur_gtdb_metadata_file,
            qc_passed_file,
            ncbi_genbank_assembly_file,
            untrustworthy_type_file,
            gtdb_type_strains_ledger,
            sp_priority_ledger,
            genus_priority_ledger,
            ncbi_env_bioproject_ledger,
            lpsn_gss_file):
        """Finalize species names based on results of manual curation."""

        # initialize species priority manager
        sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger,
                                                  genus_priority_ledger,
                                                  lpsn_gss_file,
                                                  self.output_dir)

        # identify species and genus names updated during manual curation
        self.logger.info('Parsing manually curated taxonomy.')
        mc_taxonomy = Taxonomy().read(manual_taxonomy, use_canonical_gid=True)
        self.logger.info(
            ' - read taxonomy for {:,} genomes.'.format(len(mc_taxonomy)))

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

        # establish appropriate species names for GTDB clusters with new representatives
        self.logger.info(
            'Identifying type species genomes with incongruent GTDB genus assignments.')
        fout = open(os.path.join(self.output_dir,
                                 'type_species_incongruencies.tsv'), 'w')
        fout.write(
            'Genome ID\tGTDB genus\tNCBI genus\tGTDB genus priority date\tNCBI genus priority date\tPriority status\tNCBI RefSeq note\n')
        num_incongruent = 0
        for rid, taxa in mc_taxonomy.items():
            if cur_genomes[rid].is_gtdb_type_species():
                gtdb_genus = taxa[Taxonomy.GENUS_INDEX]
                ncbi_genus = cur_genomes[rid].ncbi_taxa.genus

                if gtdb_genus != ncbi_genus:
                    priority_genus = sp_priority_mngr.genus_priority(
                        gtdb_genus, ncbi_genus)

                    if priority_genus != gtdb_genus:
                        num_incongruent += 1

                        if priority_genus == ncbi_genus:
                            priority_status = 'NCBI genus name has priority'
                        else:
                            priority_status = 'Genus with priority must be manually established'

                        fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            rid,
                            gtdb_genus,
                            ncbi_genus,
                            sp_priority_mngr.genus_priority_year(gtdb_genus),
                            sp_priority_mngr.genus_priority_year(ncbi_genus),
                            priority_status,
                            cur_genomes[rid].excluded_from_refseq_note))

        self.logger.info(
            ' - identified {:,} genomes with incongruent genus assignments.'.format(num_incongruent))
        fout.close()
