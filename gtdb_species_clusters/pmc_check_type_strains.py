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

from gtdb_species_clusters.taxon_utils import (generic_name,
                                               specific_epithet,
                                               canonical_species,
                                               is_placeholder_taxon,
                                               test_same_epithet)


class PMC_CheckTypeStrains(object):
    """Check for agreement between GTDB species and genomes assembled from type strain of species."""

    def __init__(self, output_dir):
        """Initialization."""

        self.output_dir = output_dir
        self.log = logging.getLogger('timestamp')

    def run(self,
            manual_taxonomy,
            cur_gtdb_metadata_file,
            qc_passed_file,
            ncbi_genbank_assembly_file,
            untrustworthy_type_file,
            gtdb_type_strains_ledger,
            ncbi_env_bioproject_ledger):
        """Finalize species names based on results of manual curation."""

        # identify species and genus names updated during manual curation
        self.log.info('Parsing manually curated taxonomy.')
        mc_taxonomy = Taxonomy().read(manual_taxonomy, use_canonical_gid=True)
        self.log.info(
            ' - read taxonomy for {:,} genomes.'.format(len(mc_taxonomy)))

        # create current GTDB genome sets
        self.log.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                            gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                            create_sp_clusters=False,
                                            qc_passed_file=qc_passed_file,
                                            ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                            untrustworthy_type_ledger=untrustworthy_type_file,
                                            ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)
        self.log.info(
            f' - current genome set contains {len(cur_genomes):,} genomes.')

        # get all GTDB species represented by a type strain:
        gtdb_type_species = set()
        for rid in mc_taxonomy:
            if cur_genomes[rid].is_effective_type_strain():
                gtdb_type_species.add(mc_taxonomy[rid][Taxonomy.SPECIES_INDEX])

        # establish appropriate species names for GTDB clusters with new representatives
        self.log.info(
            'Identifying type strain genomes with incongruent GTDB species assignments.')
        fout = open(os.path.join(self.output_dir,
                                 'type_strains_incongruencies.tsv'), 'w')
        fout.write(
            'Genome ID\tGTDB species\tNCBI species\tGTDB type strain\tNCBI type strain\tNCBI RefSeq note\n')
        num_incongruent = 0
        for rid, taxa in mc_taxonomy.items():
            if cur_genomes[rid].is_effective_type_strain():
                gtdb_sp = taxa[Taxonomy.SPECIES_INDEX]
                gtdb_generic = generic_name(gtdb_sp)

                ncbi_sp = cur_genomes[rid].ncbi_taxa.species
                ncbi_generic = generic_name(ncbi_sp)

                if ncbi_sp == 's__':
                    # NCBI taxonomy is sometimes behind the genome annotation pages,
                    # and do not have a species assignment even for type strain genome
                    continue

                # check if genome is a valid genus transfer into a genus
                # that already contains a species with the specific
                # name which results in a polyphyletic suffix being required
                # e.g. G002240355 is Prauserella marina at NCBI and is
                # transferred into Saccharomonospora under the GTDB. However,
                # Saccharomonospora marina already exists so this genome
                # needs to be S. marina_A.
                if (is_placeholder_taxon(gtdb_sp)
                        and gtdb_generic != ncbi_generic
                        and canonical_species(gtdb_sp) in gtdb_type_species):
                    continue

                if not test_same_epithet(specific_epithet(gtdb_sp), specific_epithet(ncbi_sp)):
                    num_incongruent += 1
                    fout.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        rid,
                        gtdb_sp,
                        ncbi_sp,
                        cur_genomes[rid].is_gtdb_type_strain(),
                        cur_genomes[rid].is_ncbi_type_strain(),
                        cur_genomes[rid].excluded_from_refseq_note))

        self.log.info(
            ' - identified {:,} genomes with incongruent species assignments.'.format(num_incongruent))
        fout.close()
