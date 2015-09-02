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
import random
from collections import defaultdict

from biolib.taxonomy import Taxonomy


class DereplicationWorkflow(object):
    """Dereplicate genomes based on taxonomy.

    Each named species is dereplicated to a fixed number of
    taxa, taking care to retain at least one genome representing
    the type species and giving preference to 'Complete' genomes.
    """

    def __init__(self, assembly_metadata_file, taxonomy_file, type_strain_file):
        """Initialization.

        Parameters
        ----------
        assembly_metadata_file : str
            File specifying assembly statistics of genomes.
        taxonomy_file : str
            File specifying 7 rank taxonomy of genomes.
        type_strain_file : str
            File specifying NCBI taxonomy ids representing type strains.
        """

        self.logger = logging.getLogger()

        self.assembly_metadata_file = assembly_metadata_file
        self.taxonomy_file = taxonomy_file
        self.type_strain_file = type_strain_file

    def _read_metadata(self, assembly_metadata_file):
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
            Set representing complete genomes.
        """

        accession_to_taxid = {}
        complete_genomes = set()
        with open(assembly_metadata_file) as f:
            header = f.readline().split('\t')
            taxid_index = header.index('Taxid')
            assembly_level_index = header.index('Assembly level')

            for line in f:
                line_split = line.split('\t')

                assembly_accession = line_split[0]
                accession_to_taxid[assembly_accession] = line_split[taxid_index]

                assembly_level = line_split[assembly_level_index]
                if assembly_level == 'Complete Genome':
                    complete_genomes.add(assembly_accession)

        return accession_to_taxid, complete_genomes

    def _read_type_strains(self, type_strain_file):
        """Parse type strain file.

        Parameters
        ----------
        type_strain_file : str
            File specifying NCBI taxonomy ids representing type strains.

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

    def _read_trusted_genomes(self, trusted_genomes_file):
        """Parse trusted genomes file.

        Parameters
        ----------
        trusted_genomes_file:
            File containing list of trusted genomes.

        Returns
        -------
        set
            Accession of trusted genomes.
        """

        trusted_accessions = set()
        for line in open(trusted_genomes_file):
            if line[0] == '#':
                continue

            line_split = line.split('\t')
            trusted_accessions.add(line_split[0])

        return trusted_accessions

    def _get_type_strains(self, genome_ids, accession_to_taxid, type_strain_taxids):
        """Get genomes representing a type strain.

        Parameters
        ----------
        genome_ids : set
            Genomes to check.
        accession_to_taxid : d[assembly_accession] -> taxonomy id
            Taxonomy id for each genomes.
        trusted_accessions : set
            Set of trusted genomes.

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

    def _dereplicate(self, max_species,
                            taxonomy,
                            complete_genomes,
                            accession_to_taxid,
                            type_strain_taxids,
                            trusted_accessions):
        """Dereplicate genomes based on taxonomy, type strain info, and an optional trusted list.

        Parameters
        ----------
        max_species : int
            Maximum number of genomes of the same species to retain.
        taxonomy : d[assembly_accession] -> [d__, ..., s__]
            Taxonomy of each genomes.
        complete_genomes : set
            Set representing complete genomes.
        accession_to_taxid : d[assembly_accession] -> taxonomy id
            Taxonomy id for each genomes.
        type_strain_taxids : set
            Set of taxonomy ids representing type strains.
        trusted_accessions : set
            Set of trusted genomes.

        Returns
        -------
        set
            Dereplicate set of assemblies.
        """

        # determine genomes belonging to each named species
        species = defaultdict(set)
        for genome_id, t in taxonomy.iteritems():
            sp = t[Taxonomy.rank_index['s__']]

            if not trusted_accessions or genome_id in trusted_accessions:
                species[sp].add(genome_id)

        # dereplicate species
        genomes_to_retain = set()
        for sp, genome_ids in species.iteritems():
            if len(genome_ids) > max_species:
                selected_genomes = set()

                type_strains = self._get_type_strains(genome_ids, accession_to_taxid, type_strain_taxids)
                complete = genome_ids.intersection(complete_genomes)

                complete_type_strains = type_strains.intersection(complete)

                # try to select a complete type strain, otherwise just take
                # any type strain if one exists
                if complete_type_strains:
                    selected_type_strain = random.sample(complete_type_strains, 1)[0]
                    selected_genomes.add(selected_type_strain)
                    genome_ids.difference_update(selected_type_strain)
                    complete.difference_update(selected_type_strain)
                elif type_strains:
                    selected_type_strain = random.sample(type_strains, 1)[0]
                    selected_genomes.add(selected_type_strain)
                    genome_ids.difference_update(selected_type_strain)

                # grab as many complete genomes as possible
                if complete:
                    genomes_to_select = min(len(complete), max_species - len(selected_genomes))
                    selected_complete_genomes = random.sample(complete, genomes_to_select)
                    selected_genomes.update(selected_complete_genomes)
                    genome_ids.difference_update(selected_complete_genomes)

                # grab incomplete genomes to get to the desired number of genomes
                if len(selected_genomes) != max_species:
                    rnd_additional_genomes = random.sample(genome_ids, max_species - len(selected_genomes))
                    selected_genomes.update(rnd_additional_genomes)

                genomes_to_retain.update(selected_genomes)
            else:
                genomes_to_retain.update(genome_ids)

        return genomes_to_retain

    def run(self, max_species,
                    trusted_genomes_file,
                    output_file):
        """Dereplicate genomes to a specific number per named species.

        Parameters
        ----------
        max_species : int
            Maximum number of genomes of the same species to retain.
        trusted_genomes_file:
            File containing list of trusted genomes.
        output_file : str
            Output file to contain list of dereplicated genomes.
        """

        accession_to_taxid, complete_genomes = self._read_metadata(self.assembly_metadata_file)
        taxonomy = Taxonomy().read(self.taxonomy_file)
        type_strain_taxids = self._read_type_strains(self.type_strain_file)

        trusted_accessions = set()
        if trusted_genomes_file:
            trusted_accessions = self._read_trusted_genomes(trusted_genomes_file)

        dereplicated = self._dereplicate(max_species,
                                         taxonomy,
                                         complete_genomes,
                                         accession_to_taxid,
                                         type_strain_taxids,
                                         trusted_accessions)

        self.logger.info('')
        self.logger.info('  Retained %d genomes after dereplication.' % len(dereplicated))

        if not trusted_genomes_file:
            trusted_genomes_file = ''

        fout = open(output_file, 'w')
        fout.write('# Selection criteria:\n')
        fout.write('# Maximum species: %d\n' % max_species)
        fout.write('# Trusted genomes file: %s\n' % trusted_genomes_file)
        fout.write('#\n')
        fout.write('# Genome Id\tComplete\tTaxonomy\n')
        for assembly_accession in dereplicated:
            complete = 'no'
            if assembly_accession in complete_genomes:
                complete = 'yes'
            fout.write('%s\t%s\t%s\n' % (assembly_accession, complete, ';'.join(taxonomy[assembly_accession])))
        fout.close()
