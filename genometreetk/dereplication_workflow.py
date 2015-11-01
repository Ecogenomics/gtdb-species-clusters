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

import sys
import logging
import random
from collections import defaultdict

from biolib.taxonomy import Taxonomy

import genometreetk.ncbi as ncbi


class DereplicationWorkflow(object):
    """Dereplicate genomes based on taxonomy.

    Each named species is dereplicated to a fixed number of
    taxa, taking care to retain at least one type strains and giving
    preference to genomes identified as 'Complete' at NCBI.
    All taxa with a fully qualified species name are retained.
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

        self.logger.info('Assembly metadata: %s' % self.assembly_metadata_file)
        self.logger.info('Taxonomy: %s' % self.taxonomy_file)
        self.logger.info('Type strains: %s' % self.type_strain_file)

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
        genomes_to_retain = set()
        for genome_id, t in taxonomy.iteritems():
            try:
                sp = t[Taxonomy.rank_index['s__']]
            except:
                self.logger.error('Genome does not have a complete taxonomy string: %s' % genome_id)
                sys.exit(0)

            if not trusted_accessions or genome_id in trusted_accessions:
                if sp == 's__' or sp.lower() == 's__unclassified':
                    genomes_to_retain.add(genome_id)
                else:
                    species[sp].add(genome_id)

        # dereplicate species
        for sp, genome_ids in species.iteritems():
            if len(genome_ids) > max_species:
                selected_genomes = set()

                type_strains = ncbi.get_type_strains(genome_ids, accession_to_taxid, type_strain_taxids)
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

        accession_to_taxid, complete_genomes = ncbi.read_metadata(self.assembly_metadata_file)
        taxonomy = Taxonomy().read(self.taxonomy_file)
        type_strain_taxids = ncbi.read_type_strains(self.type_strain_file)

        trusted_accessions = set()
        if trusted_genomes_file:
            trusted_accessions = self._read_trusted_genomes(trusted_genomes_file)

        dereplicated = self._dereplicate(max_species,
                                         taxonomy,
                                         complete_genomes,
                                         accession_to_taxid,
                                         type_strain_taxids,
                                         trusted_accessions)

        self.logger.info('Retained %d genomes after dereplication.' % len(dereplicated))

        type_strains = ncbi.get_type_strains(dereplicated, accession_to_taxid, type_strain_taxids)

        if not trusted_genomes_file:
            trusted_genomes_file = ''

        fout = open(output_file, 'w')
        fout.write('# Selection criteria:\n')
        fout.write('# Maximum species: %d\n' % max_species)
        fout.write('# Trusted genomes file: %s\n' % trusted_genomes_file)
        fout.write('# Assembly metadata: %s\n' % self.assembly_metadata_file)
        fout.write('# Taxonomy: %s\n' % self.taxonomy_file)
        fout.write('# Type strains: %s\n' % self.type_strain_file)
        fout.write('#\n')
        fout.write('# Genome Id\tTaxonomy\tType strain\tComplete\n')
        for assembly_accession in dereplicated:
            complete = 'yes' if assembly_accession in complete_genomes else 'no'
            ts = 'yes' if assembly_accession in type_strains else 'no'
            fout.write('%s\t%s\t%s\t%s\n' % (assembly_accession, ';'.join(taxonomy[assembly_accession]), ts, complete))
        fout.close()
