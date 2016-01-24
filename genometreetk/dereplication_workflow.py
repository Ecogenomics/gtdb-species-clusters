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

from genometreetk.common import (read_gtdb_genome_quality,
                                 read_gtdb_ncbi_taxonomy,
                                 read_gtdb_ncbi_type_strain)
import genometreetk.ncbi as ncbi


class DereplicationWorkflow(object):
    """Dereplicate genomes based on taxonomy.

    Each named species is dereplicated to a fixed number of
    taxa, taking care to retain all genomes marked as a
    'reference' or 'representative' at NCBI. Preference
    is then given to genomes marked as type strains at
    NCBI.

    All taxa without a fully qualified species name are
    retained.
    """

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger()
        
    def _dereplicate(self, assemble_accessions,
                            max_species,
                            taxonomy,
                            representative_genomes,
                            complete_genomes,
                            ncbi_type_strains,
                            trusted_accessions):
        """Dereplicate genomes based on taxonomy.

        Parameters
        ----------
        assemble_accessions : set
            All assemble accessions to consider.
        max_species : int
            Maximum number of genomes of the same species to retain.
        taxonomy : d[assembly_accession] -> [d__, ..., s__]
            Taxonomy of each genomes.
        representative_genomes : set
            Set of representative genomes to prioritize.
        complete_genomes : set
            Set of complete genomes  to prioritize.
        ncbi_type_strains : set
            Set of genomes marked as type strains at NCBI.
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
        for genome_id in assemble_accessions:
            t = taxonomy[genome_id]
            try:
                sp = t[Taxonomy.rank_index['s__']]
            except:
                # NCBI species not specified in GTDB
                sp = 's__unclassified'

            if not trusted_accessions or genome_id in trusted_accessions:
                if sp == 's__' or sp.lower() == 's__unclassified':
                    genomes_to_retain.add(genome_id)
                else:
                    species[sp].add(genome_id)

        self.logger.info('Retained %d genomes without a species designation.' % len(genomes_to_retain))

        # dereplicate species
        for sp, genome_ids in species.iteritems():
            if len(genome_ids) > max_species:
                selected_genomes = set()

                representatives = genome_ids.intersection(representative_genomes)
                complete = genome_ids.intersection(complete_genomes)
                type_strains = genome_ids.intersection(ncbi_type_strains)
                complete_type_strains = type_strains.intersection(complete)

                # select all representative genomes
                if representatives:
                    selected_genomes.update(representatives)
                    genome_ids.difference_update(representatives)
                    type_strains.difference_update(representatives)
                    complete_type_strains.difference_update(representatives)

                # try to select a complete type strain, otherwise just take
                # any type strain if one exists
                if len(selected_genomes) < max_species:
                    if complete_type_strains:
                        selected_type_strain = random.sample(complete_type_strains, 1)
                        selected_genomes.update(selected_type_strain)
                        genome_ids.difference_update(selected_type_strain)
                        complete.difference_update(selected_type_strain)
                    elif type_strains:
                        selected_type_strain = random.sample(type_strains, 1)
                        selected_genomes.update(selected_type_strain)
                        genome_ids.difference_update(selected_type_strain)

                # grab as many complete genomes as possible
                if len(selected_genomes) < max_species and complete:
                    genomes_to_select = min(len(complete), max_species - len(selected_genomes))
                    selected_complete_genomes = random.sample(complete, genomes_to_select)
                    selected_genomes.update(selected_complete_genomes)
                    genome_ids.difference_update(selected_complete_genomes)

                # grab incomplete genomes to get to the desired number of genomes
                if len(selected_genomes) < max_species and genome_ids:
                    genomes_to_select = min(len(genome_ids), max_species - len(selected_genomes))
                    rnd_additional_genomes = random.sample(genome_ids, genomes_to_select)
                    selected_genomes.update(rnd_additional_genomes)

                genomes_to_retain.update(selected_genomes)
            else:
                genomes_to_retain.update(genome_ids)

        return genomes_to_retain

    def run(self, max_species,
                    trusted_genomes_file,
                    metadata_file,
                    min_rep_quality,
                    output_file):
        """Dereplicate genomes to a specific number per named species.

        Parameters
        ----------
        max_species : int
            Maximum number of genomes of the same species to retain.
        trusted_genomes_file:
            File containing list of trusted genomes.
        metadata_file : str
            Metadata, including CheckM estimates, for all genomes.
        min_rep_quality : float [0, 100]
            Minimum genome quality for a genome to be a representative.
        output_file : str
            Output file to contain list of dereplicated genomes.
        """

        accession_to_taxid, complete_genomes, representative_genomes = ncbi.read_refseq_metadata(metadata_file, keep_db_prefix=True)
        self.logger.info('Identified %d RefSeq genomes.' % len(accession_to_taxid))
        self.logger.info('Identified %d representative or reference genomes.' % len(representative_genomes))
        self.logger.info('Identified %d complete genomes.' % len(complete_genomes))

        taxonomy = read_gtdb_ncbi_taxonomy(metadata_file, keep_db_prefix=True)
        self.logger.info('Identified %d genomes with taxonomy information.' % len(taxonomy))

        trusted_accessions = set()
        if trusted_genomes_file:
            trusted_accessions = self._read_trusted_genomes(trusted_genomes_file)

        # get genome quality
        genomes_to_consider = accession_to_taxid.keys()
        if metadata_file:
            genome_quality = read_gtdb_genome_quality(metadata_file, keep_db_prefix=True)
            missing_quality = set(accession_to_taxid.keys()) - set(genome_quality.keys())
            if missing_quality:
                self.logger.warning('There are %d genomes without metadata information.' % len(missing_quality))

            new_genomes_to_consider = []
            for genome_id in accession_to_taxid.keys():
                comp, cont, qual = genome_quality.get(genome_id, [-1, -1, -1])
                if qual >= min_rep_quality:
                    new_genomes_to_consider.append(genome_id)
                else:
                    # check if genome is marked as a representative at NCBI
                    if genome_id in representative_genomes:
                        self.logger.warning('Filtered RefSeq representative %s with comp = %.2f, cont = %.2f' % (genome_id, comp, cont))

            genomes_to_consider = new_genomes_to_consider
            self.logger.info('Considering %d genomes after filtering at a quality of %.2f.' % (len(genomes_to_consider), min_rep_quality))

        ncbi_type_strains = read_gtdb_ncbi_type_strain(metadata_file, keep_db_prefix=True)
        self.logger.info('Identified %d genomes marked as type strains at NCBI.' % len(ncbi_type_strains))
        
        genomes_to_retain = self._dereplicate(genomes_to_consider,
                                                max_species,
                                                taxonomy,
                                                representative_genomes,
                                                complete_genomes,
                                                ncbi_type_strains,
                                                trusted_accessions)

        self.logger.info('Retained %d genomes.' % len(genomes_to_retain))

        if not trusted_genomes_file:
            trusted_genomes_file = ''

        fout = open(output_file, 'w')
        fout.write('# Selection criteria:\n')
        fout.write('# Maximum species: %d\n' % max_species)
        fout.write('# Trusted genomes file: %s\n' % trusted_genomes_file)
        fout.write('# Genome quality metadata file: %s\n' % str(metadata_file))
        fout.write('# Min. representative quality: %s\n' % str(min_rep_quality))
        fout.write('#\n')
        fout.write('# Genome Id\tTaxonomy\tType strain\tComplete\tRepresentative\n')
        for assembly_accession in genomes_to_retain:
            representative = 'yes' if assembly_accession in representative_genomes else 'no'
            complete = 'yes' if assembly_accession in complete_genomes else 'no'
            ts = 'yes' if assembly_accession in ncbi_type_strains else 'no'
            taxa_str = ';'.join(taxonomy.get(assembly_accession, Taxonomy.rank_prefixes))

            if assembly_accession.startswith('GCF_'):
                assembly_accession = 'RS_' + assembly_accession
            elif assembly_accession.startswith('GCA_'):
                assembly_accession = 'GB_' + assembly_accession

            fout.write('%s\t%s\t%s\t%s\t%s\n' % (assembly_accession, taxa_str, ts, complete, representative))
        fout.close()
