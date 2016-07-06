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
import operator
from collections import defaultdict

from biolib.taxonomy import Taxonomy

from genometreetk.common import (read_gtdb_metadata,
                                 read_gtdb_ncbi_taxonomy,
                                 read_gtdb_ncbi_organism_name,
                                 read_gtdb_taxonomy,
                                 read_gtdb_ncbi_type_strain,
                                 species_label)
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
                            species_labels,
                            representative_genomes,
                            complete_genomes,
                            ncbi_type_strains,
                            lpsn_type_strains,
                            trusted_accessions,
                            genome_quality):
        """Dereplicate genomes based on taxonomy.

        Parameters
        ----------
        assemble_accessions : set
            All assemble accessions to consider.
        max_species : int
            Maximum number of genomes of the same species to retain.
        species_labels : d[assembly_accession] -> species
            Species label for each genome.
        representative_genomes : set
            Set of representative genomes to prioritize.
        complete_genomes : set
            Set of complete genomes to prioritize.
        ncbi_type_strains : set
            Set of genomes marked as type strains at NCBI.
        lpsn_type_strains : d[ncbi_species] -> set of genome IDs
            Genomes marked as type strains at LPSN.
        trusted_accessions : set
            Set of trusted genomes.
        genome_quality : d[genome_id] -> comp - cont
            Quality of each genome.

        Returns
        -------
        set
            Dereplicate set of assemblies.
        """

        # determine genomes belonging to each named species
        retained_reps = 0
        species = defaultdict(set)
        genomes_to_retain = set()
        for genome_id in assemble_accessions:
            sp = species_labels.get(genome_id, None)

            if not trusted_accessions or genome_id in trusted_accessions:
                if not sp:
                    if genome_id in representative_genomes:
                        retained_reps += 1
                    genomes_to_retain.add(genome_id)
                else:
                    species[sp].add(genome_id)

        self.logger.info('Retained %d genomes without a species designation.' % (len(genomes_to_retain) - retained_reps))

        # dereplicate species
        additional_reps = 0
        for sp, genome_ids in species.iteritems():
            representatives = genome_ids.intersection(representative_genomes)

            if len(genome_ids) > max_species:
                selected_genomes = set()

                # get different classes of genomes
                complete = genome_ids.intersection(complete_genomes)
                type_strains = genome_ids.intersection(ncbi_type_strains)
                complete_type_strains = type_strains.intersection(complete)
                complete = complete - complete_type_strains

                # select all representative genomes
                if representatives:
                    selected_genomes.update(representatives)
                    retained_reps += len(representatives)
                    genome_ids.difference_update(representatives)
                    type_strains.difference_update(representatives)
                    complete_type_strains.difference_update(representatives)
                    complete.difference_update(representatives)

                # select genomes based on estimated genome quality
                complete_quality = {k:genome_quality[k] for k in complete}
                complete_quality_sorted = sorted(complete_quality.items(), key=operator.itemgetter(1), reverse=True)

                type_strains_quality = {k:genome_quality[k] for k in type_strains}
                type_strains_quality_sorted = sorted(type_strains_quality.items(), key=operator.itemgetter(1), reverse=True)

                complete_type_strains_quality = {k:genome_quality[k] for k in complete_type_strains}
                complete_type_strains_quality_sorted = sorted(complete_type_strains_quality.items(), key=operator.itemgetter(1), reverse=True)

                # try to select a complete type strain, otherwise just take
                # any type strain if one exists
                if len(selected_genomes) < max_species:
                    if complete_type_strains:
                        selected_type_strain = [x[0] for x in complete_type_strains_quality_sorted[0:1]] # random.sample(complete_type_strains, 1)
                        selected_genomes.update(selected_type_strain)
                    elif type_strains:
                        selected_type_strain = [x[0] for x in type_strains_quality_sorted[0:1]]# random.sample(type_strains, 1)
                        selected_genomes.update(selected_type_strain)

                # grab as many complete genomes as possible
                if len(selected_genomes) < max_species and complete:
                    genomes_to_select = min(len(complete), max_species - len(selected_genomes))
                    selected_complete_genomes = [x[0] for x in complete_quality_sorted[0:genomes_to_select]] #random.sample(complete, genomes_to_select)
                    selected_genomes.update(selected_complete_genomes)

                # grab incomplete genomes to get to the desired number of genomes
                if len(selected_genomes) < max_species and genome_ids:
                    genome_ids.difference_update(selected_genomes)
                    genome_ids_quality = {k:genome_quality[k] for k in genome_ids}
                    genome_ids_quality_sorted = sorted(genome_ids_quality.items(), key=operator.itemgetter(1), reverse=True)
                    
                    genomes_to_select = min(len(genome_ids), max_species - len(selected_genomes))
                    rnd_additional_genomes = [x[0] for x in genome_ids_quality_sorted[0:genomes_to_select]] #random.sample(genome_ids, genomes_to_select)
                    selected_genomes.update(rnd_additional_genomes)

                genomes_to_retain.update(selected_genomes)
                additional_reps += len(selected_genomes) - len(representatives)
            else:
                genomes_to_retain.update(genome_ids)
                retained_reps += len(representatives)
                additional_reps += len(genome_ids) - len(representatives)
                
        self.logger.info('Retained %d RefSeq representatives.' % retained_reps)
        self.logger.info('Retained %d additional genomes in order to establish %s genomes per species.' % (additional_reps, max_species))
      
        # make sure to select at least one LPSN type strain for each species
        lpsn_genomes = 0
        for sp, genome_ids in lpsn_type_strains.iteritems():
            if len(genome_ids.intersection(genomes_to_retain)) >= 1:
                # an LPSN type strain has already been selected for this species
                continue
                
            # select genome with the highest quality
            genome_ids_quality = {k:genome_quality[k] for k in genome_ids}
            genome_ids_quality_sorted = sorted(genome_ids_quality.items(), key=operator.itemgetter(1), reverse=True)
            genomes_to_retain.add(genome_ids_quality_sorted[0][0])
            lpsn_genomes += 1

        self.logger.info('Retained %d additional genomes to ensure a LPSN type strain in each species.' % lpsn_genomes)

        return genomes_to_retain

    def run(self, max_species,
                    trusted_genomes_file,
                    metadata_file,
                    min_rep_comp,
                    max_rep_cont,
                    min_quality,
                    max_contigs,
                    min_N50,
                    max_ambiguous,
                    strict_filtering,
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
        min_rep_comp : float [0, 100]
            Minimum completeness for a genome to be a representative.
        max_rep_cont : float [0, 100]
            Maximum contamination for a genome to be a representative.
        min_quality : float [0, 100]
            Minimum genome quality (comp-4*cont) for a genome to be a representative.
        max_contigs : int
            Maximum number of contigs for a genome to be a representative.
        min_N50 : int
            Minimum N50 for a genome to be a representative.
        max_ambiguous : int
            Maximum number of ambiguous bases for a genome to be a representative.
        strict_filtering : boolean
            If True apply filtering to all genomes, otherise apply lenient 
            filtering to genomes where the chromosome and plasmids are reported 
            as complete.
        output_file : str
            Output file to contain list of dereplicated genomes.
        """

        accession_to_taxid, complete_genomes, representative_genomes = ncbi.read_refseq_metadata(metadata_file, keep_db_prefix=True)
        self.logger.info('Identified %d RefSeq genomes.' % len(accession_to_taxid))
        self.logger.info('Identified %d representative or reference genomes.' % len(representative_genomes))
        self.logger.info('Identified %d complete genomes.' % len(complete_genomes))

        gtdb_taxonomy = read_gtdb_taxonomy(metadata_file)
        ncbi_taxonomy = read_gtdb_ncbi_taxonomy(metadata_file)
        ncbi_organism_names = read_gtdb_ncbi_organism_name(metadata_file)
        species = species_label(gtdb_taxonomy, ncbi_taxonomy, ncbi_organism_names)

        self.logger.info('Identified %d genomes with a GTDB or NCBI species names.' % len(species))

        trusted_accessions = set()
        if trusted_genomes_file:
            trusted_accessions = self._read_trusted_genomes(trusted_genomes_file)

        # get genome quality
        genomes_to_consider = accession_to_taxid.keys()
        filtered_reps = 0

        genome_stats = read_gtdb_metadata(metadata_file, ['checkm_completeness',
                                                            'checkm_contamination',
                                                            'contig_count',
                                                            'n50_contigs',
                                                            'ambiguous_bases',
                                                            'scaffold_count',
                                                            'ssu_count',
                                                            'ncbi_molecule_count',
                                                            'ncbi_unspanned_gaps',
                                                            'ncbi_genome_representation',
                                                            'ncbi_spanned_gaps',
                                                            'ncbi_assembly_level',
                                                            'ncbi_taxonomy',
                                                            'lpsn_strain'])
        missing_quality = set(accession_to_taxid.keys()) - set(genome_stats.keys())
        if missing_quality:
            self.logger.error('There are %d genomes without metadata information.' % len(missing_quality))
            self.exit(-1)

        lpsn_type_strains = defaultdict(set)
        new_genomes_to_consider = []
        genome_quality = {}
        for genome_id in accession_to_taxid.keys():
            stats = genome_stats[genome_id]
            comp = stats.checkm_completeness
            cont = stats.checkm_contamination
            
            if (comp >= min_rep_comp
                    and cont <= max_rep_cont
                    and (comp - 4*cont) >= min_quality
                    and stats.contig_count <= max_contigs
                    and stats.n50_contigs >= min_N50
                    and stats.ambiguous_bases <= max_ambiguous):
                        new_genomes_to_consider.append(genome_id)
                        genome_quality[genome_id] = comp - cont
                        if stats.lpsn_strain:
                            ncbi_species = stats.ncbi_taxonomy.split(';')[6].strip()
                            lpsn_type_strains[ncbi_species].add(genome_id)
            elif not strict_filtering:
                # check if genome appears to consist of only an unspanned
                # chromosome and unspanned plasmids and thus can be 
                # subjected to a more lenient quality check
                if (stats.ncbi_assembly_level in ['Complete Genome', 'Chromosome']
                    and stats.ncbi_genome_representation == 'full'
                    and stats.scaffold_count == stats.ncbi_molecule_count
                    and stats.ncbi_unspanned_gaps == 0
                    and stats.ncbi_spanned_gaps <= 10
                    and stats.ambiguous_bases <= 100
                    and stats.ssu_count >= 1):
                    
                    # apply lenient quality check 
                    if comp >= 50 and cont <= 15:
                        new_genomes_to_consider.append(genome_id)
                        genome_quality[genome_id] = comp - cont
                        if stats.lpsn_strain:
                            ncbi_species = stats.ncbi_taxonomy.split(';')[6].strip()
                            lpsn_type_strains[ncbi_species].add(genome_id)
            
            # check if a representative at NCBI is being filtered
            if genome_id in representative_genomes and genome_id not in new_genomes_to_consider:
                self.logger.warning('Filtered RefSeq representative %s with comp=%.2f, cont=%.2f, contigs=%d, N50=%d' % (genome_id, 
                                                                                                                            comp, 
                                                                                                                            cont, 
                                                                                                                            stats.contig_count, 
                                                                                                                            stats.n50_contigs))
                filtered_reps += 1

        genomes_to_consider = new_genomes_to_consider
        self.logger.info('Filtered %d representative or reference genome based on genome quality.' % filtered_reps)
        self.logger.info('Considering %d genomes after filtering for genome quality.' % (len(genomes_to_consider)))

        ncbi_type_strains = read_gtdb_ncbi_type_strain(metadata_file)
        self.logger.info('Identified %d genomes marked as type strains at NCBI.' % len(ncbi_type_strains))
        self.logger.info('Identified %d genomes marked as type strains at LPSN.' % sum([len(x) for x in lpsn_type_strains.values()]))

        genomes_to_retain = self._dereplicate(genomes_to_consider,
                                                max_species,
                                                species,
                                                representative_genomes,
                                                complete_genomes,
                                                ncbi_type_strains,
                                                lpsn_type_strains,
                                                trusted_accessions,
                                                genome_quality)

        self.logger.info('Retained %d genomes.' % len(genomes_to_retain))

        if not trusted_genomes_file:
            trusted_genomes_file = ''

        fout = open(output_file, 'w')
        fout.write('# Selection criteria:\n')
        fout.write('# Maximum species: %d\n' % max_species)
        fout.write('# Trusted genomes file: %s\n' % trusted_genomes_file)
        fout.write('# Genome quality metadata file: %s\n' % str(metadata_file))
        fout.write('# Min. representative completeness: %s\n' % str(min_rep_comp))
        fout.write('# Max. representative contamination: %s\n' % str(max_rep_cont))
        fout.write('#\n')
        fout.write('# Genome Id\tGTDB Taxonomy\tNCBI Taxonomy\tType strain\tComplete\tRepresentative\n')
        for assembly_accession in genomes_to_retain:
            representative = 'yes' if assembly_accession in representative_genomes else 'no'
            complete = 'yes' if assembly_accession in complete_genomes else 'no'
            ts = 'yes' if assembly_accession in ncbi_type_strains else 'no'
            gtdb_taxa_str = ';'.join(gtdb_taxonomy.get(assembly_accession, Taxonomy.rank_prefixes))
            ncbi_taxa_str = ';'.join(ncbi_taxonomy.get(assembly_accession, Taxonomy.rank_prefixes))

            if assembly_accession.startswith('GCF_'):
                assembly_accession = 'RS_' + assembly_accession
            elif assembly_accession.startswith('GCA_'):
                assembly_accession = 'GB_' + assembly_accession

            fout.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (assembly_accession,
                                                         gtdb_taxa_str,
                                                         ncbi_taxa_str,
                                                         ts,
                                                         complete,
                                                         representative))
        fout.close()
