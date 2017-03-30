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
import sys
import logging
import operator
import tempfile
from collections import defaultdict

import biolib.seq_io as seq_io
from biolib.taxonomy import Taxonomy
from biolib.external.execute import check_on_path

from genometreetk.common import (read_gtdb_metadata,
                                    read_gtdb_taxonomy,
                                    read_gtdb_ncbi_taxonomy,
                                    read_gtdb_ncbi_type_strain)
import genometreetk.ncbi as ncbi

class Representatives(object):
    """Identify representative genomes.

    Representative genomes are identified in a
    greedy manner using an amino acid identity (AAI)
    criteria.

    To ensure good representatives are selected, genomes
    are order before processing. Genomes are order first
    based on their source: RefSeq, GenBank, User. Within
    each source, genomes are order by genome quality
    (completeness - contamination). A threshold is used
    to limit representative to genomes of sufficient
    quality. Furthermore, a genome is not clustered
    to an existing representative if they have different
    species names. NCBI genomes are also never assigned
    to User representatives.
    """

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')
        
        self.prev_rep_quality_boost = 5.0
        
        # strict threshold for clustering a 
        # query genome that may not have an
        # assigned species name
        self.mash_strict_threshold = 0.035         
        
        # clustering for genomes marked 
        # as being from the same GTDB species
        self.mash_gtdb_species_threshold = 0.05
        
        # clustering for genomes marked 
        # as being from the same GTDB species
        # and the same NCBI species
        self.mash_ncbi_species_threshold = 0.1

    def _canonical_species_name(self, gtdb_species_name):
        """Get canonical species name from GTDB species name."""
        
        if gtdb_species_name == 's__':
            return gtdb_species_name
        
        prefix = gtdb_species_name[0:3]
        full_name = gtdb_species_name[3:]
        genus, species = full_name.split(' ')
        
        underscore_pos = genus.rfind('_')
        if underscore_pos != -1:
            genus = genus[0:underscore_pos]
            
        underscore_pos = species.rfind('_')
        if underscore_pos != -1:
            species = species[0:underscore_pos]
            
        return prefix + genus + ' ' + species
        
    def _read_genome_list(self, genome_list_file):
        """Read genomes IDs in file."""
        
        genome_ids = set()
        if genome_list_file:
            for line in open(genome_list_file):
                if line[0] == '#':
                    continue
                    
                line_split = line.strip().split('\t')
                genome_ids.add(line_split[0])
           
        return genome_ids
        
    def _genome_stats(self, metadata_file):
        """Genome genome and assembly quality metadata."""
        
        stats = read_gtdb_metadata(metadata_file, ['checkm_completeness',
                                                    'checkm_contamination',
                                                    'contig_count',
                                                    'n50_scaffolds',
                                                    'ambiguous_bases',
                                                    'total_gap_length',
                                                    'scaffold_count',
                                                    'ssu_count',
                                                    'gtdb_taxonomy',
                                                    'ncbi_molecule_count',
                                                    'ncbi_unspanned_gaps',
                                                    'ncbi_genome_representation',
                                                    'ncbi_spanned_gaps',
                                                    'ncbi_assembly_level',
                                                    'ncbi_taxonomy',
                                                    'ncbi_organism_name',
                                                    'lpsn_strain'])
 
        return stats

    def _select_lpsn_type_strains(self, 
                                    lpsn_type_strains, 
                                    genomes_to_retain,
                                    genome_quality,
                                    prev_gtdb_reps):
        """Ensure a LPSN type strain has been selected for each species."""
        
        lpsn_genomes = 0
        for sp, genome_ids in lpsn_type_strains.iteritems():
            if len(genome_ids.intersection(genomes_to_retain)) >= 1:
                # an LPSN type strain has already been selected for this species
                continue
                
            # select genome with the highest quality
            q = {k:genome_quality[k]+ self.prev_rep_quality_boost*(k in prev_gtdb_reps) 
                    for k in genome_ids}
            q_sorted = sorted(q.items(), key=operator.itemgetter(1), reverse=True)
            genomes_to_retain.add(q_sorted[0][0])
            lpsn_genomes += 1

        self.logger.info('Retained %d additional genomes to ensure a LPSN type strain in each species.' % lpsn_genomes)
   
    def _dereplicate_species(self, 
                                genomes_to_consider,
                                max_species,
                                species_labels,
                                representative_genomes,
                                complete_genomes,
                                ncbi_type_strains,
                                lpsn_type_strains,
                                prev_gtdb_reps,
                                genome_quality):
        """Dereplicate named species.

        Parameters
        ----------
        genomes_to_consider : set
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
        prev_gtdb_reps : set
            Previous GTDB representative.
        genome_quality : d[genome_id] -> comp - cont
            Quality of each genome.

        Returns
        -------
        set
            Dereplicate set of assemblies.
        """
        
        self.logger.info('Boosting quality of previous representatives by %.1f%%.' % self.prev_rep_quality_boost)
        
        # determine genomes belonging to each named species
        species = defaultdict(set)
        genomes_without_species = set()
        genomes_to_retain = set()
        for genome_id in genomes_to_consider:
            sp = species_labels.get(genome_id, None)

            if not sp:
                if genome_id in representative_genomes:
                    genomes_to_retain.add(genome_id)
                else:
                    genomes_without_species.add(genome_id)
            else:
                species[sp].add(genome_id)

        self.logger.info('Retaining %d RefSeq representative genomes without a species designation.' % len(genomes_to_retain))
        self.logger.info('Identified %d non-representative genomes without a species designation.' % len(genomes_without_species))

        # dereplicate species
        retained_reps = len(genomes_to_retain)
        additional_reps = 0
        for sp, genome_ids in species.iteritems():
            representatives = genome_ids.intersection(representative_genomes)

            if len(genome_ids) > max_species:
                selected_genomes = set()
                
                # select all representative genomes
                if representatives:
                    selected_genomes.update(representatives)
                    genome_ids.difference_update(representatives)

                # get different classes of genomes
                complete = genome_ids.intersection(complete_genomes)
                type_strains = genome_ids.intersection(ncbi_type_strains)
                complete_type_strains = type_strains.intersection(complete)
                complete = complete - complete_type_strains

                # try to select a complete type strain, otherwise just take
                # any type strain if one exists
                if len(selected_genomes) < max_species:
                    if complete_type_strains:
                        q = {k:genome_quality[k] + self.prev_rep_quality_boost*(k in prev_gtdb_reps) 
                                for k in complete_type_strains}
                        q_sorted = sorted(q.items(), key=operator.itemgetter(1), reverse=True)
                        selected_genomes.add(q_sorted[0][0])
                    elif type_strains:
                        q = {k:genome_quality[k] + self.prev_rep_quality_boost*(k in prev_gtdb_reps) 
                                for k in type_strains}
                        q_sorted = sorted(q.items(), key=operator.itemgetter(1), reverse=True)
                        selected_genomes.add(q_sorted[0][0])

                # grab as many complete genomes as possible
                if len(selected_genomes) < max_species and complete:
                    q = {k:genome_quality[k] + self.prev_rep_quality_boost*(k in prev_gtdb_reps) 
                            for k in complete}
                    q_sorted = sorted(q.items(), key=operator.itemgetter(1), reverse=True)

                    genomes_to_select = min(len(complete), max_species - len(selected_genomes))
                    selected_complete_genomes = [x[0] for x in q_sorted[0:genomes_to_select]] 
                    selected_genomes.update(selected_complete_genomes)

                # grab incomplete genomes to get to the desired number of genomes
                if len(selected_genomes) < max_species and genome_ids:
                    genome_ids.difference_update(selected_genomes)
                    genome_ids_quality = {k:genome_quality[k] + self.prev_rep_quality_boost*(k in prev_gtdb_reps) for k in genome_ids}
                    genome_ids_quality_sorted = sorted(genome_ids_quality.items(), key=operator.itemgetter(1), reverse=True)
                    
                    genomes_to_select = min(len(genome_ids), max_species - len(selected_genomes))
                    additional_genomes = [x[0] for x in genome_ids_quality_sorted[0:genomes_to_select]] 
                    selected_genomes.update(additional_genomes)
            else:
                selected_genomes = genome_ids
                
            genomes_to_retain.update(selected_genomes)
            retained_reps += len(representatives)
            additional_reps += len(selected_genomes) - len(representatives)
                 
        self.logger.info('Retained %d RefSeq representatives.' % retained_reps)
        self.logger.info('Retained %d additional genomes in order to establish %s genomes per species.' % (additional_reps, max_species))
      
        # make sure to select at least one LPSN type strain for each species
        self._select_lpsn_type_strains(lpsn_type_strains, 
                                        genomes_to_retain,
                                        genome_quality,
                                        prev_gtdb_reps)
        
        return genomes_to_retain

    def dereplicate(self, metadata_file,
                    prev_rep_file,
                    exceptions_file,
                    trusted_user_file,
                    max_species,
                    min_rep_comp,
                    max_rep_cont,
                    min_quality,
                    max_contigs,
                    min_N50,
                    max_ambiguous,
                    max_gap_length,
                    strict_filtering,
                    output_file):
        """Select representative genomes from named species.
        
        Each named species is dereplicated to a fixed number of
        reprsentatives, taking care to retain all genomes marked as a
        'reference' or 'representative' at NCBI. Preference
        is then given to genomes marked as type strains at
        NCBI. Finally, genomes are selected based on estimated quality.

        Parameters
        ----------
        max_species : int
            Maximum number of genomes of the same species to retain.
        prev_rep_file : str
            File indicating previous representatives to favour during selection.
        trusted_genomes_file:
            File containing list of genomes to retain regardless of filtering criteria.
        metadata_file : str
            Metadata, including CheckM estimates, for all genomes.
        min_rep_comp : float [0, 100]
            Minimum completeness for a genome to be a representative.
        max_rep_cont : float [0, 100]
            Maximum contamination for a genome to be a representative.
        min_quality : float [0, 100]
            Minimum genome quality (comp-5*cont) for a genome to be a representative.
        max_contigs : int
            Maximum number of contigs for a genome to be a representative.
        min_N50 : int
            Minimum N50 of scaffolds for a genome to be a representative.
        max_ambiguous : int
            Maximum number of ambiguous bases for a genome to be a representative.
        max_gap_length : int
            Maximum number of ambiguous bases between contigs for a genome to be a representative.
        strict_filtering : boolean
            If True apply filtering to all genomes, otherise apply lenient 
            filtering to genomes where the chromosome and plasmids are reported 
            as complete.
        output_file : str
            Output file to contain list of dereplicated genomes.
        """
        
        # identify previous reps, genomes to treat as exceptions, 
        # and user genomes to process
        prev_gtdb_reps = self._read_genome_list(prev_rep_file)
        exception_genomes = self._read_genome_list(exceptions_file)
        trusted_user_genomes = self._read_genome_list(trusted_user_file)

        (refseq_genomes, 
            complete_genomes, 
            representative_genomes) = ncbi.read_refseq_metadata(metadata_file)
        self.logger.info('Identified %d RefSeq genomes.' % len(refseq_genomes))
        self.logger.info('Identified %d representative or reference genomes.' % len(representative_genomes))
        self.logger.info('Identified %d complete genomes.' % len(complete_genomes))
        self.logger.info('Identified %d genomes in exception list.' % len(exception_genomes))
        self.logger.info('Identified %d trusted user genomes.' % len(trusted_user_genomes))
        self.logger.info('Identified %d previous GTDB representatives.' % len(prev_gtdb_reps))
        
        # get genome and assembly quality
        genome_stats = self._genome_stats(metadata_file)
        
        # get genomes in each named GTDB species
        gtdb_taxonomy = read_gtdb_taxonomy(metadata_file)
        ncbi_taxonomy = read_gtdb_ncbi_taxonomy(metadata_file)
        
        species = {}
        species_index = Taxonomy.rank_index['s__']
        for genome_id, taxa in gtdb_taxonomy.iteritems():
            sp = taxa[species_index]
            if sp != 's__':
                species[genome_id] = sp

        self.logger.info('Identified %d genomes with a GTDB species names.' % len(species))

        # identify genomes passing filtering criteria
        filtered_reps_file = output_file + '.filtered_reps'
        fout = open(filtered_reps_file, 'w')
        fout.write('Genome ID\tCompleteness\tContamination')
        fout.write('\tContig Count\tN50\tAmbiguous Bases\tTotal Gap Length')
        fout.write('\tNote\tNCBI Organism Name\n')

        lpsn_type_strains = defaultdict(set)
        genomes_to_consider = []
        genome_quality = {}
        filtered_reps = 0
        lack_ncbi_taxonomy = 0
        for genome_id in genome_stats.keys():
            if genome_id.startswith('U_') and genome_id not in trusted_user_genomes:
                continue
                
            stats = genome_stats[genome_id]
            comp = stats.checkm_completeness
            cont = stats.checkm_contamination
            
            keep = False
            if genome_id in exception_genomes:
                keep = True
            elif (comp >= min_rep_comp
                    and cont <= max_rep_cont
                    and (comp - 5*cont) >= min_quality
                    and stats.contig_count <= max_contigs
                    and stats.n50_scaffolds >= min_N50
                    and stats.ambiguous_bases <= max_ambiguous
                    and stats.total_gap_length <= max_gap_length):
                        keep = True
            elif not strict_filtering:
                # check if genome appears to consist of only an unspanned
                # chromosome and unspanned plasmids and thus can be 
                # subjected to a more lenient quality check
                if (stats.ncbi_assembly_level in ['Complete Genome', 'Chromosome']
                    and stats.ncbi_genome_representation == 'full'
                    and stats.scaffold_count == stats.ncbi_molecule_count
                    and stats.ncbi_unspanned_gaps == 0
                    and stats.ncbi_spanned_gaps <= 10
                    and stats.ambiguous_bases <= max_ambiguous
                    and stats.total_gap_length <= max_gap_length
                    and stats.ssu_count >= 1):
                    
                    # apply lenient quality check that should pick
                    # up the vast majority (if not all) even highly
                    # reduced genomes and those with substantial genome
                    # duplication leading to high estimated contamination
                    if comp >= 40 and cont <= 15:
                        keep = True
                        
            if keep:
                genomes_to_consider.append(genome_id)
                genome_quality[genome_id] = comp - 5*cont
                if stats.lpsn_strain:
                    gtdb_species = gtdb_taxonomy[genome_id][species_index]
                    if gtdb_species != 's__':
                        lpsn_type_strains[gtdb_species].add(genome_id)
            
            # check if a representative at NCBI is being filtered
            if genome_id in representative_genomes:
                if genome_id not in genomes_to_consider:
                    if comp < min_rep_comp:
                        note = 'failed completeness criteria'
                    elif cont > max_rep_cont:
                        note = 'failed contamination criteria'
                    elif (comp - 5*cont) < min_quality:
                        note = 'failed genome quality criteria'
                    elif stats.contig_count > max_contigs:
                        note = 'failed contig count criteria'
                    elif stats.n50_scaffolds < min_N50:
                        note = 'failed scaffold N50 criteria'
                    elif stats.ambiguous_bases > max_ambiguous:
                        note = 'failed ambiguous bases criteria'
                    elif stats.total_gap_length > max_gap_length:
                        note = 'failed total gap length criteria'
                        
                    fout.write('%s\t%.2f\t%.2f\t%d\t%d\t%d\t%d\t%s\t%s\n' % (
                                genome_id, 
                                comp, 
                                cont, 
                                stats.contig_count, 
                                stats.n50_scaffolds, 
                                stats.ambiguous_bases,
                                stats.total_gap_length,
                                note,
                                stats.ncbi_organism_name
                                ))

                    warning = ('Filtered RefSeq rep %s with comp=%.2f, cont=%.2f, contigs=%d, N50=%d'
                                    % (genome_id, comp, cont, stats.contig_count, stats.n50_scaffolds))
                    self.logger.warning(warning)
                    
                    filtered_reps += 1
                    
            if genome_id in refseq_genomes and not stats.ncbi_taxonomy:
                # this should never happen, but sometimes the NCBI taxonomy
                # is missing information for some genomes probably due to when NCBI
                # updates the taxonomy database relative to RefSeq
                lack_ncbi_taxonomy += 1
                self.logger.warning('RefSeq representative %s has no assigned NCBI taxonomy.' % genome_id)

        fout.close()

        self.logger.info('Identified %d RefSeq representatives without an assigned NCBI taxonomy.' % lack_ncbi_taxonomy)
        self.logger.info('Filtered %d RefSeq representatives based on genome or assembly quality.' % filtered_reps)
        self.logger.info('Filtered RefSeq representatives written to %s' % filtered_reps_file)
        self.logger.info('Considering %d genomes after filtering for genome quality.' % (len(genomes_to_consider)))

        ncbi_type_strains = read_gtdb_ncbi_type_strain(metadata_file)
        self.logger.info('Identified %d genomes marked as type strains at NCBI.' % len(ncbi_type_strains))
        self.logger.info('Identified %d genomes marked as type strains at LPSN.' % sum([len(x) for x in lpsn_type_strains.values()]))

        # dereplicate named species
        genomes_to_retain = self._dereplicate_species(genomes_to_consider,
                                                        max_species,
                                                        species,
                                                        representative_genomes,
                                                        complete_genomes,
                                                        ncbi_type_strains,
                                                        lpsn_type_strains,
                                                        prev_gtdb_reps,
                                                        genome_quality)

        self.logger.info('Retained %d genomes.' % len(genomes_to_retain))

        # write results
        if not exceptions_file:
            exceptions_file = ''

        fout = open(output_file, 'w')
        fout.write('# Selection criteria:\n')
        fout.write('# Maximum species: %d\n' % max_species)
        fout.write('# Exception file: %s\n' % exceptions_file)
        fout.write('# Trusted user genomes file: %s\n' % trusted_user_file)
        fout.write('# Genome quality metadata file: %s\n' % str(metadata_file))
        fout.write('# Min. representative completeness: %s\n' % str(min_rep_comp))
        fout.write('# Max. representative contamination: %s\n' % str(max_rep_cont))
        fout.write('#\n')
        fout.write('# Genome Id\tGTDB Taxonomy\tNCBI Taxonomy\tNCBI Organism Name\tNCBI Type strain\tComplete\tRepresentative\n')
        for genome_id in genomes_to_retain:
            representative = 'yes' if genome_id in representative_genomes else 'no'
            complete = 'yes' if genome_id in complete_genomes else 'no'
            ts = 'yes' if genome_id in ncbi_type_strains else 'no'
            gtdb_taxa_str = ';'.join(gtdb_taxonomy.get(genome_id, Taxonomy.rank_prefixes))
            ncbi_taxa_str = ';'.join(ncbi_taxonomy.get(genome_id, Taxonomy.rank_prefixes))

            fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (genome_id,
                                                            gtdb_taxa_str,
                                                            ncbi_taxa_str,
                                                            genome_stats[genome_id].ncbi_organism_name,
                                                            ts,
                                                            complete,
                                                            representative))
        fout.close()

    def _order_genomes(self, 
                        genomes_to_consider, 
                        genome_quality, 
                        trusted_user_genomes, 
                        prev_gtdb_reps):
        """Order genomes by source and genome quality.

        Parameters
        ----------
        genomes_to_consider : iterable
          Genomes to order.
        genome_quality : d[genome_id] -> genome quality
          Estimate quality (completeness - 5*contamination) of each genome.
        trusted_user_genomes : set
          Trusted User genomes to treat as if they were in GenBank.
        prev_gtdb_reps : set
          Previous GTDB representative.

        Returns
        -------
        list
            Genomes order by source and quality.
        """
        
        # sort genomes by source repository followed by genome quality
        # giving a slight boost to genomes that were previously a representative
        
        self.logger.info('Boosting quality of previous representatives by %.1f%%.' % self.prev_rep_quality_boost)
        for genome_id, quality in  genome_quality.iteritems():
            if genome_id in prev_gtdb_reps:
                genome_quality[genome_id] = quality + self.prev_rep_quality_boost

        sorted_refseq_rep_genomes = []
        sorted_genbank_rep_genomes = []
        sorted_trusted_user_rep_genomes = []
        sorted_by_quality = sorted(genome_quality.items(), 
                                    key=operator.itemgetter(1), 
                                    reverse=True)
        for genome_id, _quality in sorted_by_quality:
            if genome_id not in genomes_to_consider:
                continue

            if genome_id.startswith('RS_'):
                sorted_refseq_rep_genomes.append(genome_id)
            elif genome_id.startswith('GB_'):
                sorted_genbank_rep_genomes.append(genome_id)
            elif genome_id in trusted_user_genomes:
                sorted_trusted_user_rep_genomes.append(genome_id)
            elif genome_id.startswith('U_'):
                # User genomes should not be selected as representatives
                # since these are not publicly available
                pass
            else:
                self.logger.error('Unrecognized genome prefix: %s' % genome_id)
                sys.exit(-1)

        return (sorted_refseq_rep_genomes
                    + sorted_genbank_rep_genomes
                    + sorted_trusted_user_rep_genomes)
                    
    def _read_mash_dists(self, mash_pairwise_file):
        """Read Mash distance file."""

        dists = defaultdict(lambda: {})
        for line in open(mash_pairwise_file):
            line_split = line.strip().split('\t')
            
            query_genome = line_split[1]
            ref_genome = line_split[0]
            if query_genome == ref_genome:
                continue
                
            dists[query_genome][ref_genome] = float(line_split[2])

        return dists

    def _greedy_representatives(self,
                                representatives,
                                ordered_genomes,
                                gtdb_taxonomy,
                                ncbi_taxonomy,
                                mash_pairwise_file):
        """Identify additional representative genomes in a greedy fashion.

        Parameters
        ----------
        representatives : set
          Initial set of representative genomes.
        ordered_genomes : list
          Order list of genomes to process for identifying new representatives.
        mash_pairwise_file : float
          File with Mash distance between genomes

        Returns
        -------
        set
            Representative genomes.
        """
        
        # read Mash distance between genomes
        self.logger.info('Reading pairwise Mash distances between genomes.')
        mash_dists = self._read_mash_dists(mash_pairwise_file)
    
        # perform greedy clustering
        self.logger.info('Preforming greedy clustering.')
        total_genomes = len(ordered_genomes)
        processed_genomes = 0
        while len(ordered_genomes):
            processed_genomes += 1
            if processed_genomes % 100 == 0:
                sys.stdout.write('==> Processed %d of %d genomes.\r' % (processed_genomes, total_genomes))
                sys.stdout.flush()

            genome_id = ordered_genomes.pop(0)
            query_dists = mash_dists[genome_id]
            
            query_gtdb_sp = gtdb_taxonomy[genome_id][6]
            query_ncbi_sp = ncbi_taxonomy[genome_id][6]
            
            assigned_rep = False
            for ref_id in representatives:
                ref_gtdb_sp = gtdb_taxonomy[ref_id][6]

                d = query_dists.get(ref_id, 1.0)
                if (d <= self.mash_strict_threshold
                        and (query_gtdb_sp == 's__' 
                        or ref_gtdb_sp == query_gtdb_sp)):
                        # genomes meet the strict threshold for
                        # clustering and don't conflict in their
                        # assigned species names
                        assigned_rep = True
                        break
                        
                if ref_gtdb_sp == 's__' or ref_gtdb_sp != query_gtdb_sp:
                    continue
                        
                if d <= self.mash_gtdb_species_threshold:
                        # genomes are from same named species and 
                        # meet the threshold for clustering
                        assigned_rep = True
                        break
                elif (d <= self.mash_ncbi_species_threshold 
                        and self._canonical_species_name(ref_gtdb_sp) == query_ncbi_sp):
                        # genomes are from same named species and 
                        # meet the threshold for clustering
                        assigned_rep = True
                        break

            if not assigned_rep:
                # genome was not assigned to an existing representative,
                # so make it a new representative genome
                representatives.add(genome_id)

        sys.stdout.write('==> Processed %d of %d genomes.\r' % (processed_genomes, total_genomes))
        sys.stdout.flush()
        sys.stdout.write('\n')

        return representatives

    def representatives(self,
                        species_derep_file,
                        metadata_file,
                        prev_rep_file,
                        mash_pairwise_file,
                        trusted_user_file,
                        min_rep_comp,
                        max_rep_cont,
                        min_quality,
                        max_contigs,
                        min_N50,
                        max_ambiguous,
                        max_gap_length,
                        output_file):
        """Identify additional representatives.
        
        Additional representatives are selected in a greedy fashion,
        by ordering genomes according to database source and estimated
        genome quality. A slight quality boost is given to genomes that 
        were previously selected as a representative in order to try and
        retain more stability between releases. Genomes only added as a new 
        representative if they cannot be clustered with an existing representative. 
        Clustering is based on a conservative Mash distance threshold that
        reflects the 95% ANI species criteria.

        Parameters
        ----------
        species_derep_file : str
            File listing selected representatives from named species.
        metadata_file : str
            Metadata, including CheckM estimates, for all genomes.
        prev_rep_file : str
            File indicating previous representatives to favour during selection.
        trusted_user_file : str
            File listing trusted User genomes that should be treated as if they are in GenBank.
        mash_pairwise_file : str
          File with pairwise Mash distances.
        min_rep_comp : float [0, 100]
            Minimum completeness for a genome to be a representative.
        max_rep_cont : float [0, 100]
            Maximum contamination for a genome to be a representative.
        min_quality : float [0, 100]
            Minimum quality (comp - 5*cont) for a genome to be a representative.
        max_contigs : int
            Maximum number of contigs for a genome to be a representative.
        min_N50 : int
            Minimum N50 of scaffolds for a genome to be a representative.
        max_ambiguous : int
            Maximum number of ambiguous bases within contigs for a genome to be a representative.
        max_gap_length : int
            Maximum number of ambiguous bases between contigs for a genome to be a representative.
        output_file : str
            Output file containing all genomes identified as representatives.
        """
        
        # read previous representatives and trusted user genomes
        prev_gtdb_reps = self._read_genome_list(prev_rep_file)
        trusted_user_genomes = self._read_genome_list(trusted_user_file)
        
        self.logger.info('Identified %d trusted User genomes.' % len(trusted_user_genomes))
        self.logger.info('Identified %d previous GTDB representatives.' % len(prev_gtdb_reps))

        # get genome and assembly quality
        genome_stats = self._genome_stats(metadata_file)

        # read initial representatives
        init_rep_genomes = set()
        for line in open(species_derep_file):
            if line[0] == '#':
                continue

            genome_id = line.strip().split('\t')[0]
            init_rep_genomes.add(genome_id)

        self.logger.info('Identified %d initial representatives.' % len(init_rep_genomes))

        # remove existing representative genomes and genomes
        # of insufficient quality to be a representative
        genome_quality = {}
        potential_reps = set()
        for genome_id, stats in genome_stats.iteritems():
            if genome_id in init_rep_genomes:
                continue
                
            if genome_id.startswith('U_') and genome_id not in trusted_user_genomes:
                continue
                
            if (stats.checkm_completeness >= min_rep_comp 
                and stats.checkm_contamination <= max_rep_cont 
                and (stats.checkm_completeness - 5*stats.checkm_contamination) >= min_quality
                and stats.contig_count <= max_contigs
                and stats.n50_scaffolds >= min_N50
                and stats.ambiguous_bases <= max_ambiguous
                and stats.total_gap_length <= max_gap_length):
                    potential_reps.add(genome_id)
                    genome_quality[genome_id] = stats.checkm_completeness - 5*stats.checkm_contamination

        # perform greedy identification of new representatives
        ordered_genomes = self._order_genomes(potential_reps, 
                                                genome_quality, 
                                                trusted_user_genomes, 
                                                prev_gtdb_reps)
        info = (('Comparing %d genomes to %d initial representatives.') % (len(ordered_genomes),
                                                                            len(init_rep_genomes)))
        self.logger.info(info)
        gtdb_taxonomy = read_gtdb_taxonomy(metadata_file)
        ncbi_taxonomy = read_gtdb_ncbi_taxonomy(metadata_file)
        representatives = self._greedy_representatives(init_rep_genomes,
                                                        ordered_genomes,
                                                        gtdb_taxonomy,
                                                        ncbi_taxonomy,
                                                        mash_pairwise_file)

        self.logger.info('Identified %d representatives.' % len(representatives))

        # read metadata for genomes
        (refseq_genomes, 
            complete_genomes, 
            representative_genomes) = ncbi.read_refseq_metadata(metadata_file)
        ncbi_type_strains = read_gtdb_ncbi_type_strain(metadata_file)
        
            
        # write out information for representative genomes
        fout = open(output_file, 'w')

        fout.write('# Selection criteria:\n')
        fout.write('# Species dereplication file: %s\n' % species_derep_file)
        fout.write('# Previous representative file: %s\n' % prev_rep_file)
        fout.write('# Trusted user genomes file: %s\n' % trusted_user_file)
        fout.write('# Genome quality metadata file: %s\n' % str(metadata_file))
        fout.write('# Min. representative completeness: %.2f\n' % min_rep_comp)
        fout.write('# Max. representative contamination: %.2f\n' % max_rep_cont)
        fout.write('# Mash strict threshold: %.3f\n' % self.mash_strict_threshold)
        fout.write('# Mash GTDB species threshold: %.3f\n' % self.mash_gtdb_species_threshold)
        fout.write('# Mash NCBI species threshold: %.3f\n' % self.mash_ncbi_species_threshold)
        fout.write('#\n')

        fout.write('# Genome Id\tGTDB Taxonomy\tNCBI Taxonomy\tNCBI Organism Name\tNCBI Type strain\tComplete\tRepresentative\n')
        for genome_id in representatives:
            representative = 'yes' if genome_id in representative_genomes else 'no'
            complete = 'yes' if genome_id in complete_genomes else 'no'
            ts = 'yes' if genome_id in ncbi_type_strains else 'no'
            gtdb_taxa_str = ';'.join(gtdb_taxonomy.get(genome_id, Taxonomy.rank_prefixes))
            ncbi_taxa_str = ';'.join(ncbi_taxonomy.get(genome_id, Taxonomy.rank_prefixes))

            fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (genome_id,
                                                            gtdb_taxa_str,
                                                            ncbi_taxa_str,
                                                            genome_stats[genome_id].ncbi_organism_name,
                                                            ts,
                                                            complete,
                                                            representative))

        fout.close()
        
    def cluster(self,
                rep_genome_file,
                metadata_file,
                mash_pairwise_file,
                output_file):
        """Cluster genomes based on Mash distances.
        
        Genomes are assigned to their closest representative,
        that is below the species cutoff. However, genomes 
        assigned to different GTDB species are never clustered
        together. This allows refinement of species to be 
        performed using alternative methods and ensures this
        will be respected.
        
        Parameters
        ----------
        rep_genome_file : str
          File indicating genome representative.
        metadata_file : str
          Metadata, including CheckM estimates, for all genomes.
        mash_pairwise_file : str
          File with pairwise Mash distances.
        output_file : str
          Output file indicating genome clusters.
        """
        
        # read previous representatives and trusted user genomes
        representatives = self._read_genome_list(rep_genome_file)
        self.logger.info('Identified %d representative genomes.' % len(representatives))
        
        # get genome and assembly quality
        genome_stats = self._genome_stats(metadata_file)
        gtdb_taxonomy = read_gtdb_taxonomy(metadata_file)
        ncbi_taxonomy = read_gtdb_ncbi_taxonomy(metadata_file)
        
        # read Mash distance between genomes
        self.logger.info('Reading pairwise Mash distances between genomes.')
        mash_dists = self._read_mash_dists(mash_pairwise_file)
        
        # cluster genomes
        self.logger.info('Clustering genomes.')
        clusters = {}
        for rep_id in representatives:
            clusters[rep_id] = []
        
        remaining_genomes = set(genome_stats) - representatives
        for i, genome_id in enumerate(remaining_genomes):
            if i % 100 == 0:
                sys.stdout.write('==> Processed %d of %d genomes.\r' % (i+1, len(remaining_genomes)))
                sys.stdout.flush()
                
            query_dists = mash_dists[genome_id]
            
            query_gtdb_sp = gtdb_taxonomy[genome_id][6]
            query_ncbi_sp = ncbi_taxonomy[genome_id][6]
            
            assigned_rep = None
            min_d = 1.0
            for ref_id in representatives:
                d = query_dists.get(ref_id, 1.0)
                if d >= min_d:
                    continue
                
                ref_gtdb_sp = gtdb_taxonomy[ref_id][6]
                
                if (d <= self.mash_strict_threshold
                        and (query_gtdb_sp == 's__' 
                        or ref_gtdb_sp == query_gtdb_sp)):
                        # genomes meet the strict threshold for
                        # clustering and don't conflict in their
                        # assigned species names
                        assigned_rep = ref_id
                        min_d = d
                        continue
                        
                if ref_gtdb_sp == 's__' or ref_gtdb_sp != query_gtdb_sp:
                    continue
                        
                if d <= self.mash_gtdb_species_threshold:
                        # genomes are from same named species and 
                        # meet the threshold for clustering
                        assigned_rep = ref_id
                        min_d = d
                elif (d <= self.mash_ncbi_species_threshold 
                        and self._canonical_species_name(ref_gtdb_sp) == query_ncbi_sp):
                        # genomes are from same named species and 
                        # meet the threshold for clustering
                        assigned_rep = ref_id
                        min_d = d
            
            if assigned_rep:
                clusters[assigned_rep].append(genome_id)
                
        sys.stdout.write('==> Processed %d of %d genomes.\r' % (len(remaining_genomes), 
                                                                len(remaining_genomes)))
        sys.stdout.flush()
        sys.stdout.write('\n')
                
        # write out clusters
        fout = open(output_file, 'w')
        clustered_genomes = 0
        for c, cluster_rep in enumerate(sorted(clusters, key=lambda x: len(clusters[x]), reverse=True)):   
            cluster_str = 'cluster_%d' % (c + 1)
            cluster = clusters[cluster_rep]
            clustered_genomes += len(cluster)
            fout.write('%s\t%s\t%d\t%s\n' % (cluster_rep, cluster_str, len(cluster) + 1, ','.join(cluster)))

        fout.close()
        
        self.logger.info('Assigned %d genomes to representatives.' % clustered_genomes)