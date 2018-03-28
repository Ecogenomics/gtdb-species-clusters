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
import shutil
import tempfile
from collections import defaultdict

import biolib.seq_io as seq_io
from biolib.taxonomy import Taxonomy
from biolib.external.execute import check_on_path

from genometreetk.common import (canonical_species_name,
                                    read_gtdb_metadata,
                                    read_gtdb_taxonomy,
                                    read_gtdb_ncbi_taxonomy)
import genometreetk.ncbi as ncbi

class Representatives(object):
    """Identify GTDB representative genomes."""

    def __init__(self):
        """Initialization."""
        
        check_on_path('ani_calculator')

        self.logger = logging.getLogger('timestamp')
        
        self.true_str = ['t', 'T', 'true', 'True']
        
        # strict threshold for clustering a 
        # query genome that may not have an
        # assigned species name
        self.strict_ani_threshold = 96.5
        
        # clustering for genomes marked 
        # as being from the same GTDB species
        self.gtdb_species_ani_threshold = 95.0
        
        # clustering for genomes marked 
        # as being from the same GTDB species
        # and the same NCBI species
        self.ncbi_species_ani_threshold = 90.0
        
        # alignment fraction threshold for 
        # assigning genomes to the same species
        self._af_threshold = 0.6
        
        # Mash distance threshold for determining
        # genome pairs from which to calculate
        # true ANI values
        self.mash_dist_threshold = 0.1
        
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
        
    def _read_ncbi_frameshift_errors(self, ncbi_assembly_file):
        """Read error status of genomes from NCBI assembly file."""
        
        ncbi_frameshift_errors = set()
             
        for line in open(ncbi_assembly_file):
            line_split = line.strip().split('\t')
            
            if line[0] == '#':
                try:
                    error_index = line_split.index('excluded_from_refseq')
                except:
                    pass
            else:
                gid = line_split[0]
                if gid.startswith('GCA_'):
                    gid = 'GB_' + gid
                else:
                    gid = 'RS_' + gid

                if len(line_split) > error_index:
                    errors = line_split[error_index]
                    if 'many frameshifted proteins' in errors:
                        ncbi_frameshift_errors.add(gid)

        return ncbi_frameshift_errors
        
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
                                                    'ssu_length',
                                                    'gtdb_taxonomy',
                                                    'ncbi_molecule_count',
                                                    'ncbi_unspanned_gaps',
                                                    'ncbi_genome_representation',
                                                    'ncbi_spanned_gaps',
                                                    'ncbi_assembly_level',
                                                    'ncbi_taxonomy',
                                                    'ncbi_organism_name',
                                                    'ncbi_refseq_category',
                                                    'gtdb_type_material',
                                                    'gtdb_type_material_sources',
                                                    'mimag_high_quality'])
 
        return stats

    def _select_highest_quality(self, genome_quality, cur_genomes, genomes_to_select):
        """Select highest quality genomes."""
                    
        q = {k:genome_quality[k] for k in cur_genomes}
        q_sorted = sorted(q.items(), key=operator.itemgetter(1), reverse=True)
        return [x[0] for x in q_sorted[0:genomes_to_select]]
        
    def _genome_quality(self, genome_metadata, prev_rep):
        """Calculate quality of genome."""
        
        q = genome_metadata.checkm_completeness - 5*genome_metadata.checkm_contamination
        q -= 5*float(genome_metadata.contig_count)/100
        q += 100*(genome_metadata.ncbi_assembly_level is not None 
                    and genome_metadata.ncbi_assembly_level.lower() == 'complete genome')
        q += 10*(genome_metadata.ncbi_refseq_category is not None 
                    and 'representative' in genome_metadata.ncbi_refseq_category.lower())
        q += 5*(1 if genome_metadata.gtdb_type_material in self.true_str else 0)
        q += 5*prev_rep
        
        # check for near-complete 16S rRNA gene
        gtdb_domain = genome_metadata.gtdb_taxonomy[0]
        min_ssu_len = 1200
        if gtdb_domain == 'd__Archaea':
            min_ssu_len = 900
            
        if genome_metadata.ssu_length and genome_metadata.ssu_length >= min_ssu_len:
            q += 5

        return q
        
    def _dereplicate_species(self, 
                                genomes_to_consider,
                                max_species,
                                genome_stats,
                                gtdb_taxonomy,
                                ncbi_taxonomy,
                                ncbi_reference_genomes,
                                ncbi_representative_genomes,
                                complete_genomes,
                                prev_gtdb_reps):
        """Dereplicate named species.

        Parameters
        ----------
        genomes_to_consider : set
            All assemble accessions to consider.
        max_species : int
            Maximum number of genomes of the same species to retain.
        gtdb_taxonomy : d[genome_id] -> taxa
            GTDB taxonomy of each genome.
        ncbi_reference_genomes : set
            Set of NCBI reference genomes to prioritize.
        ncbi_representative_genomes : set
            Set of NCBI representative genomes to prioritize.
        complete_genomes : set
            Set of complete genomes to prioritize.
        prev_gtdb_reps : set
            Previous GTDB representative.
            
        Returns
        -------
        set
            Dereplicate set of assemblies.
        """
        
        # determine type material and genome quality
        type_material = defaultdict(set)
        ssu_near_complete = defaultdict(set)
        genome_quality = {}
        mimag_hq_genomes = set()
        for gid, stats in genome_stats.iteritems():
            gtdb_species = gtdb_taxonomy[gid][6]
            if stats.gtdb_type_material in self.true_str and gtdb_species != 's__':
                type_material[gtdb_species].add(gid)

            if stats.ssu_length and stats.ssu_length >= 900:
                ssu_near_complete[gtdb_species].add(gid)
                
            if stats.mimag_high_quality in self.true_str:
                mimag_hq_genomes.add(gid)
                
            genome_quality[gid] = self._genome_quality(stats, gid in prev_gtdb_reps)

        self.logger.info('Identified %d genomes with a GTDB species designation identified as type material.' % len(type_material))
        self.logger.info('Identified %d MIMAG high-quality genomes.' % len(mimag_hq_genomes))

        # determine genomes belonging to each named species
        species = defaultdict(set)
        genomes_without_species = set()
        genomes_to_retain = set()
        for genome_id in genomes_to_consider:
            sp = gtdb_taxonomy[genome_id][6]

            if sp == 's__':
                if (genome_id in ncbi_reference_genomes 
                    or genome_id in ncbi_representative_genomes
                    or genome_stats[genome_id].gtdb_type_material in self.true_str):
                    genomes_to_retain.add(genome_id)
                else:
                    genomes_without_species.add(genome_id)
            else:
                species[sp].add(genome_id)

        self.logger.info('Retaining %d RefSeq reference, representative, or type material genomes without a GTDB species designation.' % len(genomes_to_retain))
        self.logger.info('Identified %d non-representative genomes without a species designation.' % len(genomes_without_species))
        
        # dereplicate species
        retained_reps = len(genomes_to_retain)
        additional_reps = 0
        for sp, genome_ids in species.iteritems():
            if len(genome_ids) > max_species:
                selected_genomes = set()
                
                # select all genomes marked as 'reference' at NCBI
                cur_reference_genomes = genome_ids.intersection(ncbi_reference_genomes)
                if cur_reference_genomes:
                    selected_genomes.update(cur_reference_genomes)
                    genome_ids.difference_update(cur_reference_genomes)
                    
                # select highest quality genomes meeting MIMAG "high-quality" criteria
                if len(selected_genomes) < max_species:
                    cur_mimag_hq_genomes = genome_ids.intersection(mimag_hq_genomes)
                    genomes_to_select = min(len(genome_ids), max_species - len(selected_genomes))
                    s = self._select_highest_quality(genome_quality,
                                                        cur_mimag_hq_genomes,
                                                        genomes_to_select)
                    selected_genomes.update(s)
                    genome_ids.difference_update(s)
                    
                # select highest quality genomes if maximum for species hasn't been reached
                if len(selected_genomes) < max_species:
                    genomes_to_select = min(len(genome_ids), max_species - len(selected_genomes))
                    s = self._select_highest_quality(genome_quality,
                                                        genome_ids,
                                                        genomes_to_select)
                    selected_genomes.update(s)
                    genome_ids.difference_update(s)
                    
                # select a single genome from avaliable type material if avaliable and one hasn't already been selected
                if sp in type_material and selected_genomes.intersection(type_material[sp]) == 0: 
                    cur_type_strains = genome_ids.intersection(type_material[sp])
                    s = self._select_highest_quality(genome_quality,
                                                        cur_type_strains,
                                                        1)
                    selected_genomes.update(s)
                    genome_ids.difference_update(s)
                    
                # select a single genome with a near complete 16S rRNA gene if avaliable, and one hasn't already been selected
                if sp in ssu_near_complete and selected_genomes.intersection(ssu_near_complete[sp]) == 0: 
                    cur_ssu_near_complete = genome_ids.intersection(ssu_near_complete[sp])
                    s = self._select_highest_quality(genome_quality,
                                                        cur_ssu_near_complete,
                                                        1)
                    selected_genomes.update(s)
                    genome_ids.difference_update(s)
            else:
                selected_genomes = genome_ids
                
            genomes_to_retain.update(selected_genomes)
            ncbi_refs_or_reps = selected_genomes.intersection(ncbi_reference_genomes.union(ncbi_representative_genomes))
            retained_reps += len(ncbi_refs_or_reps)
            additional_reps += len(selected_genomes) - len(ncbi_refs_or_reps)
                 
        self.logger.info('Retained %d RefSeq representatives.' % retained_reps)
        self.logger.info('Retained %d additional genomes in order to establish %s genomes per species.' % (additional_reps, max_species))
        
        return genomes_to_retain

    def dereplicate(self, metadata_file,
                    prev_rep_file,
                    exceptions_file,
                    trusted_user_file,
                    #ncbi_assembly_file,
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
        reprsentatives. Genomes are selected as follows:
        
        - select all genomes marked as 'reference' at NCBI
        - select highest quality genomes meeting MIMAG "high-quality" criteria
          - 16S, 23S, and 5S
          - 18 of 20 amino acids
        - select highest quality genomes
        - select 1 (and only 1) genome from the type strain if avaliable, and one hasn't already been selected
        - select 1 (and only 1) genome with a near complete 16S rRNA gene if avaliable, and one hasn't alraedy been selected

        quality = comp - 5*cont
                    - 5*(number of contigs/100) 
                    + 100*('NCBI complete genome') 
                    + 10*('NCBI representative') 
                    + 5*('type strain') 
                    + 5*(near-complete 16S rRNA gene)
                    + 5*(GTDB representative in previous release)

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
        #ncbi_frameshift_errors = self._read_ncbi_frameshift_errors(ncbi_assembly_file)
        #gtdb_ssu_errors = self.... # should excluded genomes identified as having erroneous 16S sequences

        (refseq_genomes, 
            complete_genomes,
            ncbi_reference_genomes,
            ncbi_representative_genomes) = ncbi.read_refseq_metadata(metadata_file)
        self.logger.info('Identified %d RefSeq genomes.' % len(refseq_genomes))
        self.logger.info('Identified %d reference genomes.' % len(ncbi_reference_genomes))
        self.logger.info('Identified %d representative genomes.' % len(ncbi_representative_genomes))
        self.logger.info('Identified %d complete genomes.' % len(complete_genomes))
        self.logger.info('Identified %d genomes in exception list.' % len(exception_genomes))
        self.logger.info('Identified %d trusted user genomes.' % len(trusted_user_genomes))
        self.logger.info('Identified %d previous GTDB representatives.' % len(prev_gtdb_reps))
        #self.logger.info('Identified %d genomes indicated as having frameshift errors.' % len(ncbi_frameshift_errors))
        
        # get genome and assembly quality
        genome_stats = self._genome_stats(metadata_file)
        
        # get genomes in each named GTDB species
        gtdb_taxonomy = read_gtdb_taxonomy(metadata_file)
        ncbi_taxonomy = read_gtdb_ncbi_taxonomy(metadata_file)

        # identify genomes passing filtering criteria
        filtered_reps_file = output_file + '.filtered_reps'
        fout = open(filtered_reps_file, 'w')
        fout.write('Genome ID\tCompleteness\tContamination')
        fout.write('\tContig Count\tN50\tAmbiguous Bases\tTotal Gap Length')
        fout.write('\tNote\tNCBI Organism Name\n')

        genomes_to_consider = []
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
                    and stats.total_gap_length <= max_gap_length): #and genome_id not in ncbi_frameshift_errors):
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
                    and stats.ssu_count >= 1): #and genome_id not in ncbi_frameshift_errors):
                    
                    # apply lenient quality check that should pick
                    # up the vast majority (if not all) even highly
                    # reduced genomes and those with substantial genome
                    # duplication leading to high estimated contamination
                    if comp >= 40 and cont <= 15:
                        keep = True
                        
            if keep:
                genomes_to_consider.append(genome_id)

            # check if a representative at NCBI is being filtered
            if genome_id in ncbi_reference_genomes or genome_id in ncbi_representative_genomes:
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
                    #self.logger.warning(warning)
                    if genome_id in ncbi_reference_genomes:
                        self.logger.warning(warning)
                        self.logger.warning('******Filtered genome is marked as a reference at NCBI.')
                    
                    filtered_reps += 1
                    
            if genome_id in refseq_genomes and not stats.ncbi_taxonomy:
                # this should never happen, but sometimes the NCBI taxonomy
                # is missing information for some genomes probably due to when NCBI
                # updates the taxonomy database relative to RefSeq
                lack_ncbi_taxonomy += 1
                self.logger.warning('RefSeq representative %s has no assigned NCBI taxonomy.' % genome_id)

        fout.close()

        self.logger.info('Identified %d RefSeq representatives without an assigned NCBI taxonomy.' % lack_ncbi_taxonomy)
        self.logger.warning('Filtered %d RefSeq representatives based on genome or assembly quality.' % filtered_reps)
        self.logger.info('Filtered RefSeq representatives written to %s' % filtered_reps_file)
        self.logger.info('Considering %d genomes after filtering for genome quality.' % (len(genomes_to_consider)))

        # dereplicate named species
        genomes_to_retain = self._dereplicate_species(genomes_to_consider,
                                                        max_species,
                                                        genome_stats,
                                                        gtdb_taxonomy,
                                                        ncbi_taxonomy,
                                                        ncbi_reference_genomes,
                                                        ncbi_representative_genomes,
                                                        complete_genomes,
                                                        prev_gtdb_reps)

        self.logger.info('Retained %d genomes.' % len(genomes_to_retain))
        
        # validate type strain information for retained representatives
        type_material_outfile = output_file + '.type_material_mismatch'
        fout = open(type_material_outfile, 'w')
        fout.write('Accession\tGTDB species\tNCBI species\tNCBI organism name\tType material sources\n')
        type_material_incongruence = 0
        for gid in genomes_to_retain:
            stats = genome_stats[gid]
            if stats.gtdb_type_material in self.true_str:
                gtdb_species = gtdb_taxonomy[gid][6]
                ncbi_species = ncbi_taxonomy[gid][6]
                
                if len(ncbi_species.split(' ')) > 2: # ***HACK to avoid corrupt NCBI species names
                    continue
                if (canonical_species_name(ncbi_species) != canonical_species_name(gtdb_species)
                        and canonical_species_name(gtdb_species)[3:] not in stats.ncbi_organism_name):
                    #self.logger.error("GTDB and NCBI disagree on species of type material: %s, %s, %s, %s" % (gid, 
                    #                                                                                            gtdb_species, 
                    #                                                                                            ncbi_species, 
                    #                                                                                            stats.ncbi_organism_name))
                    
                    type_material_incongruence += 1
                    fout.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (gid, 
                                                                gtdb_species, 
                                                                ncbi_species, 
                                                                stats.ncbi_organism_name, 
                                                                stats.ncbi_organism_name,
                                                                stats.gtdb_type_material_sources))
        fout.close()
        
        if type_material_incongruence > 0:
            self.logger.warning('Identified %d type material genomes with different GTDB and NCBI species designations.' % type_material_incongruence)
            self.logger.info('Type material with incongruent species assignments written to %s' % type_material_outfile)
            
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
        fout.write('# Genome Id\tGTDB Taxonomy\tNCBI Taxonomy\tNCBI Organism Name\tType strain\tType sources\tComplete\tNCBI reference genome\tNCBI representative genome\n')
        for gid in genomes_to_retain:
            reference = 'yes' if gid in ncbi_reference_genomes else 'no'
            representative = 'yes' if gid in ncbi_representative_genomes else 'no'
            complete = 'yes' if gid in complete_genomes else 'no'
            ts = 'yes' if genome_stats[gid].gtdb_type_material in self.true_str else 'no'

            gtdb_taxa_str = ';'.join(gtdb_taxonomy.get(gid, Taxonomy.rank_prefixes))
            ncbi_taxa_str = ';'.join(ncbi_taxonomy.get(gid, Taxonomy.rank_prefixes))

            fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (gid,
                                                            gtdb_taxa_str,
                                                            ncbi_taxa_str,
                                                            genome_stats[gid].ncbi_organism_name,
                                                            ts,
                                                            genome_stats[gid].gtdb_type_material_sources,
                                                            complete,
                                                            reference,
                                                            representative))
        fout.close()

    def _order_genomes(self, genome_quality):
        """Order genomes by source and genome quality.

        Parameters
        ----------
        genome_quality : d[genome_id] -> genome quality
          Estimate quality of each genome.

        Returns
        -------
        list
            Genomes order by source and quality.
        """
        
        # sort genomes by source repository followed by genome quality
        sorted_refseq_rep_genomes = []
        sorted_genbank_rep_genomes = []
        sorted_trusted_user_rep_genomes = []
        sorted_by_quality = sorted(genome_quality.items(), 
                                    key=operator.itemgetter(1), 
                                    reverse=True)
        for genome_id, _quality in sorted_by_quality:
            if genome_id.startswith('RS_'):
                sorted_refseq_rep_genomes.append(genome_id)
            elif genome_id.startswith('GB_'):
                sorted_genbank_rep_genomes.append(genome_id)
            elif genome_id.startswith('U_'):
                sorted_trusted_user_rep_genomes.append(genome_id)
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
        
    def _calculate_ani(self, gene_file1, gene_file2):
        """Calculate ANI between genomes."""
        
        outdir = tempfile.mkdtemp()
        outfile = os.path.join(outdir, 'ani_calculator')
        cmd = 'ani_calculator -genome1fna %s -genome2fna %s -outfile %s -outdir %s > /dev/null 2>&1' % (gene_file1, 
                                                                                                        gene_file2, 
                                                                                                        outfile, 
                                                                                                        outdir)
        os.system(cmd)
        
        with open(outfile) as f:
            f.readline()
            results = f.readline().strip().split('\t')
            ani = 0.5*(float(results[2]) + float(results[3]))
            af = 0.5*(float(results[4]) + float(results[5]))
        shutil.rmtree(outdir)
        
        return ani, af
        
    def _cluster_species(self, gene_file_rep, gene_file_query, rep_gtdb_sp, query_gtdb_sp, query_ncbi_sp):
        """Determine if genomes should be clustered based on their ANI and species assignments."""

        if query_gtdb_sp != 's__' and rep_gtdb_sp != query_gtdb_sp:
            # genomes belong to different GTDB species so
            # should not be clustered together
            return False
            
        ani, af = self._calculate_ani(gene_file_rep, gene_file_query)
        
        if af < self._af_threshold:
            # insufficient alignment fraction for genomes to
            # be considered from the same species
            return False

        if ani >= self.strict_ani_threshold and query_gtdb_sp == 's__' and query_ncbi_sp == 's__':
                # genome meets the strict threshold for
                # clustering and has no assigned species
                return True
                
        if rep_gtdb_sp == 's__' or rep_gtdb_sp != query_gtdb_sp:
            return False
                
        if ani >= self.gtdb_species_ani_threshold:
                # genomes are from same named GTDB species and 
                # meet the threshold for clustering
                return True
        elif (ani >= self.ncbi_species_ani_threshold 
                and canonical_species_name(rep_gtdb_sp) == query_ncbi_sp):
                # NCBI species of query genome is the same as the GTDB representative
                # and they meet the threshold for clustering
                return True
                
        return False

    def _greedy_representatives(self,
                                representatives,
                                ordered_genomes,
                                gtdb_taxonomy,
                                ncbi_taxonomy,
                                mash_pairwise_file,
                                gene_files):
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
        self.logger.info('Performing greedy clustering.')
        total_genomes = len(ordered_genomes)
        processed_genomes = 0
        while len(ordered_genomes):
            processed_genomes += 1
            if processed_genomes % 100 == 0:
                sys.stdout.write('==> Processed %d of %d genomes.\r' % (processed_genomes, total_genomes))
                sys.stdout.flush()

            gid = ordered_genomes.pop(0)
            query_dists = mash_dists[gid]
            
            query_gtdb_sp = gtdb_taxonomy[gid][6]
            query_ncbi_sp = ncbi_taxonomy[gid][6]
            
            query_rep_dists = []
            for rep_id in representatives:
                mash_dist = query_dists.get(rep_id, 1.0)
                if mash_dist <= self.mash_dist_threshold:
                    # estimated ANI value is sufficently small that it should be verified
                    query_rep_dists.append([rep_id, mash_dist])

            assigned_rep = False
            for rep_id, mash_dist in sorted(query_rep_dists, key = lambda x: x[1]):
                rep_gtdb_sp = gtdb_taxonomy[rep_id][6]
                if query_gtdb_sp != 's__' and rep_gtdb_sp != query_gtdb_sp:
                    continue
  
                if mash_dist < 0.035 and query_gtdb_sp == 's__':
                        # genome meets the strict threshold for
                        # clustering and has no assigned species
                        assigned_rep = True
                        break
                        
                if rep_gtdb_sp == 's__' or rep_gtdb_sp != query_gtdb_sp:
                    continue
                        
                if mash_dist < 0.05:
                        # genomes are from same named GTDB species and 
                        # meet the threshold for clustering
                        assigned_rep = True
                        break
                elif (mash_dist < 0.1
                        and canonical_species_name(rep_gtdb_sp) == canonical_species_name(query_ncbi_sp)):
                        # NCBI species of query genome is the same as the GTDB representative
                        # and they meet the threshold for clustering
                        assigned_rep = True
                        break
                
                if False:
                    rep_gtdb_sp = gtdb_taxonomy[rep_id][6]
                    cluster = self._cluster_species(gene_files[rep_id], 
                                                    gene_files[gid], 
                                                    rep_gtdb_sp, 
                                                    query_gtdb_sp, 
                                                    query_ncbi_sp)
                    if cluster:
                        assigned_rep = True
                        break

            if not assigned_rep:
                # genome was not assigned to an existing representative,
                # so make it a new representative genome
                representatives.add(gid)

        sys.stdout.write('==> Processed %d of %d genomes.\r' % (processed_genomes, total_genomes))
        sys.stdout.flush()
        sys.stdout.write('\n')

        return representatives

    def representatives(self,
                        species_derep_file,
                        metadata_file,
                        prev_rep_file,
                        mash_pairwise_file,
                        genome_dir_file,
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
        Clustering is based on conservative Mash distance threshold followed by
        ANI values calculated with ANI Calculator (Varghese et al, 2015).

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
        
        # get path to all nucleotide gene files of all genomes
        gene_files = {}
        for line in open(genome_dir_file):
            line_split = line.strip().split('\t')
            
            gtdb_genome_id = line_split[0]
            genome_id = gtdb_genome_id.replace('GB_', '').replace('RS_', '')
            genome_path = line_split[1]
            gene_files[gtdb_genome_id] = os.path.join(genome_path, 'prodigal', genome_id + '_protein.fna')
            
        self.logger.info('Read genome path for %d genomes.' % len(gene_files))
        
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
        for gid, stats in genome_stats.iteritems():
            if gid in init_rep_genomes:
                continue
                
            if gid.startswith('U_') and gid not in trusted_user_genomes:
                continue
                
            if (stats.checkm_completeness >= min_rep_comp 
                and stats.checkm_contamination <= max_rep_cont 
                and (stats.checkm_completeness - 5*stats.checkm_contamination) >= min_quality
                and stats.contig_count <= max_contigs
                and stats.n50_scaffolds >= min_N50
                and stats.ambiguous_bases <= max_ambiguous
                and stats.total_gap_length <= max_gap_length):
                    genome_quality[gid] = self._genome_quality(stats, gid in prev_gtdb_reps)

        # perform greedy identification of new representatives
        ordered_genomes = self._order_genomes(genome_quality)
        info = 'Comparing %d genomes to %d initial representatives.' % (len(ordered_genomes),
                                                                            len(init_rep_genomes))
        self.logger.info(info)
        gtdb_taxonomy = read_gtdb_taxonomy(metadata_file)
        ncbi_taxonomy = read_gtdb_ncbi_taxonomy(metadata_file)
        representatives = self._greedy_representatives(init_rep_genomes,
                                                        ordered_genomes,
                                                        gtdb_taxonomy,
                                                        ncbi_taxonomy,
                                                        mash_pairwise_file,
                                                        gene_files)

        self.logger.info('Identified %d representatives.' % len(representatives))

        # read metadata for genomes
        (refseq_genomes, 
            complete_genomes, 
            ncbi_reference_genomes,
            ncbi_representative_genomes) = ncbi.read_refseq_metadata(metadata_file)

        # write out information for representative genomes
        fout = open(output_file, 'w')

        fout.write('# Selection criteria:\n')
        fout.write('# Species dereplication file: %s\n' % species_derep_file)
        fout.write('# Previous representative file: %s\n' % prev_rep_file)
        fout.write('# Trusted user genomes file: %s\n' % trusted_user_file)
        fout.write('# Genome quality metadata file: %s\n' % str(metadata_file))
        fout.write('# Min. representative completeness: %.2f\n' % min_rep_comp)
        fout.write('# Max. representative contamination: %.2f\n' % max_rep_cont)
        fout.write('# Strict ANI threshold: %.3f\n' % self.strict_ani_threshold)
        fout.write('# GTDB species ANI threshold: %.3f\n' % self.gtdb_species_ani_threshold)
        fout.write('# NCBI species ANI threshold: %.3f\n' % self.ncbi_species_ani_threshold)
        fout.write('#\n')

        fout.write('# Genome Id\tGTDB Taxonomy\tNCBI Taxonomy\tNCBI Organism Name\tType strain\tType sources\tComplete\tNCBI reference genome\tNCBI representative genome\n')
        for gid in representatives:
            reference = 'yes' if gid in ncbi_reference_genomes else 'no'
            representative = 'yes' if gid in ncbi_representative_genomes else 'no'
            complete = 'yes' if gid in complete_genomes else 'no'
            ts = 'yes' if genome_stats[gid].gtdb_type_material in self.true_str else 'no'
            gtdb_taxa_str = ';'.join(gtdb_taxonomy.get(gid, Taxonomy.rank_prefixes))
            ncbi_taxa_str = ';'.join(ncbi_taxonomy.get(gid, Taxonomy.rank_prefixes))

            fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (gid,
                                                                gtdb_taxa_str,
                                                                ncbi_taxa_str,
                                                                genome_stats[gid].ncbi_organism_name,
                                                                ts,
                                                                genome_stats[gid].gtdb_type_material_sources,
                                                                complete,
                                                                reference,
                                                                representative))

        fout.close()
        
    def cluster(self,
                rep_genome_file,
                metadata_file,
                mash_pairwise_file,
                genome_dir_file,
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
        
        # get path to all nucleotide gene files of all genomes
        gene_files = {}
        for line in open(genome_dir_file):
            line_split = line.strip().split('\t')
            
            gtdb_genome_id = line_split[0]
            genome_id = gtdb_genome_id.replace('GB_', '').replace('RS_', '')
            genome_path = line_split[1]
            gene_files[gtdb_genome_id] = os.path.join(genome_path, 'prodigal', genome_id + '_protein.fna')
            
        self.logger.info('Read genome path for %d genomes.' % len(gene_files))
        
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
        for i, gid in enumerate(remaining_genomes):
            if i % 100 == 0:
                sys.stdout.write('==> Processed %d of %d genomes.\r' % (i+1, len(remaining_genomes)))
                sys.stdout.flush()
                
            query_dists = mash_dists[gid]
            
            query_gtdb_sp = gtdb_taxonomy[gid][6]
            query_ncbi_sp = ncbi_taxonomy[gid][6]
            
            assigned_rep = False
            for rep_id, mash_dist in sorted(query_dists.items(), key=operator.itemgetter(1)):
                if mash_dist > self.mash_dist_threshold:
                    # estimated ANI value is sufficently large
                    # that it is almost certainly not of interest
                    break
                    
                rep_gtdb_sp = gtdb_taxonomy[rep_id][6]
                if query_gtdb_sp != 's__' and rep_gtdb_sp != query_gtdb_sp:
                    continue
  
                if mash_dist < 0.035 and query_gtdb_sp == 's__':
                        # genome meets the strict threshold for
                        # clustering and has no assigned species
                        assigned_rep = True
                        break
                        
                if rep_gtdb_sp == 's__' or rep_gtdb_sp != query_gtdb_sp:
                    continue
                        
                if mash_dist < 0.05:
                        # genomes are from same named GTDB species and 
                        # meet the threshold for clustering
                        assigned_rep = True
                        break
                elif (mash_dist < 0.1
                        and canonical_species_name(rep_gtdb_sp) == canonical_species_name(query_ncbi_sp)):
                        # NCBI species of query genome is the same as the GTDB representative
                        # and they meet the threshold for clustering
                        assigned_rep = True
                        break
                
                if False:
                    cluster = self._cluster_species(gene_files[rep_id], 
                                                    gene_files[gid], 
                                                    rep_gtdb_sp, 
                                                    query_gtdb_sp, 
                                                    query_ncbi_sp)
                    if cluster:
                        assigned_rep = rep_id
           
            if assigned_rep:
                clusters[assigned_rep].append(gid)
                
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