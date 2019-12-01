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
import ntpath
import pickle
from itertools import combinations, product
from collections import defaultdict, namedtuple

from biolib.taxonomy import Taxonomy
from biolib.external.execute import check_dependencies

from numpy import (mean as np_mean,
                    std as np_std)

from gtdb_species_clusters.common import read_gtdb_metadata
                                    
from gtdb_species_clusters.genome_utils import (read_genome_path,
                                                exclude_from_refseq,
                                                read_qc_file)
                                    
from gtdb_species_clusters.taxon_utils import (binomial_species,
                                                genome_species_assignments,
                                                read_gtdb_taxonomy,
                                                read_gtdb_ncbi_taxonomy)
                                    
from gtdb_species_clusters.type_genome_utils import (NCBI_TYPE_SPECIES,
                                            NCBI_PROXYTYPE,
                                            NCBI_TYPE_SUBSP,
                                            GTDB_TYPE_SPECIES,
                                            GTDB_TYPE_SUBSPECIES,
                                            GTDB_NOT_TYPE_MATERIAL,
                                            check_ncbi_subsp,
                                            gtdb_type_strain_of_species,
                                            symmetric_ani,
                                            quality_score)
                                    
from gtdb_species_clusters.fastani import FastANI
from gtdb_species_clusters.mash import Mash

class SelectTypeGenomes(object):
    """Select GTDB type genomes for named species."""

    def __init__(self, ani_cache_file, cpus, output_dir):
        """Initialization."""
        
        check_dependencies(['fastANI', 'mash'])
        
        self.cpus = cpus
        self.output_dir = output_dir

        self.logger = logging.getLogger('timestamp')

        self.min_intra_strain_ani = 99.0
        self.min_mash_ani = 90.0
        
        self.max_ani_neighbour = 97.0
        
        self.fastani = FastANI(ani_cache_file, cpus)
        
        self.BlastHit = namedtuple('BlastHit', ['ltp_species', 'ssu_len', 'align_len', 'perc_identity', 'bitscore', 'evalue'])
        
    def  _type_metadata(self, metadata_file):
        """Read and parse type material metadata."""
        
        type_metadata = read_gtdb_metadata(metadata_file, ['ncbi_taxonomy_unfiltered',
                                                                'ncbi_refseq_category',
                                                                'ncbi_strain_identifiers',
                                                                'ncbi_type_material_designation',
                                                                'gtdb_type_designation',
                                                                'gtdb_type_designation_sources',
                                                                'gtdb_type_species_of_genus'])
                                                                
        for gid in type_metadata:
            strain_ids = set()
            if type_metadata[gid].ncbi_strain_identifiers:
                s = [s.strip() for s in str(type_metadata[gid].ncbi_strain_identifiers).split(';')]
                for strain_id in s:
                    if strain_id:
                        strain_ids.add(strain_id)
                    
            type_metadata[gid] = type_metadata[gid]._replace(ncbi_strain_identifiers = strain_ids)
            
            if type_metadata[gid].gtdb_type_designation_sources:
                sources = set([t.strip() for t in type_metadata[gid].gtdb_type_designation_sources.split(';')])
                type_metadata[gid] = type_metadata[gid]._replace(gtdb_type_designation_sources = sources)
            else:
                type_metadata[gid] = type_metadata[gid]._replace(gtdb_type_designation_sources = set())
                
        return type_metadata
        
    def _parse_ltp_blast_table(self, ltp_blast_file, ncbi_taxonomy):
        """Parse Living Tree Project (LTP) BLAST result to identify genomes that are likely type material."""
        
        ltp_type_species_of_genus = defaultdict(set)
        ltp_type_strain_of_species = defaultdict(set)
        ltp_type_strain_of_subspecies = defaultdict(set)
        ltp_top_blast_hit = {}

        with open(ltp_blast_file) as f:
            header = f.readline().strip().split('\t')
            
            accession_index = header.index('accession')
            taxonomy_index = header.index('taxonomy')
            length_index = header.index('length')
            align_len_index = header.index('align_len')
            perc_identity_index = header.index('perc_identity')
            bitscore_index = header.index('bitscore')
            evalue_index = header.index('evalue')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gid = line_split[accession_index]
                if gid.startswith('GCF_'):
                    gid = 'RS_' + gid
                elif gid.startswith('GCA_'):
                    gid = 'GB_' + gid
                    
                taxonomy_info = line_split[taxonomy_index]
                taxonomy, strain_info = taxonomy_info.split('|')
                taxa = [t.strip() for t in taxonomy.split(';')]
                type_info = taxa[-1]
                ltp_sp = 's__' + taxa[-2]
                ltp_sp_canonical = ltp_sp
                if 'subsp.' in ltp_sp_canonical:
                    ltp_sp_canonical = ltp_sp_canonical[0:ltp_sp_canonical.find('subsp.')].strip()
                
                length = int(line_split[length_index])
                align_len = int(line_split[align_len_index])
                perc_identity = float(line_split[perc_identity_index])
                bitscore = float(line_split[bitscore_index])
                evalue = float(line_split[evalue_index])

                if gid not in ltp_top_blast_hit or bitscore > ltp_top_blast_hit[gid].bitscore:
                    ltp_top_blast_hit[gid] = self.BlastHit(ltp_species=ltp_sp,
                                                            ssu_len=length,
                                                            align_len=align_len,
                                                            perc_identity=perc_identity,
                                                            bitscore=bitscore,
                                                            evalue=evalue)
                
                if perc_identity >= 99.8 and align_len*100.0/length >= 99.8:
                    if gid not in ncbi_taxonomy:
                        self.logger.warning('Genome %s is missing NCBI taxonomy information.' % gid)
                        continue
                        
                    ncbi_sp = ncbi_taxonomy[gid][6]
                    if ltp_sp_canonical == ncbi_sp:
                        if "l[T]" in strain_info: # type strain according to LTP
                            if type_info == "type sp.":
                                ltp_type_species_of_genus[ncbi_sp].add(gid)
                                
                            if "subsp." in ltp_sp:
                                ltp_type_strain_of_subspecies[ncbi_sp].add(gid)
                            else:
                                ltp_type_strain_of_species[ncbi_sp].add(gid)

        return ltp_type_species_of_genus, ltp_type_strain_of_species, ltp_type_strain_of_subspecies, ltp_top_blast_hit
        
    def _genome_quality(self, metadata_file):
        """Calculate quality of genomes."""

        quality_metadata = read_gtdb_metadata(metadata_file, ['checkm_completeness',
                                                                'checkm_contamination',
                                                                'checkm_strain_heterogeneity_100',
                                                                'contig_count',
                                                                'scaffold_count',
                                                                'n50_contigs',
                                                                'ambiguous_bases',
                                                                'gtdb_taxonomy',
                                                                'ssu_count',
                                                                'ssu_length',
                                                                'total_gap_length',
                                                                'ncbi_assembly_level',
                                                                'ncbi_genome_representation',
                                                                'ncbi_molecule_count',
                                                                'ncbi_unspanned_gaps',
                                                                'ncbi_spanned_gaps',
                                                                'ncbi_refseq_category',
                                                                'ncbi_type_material_designation',
                                                                'mimag_high_quality',
                                                                'genome_size',
                                                                'ncbi_genome_category'])
                                                    
        quality = quality_score(quality_metadata.keys(), quality_metadata)
 
        return quality, quality_metadata
        
    def _select_highest_quality(self, gids, genome_quality):
        """Select highest quality genome."""
                    
        q = {k:genome_quality[k] for k in gids}
        q_sorted = sorted(q.items(), key=operator.itemgetter(1), reverse=True)
        return q_sorted[0][0]
        
    def _select_ani_neighbours(self, species, gids, genome_quality, ani_af):
        """Select highest-quality genome with sufficient number of ANI neighbours."""
        
        # calculate mean ANI of all genome pairs
        anis = []
        for gid1, gid2 in combinations(gids, 2):
            ani, af = symmetric_ani(ani_af, gid1, gid2)
            if ani > 0:
                anis.append(ani)
                
        if not anis:
            self.logger.warning('Could not calculate ANI between %d genomes in %s.' % (len(gids), species))
            self.logger.warning('Selecting highest-quality genome.')
            return self._select_highest_quality(gids, genome_quality)
                    
        # calculate number of ANI neighbours for each genome
        mean_ani = np_mean(anis)
        std_ani = np_std(anis)
        ani_neighbours = defaultdict(int)
        for gid1, gid2 in product(gids, repeat=2):
            if gid1 == gid2:
                ani_neighbours[gid1] += 1
            else:
                ani, af = symmetric_ani(ani_af, gid1, gid2)
                if ani >= mean_ani - std_ani:
                    ani_neighbours[gid1] += 1
                    
        # get all genomes that are neighbours with at least half the other genomes
        neighbour_gids = []
        for gid in ani_neighbours:
            if ani_neighbours[gid] >= 0.5*len(gids):
                neighbour_gids.append(gid)
                
        if len(neighbour_gids) == 0:
            self.logger.error('No ANI neighbours identified in _select_ani_neighbours')
            sys.exit(-1)

        return self._select_highest_quality(neighbour_gids, genome_quality)

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

    def _get_type_designations(self, ncbi_taxonomy, type_metadata, passed_qc):
        """Get type material designations for each genome."""

        gtdb_type_genus = defaultdict(set)  # type species of genus
        gtdb_type_sp = defaultdict(set)     # type strain of species
        gtdb_type_subsp = defaultdict(set)  # type strain of subspecies
        ncbi_type_sp = defaultdict(set)
        ncbi_proxy = defaultdict(set)
        ncbi_type_subsp = defaultdict(set)
        ncbi_reps = defaultdict(set)
        
        for gid, metadata in type_metadata.items():
            if gid not in passed_qc:
                continue
                
            ncbi_taxa = ncbi_taxonomy[gid]
            ncbi_species = ncbi_taxa[6]
            
            # sanity check NCBI species assignments for all designated GTDB type material
            if (metadata.gtdb_type_designation 
                and metadata.gtdb_type_designation not in GTDB_NOT_TYPE_MATERIAL
                and ncbi_species == 's__'):
                self.logger.error('Missing NCBI species name for GTDB designated type material: %s (%s)' % (
                                    gid, 
                                    metadata.gtdb_type_designation))
                sys.exit(-1)

            # determine type status of genome
            if metadata.gtdb_type_species_of_genus:
                gtdb_type_genus[ncbi_species].add(gid)
                
            if metadata.gtdb_type_designation:
                if metadata.gtdb_type_designation in GTDB_TYPE_SPECIES:
                    gtdb_type_sp[ncbi_species].add(gid)
                elif metadata.gtdb_type_designation in GTDB_TYPE_SUBSPECIES:
                    gtdb_type_subsp[ncbi_species].add(gid)
                elif metadata.gtdb_type_designation in GTDB_NOT_TYPE_MATERIAL:
                    pass # not type material
                else:
                    self.logger.error('Unrecognized GTDB type designation for %s: %s' % (
                                        gid, 
                                        metadata.gtdb_type_designation))
                    sys.exit(-1)
                    
            if metadata.ncbi_type_material_designation:
                if metadata.ncbi_type_material_designation in NCBI_TYPE_SPECIES:
                    if check_ncbi_subsp(metadata.ncbi_taxonomy_unfiltered):
                        # genome is marked as 'assembled from type material', but
                        # is a subspecies according to the NCBI taxonomy
                        ncbi_type_subsp[ncbi_species].add(gid)
                    else:
                        ncbi_type_sp[ncbi_species].add(gid)
                elif metadata.ncbi_type_material_designation in NCBI_PROXYTYPE:
                    ncbi_proxy[ncbi_species].add(gid)
                elif metadata.ncbi_type_material_designation in NCBI_TYPE_SUBSP:
                    ncbi_type_subsp[ncbi_species].add(gid)
                else:
                    self.logger.error('Unrecognized NCBI type designation for %s: %s' % (
                                        gid,
                                        metadata.ncbi_type_material_designation))
                    sys.exit(-1)

            if metadata.ncbi_refseq_category and metadata.ncbi_refseq_category != 'na':
                ncbi_reps[ncbi_species].add(gid)
                assert('reference' in metadata.ncbi_refseq_category or 'representative' in metadata.ncbi_refseq_category)
                
        # clean up NCBI representative genomes so only the representative of 
        # a species and not subspecies is retained in multiple representatives exist
        for sp, gids in ncbi_reps.items():
            sp_gids = set()
            for gid in gids:
                if not check_ncbi_subsp(type_metadata[gid].ncbi_taxonomy_unfiltered):
                    sp_gids.add(gid)
                    
            if len(sp_gids) > 0:
                ncbi_reps[sp] = sp_gids
            
        # sanity check
        if 's__' in gtdb_type_sp:
            self.logger.error('Type strain of species has no NCBI species assignment:')
            self.logger.error('%s' % ','.join([t[0] for t in gtdb_type_sp['s__']]))
        if 's__' in gtdb_type_subsp:
            self.logger.error('Type strain of subspecies has no NCBI species assignment:')
            self.logger.error('%s' % ','.join([t[0] for t in gtdb_type_subsp['s__']]))
            
        return (gtdb_type_genus,
                gtdb_type_sp, 
                gtdb_type_subsp, 
                ncbi_type_sp,
                ncbi_proxy,
                ncbi_type_subsp,
                ncbi_reps)

    def _type_sources(self, type_metadata, gids):
        """Count number of genomes indicated as type material from each source."""
        
        type_sources = defaultdict(int)
        strain_ids = set()
        for gid in gids:
            if type_metadata[gid].gtdb_type_designation_sources:
                for source in type_metadata[gid].gtdb_type_designation_sources:
                    type_sources[source] += 1

            if type_metadata[gid].ncbi_strain_identifiers:
                strain_ids.update(type_metadata[gid].ncbi_strain_identifiers)
                
        return type_sources, strain_ids
        
    def _validate_type_designations(self,
                                    type_metadata,
                                    ncbi_taxonomy, 
                                    gtdb_type_sp, 
                                    gtdb_type_subsp,
                                    ncbi_type_sp,
                                    ncbi_proxy,
                                    ncbi_type_subsp,
                                    ncbi_reps):
        """Produce files indicating potentially missing type material."""
        
        # determine all binomial NCBI species names
        ncbi_species = binomial_species(ncbi_taxonomy)
        
        # identify species without a type strain
        fout = open(os.path.join(self.output_dir, 'validate_type_strain_of_species.tsv'), 'w')
        fout.write('NCBI species\tType designation\tNo. genomes with designation\tGTDB type sources\tNCBI type designations\tAccession(s)\n')
        missing_type_strain = 0
        missing_type_strain_ncbi_type = 0
        not_ncbi_type_sp = 0
        for sp in ncbi_species:
            ncbi_types = defaultdict(int)
            for gid in ncbi_type_sp.get(sp, []):
                ncbi_types[type_metadata[gid].ncbi_type_material_designation] += 1

            type_sources, strain_ids = self._type_sources(type_metadata, gtdb_type_sp.get(sp, []))
            
            if sp in gtdb_type_sp:
                fout.write('%s\t%s\t%d\t%s\t%s\t%s\n' % (sp,
                                                        'type strain of species',
                                                        len(gtdb_type_sp[sp]),
                                                        ', '.join("%s=%r" % (key,val) for (key,val) in type_sources.iteritems()),
                                                        ', '.join("%s=%r" % (key,val) for (key,val) in ncbi_types.iteritems()),
                                                        ', '.join([gid for d in gtdb_type_sp[sp]])))
                if sp not in ncbi_type_sp:
                    not_ncbi_type_sp += 1
            elif sp in ncbi_type_sp:
                fout.write('%s\t%s\t%d\t%s\t%s\t%s\n' % (sp,
                                                        'assembly from type material',
                                                        len(ncbi_type_sp[sp]),
                                                        ', '.join("%s=%r" % (key,val) for (key,val) in type_sources.iteritems()),
                                                        ', '.join("%s=%r" % (key,val) for (key,val) in ncbi_types.iteritems()),
                                                        ', '.join(ncbi_type_sp[sp])))
                missing_type_strain_ncbi_type += 1
            elif sp in ncbi_proxy:
                fout.write('%s\t%s\t%d\t%s\t%s\t%s\n' % (sp,
                                                        'assembly from proxytype material',
                                                        len(ncbi_proxy[sp]),
                                                        '',
                                                        '',
                                                        ', '.join(ncbi_proxy[sp])))
            elif sp in ncbi_reps:
                fout.write('%s\t%s\t%d\t%s\t%s\t%s\n' % (sp,
                                                        'representative at NCBI',
                                                        len(ncbi_reps[sp]),
                                                        '',
                                                        '',
                                                        ', '.join(ncbi_reps[sp])))
            elif sp in gtdb_type_subsp or sp in ncbi_type_subsp:
                gids = gtdb_type_subsp[sp].union(ncbi_type_subsp[sp])
                fout.write('%s\t%s\t%d\t%s\t%s\t%s\n' % (sp,
                                                        'type strain of subspecies or assembly from synonym type material',
                                                        len(gids),
                                                        '',
                                                        '',
                                                        ', '.join(gids)))
            else:
                fout.write('%s\t%s\t%d\t%s\t%s\t%s\n' % (sp,
                                                        'no type genome',
                                                        len(ncbi_species[sp]),
                                                        '',
                                                        '',
                                                        ', '.join(ncbi_species[sp])))
                                                        
                missing_type_strain += 1
                    
        fout.close()
        
        self.logger.info('Identified %d NCBI species.' % len(ncbi_species))
        self.logger.info('Identified %d species without a type strain of species, type strain of subspecies, or representative genome.' % missing_type_strain)
        self.logger.info('Identified %d species without a GTDB designated type strain of species that have genome(s) designated as assembled from type material at NCBI.' % missing_type_strain_ncbi_type)
        self.logger.info('Identified %d species with a GTDB designated type strain of species that are not designated as assembled from type material at NCBI.' % not_ncbi_type_sp)
        
    def _select_type_genomes(self,
                                passed_qc,
                                genome_files,
                                genome_quality,
                                quality_metadata,
                                type_metadata,
                                gtdb_type_sp, 
                                gtdb_type_subsp,
                                ncbi_type_sp,
                                ncbi_proxy,
                                ncbi_type_subsp,
                                ncbi_reps,
                                ltp_top_blast_hit,
                                ncbi_taxonomy,
                                gtdb_taxonomy,
                                excluded_from_refseq_note,
                                manual_type_genomes):
        """Select type genome for each species."""

        # determine all NCBI species names
        ncbi_species = binomial_species(ncbi_taxonomy)

        # identify species without a type strain
        self.logger.info('Selecting type genome for each of the %d species.' % len(ncbi_species))
        
        fout = open(os.path.join(self.output_dir, 'gtdb_type_genomes_initial.tsv'), 'w')
        fout.write('NCBI species\tType genome\tStrain IDs\tType status\tType sources\tNCBI assembly types\tNCBI representative\tNCBI assembly level')
        fout.write('\tNCBI genome category\tGenome size (bp)\tQuality score\tCompleteness (%)\tContamination (%)\tNo. scaffolds\tNo. contigs\tN50 contigs\tAmbiguous bases\tSSU count\tSSU length (bp)')
        fout.write('\tLTP species\tBLAST alignment length (bp)\tBLAST percent identity\tBLAST bitscore\tBLAST e-value')
        fout.write('\tNo. type genomes\tNo. species genomes\tMean ANI\tMean AF\tMin ANI\tMin AF\tNCBI exclude from RefSeq\tSelection note\n')
        
        fout_manual = open(os.path.join(self.output_dir, 'gtdb_type_genomes.manual.tsv'), 'w')
        fout_manual.write('NCBI species\tAccession\tGTDB taxonomy\tNCBI taxonomy\tGTDB type genome\tType status\tANI to type genome\tAF to type genome\tSelection note')
        fout_manual.write('\tNCBI genome category\tGenome size (bp)\tQuality score\tCompleteness (%)\tContamination (%)\tNo. scaffolds\tNo. contigs\tN50 contigs\tAmbiguous bases\tSSU count\tSSU length (bp)')
        fout_manual.write('\tLTP species\tBLAST alignment length (bp)\tBLAST percent identity\tBLAST bitscore\tBLAST e-value')
        fout_manual.write('\tNo. type genomes\tNo. species genomes\tMean ANI\tMean AF\tMin ANI\tMin AF\tNCBI exclude from RefSeq\tType accessions\tSpecies accessions\n')

        num_type_strain_of_species = 0
        num_type_strain_of_subspecies = 0
        num_ncbi_assembled_from_type = 0
        num_ncbi_assembled_from_proxytype = 0
        num_ncbi_rep = 0
        num_de_novo = 0
        
        sp_manual_curation = []
        num_type_strain_of_species_manual = 0
        num_type_strain_of_subspecies_manual = 0
        num_ncbi_assembled_from_type_manual = 0
        num_ncbi_assembled_from_proxytype_manual = 0
        num_ncbi_rep_manual = 0
        num_de_novo_manual = 0
        type_genomes = {}
        multi_gids = 0
        for idx, sp in enumerate(ncbi_species):
            species_gids = ncbi_species[sp].intersection(passed_qc)
            if len(species_gids) == 0:
                continue
            
            statusStr = '-> Processing %d of %d (%.2f%%) species [%s: %d].'.ljust(86) % (idx+1, 
                                                                                len(ncbi_species), 
                                                                                float(idx+1)*100/len(ncbi_species),
                                                                                sp,
                                                                                len(species_gids))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            if sp in gtdb_type_sp:
                gid, manual_inspection = self._select_type_genome(sp,
                                                'type strain of species',
                                                gtdb_type_sp[sp],
                                                ncbi_type_sp,
                                                ncbi_reps,
                                                genome_quality,
                                                quality_metadata,
                                                type_metadata,
                                                ltp_top_blast_hit,
                                                excluded_from_refseq_note,
                                                genome_files,
                                                species_gids,
                                                ncbi_taxonomy,
                                                gtdb_taxonomy,
                                                manual_type_genomes,
                                                fout,
                                                fout_manual)
                                                
                if len(gtdb_type_sp[sp]) > 1:
                    multi_gids += 1
                    
                if manual_inspection:
                    num_type_strain_of_species_manual += 1

                num_type_strain_of_species += 1
            elif sp in ncbi_type_sp:
                gid, manual_inspection = self._select_type_genome(sp,
                                                'NCBI assembled from type material',
                                                ncbi_type_sp[sp],
                                                ncbi_type_sp,
                                                ncbi_reps,
                                                genome_quality,
                                                quality_metadata,
                                                type_metadata,
                                                ltp_top_blast_hit,
                                                excluded_from_refseq_note,
                                                genome_files,
                                                species_gids,
                                                ncbi_taxonomy,
                                                gtdb_taxonomy,
                                                manual_type_genomes,
                                                fout,
                                                fout_manual)
                                                
                if len(ncbi_type_sp[sp]) > 1:
                    multi_gids += 1
                    
                if manual_inspection:
                    num_ncbi_assembled_from_type_manual += 1
                
                num_ncbi_assembled_from_type += 1
            elif sp in ncbi_proxy:
                gid, manual_inspection = self._select_type_genome(sp,
                                                'NCBI assembled from proxytype material',
                                                ncbi_proxy[sp],
                                                ncbi_type_sp,
                                                ncbi_reps,
                                                genome_quality,
                                                quality_metadata,
                                                type_metadata,
                                                ltp_top_blast_hit,
                                                excluded_from_refseq_note,
                                                genome_files,
                                                species_gids,
                                                ncbi_taxonomy,
                                                gtdb_taxonomy,
                                                manual_type_genomes,
                                                fout,
                                                fout_manual)
                                                
                if len(ncbi_proxy[sp]) > 1:
                    multi_gids += 1
                    
                if manual_inspection:
                    num_ncbi_assembled_from_proxytype_manual += 1
                
                num_ncbi_assembled_from_proxytype += 1
            elif sp in ncbi_reps:
                gid, manual_inspection = self._select_type_genome(sp,
                                                'NCBI representative genome',
                                                ncbi_reps[sp],
                                                ncbi_type_sp,
                                                ncbi_reps,
                                                genome_quality,
                                                quality_metadata,
                                                type_metadata,
                                                ltp_top_blast_hit,
                                                excluded_from_refseq_note,
                                                genome_files,
                                                species_gids,
                                                ncbi_taxonomy,
                                                gtdb_taxonomy,
                                                manual_type_genomes,
                                                fout,
                                                fout_manual)
                
                if len(ncbi_reps[sp]) > 1:
                    multi_gids += 1
                    
                if manual_inspection:
                    num_ncbi_rep_manual += 1
                    
                num_ncbi_rep += 1
            elif sp in gtdb_type_subsp or sp in ncbi_type_subsp:
                gids = gtdb_type_subsp[sp].union(ncbi_type_subsp[sp])
                gid, manual_inspection = self._select_type_genome(sp,
                                                'type strain of subspecies',
                                                gids,
                                                ncbi_type_sp,
                                                ncbi_reps,
                                                genome_quality,
                                                quality_metadata,
                                                type_metadata,
                                                ltp_top_blast_hit,
                                                excluded_from_refseq_note,
                                                genome_files,
                                                species_gids,
                                                ncbi_taxonomy,
                                                gtdb_taxonomy,
                                                manual_type_genomes,
                                                fout,
                                                fout_manual)
                                                
                if len(gids) > 1:
                    multi_gids += 1
                    
                if manual_inspection:
                    num_type_strain_of_subspecies_manual += 1
                
                num_type_strain_of_subspecies += 1
            else:
                gid, manual_inspection = self._select_type_genome(sp,
                                                'no type genome',
                                                species_gids,
                                                ncbi_type_sp,
                                                ncbi_reps,
                                                genome_quality,
                                                quality_metadata,
                                                type_metadata,
                                                ltp_top_blast_hit,
                                                excluded_from_refseq_note,
                                                genome_files,
                                                species_gids,
                                                ncbi_taxonomy,
                                                gtdb_taxonomy,
                                                manual_type_genomes,
                                                fout,
                                                fout_manual)
                if len(species_gids) > 1:
                    multi_gids += 1
                    
                if manual_inspection:
                    num_de_novo_manual += 1
                    
                num_de_novo += 1
                
            type_genomes[sp] = gid

        sys.stdout.write('\n')
        fout.close()
        fout_manual.close()
        
        self.logger.info('GTDB type genome is type strain of species: %d (%.1f%%)' % (num_type_strain_of_species, 
                                                                                        num_type_strain_of_species*100.0/len(ncbi_species)))
        self.logger.info('GTDB type genome is assembled from type material according to NCBI: %d (%.1f%%)' % (num_ncbi_assembled_from_type, 
                                                                                        num_ncbi_assembled_from_type*100.0/len(ncbi_species)))
        self.logger.info('GTDB type genome is assembled from proxytype material according to NCBI: %d (%.1f%%)' % (num_ncbi_assembled_from_proxytype, 
                                                                                        num_ncbi_assembled_from_proxytype*100.0/len(ncbi_species)))
        self.logger.info('GTDB type genome is a representative genome at NCBI: %d (%.1f%%)' % (num_ncbi_rep, 
                                                                                        num_ncbi_rep*100.0/len(ncbi_species)))
        self.logger.info('GTDB type genome is type strain of subspecies: %d (%.1f%%)' % (num_type_strain_of_subspecies, 
                                                                                        num_type_strain_of_subspecies*100.0/len(ncbi_species)))
        self.logger.info('Species with de novo selected type genome: %d (%.1f%%)' % (num_de_novo, num_de_novo*100.0/len(ncbi_species)))

        self.logger.info('Identified %d species where multiple potential type genomes exist.' % multi_gids)
        self.logger.info('Identified species requiring manual inspection of selected type genome:')
        self.logger.info('  TS = %d; NTS = %d; NP %d; NR = %d; TSS = %d; DN = %d' % (
                            num_type_strain_of_species_manual,
                            num_ncbi_assembled_from_type_manual,
                            num_ncbi_assembled_from_proxytype_manual,
                            num_ncbi_rep_manual,
                            num_type_strain_of_subspecies_manual,
                            num_de_novo_manual))
                                                                                
        return type_genomes

    def _select_type_genome(self,
                            species,
                            type_status,
                            gids,
                            ncbi_type_sp,
                            ncbi_reps,
                            genome_quality,
                            quality_metadata,
                            type_metadata,
                            ltp_top_blast_hit,
                            excluded_from_refseq_note,
                            genome_files,
                            species_gids,
                            ncbi_taxonomy,
                            gtdb_taxonomy,
                            manual_type_genomes,
                            fout, 
                            fout_manual):
        """Select type genome."""
        
        ncbi_types = defaultdict(int)
        ncbi_type_strain_ids = set()
        for gid in ncbi_type_sp.get(species, []):
            ncbi_types[type_metadata[gid].ncbi_type_material_designation] += 1
            ncbi_type_strain_ids.update(type_metadata[gid].ncbi_strain_identifiers)
                
        ncbi_rep_categories = defaultdict(int)
        ncbi_rep_strain_ids = set()
        for gid in ncbi_reps.get(species, []):
            ncbi_rep_categories[type_metadata[gid].ncbi_refseq_category] += 1
            ncbi_rep_strain_ids.update(type_metadata[gid].ncbi_strain_identifiers)

        # select type genome
        mean_ani = mean_af = min_ani = min_af = 'n/a'
        note = ''
        require_manual_inspection = False
        if len(gids) == 1:
            gid = next(iter(gids))
            note = 'select single genome'
        else:
            # check if single 'assembled from type material' genome selected by NCBI
            ncbi_type_gids = gids.intersection(ncbi_type_sp[species])
            if len(ncbi_type_gids) == 1:
                gid = ncbi_type_gids.pop()
                note = 'select single genome annotated as assembled from type material at NCBI'
            else:
                # check for manually selected type genomes
                if species in manual_type_genomes:
                    gid = manual_type_genomes[species]
                    if gid not in gids:
                        self.logger.error('Manually selected type genome specified for %s, but genome %s not in genome list.' % (species, gid))
                    note = 'manually selected type genome'
                else:
                    # calculate ANI between genomes
                    ani_af = self.fastani.pairwise(gids, genome_files)
                    anis = []
                    afs = []
                    for q in ani_af:
                        anis += [d[0] for d in ani_af[q].values()]
                        afs += [d[1] for d in ani_af[q].values()]
                    
                    if anis:
                        mean_ani = '%.1f' % np_mean(anis)
                        mean_af = '%.2f' % np_mean(afs)
                        min_ani = '%.1f' % min(anis)
                        min_af = '%.2f' % min(afs)
                    else:
                        mean_ani = mean_af = min_ani = min_af = 'n/a'

                    gid = self._select_highest_quality(gids, genome_quality)
                    note = 'selected highest-quality genome'
                    
                    if float(min_ani) < self.min_intra_strain_ani:
                        # check if NCBI has designated a reference or representative genome
                        ncbi_rep_gids = gids.intersection(ncbi_reps[species])
                        if len(ncbi_rep_gids) == 1:
                            gid = ncbi_rep_gids.pop()
                            note = 'selected single NCBI representative genome'
                        else:
                            require_manual_inspection = True
                            gid = self._select_ani_neighbours(species, gids, genome_quality, ani_af)
                            note = 'selected highest-quality genome with sufficient ANI neighbours'
                            
                            for cur_gid in [gid] + list(gids.difference([gid])): # write results with selected genome first
                                fout_manual.write('%s\t%s\t%s\t%s\t%s\t%s' % (
                                                    species, 
                                                    cur_gid,
                                                    '; '.join(gtdb_taxonomy[cur_gid]),
                                                    '; '.join(ncbi_taxonomy[cur_gid]),
                                                    cur_gid==gid, 
                                                    type_status))
                                if cur_gid != gid:
                                    if cur_gid in ani_af and gid in ani_af[cur_gid]:
                                        cur_ani, cur_af = ani_af[cur_gid][gid]
                                    else:
                                        cur_ani, cur_af = 0.0, 0.0
                                    fout_manual.write('\t%.1f\t%.2f' % (cur_ani, cur_af))
                                else:
                                    fout_manual.write('\t%.1f\t%.2f' % (100.0, 1.0))
                                fout_manual.write('\t%s\t%s\t%d\t%.2f\t%.2f\t%.2f\t%d\t%d\t%.1f\t%d\t%d\t%d' % (
                                                        note,
                                                        quality_metadata[cur_gid].ncbi_genome_category,
                                                        quality_metadata[cur_gid].genome_size,
                                                        genome_quality[cur_gid], 
                                                        quality_metadata[cur_gid].checkm_completeness,
                                                        quality_metadata[cur_gid].checkm_contamination,
                                                        quality_metadata[cur_gid].scaffold_count,
                                                        quality_metadata[cur_gid].contig_count,
                                                        quality_metadata[cur_gid].n50_contigs,
                                                        quality_metadata[cur_gid].ambiguous_bases,
                                                        quality_metadata[cur_gid].ssu_count,
                                                        quality_metadata[cur_gid].ssu_length if quality_metadata[cur_gid].ssu_length else 0))
                                if cur_gid in ltp_top_blast_hit:
                                    fout_manual.write('\t%s\t%d\t%.2f\t%.2g\t%.2g' % (
                                                            ltp_top_blast_hit[cur_gid].ltp_species,
                                                            ltp_top_blast_hit[cur_gid].align_len,
                                                            ltp_top_blast_hit[cur_gid].perc_identity,
                                                            ltp_top_blast_hit[cur_gid].bitscore,
                                                            ltp_top_blast_hit[cur_gid].evalue))
                                else:
                                    fout_manual.write('\t%s\t%d\t%.2f\t%.2g\t%s' % ('N/A', 0, 0, 0, 'N/A'))
                                fout_manual.write('\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                                                        len(gids),
                                                        len(species_gids),
                                                        mean_ani, mean_af,
                                                        min_ani, min_af,
                                                        excluded_from_refseq_note[cur_gid],
                                                        ','.join(gids),
                                                        ','.join(species_gids)))

        # report selection
        type_sources, strain_ids = self._type_sources(type_metadata, [gid])
        fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
                        species, 
                        gid, 
                        ', '.join(sorted(strain_ids)),
                        type_status,
                        ', '.join("%s=%r" % (key,val) for (key,val) in type_sources.iteritems()),
                        ', '.join("%s=%r" % (key,val) for (key,val) in ncbi_types.iteritems()),
                        ', '.join("%s=%r" % (key,val) for (key,val) in ncbi_rep_categories.iteritems()),
                        quality_metadata[gid].ncbi_assembly_level if quality_metadata[gid].ncbi_assembly_level else ''))

        fout.write('\t%s\t%d\t%.2f\t%.2f\t%.2f\t%d\t%d\t%.1f\t%d\t%d\t%d' % (
                        quality_metadata[gid].ncbi_genome_category,
                        quality_metadata[gid].genome_size,
                        genome_quality[gid], 
                        quality_metadata[gid].checkm_completeness,
                        quality_metadata[gid].checkm_contamination,
                        quality_metadata[gid].scaffold_count,
                        quality_metadata[gid].contig_count,
                        quality_metadata[gid].n50_contigs,
                        quality_metadata[gid].ambiguous_bases,
                        quality_metadata[gid].ssu_count,
                        quality_metadata[gid].ssu_length if quality_metadata[gid].ssu_length else 0))
                                                                    
        if gid in ltp_top_blast_hit:
            fout.write('\t%s\t%d\t%.2f\t%.2g\t%.2g' % (
                            ltp_top_blast_hit[gid].ltp_species,
                            ltp_top_blast_hit[gid].align_len,
                            ltp_top_blast_hit[gid].perc_identity,
                            ltp_top_blast_hit[gid].bitscore,
                            ltp_top_blast_hit[gid].evalue))
        else:
            fout.write('\t%s\t%d\t%.2f\t%.2g\t%s' % ('N/A', 0, 0, 0, 'N/A'))
                                                            
        fout.write('\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n' % (len(gids),
                                                            len(species_gids),
                                                            mean_ani, mean_af,
                                                            min_ani, min_af,
                                                            excluded_from_refseq_note[gid],
                                                            note))

        return gid, require_manual_inspection

    def _ani_type_genomes(self, genome_files, type_genomes, ncbi_taxonomy):
        """Calculate ANI between type genomes."""
        
        mash = Mash(self.cpus)
        
        # create Mash sketch for potential representative genomes
        genome_list_file = os.path.join(self.output_dir, 'gtdb_type_genomes.lst')
        sketch = os.path.join(self.output_dir, 'gtdb_type_genomes.msh')
        mash.sketch(type_genomes.values(), genome_files, genome_list_file, sketch)

        # get Mash distances
        mash_dist_file = os.path.join(self.output_dir, 'gtdb_type_genomes.dst')
        mash.dist_pairwise(float(100 - self.min_mash_ani)/100, sketch, mash_dist_file)

        # read Mash distances
        mash_ani = mash.read_ani(mash_dist_file)

        # get pairs above Mash threshold
        mash_ani_pairs = []
        for qid in mash_ani:
            for rid in mash_ani[qid]:
                if mash_ani[qid][rid] >= self.min_mash_ani:
                    if qid != rid:
                        mash_ani_pairs.append((qid, rid))
                        mash_ani_pairs.append((rid, qid))
                
        self.logger.info('Identified %d genome pairs with a Mash ANI >= %.1f%%.' % (len(mash_ani_pairs), self.min_mash_ani))

        # compare genomes in the same genus
        genus_ani_pairs = []
        for rep_idA, rep_idB in combinations(type_genomes.values(), 2):
            ncbi_genusA = ncbi_taxonomy[rep_idA][5]
            ncbi_genusB = ncbi_taxonomy[rep_idB][5]
            if ncbi_genusA == ncbi_genusB:
                genus_ani_pairs.append((rep_idA, rep_idB))
                genus_ani_pairs.append((rep_idB, rep_idA))
        
        self.logger.info('Identified %d genome pairs within the same genus.' % len(genus_ani_pairs))
        
        # calculate ANI between pairs
        gid_pairs = mash_ani_pairs + genus_ani_pairs
        self.logger.info('Calculating ANI between %d genome pairs:' % len(gid_pairs))
        if True: #***
            ani_af = self.fastani.pairs(gid_pairs, genome_files)
            pickle.dump(ani_af, open(os.path.join(self.output_dir, 'type_genomes_ani_af.pkl'), 'wb'))
        else:
            ani_af = pickle.load(open(os.path.join(self.output_dir, 'type_genomes_ani_af.pkl'), 'rb'))
            
        return ani_af

    def _ani_neighbours(self, ani_af, type_genomes, ncbi_taxonomy):
        """Find all ANI neighbours."""

        # find nearest ANI neighbours
        ani_neighbours = defaultdict(lambda: set())
        fout = open(os.path.join(self.output_dir, 'gtdb_type_genome_pairwise_ani.tsv'), 'w')
        fout.write('Species 1\tType genome 1\tSpecies 2\tType genome 2\tANI\tAF\tANI12\tAF12\tANI21\tAF21\n')
        for sp1, gid1 in type_genomes.items():
            for gid2 in ani_af.get(gid1, []):
                if gid1 == gid2:
                    continue
                    
                cur_ani, cur_af = ani_af[gid1][gid2]
                rev_ani, rev_af = ani_af[gid2][gid1]
                
                # ANI should be the larger of the two values as this
                # is the most conservative circumscription and reduces the
                # change of creating polyphyletic species clusters
                ani = max(rev_ani, cur_ani)
                
                # AF should be the larger of the two values in order to 
                # accomodate incomplete and contaminated genomes
                af = max(rev_af, cur_af)
                
                # sanity check
                test_ani, test_af = symmetric_ani(ani_af, gid1, gid2)
                assert(test_ani == ani)
                assert(test_af == af)
                
                sp2 = ncbi_taxonomy[gid2][6]
                fout.write('%s\t%s\t%s\t%s' % (sp1, gid1, sp2, gid2))
                fout.write('\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                            ani, af,
                            cur_ani, cur_af, 
                            rev_ani, rev_af))

                if ani >= self.max_ani_neighbour:
                    ani_neighbours[gid1].add(gid2)
                
        fout.close()
        
        return ani_neighbours
        
    def _year_of_priority(self, metadata_file, ncbi_taxonomy):
        """Get year of priority for type strains of species."""
        
        NO_PRIORITY_YEAR = 1e6
        PriorityYear = namedtuple('PriorityYear', 'lpsn dsmz straininfo')
        
        
        gid_to_species = genome_species_assignments(ncbi_taxonomy)
        priority = read_gtdb_metadata(metadata_file, ['gtdb_type_designation',
                                                        'lpsn_priority_year', 
                                                        'dsmz_priority_year', 
                                                        'straininfo_priority_year'])
                                                        
        gtdb_type_species = gtdb_type_strain_of_species(priority)
        
        priority_years = {}
        sp_gids = defaultdict(lambda: [])
        reported_sp = set()
        for gid in gid_to_species:
            if not gid in gtdb_type_species:
                continue
                
            lpsn_year = priority[gid].lpsn_priority_year
            if not lpsn_year:
                lpsn_year = NO_PRIORITY_YEAR
                
            dsmz_year = priority[gid].dsmz_priority_year
            if not dsmz_year:
                dsmz_year = NO_PRIORITY_YEAR
                
            straininfo_year = priority[gid].straininfo_priority_year
            if not straininfo_year:
                straininfo_year = NO_PRIORITY_YEAR
                
            sp = gid_to_species[gid]
            if sp in priority_years:
                lpsn_year = min(lpsn_year, priority_years[sp].lpsn)
                dsmz_year = min(dsmz_year, priority_years[sp].dsmz)
                straininfo_year = min(straininfo_year, priority_years[sp].straininfo)
                
            priority_years[sp] = PriorityYear(lpsn=lpsn_year,
                                                dsmz=dsmz_year,
                                                straininfo=straininfo_year)

        year_of_priority = {}
        for sp in priority_years:
            year_of_priority[sp] = min(priority_years[sp].lpsn, priority_years[sp].dsmz)
            if year_of_priority[sp] == NO_PRIORITY_YEAR:
                year_of_priority[sp] = priority_years[sp].straininfo
                
        return year_of_priority
        
    def _resolve_close_ani_neighbours(self, 
                                        ani_neighbours,
                                        gtdb_type_genus,
                                        gtdb_type_sp, 
                                        gtdb_type_subsp,
                                        ncbi_type_sp,
                                        ncbi_proxy,
                                        ncbi_type_subsp,
                                        ncbi_reps,
                                        genome_quality,
                                        ncbi_taxonomy,
                                        metadata_file):
        """Resolve type genomes that have ANI neighbours deemed to be too close."""

        self.logger.info('Resolving %d type genomes with one or more neighbours within a %.1f%% ANI radius.' % (len(ani_neighbours), self.max_ani_neighbour))
        
        # get priority dates
        year_of_priority = self._year_of_priority(metadata_file, ncbi_taxonomy)
        
        # sanity check that unspecified years are set to a large number (i.e. low priority)
        for sp, year in year_of_priority.items():
            if year is None: 
                self.logger.error('Year of priority found to be None: %s' % sp)
                sys.exit(-1)
        
        # sanity check ANI neighbours
        for cur_gid in ani_neighbours:
            for neighbour_gid in ani_neighbours[cur_gid]:
                if cur_gid not in ani_neighbours[neighbour_gid]:
                    self.logger.info('ANI neighbours is not symmetrical.')
                    sys.exit(-1)
        
        # get species for each genome
        gid_to_species = genome_species_assignments(ncbi_taxonomy)
        
        # get genomes that are the type species of genus
        type_species_of_genus = set()
        for gids in gtdb_type_genus.values():
            type_species_of_genus.update(gids)

        # get type status of each genome
        type_status = defaultdict(lambda: [])
        gid_type_status = {}
        for gid in ani_neighbours:
            sp = gid_to_species[gid]
            
            if gid in gtdb_type_sp[sp]:
                type_status['TS'].append(gid)
                gid_type_status[gid] = 'TS'
            elif gid in ncbi_type_sp[sp]:
                type_status['NT'].append(gid)
                gid_type_status[gid] = 'NT'
            elif gid in ncbi_proxy[sp]:
                type_status['NP'].append(gid)
                gid_type_status[gid] = 'NP'
            elif gid in ncbi_reps[sp]:
                type_status['NR'].append(gid)
                gid_type_status[gid] = 'NR'
            elif gid in gtdb_type_subsp[sp] or gid in ncbi_type_subsp[sp]:
                type_status['TSS'].append(gid)
                gid_type_status[gid] = 'TSS'
            else:
                type_status['DN'].append(gid)
                gid_type_status[gid] = 'DN'
                
        self.logger.info('  TS = %d; NTS = %d; NP = %d; NR = %d; TSS = %d; DN = %d' % (len(type_status['TS']),
                                                                                len(type_status['NT']),
                                                                                len(type_status['NP']),
                                                                                len(type_status['NR']),
                                                                                len(type_status['TSS']),
                                                                                len(type_status['DN'])))

        # select type genomes to exclude, processing conflicts in 
        # order from least to most official in terms of type status
        excluded_gids = []
        for cur_type_status in ['DN', 'TSS', 'NR', 'NP', 'NT', 'TS']:
            if cur_type_status == 'TS':
                # greedily exclude genomes by sorting by inversely by type species of genus status,
                # year of priority, and inversely by genome quality
                cur_gids = []
                for gid in type_status[cur_type_status]:
                    type_genus = 1 if gid in type_species_of_genus else 0
                    cur_gids.append((gid, 
                                        type_genus,
                                        year_of_priority[gid_to_species[gid]], 
                                        genome_quality[gid]))
                sorted_gids = sorted(cur_gids, key=lambda x: (-x[1], x[2], -x[3]), reverse=True)
            else:
                # greedily exclude genomes by sorting by number of neighbours
                # and inversely by genome quality
                cur_gids = [(gid, len(ani_neighbours[gid]), genome_quality[gid]) for gid in type_status[cur_type_status]]
                sorted_gids = sorted(cur_gids, key=lambda x: (x[1], -x[2]), reverse=True)
                
            for d in sorted_gids:
                cur_gid = d[0]
                if len(ani_neighbours[cur_gid] - set(excluded_gids)) == 0:
                    # all ANI neighbours already added to exclusion list
                    # so no need to also exclude this type genome
                    continue
                    
                # add type genome to exclusion list
                excluded_gids.append(cur_gid)
                
        # reinstate excluded genomes from most to least official in terms of type status
        # (this resolves transitive cases: A is synonym of B which is a synonym of C)
        final_excluded_gids = set(excluded_gids)
        for gid in reversed(excluded_gids):
            if len(ani_neighbours[gid] - final_excluded_gids) == 0:
                # genome can be reinstated
                final_excluded_gids.remove(gid)
                
        self.logger.info('Identified %d type genomes for exclusion.', len(final_excluded_gids))
                
        # write out details about excluded genomes
        fout = open(os.path.join(self.output_dir, 'gtdb_excluded_ani_neighbours.tsv'), 'w')
        fout.write('Species\tType genome\tType status\tGenome quality')
        fout.write('\tPriority year\tNo. ANI neighbours\tNeighbour species\tNeighbour priority years\tNeighbour accessions\n')

        for gid in final_excluded_gids:
            sp = gid_to_species[gid]
            fout.write('%s\t%s\t%s\t%.1f\t%s\t%d\t%s\t%s\t%s\n' % (
                        sp, 
                        gid, 
                        gid_type_status[gid], 
                        genome_quality[gid],
                        str(year_of_priority.get(sp, 'N/A')),
                        len(ani_neighbours[gid]),
                        ', '.join([gid_to_species[gid] for gid in ani_neighbours[gid]]),
                        ', '.join([str(year_of_priority.get(gid_to_species[gid], 'N/A')) for gid in ani_neighbours[gid]]),
                        ', '.join([gid for gid in ani_neighbours[gid]])))

        fout.close()

        # sanity check results
        for gid in ani_neighbours:
            if gid in final_excluded_gids:
                if len(ani_neighbours[gid] - final_excluded_gids) == 0:
                    self.logger.info('Excluded type genomes %s has no remaining ANI neighbours.' % gid)
                    sys.exit(-1)
            else:
                if len(ani_neighbours[gid] - final_excluded_gids) != 0:
                    self.logger.info('Type genomes %s still has ANI neighbours.' % gid)
                    sys.exit(-1)
        
        return final_excluded_gids
        
    def write_final_type_genomes(self, initial_type_genomes_file, excluded_gids):
        """Write out final set of selected type genomes."""
        
        out_file = os.path.join(self.output_dir, 'gtdb_type_genomes_final.tsv')
        self.logger.info('Writing final type genomes to: %s' % out_file)
        fout = open(out_file, 'w')
        fout = open(out_file, 'w')
        
        with open(initial_type_genomes_file) as f:
            header = f.readline()
            fout.write(header)
            
            headers = header.strip().split('\t')
            type_genome_index = headers.index('Type genome')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gid = line_split[type_genome_index]
                
                if gid not in excluded_gids:
                    fout.write(line)
                    
        fout.close()

    def write_synonym_table(self,
                            ani_af,
                            ani_neighbours,
                            excluded_gids,
                            gtdb_type_sp,
                            ncbi_type_sp,
                            type_metadata,
                            ncbi_taxonomy,
                            metadata_file):
        """Create table indicating species names that should be considered synonyms based on ANI."""
        
        year_of_priority = self._year_of_priority(metadata_file, ncbi_taxonomy)
        
        gid_to_species = genome_species_assignments(ncbi_taxonomy)
        
        out_file = os.path.join(self.output_dir, 'synonyms.tsv')
        self.logger.info('Writing synonyms to: %s' % out_file)
        fout = open(out_file, 'w')
        fout.write('NCBI species\tType genome\tStrain IDs\tType sources\tPriority year\tNCBI assembly type')
        fout.write('\tSynonym\tSynonym type genome\tSynonym strain IDs\tSynonym type sources\tPriority year\tSynonym NCBI assembly type')
        fout.write('\tANI\tAF\n')
        
        # find closest neighbour for each excluded genome ID and 
        # report this as a synonym 
        for ex_gid in excluded_gids:
            closest_ani = 0
            closest_gid = None
            for n_gid in ani_neighbours[ex_gid]:
                if n_gid in excluded_gids:
                    continue
                    
                ani, af = symmetric_ani(ani_af, ex_gid, n_gid)
                if ani > closest_ani:
                    closest_ani = ani
                    closest_gid = n_gid

            closest_sp = gid_to_species[closest_gid]
            ex_sp = gid_to_species[ex_gid]

            ani, af = symmetric_ani(ani_af, ex_gid, closest_gid)

            fout.write('%s\t%s\t%s\t%s\t%s\t%s' % (
                        closest_sp,
                        closest_gid,
                        ','.join(sorted(type_metadata[closest_gid].ncbi_strain_identifiers)),
                        ','.join(sorted(type_metadata[closest_gid].gtdb_type_designation_sources)).upper().replace('STRAININFO', 'StrainInfo'),
                        str(year_of_priority.get(closest_sp, 'n/a')),
                        type_metadata[closest_gid].ncbi_type_material_designation))
            fout.write('\t%s\t%s\t%s\t%s\t%s\t%s' % (
                        ex_sp,
                        ex_gid,
                        ','.join(sorted(type_metadata[ex_gid].ncbi_strain_identifiers)),
                        ','.join(sorted(type_metadata[ex_gid].gtdb_type_designation_sources)).upper().replace('STRAININFO', 'StrainInfo'),
                        str(year_of_priority.get(ex_sp, 'n/a')),
                        type_metadata[ex_gid].ncbi_type_material_designation))
            fout.write('\t%.2f\t%.2f\n' % (ani, af))
            
    def verify_domain_assignments(self, gtdb_domain_report, gtdb_taxonomy):
        """Parse percentage of marker genes for each genome."""
        
        marker_perc = {}
        with open(gtdb_domain_report, encoding='utf-8') as f:
            header = f.readline().rstrip().split('\t')
            
            domain_index = header.index('Predicted domain')
            bac_marker_perc_index = header.index('Bacterial Marker Percentage')
            ar_marker_perc_index = header.index('Archaeal Marker Percentage')
            ncbi_taxonomy_index = header.index('NCBI taxonomy')
            gtdb_taxonomy_index = header.index('GTDB taxonomy')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gid = canonical_gid(line_split[0])
                domain = line_split[domain_index]
                bac_perc = float(line_split[bac_marker_perc_index])
                ar_perc = float(line_split[ar_marker_perc_index])
                ncbi_domain = [t.strip() for t in line_split[ncbi_taxonomy_index].split(';')][0]
                gtdb_domain = [t.strip() for t in line_split[gtdb_taxonomy_index].split(';')][0]
                assert(gtdb_domain == gtdb_taxonomy[gid][0])

                marker_perc[gid] = max(bac_perc, ar_perc)
                if not gid.startswith('U'):
                    if marker_perc[gid] > 10:
                        if ncbi_domain != gtdb_domain and ncbi_domain != 'None':
                            print(f'[WARNING] NCBI and GTDB domains disagree in domain report: {gid}')
                            
                            if ncbi_domain != domain and domain != 'None':
                                print(f' ... NCBI domain {ncbi_domain} also disagrees with predicted domain {domain}.')
                            
                    if marker_perc[gid] > 25 and abs(bac_perc - ar_perc) > 5:
                        if domain != gtdb_domain and domain != 'None':
                            print(f'[WARNING] GTDB and predicted domain (Bac = {bac_perc:.1f}%; Ar = {ar_perc:.1f}%) disagree in domain report: {gid}')

        return marker_perc
 
    def run(self, qc_file,
                metadata_file,
                ltp_blast_file,
                genome_path_file,
                prev_rep_file,
                ncbi_refseq_assembly_file,
                ncbi_genbank_assembly_file,
                gtdb_domain_report,
                species_exception_file,
                gtdb_type_genome_file):
        """Select GTDB type genomes for named species."""
        
        # identify genomes failing quality criteria
        self.logger.info('Reading QC file.')
        passed_qc = read_qc_file(qc_file)
        self.logger.info('Identified %d genomes passing QC.' % len(passed_qc))

        # get GTDB and NCBI taxonomy strings for each genome
        self.logger.info('Reading NCBI taxonomy from GTDB metadata file.')
        ncbi_taxonomy, ncbi_update_count = read_gtdb_ncbi_taxonomy(metadata_file, species_exception_file)
        gtdb_taxonomy = read_gtdb_taxonomy(metadata_file)
        self.logger.info('Read NCBI taxonomy for %d genomes with %d manually defined updates.' % (len(ncbi_taxonomy), ncbi_update_count))
        self.logger.info('Read GTDB taxonomy for %d genomes.' % len(gtdb_taxonomy))
        
        # verify domain assignments
        self.logger.info('Verifying domain assignment for each genome.')
        self.verify_domain_assignments(gtdb_domain_report, gtdb_taxonomy)
        
        # read manually selected type genomes
        manual_type_genomes = {}
        with open(gtdb_type_genome_file) as f:
            header = f.readline().strip().split('\t')
            
            sp_index = header.index('NCBI species')
            gid_index = header.index('Accession')
            manual_selection_index = header.index('Manual selection')
            
            for line in f:
                line_split = line.strip().split('\t')

                manual_selection = line_split[manual_selection_index].lower()
                if manual_selection.startswith('true') or manual_selection.startswith('ok'):
                    sp = line_split[sp_index].replace('Candidatus ', '')
                    gid = line_split[gid_index]
                    
                    ncbi_sp = [t.strip() for t in ncbi_taxonomy[gid]][6]
                    if ncbi_sp != sp:
                        self.logger.warning('NCBI species %s does not match species %s in manual type genome file: %s' % (ncbi_sp, sp, gid))
                    manual_type_genomes[ncbi_sp] = gid
        self.logger.info('Identified %d species with manually selected type genomes.' % len(manual_type_genomes))

        # get path to genome FASTA files
        self.logger.info('Reading path to genome FASTA files.')
        genome_files = read_genome_path(genome_path_file)
        self.logger.info('Read path for %d genomes.' % len(genome_files))
        
        # parse NCBI assembly files
        self.logger.info('Parsing NCBI assembly files.')
        excluded_from_refseq_note = exclude_from_refseq(ncbi_refseq_assembly_file, ncbi_genbank_assembly_file)
        
        # calculate quality score for genomes
        self.logger.info('Calculate quality score for all genomes.')
        genome_quality, quality_metadata = self._genome_quality(metadata_file)
        self.logger.info('Read genome quality for %d genomes.' % len(genome_quality))
        
        # parse LTP blast results to identify potential type genomes based on their 16S rRNA
        fout = open(os.path.join(self.output_dir, 'type_material_stats.tsv'), 'w')
        fout.write('Type\tNo. taxa\tNo. genomes\n')
        
        self.logger.info('Reading Living Tree Project (LTP) BLAST results to identify putative type genomes.')
        (ltp_type_species_of_genus, 
            ltp_type_strain_of_species, 
            ltp_type_strain_of_subspecies,
            ltp_top_blast_hit) = self._parse_ltp_blast_table(ltp_blast_file, ncbi_taxonomy)
            
        self.logger.info('Identified %d LTP types species of genus spanning %d genomes.' % (len(ltp_type_species_of_genus), 
                                                                                        sum([len(gids) for gids in ltp_type_species_of_genus.values()])))
        self.logger.info('Identified %d LTP types strains of species spanning %d genomes.' % (len(ltp_type_strain_of_species),
                                                                                        sum([len(gids) for gids in ltp_type_strain_of_species.values()])))
        self.logger.info('Identified %d LTP types strains of subspecies spanning %d genomes.' % (len(ltp_type_strain_of_subspecies),
                                                                                        sum([len(gids) for gids in ltp_type_strain_of_subspecies.values()])))
                                                                                        
        fout.write('%s\t%d\t%d\n' % ('LTP type species of genus', len(ltp_type_species_of_genus), sum([len(gids) for gids in ltp_type_species_of_genus.values()])))
        fout.write('%s\t%d\t%d\n' % ('LTP type strain of species', len(ltp_type_strain_of_species), sum([len(gids) for gids in ltp_type_strain_of_species.values()])))
        fout.write('%s\t%d\t%d\n' % ('LTP type strain of subspecies', len(ltp_type_strain_of_subspecies), sum([len(gids) for gids in ltp_type_strain_of_subspecies.values()])))

        # get type material designations for each genome
        self.logger.info('Reading type material designations for genomes from GTDB metadata file.')
        type_metadata = self._type_metadata(metadata_file)
                                                                
        d = self._get_type_designations(ncbi_taxonomy, type_metadata, passed_qc)
        gtdb_type_genus, gtdb_type_sp, gtdb_type_subsp, ncbi_type_sp, ncbi_proxy, ncbi_type_subsp, ncbi_reps = d

        self.logger.info('Identified %d species spanning %d genomes designated as type species of genus by GTDB.' % (
                            len(gtdb_type_genus),
                            sum([len(gids) for gids in gtdb_type_genus.values()])))
        self.logger.info('Identified %d species spanning %d genomes designated as type strain of species by GTDB.' % (
                            len(gtdb_type_sp),
                            sum([len(gids) for gids in gtdb_type_sp.values()])))
        self.logger.info('Identified %d species spanning %d genomes designated as type strain of subspecies or synonym by GTDB.' % (
                            len(gtdb_type_subsp),
                            sum([len(gids) for gids in gtdb_type_subsp.values()])))
        self.logger.info('Identified %d species spanning %d genomes designated as assembled from type material by NCBI.' % (
                            len(ncbi_type_sp),
                            sum([len(gids) for gids in ncbi_type_sp.values()])))
        self.logger.info('Identified %d species spanning %d genomes designated as assembled from proxytype material by NCBI.' % (
                            len(ncbi_proxy),
                            sum([len(gids) for gids in ncbi_proxy.values()])))
        self.logger.info('Identified %d species spanning %d genomes designated as assembly from synonym type material by NCBI.' % (
                            len(ncbi_type_subsp),
                            sum([len(gids) for gids in ncbi_type_subsp.values()])))
        self.logger.info('Identified %d species spanning %d genomes designated as a reference or representative by NCBI.' % (
                            len(ncbi_reps),
                            sum([len(gids) for gids in ncbi_reps.values()])))
        
        fout.write('%s\t%d\t%d\n' % ('GTDB type species of genus', len(gtdb_type_genus), sum([len(gids) for gids in gtdb_type_genus.values()])))
        fout.write('%s\t%d\t%d\n' % ('GTDB type strain of species', len(gtdb_type_sp), sum([len(gids) for gids in gtdb_type_sp.values()])))
        fout.write('%s\t%d\t%d\n' % ('GTDB type strain of subspecies', len(gtdb_type_subsp), sum([len(gids) for gids in gtdb_type_subsp.values()])))
        fout.write('%s\t%d\t%d\n' % ('NCBI assembled from type material', len(ncbi_type_sp), sum([len(gids) for gids in ncbi_type_sp.values()])))
        fout.write('%s\t%d\t%d\n' % ('NCBI assembled from proxytype material', len(ncbi_proxy), sum([len(gids) for gids in ncbi_proxy.values()])))
        fout.write('%s\t%d\t%d\n' % ('NCBI assembly from synonym type material', len(ncbi_type_subsp), sum([len(gids) for gids in ncbi_type_subsp.values()])))
        fout.write('%s\t%d\t%d\n' % ('NCBI reference or representative', len(ncbi_reps), sum([len(gids) for gids in ncbi_reps.values()])))
        fout.close()

        self._validate_type_designations(type_metadata,
                                            ncbi_taxonomy, 
                                            gtdb_type_sp, 
                                            gtdb_type_subsp,
                                            ncbi_type_sp,
                                            ncbi_proxy,
                                            ncbi_type_subsp,
                                            ncbi_reps)
        
        if True: #***
            type_genomes = self._select_type_genomes(passed_qc,
                                                        genome_files,
                                                        genome_quality,
                                                        quality_metadata,
                                                        type_metadata,
                                                        gtdb_type_sp, 
                                                        gtdb_type_subsp,
                                                        ncbi_type_sp,
                                                        ncbi_proxy,
                                                        ncbi_type_subsp,
                                                        ncbi_reps,
                                                        ltp_top_blast_hit,
                                                        ncbi_taxonomy,
                                                        gtdb_taxonomy,
                                                        excluded_from_refseq_note,
                                                        manual_type_genomes)
                                                        
            pickle.dump(type_genomes, open(os.path.join(self.output_dir, 'type_genomes.pkl'), 'wb'))
        else:
            type_genomes = pickle.load(open(os.path.join(self.output_dir, 'type_genomes.pkl'), 'rb'))
            
        # calculate ANI between type genomes and resolve cases where type genomes have close ANI neighbours
        ani_af = self._ani_type_genomes(genome_files, type_genomes, ncbi_taxonomy)
        ani_neighbours = self._ani_neighbours(ani_af, type_genomes, ncbi_taxonomy)
        excluded_gids = self._resolve_close_ani_neighbours(ani_neighbours,
                                                            gtdb_type_genus,
                                                            gtdb_type_sp, 
                                                            gtdb_type_subsp,
                                                            ncbi_type_sp,
                                                            ncbi_proxy,
                                                            ncbi_type_subsp,
                                                            ncbi_reps,
                                                            genome_quality,
                                                            ncbi_taxonomy,
                                                            metadata_file)
                                                            
        self.write_final_type_genomes(os.path.join(self.output_dir, 'gtdb_type_genomes_initial.tsv'), excluded_gids)
        
        self.write_synonym_table(ani_af,
                                    ani_neighbours, 
                                    excluded_gids, 
                                    gtdb_type_sp, 
                                    ncbi_type_sp,
                                    type_metadata,
                                    ncbi_taxonomy,
                                    metadata_file)
