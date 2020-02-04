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

from gtdb_species_clusters.mash import Mash
from gtdb_species_clusters.fastani import FastANI

from gtdb_species_clusters.genome import Genome
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_clusters import SpeciesClusters
                                    
from gtdb_species_clusters.genome_utils import (canonical_gid,
                                                read_qc_file,
                                                exclude_from_refseq)

from gtdb_species_clusters.type_genome_utils import symmetric_ani
                                    


class UpdateSelectRepresentatives(object):
    """Select GTDB representatives for named species."""

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
        
    def _select_highest_quality(self, gids, cur_genomes):
        """Select highest quality genome."""
                    
        q = {k:cur_genomes[k].score_type_strain() for k in gids}
        q_sorted = sorted(q.items(), key=operator.itemgetter(1), reverse=True)
        return q_sorted[0][0]
        
    def _select_ani_neighbours(self, species, gids, cur_genomes, ani_af):
        """Select highest-quality genome with sufficient number of ANI neighbours."""
        
        # calculate mean ANI of all genome pairs
        anis = []
        for gid1, gid2 in combinations(gids, 2):
            ani, af = symmetric_ani(ani_af, gid1, gid2)
            if ani > 0:
                anis.append(ani)
                
        if not anis:
            self.logger.warning('Could not calculate ANI between {:,} genomes in {}.'.format(len(gids), species))
            self.logger.warning('Selecting highest-quality genome.')
            return self._select_highest_quality(gids, cur_genomes)
                    
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

        return self._select_highest_quality(neighbour_gids, cur_genomes)

    def _type_sources(self, cur_genomes, gids):
        """Count number of genomes indicated as type material from each source."""
        
        type_sources = defaultdict(int)
        strain_ids = set()
        for gid in gids:
            for source in cur_genomes[gid].gtdb_type_sources():
                type_sources[source] += 1

            if cur_genomes[gid].strain_ids():
                strain_ids.update(cur_genomes[gid].strain_ids())
                
        return type_sources, strain_ids
        
    def _validate_type_designations(self,
                                    cur_genomes, 
                                    gtdb_type_sp, 
                                    gtdb_type_subsp,
                                    ncbi_type_sp,
                                    ncbi_proxy,
                                    ncbi_type_subsp,
                                    ncbi_reps):
        """Produce files indicating potentially missing type material."""
        
        # determine all binomial NCBI species names
        ncbi_species = cur_genomes.named_ncbi_species()
        
        # identify species without a type strain
        fout = open(os.path.join(self.output_dir, 'validate_type_strain_of_species.tsv'), 'w')
        fout.write('NCBI species\tType designation\tNo. genomes with designation\tGTDB type sources\tNCBI type designations\tAccession(s)\n')
        missing_type_strain = 0
        missing_type_strain_ncbi_type = 0
        not_ncbi_type_sp = 0
        for sp in ncbi_species:
            ncbi_types = defaultdict(int)
            for gid in ncbi_type_sp.get(sp, []):
                ncbi_types[cur_genomes[gid].ncbi_type_material] += 1

            type_sources, strain_ids = self._type_sources(cur_genomes, gtdb_type_sp.get(sp, []))
            
            if sp in gtdb_type_sp:
                fout.write('%s\t%s\t%d\t%s\t%s\t%s\n' % (sp,
                                                        'type strain of species',
                                                        len(gtdb_type_sp[sp]),
                                                        ', '.join("%s=%r" % (key,val) for (key,val) in type_sources.items()),
                                                        ', '.join("%s=%r" % (key,val) for (key,val) in ncbi_types.items()),
                                                        ', '.join([gid for d in gtdb_type_sp[sp]])))
                if sp not in ncbi_type_sp:
                    not_ncbi_type_sp += 1
            elif sp in ncbi_type_sp:
                fout.write('%s\t%s\t%d\t%s\t%s\t%s\n' % (sp,
                                                        'assembly from type material',
                                                        len(ncbi_type_sp[sp]),
                                                        ', '.join("%s=%r" % (key,val) for (key,val) in type_sources.items()),
                                                        ', '.join("%s=%r" % (key,val) for (key,val) in ncbi_types.items()),
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
        
        self.logger.info('Identified {:,} NCBI species.'.format(len(ncbi_species)))
        self.logger.info('Identified {:,} species without a type strain of species, type strain of subspecies, or representative genome.'.format(missing_type_strain))
        self.logger.info('Identified {:,} species without a GTDB designated type strain of species that have genome(s) designated as assembled from type material at NCBI.'.format(missing_type_strain_ncbi_type))
        self.logger.info('Identified {:,} species with a GTDB designated type strain of species that are not designated as assembled from type material at NCBI.'.format(not_ncbi_type_sp))
        
    def _select_rep_genomes(self,
                                cur_genomes,
                                existing_sp_reps,
                                gtdb_type_sp, 
                                gtdb_type_subsp,
                                ncbi_type_sp,
                                ncbi_proxy,
                                ncbi_type_subsp,
                                ncbi_reps,
                                excluded_from_refseq_note):
        """Select representative genome for each species."""

        # determine all NCBI species names
        ncbi_species = cur_genomes.named_ncbi_species()
        
        # identify species without a type strain
        self.logger.info('Selecting representative for each of the {:,} named NCBI species.'.format(len(ncbi_species)))
        
        init_rep_out_file = os.path.join(self.output_dir, 'gtdb_rep_genomes_initial.tsv')
        fout = open(init_rep_out_file, 'w')
        fout.write('NCBI species\tRepresentative\tStrain IDs\tType status\tType sources\tNCBI assembly types\tNCBI representative\tNCBI assembly level')
        fout.write('\tNCBI genome category\tGenome size (bp)\tQuality score\tCompleteness (%)\tContamination (%)\tNo. scaffolds\tNo. contigs\tN50 contigs\tAmbiguous bases\tSSU count\tSSU length (bp)')
        fout.write('\tNo. type genomes\tNo. species genomes\tMean ANI\tMean AF\tMin ANI\tMin AF\tNCBI exclude from RefSeq\tSelection note\n')
        
        fout_manual = open(os.path.join(self.output_dir, 'gtdb_rep_genomes.manual.tsv'), 'w')
        fout_manual.write('NCBI species\tAccession\tGTDB taxonomy\tNCBI taxonomy\tGTDB type genome\tType status\tANI to rep\tAF to rep\tSelection note')
        fout_manual.write('\tNCBI genome category\tGenome size (bp)\tQuality score\tCompleteness (%)\tContamination (%)\tNo. scaffolds\tNo. contigs\tN50 contigs\tAmbiguous bases\tSSU count\tSSU length (bp)')
        fout_manual.write('\tNo. type genomes\tNo. species genomes\tMean ANI\tMean AF\tMin ANI\tMin AF\tNCBI exclude from RefSeq\tType accessions\tSpecies accessions\n')

        num_existing_reps = 0
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
        rep_genomes = {}
        multi_gids = 0
        for idx, ncbi_sp in enumerate(ncbi_species):
            species_gids = ncbi_species[ncbi_sp]

            statusStr = '-> Processing {:,} of {:,} ({:.2f}%) species [{}: {:,}].'.format(idx+1, 
                                                                                len(ncbi_species), 
                                                                                float(idx+1)*100/len(ncbi_species),
                                                                                ncbi_sp,
                                                                                len(species_gids)).ljust(86)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            existing_reps = set(existing_sp_reps.keys()).intersection(species_gids)
            if len(existing_reps):
                num_existing_reps += 1
                continue

            if ncbi_sp in gtdb_type_sp:
                gid, manual_inspection = self._select_rep(cur_genomes,
                                                ncbi_sp,
                                                'type strain of species',
                                                gtdb_type_sp[ncbi_sp],
                                                ncbi_type_sp,
                                                ncbi_reps,
                                                excluded_from_refseq_note,
                                                species_gids,
                                                fout,
                                                fout_manual)
                                                
                if len(gtdb_type_sp[ncbi_sp]) > 1:
                    multi_gids += 1
                    
                if manual_inspection:
                    num_type_strain_of_species_manual += 1

                num_type_strain_of_species += 1
            elif ncbi_sp in ncbi_type_sp:
                gid, manual_inspection = self._select_rep(cur_genomes,
                                                ncbi_sp,
                                                'NCBI assembled from type material',
                                                ncbi_type_sp[ncbi_sp],
                                                ncbi_type_sp,
                                                ncbi_reps,
                                                excluded_from_refseq_note,
                                                species_gids,
                                                fout,
                                                fout_manual)
                                                
                if len(ncbi_type_sp[ncbi_sp]) > 1:
                    multi_gids += 1
                    
                if manual_inspection:
                    num_ncbi_assembled_from_type_manual += 1
                
                num_ncbi_assembled_from_type += 1
            elif ncbi_sp in ncbi_proxy:
                gid, manual_inspection = self._select_rep(cur_genomes,
                                                ncbi_sp,
                                                'NCBI assembled from proxytype material',
                                                ncbi_proxy[ncbi_sp],
                                                ncbi_type_sp,
                                                ncbi_reps,
                                                excluded_from_refseq_note,
                                                species_gids,
                                                fout,
                                                fout_manual)
                                                
                if len(ncbi_proxy[ncbi_sp]) > 1:
                    multi_gids += 1
                    
                if manual_inspection:
                    num_ncbi_assembled_from_proxytype_manual += 1
                
                num_ncbi_assembled_from_proxytype += 1
            elif ncbi_sp in ncbi_reps:
                gid, manual_inspection = self._select_rep(cur_genomes,
                                                ncbi_sp,
                                                'NCBI representative genome',
                                                ncbi_reps[ncbi_sp],
                                                ncbi_type_sp,
                                                ncbi_reps,
                                                excluded_from_refseq_note,
                                                species_gids,
                                                fout,
                                                fout_manual)
                
                if len(ncbi_reps[ncbi_sp]) > 1:
                    multi_gids += 1
                    
                if manual_inspection:
                    num_ncbi_rep_manual += 1
                    
                num_ncbi_rep += 1
            elif ncbi_sp in gtdb_type_subsp or ncbi_sp in ncbi_type_subsp:
                gids = gtdb_type_subsp[ncbi_sp].union(ncbi_type_subsp[ncbi_sp])
                gid, manual_inspection = self._select_rep(cur_genomes,
                                                ncbi_sp,
                                                'type strain of subspecies',
                                                gids,
                                                ncbi_type_sp,
                                                ncbi_reps,
                                                excluded_from_refseq_note,
                                                species_gids,
                                                fout,
                                                fout_manual)
                                                
                if len(gids) > 1:
                    multi_gids += 1
                    
                if manual_inspection:
                    num_type_strain_of_subspecies_manual += 1
                
                num_type_strain_of_subspecies += 1
            else:
                gid, manual_inspection = self._select_rep(cur_genomes,
                                                ncbi_sp,
                                                'no type genome',
                                                species_gids,
                                                ncbi_type_sp,
                                                ncbi_reps,
                                                excluded_from_refseq_note,
                                                species_gids,
                                                fout,
                                                fout_manual)
                if len(species_gids) > 1:
                    multi_gids += 1
                    
                if manual_inspection:
                    num_de_novo_manual += 1
                    
                num_de_novo += 1
            
            if gid in rep_genomes:
                self.logger.error('Representative genome selected for multiple species: {} {} {}'.format(gid, ncbi_sp, rep_genomes[gid]))
                sys.exit(-1)
                
            rep_genomes[gid] = ncbi_sp

        sys.stdout.write('\n')
        fout.close()
        fout_manual.close()
        
        self.logger.info('GTDB representative is existing representative: {:,} ({:.1f}%)'.format(num_existing_reps, 
                                                                                        num_existing_reps*100.0/len(ncbi_species)))
        self.logger.info('GTDB representative is type strain of species: {:,} ({:.1f}%)'.format(num_type_strain_of_species, 
                                                                                        num_type_strain_of_species*100.0/len(ncbi_species)))
        self.logger.info('GTDB representative is assembled from type material according to NCBI: {:,} ({:.1f}%)'.format(num_ncbi_assembled_from_type, 
                                                                                        num_ncbi_assembled_from_type*100.0/len(ncbi_species)))
        self.logger.info('GTDB representative is assembled from proxytype material according to NCBI: {:,} ({:.1f}%)'.format(num_ncbi_assembled_from_proxytype, 
                                                                                        num_ncbi_assembled_from_proxytype*100.0/len(ncbi_species)))
        self.logger.info('GTDB representative is a representative genome at NCBI: {:,} ({:.1f}%)'.format(num_ncbi_rep, 
                                                                                        num_ncbi_rep*100.0/len(ncbi_species)))
        self.logger.info('GTDB representative is type strain of subspecies: {:,} ({:.1f}%)'.format(num_type_strain_of_subspecies, 
                                                                                        num_type_strain_of_subspecies*100.0/len(ncbi_species)))
        self.logger.info('Species with de novo selected representative: {:,} ({:.1f}%)'.format(num_de_novo, num_de_novo*100.0/len(ncbi_species)))

        self.logger.info('Identified {:,} species where multiple potential representatives exist.'.format(multi_gids))
        self.logger.info('Identified species requiring manual inspection of selected representative:')
        self.logger.info('  TS = {:,}; NTS = {:,}; NP = {:,}; NR = {:,}; TSS = {:,}; DN = {:,}'.format(
                            num_type_strain_of_species_manual,
                            num_ncbi_assembled_from_type_manual,
                            num_ncbi_assembled_from_proxytype_manual,
                            num_ncbi_rep_manual,
                            num_type_strain_of_subspecies_manual,
                            num_de_novo_manual))
                                                                                
        return rep_genomes

    def _select_rep(self,
                    cur_genomes,
                    ncbi_sp,
                    type_status,
                    gids,
                    ncbi_type_sp,
                    ncbi_reps,
                    excluded_from_refseq_note,
                    species_gids,
                    fout,
                    fout_manual):
        """Select type genome."""
        
        ncbi_types = defaultdict(int)
        ncbi_type_strain_ids = set()
        for gid in ncbi_type_sp.get(ncbi_sp, []):
            ncbi_types[cur_genomes[gid].ncbi_type_material] += 1
            ncbi_type_strain_ids.update(cur_genomes[gid].strain_ids())
                
        ncbi_rep_categories = defaultdict(int)
        ncbi_rep_strain_ids = set()
        for gid in ncbi_reps.get(ncbi_sp, []):
            ncbi_rep_categories[cur_genomes[gid].ncbi_refseq_category] += 1
            ncbi_rep_strain_ids.update(cur_genomes[gid].strain_ids())

        # select representative genome
        mean_ani = mean_af = min_ani = min_af = 'n/a'
        note = ''
        require_manual_inspection = False
        if len(gids) == 1:
            rep_gid = next(iter(gids))
            note = 'select single genome'
        else:
            # check if single 'assembled from type material' genome selected by NCBI
            ncbi_type_gids = gids.intersection(ncbi_type_sp[ncbi_sp])
            if len(ncbi_type_gids) == 1:
                rep_gid = ncbi_type_gids.pop()
                note = 'select single genome annotated as assembled from type material at NCBI'
            else:
                # calculate ANI between genomes
                ani_af = self.fastani.pairwise(gids, cur_genomes.genomic_files)
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

                rep_gid = self._select_highest_quality(gids, cur_genomes)
                note = 'selected highest-quality genome'
                
                if float(min_ani) < self.min_intra_strain_ani:
                    # check if NCBI has designated a reference or representative genome
                    ncbi_rep_gids = gids.intersection(ncbi_reps[ncbi_sp])
                    if len(ncbi_rep_gids) == 1:
                        rep_gid = ncbi_rep_gids.pop()
                        note = 'selected single NCBI representative genome'
                    else:
                        require_manual_inspection = True
                        rep_gid = self._select_ani_neighbours(ncbi_sp, gids, cur_genomes, ani_af)
                        note = 'selected highest-quality genome with sufficient ANI neighbours'
                        
                        for cur_gid in [rep_gid] + list(gids.difference([rep_gid])): # write results with selected genome first
                            fout_manual.write('%s\t%s\t%s\t%s\t%s\t%s' % (
                                                ncbi_sp, 
                                                cur_gid,
                                                cur_genomes[cur_gid].gtdb_taxa,
                                                cur_genomes[cur_gid].ncbi_taxa,
                                                cur_gid==rep_gid, 
                                                type_status))
                            if cur_gid != rep_gid:
                                if cur_gid in ani_af and rep_gid in ani_af[cur_gid]:
                                    cur_ani, cur_af = ani_af[cur_gid][rep_gid]
                                else:
                                    cur_ani, cur_af = 0.0, 0.0
                                fout_manual.write('\t%.1f\t%.2f' % (cur_ani, cur_af))
                            else:
                                fout_manual.write('\t%.1f\t%.2f' % (100.0, 1.0))
                            fout_manual.write('\t%s\t%s\t%d\t%.2f\t%.2f\t%.2f\t%d\t%d\t%.1f\t%d\t%d\t%d' % (
                                                    note,
                                                    cur_genomes[cur_gid].ncbi_genome_category,
                                                    cur_genomes[cur_gid].length,
                                                    cur_genomes[cur_gid].score_assembly(), 
                                                    cur_genomes[cur_gid].comp,
                                                    cur_genomes[cur_gid].cont,
                                                    cur_genomes[cur_gid].scaffold_count,
                                                    cur_genomes[cur_gid].contig_count,
                                                    cur_genomes[cur_gid].contig_n50,
                                                    cur_genomes[cur_gid].ambiguous_bases,
                                                    cur_genomes[cur_gid].ssu_count,
                                                    cur_genomes[cur_gid].ssu_length))

                            fout_manual.write('\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                                                    len(gids),
                                                    len(species_gids),
                                                    mean_ani, mean_af,
                                                    min_ani, min_af,
                                                    excluded_from_refseq_note[cur_gid],
                                                    ','.join(gids),
                                                    ','.join(species_gids)))

        # report selection
        type_sources, strain_ids = self._type_sources(cur_genomes, [rep_gid])
        fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
                        ncbi_sp, 
                        rep_gid, 
                        ', '.join(sorted(strain_ids)),
                        type_status,
                        ', '.join("%s=%r" % (key,val) for (key,val) in type_sources.items()),
                        ', '.join("%s=%r" % (key,val) for (key,val) in ncbi_types.items()),
                        ', '.join("%s=%r" % (key,val) for (key,val) in ncbi_rep_categories.items()),
                        cur_genomes[rep_gid].ncbi_assembly_level if cur_genomes[rep_gid].ncbi_assembly_level else ''))

        fout.write('\t%s\t%d\t%.2f\t%.2f\t%.2f\t%d\t%d\t%.1f\t%d\t%d\t%d' % (
                        cur_genomes[rep_gid].ncbi_genome_category,
                        cur_genomes[rep_gid].length,
                        cur_genomes[rep_gid].score_assembly(), 
                        cur_genomes[rep_gid].comp,
                        cur_genomes[rep_gid].cont,
                        cur_genomes[rep_gid].scaffold_count,
                        cur_genomes[rep_gid].contig_count,
                        cur_genomes[rep_gid].contig_n50,
                        cur_genomes[rep_gid].ambiguous_bases,
                        cur_genomes[rep_gid].ssu_count,
                        cur_genomes[rep_gid].ssu_length))
                                                            
        fout.write('\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n' % (len(gids),
                                                            len(species_gids),
                                                            mean_ani, mean_af,
                                                            min_ani, min_af,
                                                            excluded_from_refseq_note[gid],
                                                            note))

        return gid, require_manual_inspection

    def _ani_reps(self, cur_genomes, all_rep_genomes):
        """Calculate ANI between representative genomes."""
        
        if False: #***
            mash = Mash(self.cpus)
            
            # sanity check
            for gid, sp in all_rep_genomes.items():
                if gid not in cur_genomes.genomic_files:
                    self.logger.error(f'Missing genomic file for {gid} from {sp}.')
                    sys.exit(-1)
            
            # create Mash sketch for potential representative genomes
            genome_list_file = os.path.join(self.output_dir, 'gtdb_reps.lst')
            sketch = os.path.join(self.output_dir, 'gtdb_reps.msh')
            mash.sketch(all_rep_genomes, cur_genomes.genomic_files, genome_list_file, sketch)

            # get Mash distances
            mash_dist_file = os.path.join(self.output_dir, 'gtdb_reps.dst')
            mash.dist_pairwise(float(100 - self.min_mash_ani)/100, sketch, mash_dist_file)

            # read Mash distances
            mash_ani = mash.read_ani(mash_dist_file)

            # get pairs above Mash threshold
            mash_ani_pairs = []
            for qid in mash_ani:
                for rid in mash_ani[qid]:
                    if mash_ani[qid][rid] >= self.min_mash_ani:
                        if qid != rid:
                            q = cur_genomes.user_uba_id_map.get(qid, qid)
                            r = cur_genomes.user_uba_id_map.get(rid, rid)
                            mash_ani_pairs.append((q, r))
                            mash_ani_pairs.append((r, q))
                    
            self.logger.info('Identified {:,} genome pairs with a Mash ANI >= {:.1f}%.'.format(len(mash_ani_pairs), self.min_mash_ani))

            # compare genomes in the same genus
            genus_ani_pairs = []
            for rep_idA, rep_idB in combinations(all_rep_genomes, 2):
                ncbi_genusA = cur_genomes[rep_idA].ncbi_taxa.genus
                ncbi_genusB = cur_genomes[rep_idB].ncbi_taxa.genus
                gtdb_genusA = cur_genomes[rep_idA].gtdb_taxa.genus
                gtdb_genusB = cur_genomes[rep_idB].gtdb_taxa.genus
                if (ncbi_genusA != 'g__' and ncbi_genusA == ncbi_genusB) or (gtdb_genusA != 'g__' and gtdb_genusA == gtdb_genusB):
                    genus_ani_pairs.append((rep_idA, rep_idB))
                    genus_ani_pairs.append((rep_idB, rep_idA))
            
            self.logger.info('Identified {:,} genome pairs within the same genus.'.format(len(genus_ani_pairs)))
            
            # calculate ANI between pairs
            gid_pairs = set(mash_ani_pairs + genus_ani_pairs)
            self.logger.info('Calculating ANI between {:,} genome pairs:'.format(len(gid_pairs)))
            ani_af = self.fastani.pairs(gid_pairs, cur_genomes.genomic_files)
            pickle.dump(ani_af, open(os.path.join(self.output_dir, 'reps_ani_af.pkl'), 'wb'))
        else:
            ani_af = pickle.load(open(os.path.join(self.output_dir, 'reps_ani_af.pkl'), 'rb'))
            
        return ani_af

    def _ani_neighbours(self, cur_genomes, ani_af, all_rep_genomes):
        """Find all ANI neighbours."""

        # find nearest ANI neighbours
        ani_neighbours = defaultdict(lambda: set())
        fout = open(os.path.join(self.output_dir, 'gtdb_rep_pairwise_ani.tsv'), 'w')
        fout.write('NCBI species 1\tRepresentative 1\tNCBI species 2\tRepresentative 2\tANI\tAF\tANI12\tAF12\tANI21\tAF21\n')
        for gid1, ncbi_sp1 in all_rep_genomes.items():
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
                # accommodate incomplete and contaminated genomes
                af = max(rev_af, cur_af)
                
                ncbi_sp2 = cur_genomes[gid2].ncbi_taxa.species
                fout.write('%s\t%s\t%s\t%s' % (ncbi_sp1, gid1, ncbi_sp2, gid2))
                fout.write('\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                            ani, af,
                            cur_ani, cur_af, 
                            rev_ani, rev_af))

                if ani >= self.max_ani_neighbour:
                    ani_neighbours[gid1].add(gid2)
                
        fout.close()
        
        # sanity check ANI neighbours
        for cur_gid in ani_neighbours:
            for neighbour_gid in ani_neighbours[cur_gid]:
                if cur_gid not in ani_neighbours[neighbour_gid]:
                    self.logger.info('ANI neighbours is not symmetrical: {} {}'.format(ani_af[cur_gid][neighbour_gid], ani_af[neighbour_gid][cur_gid]))
                    print(cur_gid, neighbour_gid)
                    print('cur_gid in all_rep_genomes', cur_gid in all_rep_genomes)
                    print('neighbour_gid in all_rep_genomes', neighbour_gid in all_rep_genomes)
                    sys.exit(-1)
        
        return ani_neighbours
        
    def _year_of_priority(self, cur_genomes):
        """Get year of priority for type strains of species."""

        PriorityYear = namedtuple('PriorityYear', 'lpsn dsmz straininfo')

        priority_years = {}
        for gid in cur_genomes:
            if not cur_genomes[gid].is_gtdb_type_strain:
                continue
                
            lpsn_year = cur_genomes[gid].lpsn_priority_year
            dsmz_year = cur_genomes[gid].dsmz_priority_year
            straininfo_year = cur_genomes[gid].straininfo_priority_year

            ncbi_sp = cur_genomes[gid].ncbi_taxa.species
            if ncbi_sp in priority_years:
                lpsn_year = min(lpsn_year, priority_years[ncbi_sp].lpsn)
                dsmz_year = min(dsmz_year, priority_years[ncbi_sp].dsmz)
                straininfo_year = min(straininfo_year, priority_years[ncbi_sp].straininfo)
                
            priority_years[ncbi_sp] = PriorityYear(lpsn=lpsn_year,
                                                dsmz=dsmz_year,
                                                straininfo=straininfo_year)

        year_of_priority = {}
        for sp in priority_years:
            year_of_priority[sp] = min(priority_years[sp].lpsn, priority_years[sp].dsmz)
            if year_of_priority[sp] == Genome.NO_PRIORITY_YEAR:
                year_of_priority[sp] = priority_years[sp].straininfo
                
        return year_of_priority
        
    def _resolve_close_ani_neighbours(self, 
                                        cur_genomes,
                                        ani_neighbours,
                                        gtdb_type_genus,
                                        gtdb_type_sp, 
                                        gtdb_type_subsp,
                                        ncbi_type_sp,
                                        ncbi_proxy,
                                        ncbi_type_subsp,
                                        ncbi_reps):
        """Resolve representatives that have ANI neighbours deemed to be too close."""

        self.logger.info('Resolving {:,} representatives with one or more neighbours within a {:.1f}% ANI radius.'.format(len(ani_neighbours), self.max_ani_neighbour))
        
        # get priority dates
        year_of_priority = self._year_of_priority(cur_genomes)
        
        # sanity check that unspecified years are set to a large number (i.e. low priority)
        for sp, year in year_of_priority.items():
            if year is None: 
                self.logger.error('Year of priority found to be None: %s' % sp)
                sys.exit(-1)

        # get genomes that are the type species of genus
        type_species_of_genus = set()
        for gids in gtdb_type_genus.values():
            type_species_of_genus.update(gids)

        # get type status of each genome
        type_status = defaultdict(lambda: [])
        gid_type_status = {}
        for gid in ani_neighbours:
            sp = cur_genomes[gid].ncbi_taxa.species
            
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

        # select representatives to exclude, processing conflicts in 
        # order from least to most official in terms of type status
        excluded_gids = []
        for cur_type_status in ['DN', 'TSS', 'NR', 'NP', 'NT', 'TS']:
            if cur_type_status == 'TS':
                # greedily exclude representatives by sorting by inversely by type species of genus status,
                # year of priority, and inversely by genome quality
                cur_gids = []
                for gid in type_status[cur_type_status]:
                    type_genus = 1 if gid in type_species_of_genus else 0
                    cur_gids.append((gid, 
                                        type_genus,
                                        year_of_priority[cur_genomes[gid].ncbi_taxa.species], 
                                        cur_genomes[gid].score_assembly()))
                sorted_gids = sorted(cur_gids, key=lambda x: (-x[1], x[2], -x[3]), reverse=True)
            else:
                # greedily exclude representatives by sorting by number of neighbours
                # and inversely by genome quality
                cur_gids = [(gid, len(ani_neighbours[gid]), cur_genomes[gid].score_assembly()) for gid in type_status[cur_type_status]]
                sorted_gids = sorted(cur_gids, key=lambda x: (x[1], -x[2]), reverse=True)
                
            for d in sorted_gids:
                cur_gid = d[0]
                if len(ani_neighbours[cur_gid] - set(excluded_gids)) == 0:
                    # all ANI neighbours already added to exclusion list
                    # so no need to also exclude this representative
                    continue
                    
                # add representative to exclusion list
                excluded_gids.append(cur_gid)
                
        # reinstate excluded genomes from most to least official in terms of type status
        # (this resolves transitive cases: A is synonym of B which is a synonym of C)
        final_excluded_gids = set(excluded_gids)
        for gid in reversed(excluded_gids):
            if len(ani_neighbours[gid] - final_excluded_gids) == 0:
                # genome can be reinstated
                final_excluded_gids.remove(gid)
                
        self.logger.info('Identified {:,} representatives for exclusion.'.format(len(final_excluded_gids)))
                
        # write out details about excluded genomes
        fout = open(os.path.join(self.output_dir, 'gtdb_excluded_ani_neighbours.tsv'), 'w')
        fout.write('Species\tRepresentative\tRepresentative status\tGenome quality')
        fout.write('\tPriority year\tNo. ANI neighbours\tNeighbour species\tNeighbour priority years\tNeighbour accessions\n')

        for gid in final_excluded_gids:
            sp = cur_genomes[gid].ncbi_taxa.species
            fout.write('%s\t%s\t%s\t%.1f\t%s\t%d\t%s\t%s\t%s\n' % (
                        sp, 
                        gid, 
                        gid_type_status[gid], 
                        cur_genomes[gid].score_assembly(),
                        str(year_of_priority.get(sp, 'N/A')),
                        len(ani_neighbours[gid]),
                        ', '.join([cur_genomes[gid].ncbi_taxa.species for gid in ani_neighbours[gid]]),
                        ', '.join([str(year_of_priority.get(cur_genomes[gid].ncbi_taxa.species, 'N/A')) for gid in ani_neighbours[gid]]),
                        ', '.join([gid for gid in ani_neighbours[gid]])))

        fout.close()

        # sanity check results
        for gid in ani_neighbours:
            if gid in final_excluded_gids:
                if len(ani_neighbours[gid] - final_excluded_gids) == 0:
                    self.logger.info('Excluded representative {} has no remaining ANI neighbours.'.format(gid))
                    sys.exit(-1)
            else:
                if len(ani_neighbours[gid] - final_excluded_gids) != 0:
                    self.logger.info('Representative genome {} still has ANI neighbours.'.format(gid))
                    sys.exit(-1)
        
        return final_excluded_gids
        
    def write_final_reps(self, cur_genomes, all_rep_genomes, excluded_gids, excluded_from_refseq_note):
        """Write out final set of selected representatives."""
        
        out_file = os.path.join(self.output_dir, 'gtdb_reps_final.tsv')
        self.logger.info('Writing final representatives to: {}'.format(out_file))
        fout = open(out_file, 'w')
        fout.write('Representative\tProposed species\tNCBI species\tGTDB species\tStrain IDs\tGTDB type designation\tType sources\tNCBI type strain\tNCBI representative\tNCBI assembly level')
        fout.write('\tNCBI genome category\tGenome size (bp)\tQuality score\tCompleteness (%)\tContamination (%)\tNo. scaffolds\tNo. contigs\tN50 contigs\tAmbiguous bases\tSSU count\tSSU length (bp)')
        fout.write('\tNCBI exclude from RefSeq\n')
        
        num_reps = 0
        for rid, proposed_species in all_rep_genomes.items():
            if rid in excluded_gids:
                continue 
           
            num_reps += 1
            type_sources, strain_ids = self._type_sources(cur_genomes, [rid])
            fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                            rid,
                            proposed_species,
                            cur_genomes[rid].ncbi_taxa.species,
                            cur_genomes[rid].gtdb_taxa.species,
                            ', '.join(sorted(strain_ids)),
                            cur_genomes[rid].gtdb_type_designation,
                            ', '.join("%s=%r" % (key,val) for (key,val) in type_sources.items()),
                            cur_genomes[rid].is_ncbi_type_strain(),
                            cur_genomes[rid].is_ncbi_representative(),
                            cur_genomes[rid].ncbi_assembly_level))

            fout.write('\t{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{}\t{:.1f}\t{}\t{}\t{}\t{}\n'.format(
                            cur_genomes[rid].ncbi_genome_category,
                            cur_genomes[rid].length,
                            cur_genomes[rid].score_assembly(), 
                            cur_genomes[rid].comp,
                            cur_genomes[rid].cont,
                            cur_genomes[rid].scaffold_count,
                            cur_genomes[rid].contig_count,
                            cur_genomes[rid].contig_n50,
                            cur_genomes[rid].ambiguous_bases,
                            cur_genomes[rid].ssu_count,
                            cur_genomes[rid].ssu_length,
                            excluded_from_refseq_note.get(rid, None)))

        fout.close()
        
        self.logger.info('Wrote {:,} named GTDB representative to file.'.format(num_reps))

    def write_synonym_table(self,
                            cur_genomes,
                            ani_af,
                            ani_neighbours, 
                            excluded_gids, 
                            gtdb_type_sp, 
                            ncbi_type_sp):
        """Create table indicating species names that should be considered synonyms based on ANI."""
        
        year_of_priority = self._year_of_priority(cur_genomes)

        out_file = os.path.join(self.output_dir, 'synonyms.tsv')
        self.logger.info('Writing synonyms to: %s' % out_file)
        fout = open(out_file, 'w')
        fout.write('NCBI species\tRepresentative\tStrain IDs\tRepresentative type sources\tPriority year\tNCBI assembly type')
        fout.write('\tSynonym\tSynonym genome\tSynonym strain IDs\tSynonym type sources\tPriority year\tSynonym NCBI assembly type')
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

            closest_sp = cur_genomes[closest_gid].ncbi_taxa.species
            ex_sp = cur_genomes[ex_gid].ncbi_taxa.species

            ani, af = symmetric_ani(ani_af, ex_gid, closest_gid)

            fout.write('%s\t%s\t%s\t%s\t%s\t%s' % (
                        closest_sp,
                        closest_gid,
                        ','.join(sorted(cur_genomes[closest_gid].strain_ids())),
                        ','.join(sorted(cur_genomes[closest_gid].gtdb_type_sources())).upper().replace('STRAININFO', 'StrainInfo'),
                        str(year_of_priority.get(closest_sp, 'n/a')),
                        cur_genomes[closest_gid].ncbi_type_material))
            fout.write('\t%s\t%s\t%s\t%s\t%s\t%s' % (
                        ex_sp,
                        ex_gid,
                        ','.join(sorted(cur_genomes[ex_gid].strain_ids())),
                        ','.join(sorted(cur_genomes[ex_gid].gtdb_type_sources())).upper().replace('STRAININFO', 'StrainInfo'),
                        str(year_of_priority.get(ex_sp, 'n/a')),
                        cur_genomes[ex_gid].ncbi_type_material))
            fout.write('\t%.2f\t%.2f\n' % (ani, af))
            
    def read_existing_sp_reps(self, existing_sp_rep_file):
        """Parse representative for existing species clusters."""
        
        existing_sp_reps = {}
        lost_prev_sp = set()
        updated_prev_sp = set()
        with open(existing_sp_rep_file) as f:
            header = f.readline().strip().split('\t')
            
            prev_sp_index = header.index('Previous species name')
            cur_rep_id_index = header.index('New representative ID')
            cur_sp_index = header.index('New species name')
            rep_status_index = header.index('Representative status')
            
            for line in f:
                tokens = line.strip().split('\t')
                
                prev_sp = tokens[prev_sp_index]
                
                rep_status = tokens[rep_status_index]
                if rep_status != 'LOST':
                    cur_sp = tokens[cur_sp_index]
                    cur_rep_id = tokens[cur_rep_id_index]
                    assert cur_rep_id != 'None'
                    
                    existing_sp_reps[cur_rep_id] = cur_sp
                    
                    if prev_sp != cur_sp:
                        updated_prev_sp.add(prev_sp)
                else:
                    lost_prev_sp.add(prev_sp)
        
        return existing_sp_reps, lost_prev_sp, updated_prev_sp
 
    def run(self, existing_sp_reps,
                    qc_passed_file,
                    cur_gtdb_metadata_file,
                    cur_genomic_path_file,
                    uba_genome_paths,
                    ncbi_genbank_assembly_file,
                    gtdb_domain_report,
                    species_exception_file,
                    genus_exception_file,
                    gtdb_type_strains_ledger):
        """Select GTDB type genomes for named species."""
        
        # read species representative for existing species clusters
        self.logger.info('Reading existing species representatives.')
        existing_sp_reps, lost_prev_sp, updated_prev_sp = self.read_existing_sp_reps(existing_sp_reps)
        self.logger.info(f' ... identified representatives for {len(existing_sp_reps):,} existing species clusters.')
        self.logger.info(f' ... identified {len(lost_prev_sp):,} lost species names.')
        self.logger.info(f' ... identified {len(updated_prev_sp):,} updated species names.')

        # create current GTDB genome sets
        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                                species_exception_file,
                                                genus_exception_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                create_sp_clusters=False,
                                                uba_genome_file=uba_genome_paths,
                                                qc_passed_file=qc_passed_file)
        self.logger.info(f' ... current genome set contains {len(cur_genomes):,} genomes.')
        
        # get path to previous and current genomic FASTA files
        self.logger.info('Reading path to current genomic FASTA files.')
        cur_genomes.load_genomic_file_paths(cur_genomic_path_file)
        cur_genomes.load_genomic_file_paths(uba_genome_paths)

        # report type material
        self.logger.info('Tabulating genomes assigned as type material.')
        gtdb_type_genus = defaultdict(set)
        gtdb_type_sp = defaultdict(set)
        gtdb_type_subsp = defaultdict(set)
        ncbi_type_sp = defaultdict(set)
        ncbi_proxy = defaultdict(set)
        ncbi_type_subsp = defaultdict(set)
        ncbi_reps = defaultdict(set)
        for gid in cur_genomes:
            ncbi_sp = cur_genomes[gid].ncbi_taxa.species
            if cur_genomes[gid].gtdb_type_species_of_genus:
                gtdb_type_genus[ncbi_sp].add(gid)
                
            if cur_genomes[gid].is_gtdb_type_strain():
                gtdb_type_sp[ncbi_sp].add(gid)
            elif cur_genomes[gid].is_gtdb_type_subspecies():
                gtdb_type_subsp[ncbi_sp].add(gid)
                
            if cur_genomes[gid].is_ncbi_type_strain():
                ncbi_type_sp[ncbi_sp].add(gid)
            elif cur_genomes[gid].is_ncbi_proxy():
                ncbi_proxy[ncbi_sp].add(gid)
            elif cur_genomes[gid].is_ncbi_type_subspecies():
                ncbi_type_subsp[ncbi_sp].add(gid)
                
            if cur_genomes[gid].is_ncbi_representative():
                ncbi_reps[ncbi_sp].add(gid)
                
        # clean up NCBI representative genomes so only the representative of 
        # a species and not subspecies is retained if multiple representatives exist
        for sp, gids in ncbi_reps.items():
            sp_gids = set()
            for gid in gids:
                if not cur_genomes[gid].is_ncbi_subspecies():
                    sp_gids.add(gid)
                    
            if len(sp_gids) > 0:
                ncbi_reps[sp] = sp_gids
                
        self.logger.info('Identified {:,} species spanning {:,} genomes designated as type species of genus by GTDB.'.format(
                            len(gtdb_type_genus),
                            sum([len(gids) for gids in gtdb_type_genus.values()])))
        self.logger.info('Identified {:,} species spanning {:,} genomes designated as type strain of species by GTDB.'.format(
                            len(gtdb_type_sp),
                            sum([len(gids) for gids in gtdb_type_sp.values()])))
        self.logger.info('Identified {:,} species spanning {:,} genomes designated as type strain of subspecies or synonym by GTDB.'.format(
                            len(gtdb_type_subsp),
                            sum([len(gids) for gids in gtdb_type_subsp.values()])))
        self.logger.info('Identified {:,} species spanning {:,} genomes designated as assembled from type material by NCBI.'.format(
                            len(ncbi_type_sp),
                            sum([len(gids) for gids in ncbi_type_sp.values()])))
        self.logger.info('Identified {:,} species spanning {:,} genomes designated as assembled from proxytype material by NCBI.'.format(
                            len(ncbi_proxy),
                            sum([len(gids) for gids in ncbi_proxy.values()])))
        self.logger.info('Identified {:,} species spanning {:,} genomes designated as assembly from synonym type material by NCBI.'.format(
                            len(ncbi_type_subsp),
                            sum([len(gids) for gids in ncbi_type_subsp.values()])))
        self.logger.info('Identified {:,} species spanning {:,} genomes designated as a reference or representative by NCBI.'.format(
                            len(ncbi_reps),
                            sum([len(gids) for gids in ncbi_reps.values()])))
        
        fout = open(os.path.join(self.output_dir, 'type_material_stats.tsv'), 'w')
        fout.write('Type\tNo. taxa\tNo. genomes\n')
        fout.write('%s\t%d\t%d\n' % ('GTDB type species of genus', len(gtdb_type_genus), sum([len(gids) for gids in gtdb_type_genus.values()])))
        fout.write('%s\t%d\t%d\n' % ('GTDB type strain of species', len(gtdb_type_sp), sum([len(gids) for gids in gtdb_type_sp.values()])))
        fout.write('%s\t%d\t%d\n' % ('GTDB type strain of subspecies', len(gtdb_type_subsp), sum([len(gids) for gids in gtdb_type_subsp.values()])))
        fout.write('%s\t%d\t%d\n' % ('NCBI assembled from type material', len(ncbi_type_sp), sum([len(gids) for gids in ncbi_type_sp.values()])))
        fout.write('%s\t%d\t%d\n' % ('NCBI assembled from proxytype material', len(ncbi_proxy), sum([len(gids) for gids in ncbi_proxy.values()])))
        fout.write('%s\t%d\t%d\n' % ('NCBI assembly from synonym type material', len(ncbi_type_subsp), sum([len(gids) for gids in ncbi_type_subsp.values()])))
        fout.write('%s\t%d\t%d\n' % ('NCBI reference or representative', len(ncbi_reps), sum([len(gids) for gids in ncbi_reps.values()])))
        fout.close()
        
        # identify species without type material
        self._validate_type_designations(cur_genomes, 
                                            gtdb_type_sp, 
                                            gtdb_type_subsp,
                                            ncbi_type_sp,
                                            ncbi_proxy,
                                            ncbi_type_subsp,
                                            ncbi_reps)
        
        # parse NCBI assembly files
        self.logger.info('Parsing NCBI assembly files.')
        excluded_from_refseq_note = exclude_from_refseq(ncbi_genbank_assembly_file)
        
        self.logger.info('Selecting representative for new named NCBI species.')
        if False: #***
            new_rep_genomes = self._select_rep_genomes(cur_genomes,
                                                        existing_sp_reps,
                                                        gtdb_type_sp, 
                                                        gtdb_type_subsp,
                                                        ncbi_type_sp,
                                                        ncbi_proxy,
                                                        ncbi_type_subsp,
                                                        ncbi_reps,
                                                        excluded_from_refseq_note)
                                                        
            pickle.dump(new_rep_genomes, open(os.path.join(self.output_dir, 'new_rep_genomes.pkl'), 'wb'))
        else:
            new_rep_genomes = pickle.load(open(os.path.join(self.output_dir, 'new_rep_genomes.pkl'), 'rb'))
            
        self.logger.info(f'Identified representatives for {len(new_rep_genomes):,} new NCBI named species.') 
            
        # calculate ANI between representatives and resolve cases where type genomes have close ANI neighbours
        dup_sp = set(existing_sp_reps.values()).intersection(new_rep_genomes.values())
        for cur_sp in dup_sp:
            print(cur_sp, [gid for gid, sp in existing_sp_reps.items() if sp == cur_sp], [gid for gid, sp in new_rep_genomes.items() if sp == cur_sp])
            
        all_rep_genomes = {}
        for gid, sp in existing_sp_reps.items():
            all_rep_genomes[gid] = sp
            
        for gid, sp in new_rep_genomes.items():
            all_rep_genomes[gid] = sp
            
        self.logger.info(f'Calculate ANI between {len(all_rep_genomes):,} representatives of GTDB species clusters.')

        ani_af = self._ani_reps(cur_genomes, all_rep_genomes)

        self.logger.info('Establishing ANI neighbours.')
        ani_neighbours = self._ani_neighbours(cur_genomes, ani_af, all_rep_genomes)
        
        self.logger.info('Resolving ANI neighbours.')
        excluded_gids = self._resolve_close_ani_neighbours(cur_genomes,
                                                            ani_neighbours,
                                                            gtdb_type_genus,
                                                            gtdb_type_sp, 
                                                            gtdb_type_subsp,
                                                            ncbi_type_sp,
                                                            ncbi_proxy,
                                                            ncbi_type_subsp,
                                                            ncbi_reps)
                                                            
        self.write_final_reps(cur_genomes,
                                all_rep_genomes, 
                                excluded_gids,
                                excluded_from_refseq_note)
        
        self.write_synonym_table(cur_genomes,
                                    ani_af,
                                    ani_neighbours, 
                                    excluded_gids, 
                                    gtdb_type_sp, 
                                    ncbi_type_sp)
