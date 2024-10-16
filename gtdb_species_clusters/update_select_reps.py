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
import pickle
from itertools import combinations, product
from collections import defaultdict

from numpy import (mean as np_mean,
                   std as np_std)

from gtdblib.util.shell.execute import check_dependencies

from gtdb_species_clusters.skani import Skani
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_clusters import SpeciesClusters
from gtdb_species_clusters.genome_utils import select_highest_quality
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters import defaults as Defaults
from gtdb_species_clusters.taxon_utils import (specific_epithet, is_placeholder_sp_epithet)


class UpdateSelectRepresentatives():
    """Select GTDB representatives for named species."""

    def __init__(self, cpus, output_dir):
        """Initialization."""

        # List of validly published NCBI species that are
        # considered synonyms under the GTDB and are ignored
        # in the sense that representatives are not selected
        # for these species. A more formal way should be
        # considered for handling this (e.g., a ledger
        # indicating GTDB synonyms), but since this is
        # a special case for Shigella at the moment it
        # is being explicitly handled here.
        self.invalid_ncbi_sp = set(['s__Shigella flexneri',
                                    's__Shigella sonnei',
                                    's__Shigella boydii',
                                    's__Shigella dysenteriae'])

        check_dependencies(['skani'])

        self.cpus = cpus
        self.output_dir = output_dir

        self.log = logging.getLogger('rich')

        self.min_intra_strain_ani = 99.0
        self.max_ani_neighbour = Defaults.ANI_SYNONYMS
        self.max_af_neighbour = Defaults.AF_SP

        self.skani = Skani(None, cpus)

    def select_ani_neighbours(self, species, gids, cur_genomes, ani_af):
        """Select highest-quality genome with sufficient number of ANI neighbours."""

        # calculate mean ANI of all genome pairs
        anis = []
        for gid1, gid2 in combinations(gids, 2):
            ani, _af = Skani.symmetric_ani_af(ani_af, gid1, gid2)
            if ani > 0:
                anis.append(ani)

        if not anis:
            self.log.warning(
                'Could not calculate ANI between {:,} genomes in {}.'.format(len(gids), species))

            hq_gid = select_highest_quality(gids, cur_genomes)
            self.log.warning(f'Selecting highest-quality genome: {hq_gid}')

            return hq_gid

        # calculate number of ANI neighbours for each genome
        mean_ani = np_mean(anis)
        std_ani = np_std(anis)
        ani_neighbours = defaultdict(int)
        for gid1, gid2 in product(gids, repeat=2):
            if gid1 == gid2:
                ani_neighbours[gid1] += 1
            else:
                ani, _af = Skani.symmetric_ani_af(ani_af, gid1, gid2)
                if ani >= mean_ani - std_ani:
                    ani_neighbours[gid1] += 1

        # get all genomes that are neighbours with at least half the other genomes,
        # or if no genomes are neighbours with at least half the genomes has the most
        # number of neighbours
        majority_neighbour_gids = []
        most_neighbours_gids = []
        most_neighbours = max(ani_neighbours.values())
        for gid in ani_neighbours:
            if ani_neighbours[gid] >= 0.5*len(gids):
                majority_neighbour_gids.append(gid)
            if ani_neighbours[gid] >= most_neighbours:
                most_neighbours_gids.append(gid)

        if len(majority_neighbour_gids) == 0:
            return select_highest_quality(most_neighbours_gids, cur_genomes)

        return select_highest_quality(majority_neighbour_gids, cur_genomes)

    def type_sources(self, cur_genomes, gids):
        """Count number of genomes indicated as type material from each source."""

        type_sources = defaultdict(int)
        strain_ids = set()
        for gid in gids:
            for source in cur_genomes[gid].gtdb_type_sources():
                type_sources[source] += 1

            if cur_genomes[gid].strain_ids():
                strain_ids.update(cur_genomes[gid].strain_ids())

        return type_sources, strain_ids

    def validate_type_designations(self,
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
        fout = open(os.path.join(self.output_dir,
                                 'validate_type_strain_of_species.tsv'), 'w')
        fout.write(
            'NCBI species\tType designation\tNo. genomes with designation\tGTDB type sources\tNCBI type designations\tAccession(s)\n')
        missing_type_strain = 0
        missing_type_strain_ncbi_type = 0
        not_ncbi_type_sp = 0
        for sp in ncbi_species:
            ncbi_types = defaultdict(int)
            for gid in ncbi_type_sp.get(sp, []):
                ncbi_types[cur_genomes[gid].ncbi_type_material] += 1

            type_sources, _strain_ids = self.type_sources(
                cur_genomes, gtdb_type_sp.get(sp, []))

            if sp in gtdb_type_sp:
                fout.write('%s\t%s\t%d\t%s\t%s\t%s\n' % (sp,
                                                         'type strain of species',
                                                         len(gtdb_type_sp[sp]),
                                                         ', '.join("%s=%r" % (key, val) for (
                                                             key, val) in type_sources.items()),
                                                         ', '.join("%s=%r" % (key, val) for (
                                                             key, val) in ncbi_types.items()),
                                                         ', '.join(gtdb_type_sp[sp])))
                if sp not in ncbi_type_sp:
                    not_ncbi_type_sp += 1
            elif sp in ncbi_type_sp:
                fout.write('%s\t%s\t%d\t%s\t%s\t%s\n' % (sp,
                                                         'assembly from type material',
                                                         len(ncbi_type_sp[sp]),
                                                         ', '.join("%s=%r" % (key, val) for (
                                                             key, val) in type_sources.items()),
                                                         ', '.join("%s=%r" % (key, val) for (
                                                             key, val) in ncbi_types.items()),
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

        self.log.info(
            'Identified {:,} NCBI species.'.format(len(ncbi_species)))
        self.log.info('Identified {:,} species without a type strain of species, type strain of subspecies, or representative genome.'.format(
            missing_type_strain))
        self.log.info('Identified {:,} species without a GTDB designated type strain of species that have genome(s) designated as assembled from type material at NCBI.'.format(
            missing_type_strain_ncbi_type))
        self.log.info(
            'Identified {:,} species with a GTDB designated type strain of species that are not designated as assembled from type material at NCBI.'.format(not_ncbi_type_sp))

    def select_rep_genomes(self,
                           unrepresented_ncbi_sp,
                           cur_genomes,
                           gtdb_type_sp,
                           gtdb_type_subsp,
                           ncbi_type_sp,
                           ncbi_proxy,
                           ncbi_type_subsp,
                           ncbi_reps):
        """Select representative genome for each species."""

        # identify species without a type strain
        init_rep_out_file = os.path.join(
            self.output_dir, 'gtdb_rep_genomes_initial.tsv')
        fout = open(init_rep_out_file, 'w')
        fout.write('NCBI species\tRepresentative\tStrain IDs\tType status\tType sources\tNCBI assembly types\tNCBI representative\tNCBI assembly level')
        fout.write('\tNCBI genome category\tGenome size (bp)\tQuality score\tCompleteness (%)\tContamination (%)\tNo. scaffolds\tNo. contigs\tN50 contigs\tAmbiguous bases\tSSU count\tSSU length (bp)')
        fout.write(
            '\tNo. type genomes\tNo. species genomes\tMean ANI\tMean AF\tMin ANI\tMin AF\tNCBI exclude from RefSeq\tSelection note\n')

        fout_manual = open(os.path.join(
            self.output_dir, 'gtdb_rep_genomes.manual.tsv'), 'w')
        fout_manual.write(
            'NCBI species\tAccession\tGTDB taxonomy\tNCBI taxonomy\tGTDB type genome\tType status\tANI to rep\tAF to rep\tSelection note')
        fout_manual.write(
            '\tNCBI genome category\tGenome size (bp)\tQuality score\tCompleteness (%)\tContamination (%)\tNo. scaffolds\tNo. contigs\tN50 contigs\tAmbiguous bases\tSSU count\tSSU length (bp)')
        fout_manual.write(
            '\tNo. type genomes\tNo. species genomes\tMean ANI\tMean AF\tMin ANI\tMin AF\tNCBI exclude from RefSeq\tType accessions\tSpecies accessions\n')

        num_type_strain_of_species = 0
        num_type_strain_of_subspecies = 0
        num_ncbi_assembled_from_type = 0
        num_ncbi_assembled_from_proxytype = 0
        num_ncbi_rep = 0
        num_de_novo = 0

        num_type_strain_of_species_manual = 0
        num_type_strain_of_subspecies_manual = 0
        num_ncbi_assembled_from_type_manual = 0
        num_ncbi_assembled_from_proxytype_manual = 0
        num_ncbi_rep_manual = 0
        num_de_novo_manual = 0
        rep_genomes = {}
        multi_gids = 0

        ncbi_species = cur_genomes.named_ncbi_species()
        for idx, ncbi_sp in enumerate(unrepresented_ncbi_sp):
            species_gids = ncbi_species[ncbi_sp]

            statusStr = '-> Processing {:,} of {:,} ({:.2f}%) species [{}: {:,}].'.format(
                idx+1,
                len(unrepresented_ncbi_sp),
                float(idx+1)*100/len(unrepresented_ncbi_sp),
                ncbi_sp,
                len(species_gids)).ljust(86)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            if ncbi_sp in gtdb_type_sp:
                gid, manual_inspection = self.select_rep(cur_genomes,
                                                         ncbi_sp,
                                                         'type strain of species',
                                                         gtdb_type_sp[ncbi_sp],
                                                         ncbi_type_sp,
                                                         ncbi_reps,
                                                         species_gids,
                                                         fout,
                                                         fout_manual)

                if len(gtdb_type_sp[ncbi_sp]) > 1:
                    multi_gids += 1

                if manual_inspection:
                    num_type_strain_of_species_manual += 1

                num_type_strain_of_species += 1
            elif ncbi_sp in ncbi_type_sp:
                gid, manual_inspection = self.select_rep(cur_genomes,
                                                         ncbi_sp,
                                                         'NCBI assembled from type material',
                                                         ncbi_type_sp[ncbi_sp],
                                                         ncbi_type_sp,
                                                         ncbi_reps,
                                                         species_gids,
                                                         fout,
                                                         fout_manual)

                if len(ncbi_type_sp[ncbi_sp]) > 1:
                    multi_gids += 1

                if manual_inspection:
                    num_ncbi_assembled_from_type_manual += 1

                num_ncbi_assembled_from_type += 1
            elif ncbi_sp in ncbi_proxy:
                gid, manual_inspection = self.select_rep(cur_genomes,
                                                         ncbi_sp,
                                                         'NCBI assembled from proxytype material',
                                                         ncbi_proxy[ncbi_sp],
                                                         ncbi_type_sp,
                                                         ncbi_reps,
                                                         species_gids,
                                                         fout,
                                                         fout_manual)

                if len(ncbi_proxy[ncbi_sp]) > 1:
                    multi_gids += 1

                if manual_inspection:
                    num_ncbi_assembled_from_proxytype_manual += 1

                num_ncbi_assembled_from_proxytype += 1
            elif ncbi_sp in ncbi_reps:
                gid, manual_inspection = self.select_rep(cur_genomes,
                                                         ncbi_sp,
                                                         'NCBI representative genome',
                                                         ncbi_reps[ncbi_sp],
                                                         ncbi_type_sp,
                                                         ncbi_reps,
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
                gid, manual_inspection = self.select_rep(cur_genomes,
                                                         ncbi_sp,
                                                         'type strain of subspecies',
                                                         gids,
                                                         ncbi_type_sp,
                                                         ncbi_reps,
                                                         species_gids,
                                                         fout,
                                                         fout_manual)

                if len(gids) > 1:
                    multi_gids += 1

                if manual_inspection:
                    num_type_strain_of_subspecies_manual += 1

                num_type_strain_of_subspecies += 1
            else:
                gid, manual_inspection = self.select_rep(cur_genomes,
                                                         ncbi_sp,
                                                         'no type genome',
                                                         species_gids,
                                                         ncbi_type_sp,
                                                         ncbi_reps,
                                                         species_gids,
                                                         fout,
                                                         fout_manual)
                if len(species_gids) > 1:
                    multi_gids += 1

                if manual_inspection:
                    num_de_novo_manual += 1

                num_de_novo += 1

            if gid in rep_genomes:
                self.log.error('Representative genome selected for multiple species: {} {} {}'.format(
                    gid, ncbi_sp, rep_genomes[gid]))
                sys.exit(-1)

            rep_genomes[gid] = ncbi_sp

        sys.stdout.write('\n')
        fout.close()
        fout_manual.close()

        self.log.info(' - GTDB representative is type strain of species: {:,} ({:.1f}%)'.format(num_type_strain_of_species,
                                                                                             num_type_strain_of_species*100.0/len(unrepresented_ncbi_sp)))
        self.log.info(' - GTDB representative is assembled from type material according to NCBI: {:,} ({:.1f}%)'.format(num_ncbi_assembled_from_type,
                                                                                                                     num_ncbi_assembled_from_type*100.0/len(unrepresented_ncbi_sp)))
        self.log.info(' - GTDB representative is assembled from proxytype material according to NCBI: {:,} ({:.1f}%)'.format(num_ncbi_assembled_from_proxytype,
                                                                                                                          num_ncbi_assembled_from_proxytype*100.0/len(unrepresented_ncbi_sp)))
        self.log.info(' - GTDB representative is a representative genome at NCBI: {:,} ({:.1f}%)'.format(num_ncbi_rep,
                                                                                                      num_ncbi_rep*100.0/len(unrepresented_ncbi_sp)))
        self.log.info(' - GTDB representative is type strain of subspecies: {:,} ({:.1f}%)'.format(num_type_strain_of_subspecies,
                                                                                                num_type_strain_of_subspecies*100.0/len(unrepresented_ncbi_sp)))
        self.log.info(' - species with de novo selected representative: {:,} ({:.1f}%)'.format(
            num_de_novo, num_de_novo*100.0/len(unrepresented_ncbi_sp)))

        self.log.info(
            ' - identified {:,} species where multiple potential representatives exist'.format(multi_gids))
        self.log.info(
            ' - identified species requiring manual inspection of selected representative:')
        self.log.info('  TS = {:,}; NTS = {:,}; NP = {:,}; NR = {:,}; TSS = {:,}; DN = {:,}'.format(
            num_type_strain_of_species_manual,
            num_ncbi_assembled_from_type_manual,
            num_ncbi_assembled_from_proxytype_manual,
            num_ncbi_rep_manual,
            num_type_strain_of_subspecies_manual,
            num_de_novo_manual))

        return rep_genomes

    def select_rep(self,
                   cur_genomes,
                   ncbi_sp,
                   type_status,
                   gids,
                   ncbi_type_sp,
                   ncbi_reps,
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
                ani_af = self.skani.pairwise(gids, 
                                             cur_genomes.genomic_files, 
                                             preset = Defaults.SKANI_PRESET,
                                             report_progress=False)
                anis = []
                afs = []
                for q in ani_af:
                    anis += [d[0] for d in ani_af[q].values()]
                    afs += [d[1] for d in ani_af[q].values()]
                    afs += [d[2] for d in ani_af[q].values()]

                if anis:
                    mean_ani = f'{np_mean(anis):.1f}'
                    mean_af = f'{np_mean(afs):.2f}'
                    min_ani = f'{min(anis):.1f}'
                    min_af = f'{min(afs):.2f}'
                else:
                    mean_ani = mean_af = min_ani = min_af = 'n/a'

                rep_gid = select_highest_quality(gids, cur_genomes)
                note = 'selected highest-quality genome'

                if float(min_ani) < self.min_intra_strain_ani:
                    # check if NCBI has designated a reference or representative genome
                    ncbi_rep_gids = gids.intersection(ncbi_reps[ncbi_sp])
                    if len(ncbi_rep_gids) == 1:
                        rep_gid = ncbi_rep_gids.pop()
                        note = 'selected single NCBI representative genome'
                    else:
                        require_manual_inspection = True
                        rep_gid = self.select_ani_neighbours(
                            ncbi_sp, gids, cur_genomes, ani_af)
                        note = 'selected highest-quality genome with sufficient ANI neighbours'

                        # write results with selected genome first
                        for cur_gid in [rep_gid] + list(gids.difference([rep_gid])):
                            fout_manual.write('%s\t%s\t%s\t%s\t%s\t%s' % (
                                ncbi_sp,
                                cur_gid,
                                cur_genomes[cur_gid].gtdb_taxa,
                                cur_genomes[cur_gid].ncbi_taxa,
                                cur_gid == rep_gid,
                                type_status))
                            if cur_gid != rep_gid:
                                cur_ani, cur_af = Skani.symmetric_ani_af(ani_af, cur_gid, rep_gid)
                                fout_manual.write(
                                    '\t%.1f\t%.2f' % (cur_ani, cur_af))
                            else:
                                fout_manual.write(
                                    '\t%.1f\t%.2f' % (100.0, 1.0))
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
                                cur_genomes[cur_gid].excluded_from_refseq_note,
                                ','.join(gids),
                                ','.join(species_gids)))

        # report selection
        type_sources, strain_ids = self.type_sources(cur_genomes, [rep_gid])
        fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
            ncbi_sp,
            rep_gid,
            ', '.join(sorted(strain_ids)),
            type_status,
            ', '.join("%s=%r" % (key, val)
                      for (key, val) in type_sources.items()),
            ', '.join("%s=%r" % (key, val)
                      for (key, val) in ncbi_types.items()),
            ', '.join("%s=%r" % (key, val)
                      for (key, val) in ncbi_rep_categories.items()),
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
                                                           cur_genomes[rep_gid].excluded_from_refseq_note,
                                                           note))

        return rep_gid, require_manual_inspection

    def ani_reps(self, cur_genomes, all_rep_genomes):
        """Calculate ANI between representative genomes."""

        if True:  # ***DEBUGGING
            # create dictionary indicating path to species representative genomes FASTA files
            rep_genomic_files = {}
            for gid, sp in all_rep_genomes.items():
                if gid not in cur_genomes.genomic_files:
                    self.log.error(
                        f'Missing genomic file for {gid} from {sp}.')
                    sys.exit(-1)

                rep_genomic_files[gid] = cur_genomes.genomic_files[gid]

            # calculate ANI between pairs
            ani_af = self.skani.triangle(rep_genomic_files, 
                                         self.output_dir, 
                                         preset = Defaults.SKANI_PRESET,
                                         min_af = Defaults.AF_SP,
                                         min_sketch_ani = Defaults.SKANI_PREFILTER_THRESHOLD)
            pickle.dump(ani_af, open(os.path.join(
                self.output_dir, 'reps_ani_af.pkl'), 'wb'))
        else:
            ani_af = pickle.load(
                open(os.path.join(self.output_dir, 'reps_ani_af.pkl'), 'rb'))

        return ani_af

    def ani_neighbours(self, updated_sp_clusters, cur_genomes, ani_af):
        """Find all ANI neighbours."""

        # find nearest ANI neighbours
        ani_neighbours = defaultdict(set)
        fout = open(os.path.join(self.output_dir,
                                 'gtdb_rep_pairwise_ani.tsv'), 'w')
        fout.write(
            'NCBI species 1\tRepresentative 1\tNCBI species 2\tRepresentative 2\tANI\tAF12\tAF21\n')
        
        fout_n = open(os.path.join(self.output_dir, 'ani_neighbour_table.init.tsv'), 'w')
        fout_n.write('Genome ID 1\tGTDB species\tNCBI species\tGenome ID 2\tGTDB species\tNCBI species\tANI\tAF\n')

        for gid1 in ani_af:
            for gid2 in ani_af.get(gid1, []):
                ani, af_r, af_q = Skani.result(ani_af, gid1, gid2)
                if ani == 0.0:
                    continue

                af = max(af_r, af_q)

                ncbi_sp1 = cur_genomes[gid1].ncbi_taxa.species
                ncbi_sp2 = cur_genomes[gid2].ncbi_taxa.species
                fout.write('{}\t{}\t{}\t{}'.format(
                    ncbi_sp1, gid1, ncbi_sp2, gid2))
                fout.write('\t{}\t{}\t{}\n'.format(
                    ani, 
                    af_r, 
                    af_q))

                prev_reps_factor = 0
                if gid1 in updated_sp_clusters and gid2 in updated_sp_clusters:
                    # HACK: representatives were distinct GTDB species clusters in
                    # the previous release so only consider them ANI neibhours
                    # if they exceed an additional 'fudge factor' for ANI
                    # similarity. This accounts for small differences in ANI
                    # calculated by different versions of FastANI / skani that can
                    # cause two GTDB clusters to just miss the neighbour
                    # cutoff in one release and then exceeed it in the next
                    prev_reps_factor = 0.1

                # TODO: in order to preserve Latin species names, these are allowed to 
                # have an ANI upto 97%. In contrast, placeholder species should be
                # clustered at 95%.
                if False:
                    specific1 = specific_epithet(cur_genomes[gid1].gtdb_taxa.species)
                    specific2 = specific_epithet(cur_genomes[gid2].gtdb_taxa.species)

                    if ((specific1 and is_placeholder_sp_epithet(specific1)) or 
                        (specific2 and is_placeholder_sp_epithet(specific2))):
                        ani_neighbour_cutoff = Defaults.ANI_SP
                    else:
                        ani_neighbour_cutoff = self.max_ani_neighbour
                else:
                    ani_neighbour_cutoff = self.max_ani_neighbour

                if (ani >= ani_neighbour_cutoff + prev_reps_factor
                        and af >= self.max_af_neighbour):
                    ani_neighbours[gid1].add(gid2)
                    ani_neighbours[gid2].add(gid1)

                    fout_n.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        gid1,
                        cur_genomes[gid1].gtdb_taxa.species,
                        cur_genomes[gid1].ncbi_taxa.species,
                        gid2,
                        cur_genomes[gid2].gtdb_taxa.species,
                        cur_genomes[gid2].ncbi_taxa.species,
                        ani,
                        af
                    ))

        fout.close()
        fout_n.close()

        return ani_neighbours

    def resolve_close_ani_neighbours(self,
                                     cur_genomes,
                                     ani_neighbours,
                                     gtdb_type_genus,
                                     gtdb_type_sp,
                                     gtdb_type_subsp,
                                     ncbi_type_sp,
                                     ncbi_proxy,
                                     ncbi_type_subsp,
                                     ncbi_reps,
                                     sp_priority_ledger,
                                     genus_priority_ledger,
                                     lpsn_gss_file):
        """Resolve representatives that have ANI neighbours deemed to be too close."""

        self.log.info('Resolving {:,} representatives with neighbours within a {:.1f}% ANI radius and >= {:.1f}% AF.'.format(
            len(ani_neighbours),
            self.max_ani_neighbour,
            self.max_af_neighbour))

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

        self.log.info('  TS = {}; NTS = {}; NP = {}; NR = {}; TSS = {}; DN = {}'.format(
            len(type_status['TS']),
            len(
                type_status['NT']),
            len(
                type_status['NP']),
            len(
                type_status['NR']),
            len(
                type_status['TSS']),
            len(type_status['DN'])))

        # select representatives to exclude, processing conflicts in
        # order from least to most official in terms of type status
        excluded_gids = []
        for cur_type_status in ['DN', 'TSS', 'NR', 'NP', 'NT', 'TS']:
            if cur_type_status == 'TS':
                # greedily exclude representatives in reverse order of priority
                sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger,
                                                          genus_priority_ledger,
                                                          lpsn_gss_file,
                                                          self.output_dir)
                sorted_gids = sp_priority_mngr.sort_by_sp_priority(cur_genomes,
                                                                   type_status[cur_type_status],
                                                                   reverse=True)
            else:
                # greedily exclude representatives by sorting by number of neighbours
                # and inversely by genome quality
                cur_gids = [(gid, len(ani_neighbours[gid]), cur_genomes[gid].score_assembly(
                )) for gid in type_status[cur_type_status]]
                sorted_gids = sorted(
                    cur_gids, key=lambda x: (x[1], -x[2]), reverse=True)
                sorted_gids = [d[0] for d in sorted_gids]

            for cur_gid in sorted_gids:
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

        self.log.info('Identified {:,} representatives for exclusion.'.format(
            len(final_excluded_gids)))

        # write out details about excluded genomes
        fout = open(os.path.join(self.output_dir,
                                 'gtdb_excluded_ani_neighbours.tsv'), 'w')
        fout.write(
            'Species\tRepresentative\tStatus\tGenome quality')
        fout.write(
            '\tPriority year\tNo. ANI neighbours\tNeighbour species\tNeighbour status\tNeighbour priority years\tNeighbour accessions\n')

        for gid in final_excluded_gids:
            sp = cur_genomes[gid].ncbi_taxa.species
            fout.write('{}\t{}\t{}\t{:.1f}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                sp,
                gid,
                gid_type_status[gid],
                cur_genomes[gid].score_assembly(),
                str(cur_genomes[gid].year_of_priority()),
                len(ani_neighbours[gid]),
                ', '.join(
                    [cur_genomes[gid].ncbi_taxa.species for gid in ani_neighbours[gid]]),
                ', '.join([gid_type_status[gid] for gid in ani_neighbours[gid]]),
                ', '.join([str(cur_genomes[gid].year_of_priority())
                           for gid in ani_neighbours[gid]]),
                ', '.join([gid for gid in ani_neighbours[gid]])))

        fout.close()

        # sanity check results
        for gid in ani_neighbours:
            if gid in final_excluded_gids:
                if len(ani_neighbours[gid] - final_excluded_gids) == 0:
                    self.log.info(
                        'Excluded representative {} has no remaining ANI neighbours.'.format(gid))
                    sys.exit(-1)
            else:
                if len(ani_neighbours[gid] - final_excluded_gids) != 0:
                    self.log.info(
                        'Representative genome {} still has ANI neighbours.'.format(gid))
                    sys.exit(-1)

        return final_excluded_gids

    def write_final_reps(self, cur_genomes, all_rep_genomes, excluded_gids):
        """Write out final set of selected representatives."""

        out_file = os.path.join(self.output_dir, 'gtdb_named_reps_final.tsv')
        self.log.info(
            'Writing final representatives to: {}'.format(out_file))
        fout = open(out_file, 'w')
        fout.write('Representative\tProposed species\tNCBI species\tGTDB species\tStrain IDs\tGTDB type designation\tType sources\tNCBI type strain\tNCBI representative\tNCBI assembly level')
        fout.write('\tNCBI genome category\tGenome size (bp)\tQuality score\tCompleteness (%)\tContamination (%)\tNo. scaffolds\tNo. contigs\tN50 contigs\tAmbiguous bases\tSSU count\tSSU length (bp)')
        fout.write('\tNCBI exclude from RefSeq\n')

        num_reps = 0
        for rid, proposed_species in all_rep_genomes.items():
            if rid in excluded_gids:
                continue

            num_reps += 1
            type_sources, strain_ids = self.type_sources(cur_genomes, [rid])
            fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                rid,
                proposed_species,
                cur_genomes[rid].ncbi_taxa.species,
                cur_genomes[rid].gtdb_taxa.species,
                ', '.join(sorted(strain_ids)),
                cur_genomes[rid].gtdb_type_designation,
                ', '.join("%s=%r" % (key, val)
                          for (key, val) in type_sources.items()),
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
                cur_genomes[rid].excluded_from_refseq_note))

        fout.close()

        self.log.info(
            'Wrote {:,} named GTDB representative to file.'.format(num_reps))

    def write_ani_neighbour_table(self,
                                  updated_sp_clusters,
                                  cur_genomes,
                                  ani_af,
                                  ani_neighbours,
                                  excluded_gids):
        """Create table indicating GTDB clusters that were merged as they had relatively high ANI."""

        out_file = os.path.join(self.output_dir, 'ani_neighbour_table.merged.tsv')
        self.log.info(f'Writing ANI neighbour table to: {out_file}')
        fout = open(out_file, 'w')
        fout.write('NCBI species\tGTDB species\tRepresentative\tStrain IDs\tRepresentative type sources\tPriority year\tGTDB type species\tGTDB type strain\tNCBI assembly type')
        fout.write('\tNCBI synonym\tGTDB synonym\tSynonym genome\tSynonym strain IDs\tSynonym type sources\tPriority year\tSynonym GTDB type species\tSynonym GTDB type strain\tSynonym NCBI assembly type')
        fout.write('\tANI\tAF\tAction\n')

        # find closest neighbour for each excluded genome ID and
        # report this as a synonym
        for ex_gid in excluded_gids:
            closest_ani = 0
            closest_gid = None
            for n_gid in ani_neighbours[ex_gid]:
                if n_gid in excluded_gids:
                    continue

                ani, af = Skani.symmetric_ani_af(ani_af, ex_gid, n_gid)
                if ani > closest_ani:
                    closest_ani = ani
                    closest_gid = n_gid

            ani, af = Skani.symmetric_ani_af(ani_af, ex_gid, closest_gid)

            fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                cur_genomes[closest_gid].ncbi_taxa.species,
                cur_genomes[closest_gid].gtdb_taxa.species,
                closest_gid,
                ','.join(sorted(cur_genomes[closest_gid].strain_ids())),
                ','.join(sorted(cur_genomes[closest_gid].gtdb_type_sources())).upper().replace(
                    'STRAININFO', 'StrainInfo'),
                cur_genomes[closest_gid].year_of_priority(),
                cur_genomes[closest_gid].is_gtdb_type_species(),
                cur_genomes[closest_gid].is_gtdb_type_strain(),
                cur_genomes[closest_gid].ncbi_type_material))
            fout.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                cur_genomes[ex_gid].ncbi_taxa.species,
                cur_genomes[ex_gid].gtdb_taxa.species,
                ex_gid,
                ','.join(sorted(cur_genomes[ex_gid].strain_ids())),
                ','.join(sorted(cur_genomes[ex_gid].gtdb_type_sources())).upper().replace(
                    'STRAININFO', 'StrainInfo'),
                cur_genomes[ex_gid].year_of_priority(),
                cur_genomes[ex_gid].is_gtdb_type_species(),
                cur_genomes[ex_gid].is_gtdb_type_strain(),
                cur_genomes[ex_gid].ncbi_type_material))
            fout.write('\t{:.3f}\t{:.4f}'.format(ani, af))

            if closest_gid in updated_sp_clusters and ex_gid in updated_sp_clusters:
                action = 'ANI_NEIGHBOURS:MERGED:PREVIOUS_REPRESENTATIVES'
            elif closest_gid not in updated_sp_clusters and ex_gid not in updated_sp_clusters:
                action = 'ANI_NEIGHBOURS:MERGED:NEW_REPRESENTATIVES'
            else:
                action = 'ANI_NEIGHBOURS:MERGED:NEW_AND_PREVIOUS_REPRESENTATIVES'

            fout.write('\t{}\n'.format(action))

    def unrepresented_ncbi_sp(self, updated_sp_clusters, cur_genomes):
        """Determining NCBI species unrepresented by GTDB species clusters."""

        ncbi_species = cur_genomes.named_ncbi_species()
        type_strain_genomes = cur_genomes.get_ncbi_type_strain_genomes()
        ncbi_sp_with_type_strains = len(
            set(ncbi_species).intersection(type_strain_genomes))
        self.log.info(
            f' - identified {len(ncbi_species):,} named NCBI species')
        self.log.info(
            f' - identified {ncbi_sp_with_type_strains:,} named NCBI species with type strain genomes')

        # count number of genomes assigned to each named NCBI species
        # and number of these genomes already in GTDB species clusters
        sp_gids = defaultdict(set)
        sp_genome_count = defaultdict(int)
        sp_cluster_count = defaultdict(int)
        for gid in cur_genomes:
            if cur_genomes[gid].ncbi_untrustworthy_sp:
                continue

            ncbi_sp = cur_genomes[gid].ncbi_taxa.species
            sp_genome_count[ncbi_sp] += 1
            sp_gids[ncbi_sp].add(gid)

            if updated_sp_clusters.get_representative(gid) is not None:
                sp_cluster_count[ncbi_sp] += 1

        # determine if NCBI species is not represented by current
        # GTDB species clusters
        unrepresented_ncbi_sp = set()
        type_count = 0
        nontype_count = 0
        for ncbi_sp in ncbi_species:
            if ncbi_sp in self.invalid_ncbi_sp:
                # skip species not recognized under the GTDB
                continue

            type_gids = type_strain_genomes[ncbi_sp]
            type_cluster_count = sum([1 for gid in type_gids
                                      if updated_sp_clusters.get_representative(gid) is not None])

            if len(type_gids) > 0:
                # GTDB is covering this species if it has one of the type strain
                # genomes as a cluster representative or all of the type strains
                # contained in GTDB clusters (i.e. it is a synonym)
                if (len(type_gids.intersection(updated_sp_clusters)) == 0
                        and len(type_gids) != type_cluster_count):
                    unrepresented_ncbi_sp.add(ncbi_sp)
                    type_count += 1

                    if type_cluster_count > 0:
                        self.log.warning('Identified {:,} of {:,} type strain genomes for {} in GTDB species clusters.'.format(
                            type_cluster_count,
                            len(type_gids),
                            ncbi_sp))
            else:
                # species has no type strain genomes, so we check if a genome with
                # this species name is a GTDB representative or if the majority of
                # genomes from this species are already in GTDB species clusters
                if (len(sp_gids[ncbi_sp].intersection(updated_sp_clusters)) == 0
                        and sp_cluster_count[ncbi_sp] <= 0.5*sp_genome_count[ncbi_sp]):
                    unrepresented_ncbi_sp.add(ncbi_sp)
                    nontype_count += 1

                    if sp_cluster_count[ncbi_sp] > 0.5*sp_genome_count[ncbi_sp]:
                        self.log.warning('Identified {:,} of {:,} {} genomes in GTDB species clusters.'.format(
                            sp_cluster_count[ncbi_sp],
                            sp_genome_count[ncbi_sp],
                            ncbi_sp))

        self.log.info(
            f' - identified {type_count:,} unrepresented NCBI species with type strain genomes')
        self.log.info(
            f' - identified {nontype_count:,} unrepresented NCBI species without type strain genomes')
        self.log.info(
            f' - identified {len(unrepresented_ncbi_sp):,} total unrepresented NCBI species')

        return unrepresented_ncbi_sp

    def run(self,
            updated_sp_cluster_file,
            cur_gtdb_metadata_file,
            cur_genomic_path_file,
            qc_passed_file,
            ncbi_genbank_assembly_file,
            untrustworthy_type_file,
            gtdb_type_strains_ledger,
            ncbi_untrustworthy_sp_ledger,
            sp_priority_ledger,
            genus_priority_ledger,
            ncbi_env_bioproject_ledger,
            lpsn_gss_file):
        """Select GTDB type genomes for named species."""

        # read updated GTDB species clusters
        self.log.info('Reading updated GTDB species clusters:')
        updated_sp_clusters = SpeciesClusters()
        updated_sp_clusters.load_from_sp_cluster_file(updated_sp_cluster_file)
        self.log.info(' - identified {:,} updated species clusters spanning {:,} genomes.'.format(
            len(updated_sp_clusters),
            updated_sp_clusters.total_num_genomes()))

        # create current GTDB genome sets
        self.log.info('Creating current GTDB genome set:')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                            gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                            create_sp_clusters=False,
                                            qc_passed_file=qc_passed_file,
                                            ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                            untrustworthy_type_ledger=untrustworthy_type_file,
                                            ncbi_untrustworthy_sp_ledger=ncbi_untrustworthy_sp_ledger,
                                            ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)

        # get path to previous and current genomic FASTA files
        self.log.info('Reading path to current genomic FASTA files.')
        cur_genomes.load_genomic_file_paths(cur_genomic_path_file)

        # determine new NCBI species requiring a GTDB representative to be selected
        self.log.info(
            'Determining NCBI species unrepresented by GTDB species clusters:')
        unrepresented_ncbi_sp = self.unrepresented_ncbi_sp(
            updated_sp_clusters, cur_genomes)

        # report type material
        self.log.info('Tabulating genomes assigned as type material:')
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

        self.log.info(' - identified {:,} species spanning {:,} genomes designated as type species of genus by GTDB'.format(
            len(gtdb_type_genus),
            sum([len(gids) for gids in gtdb_type_genus.values()])))
        self.log.info(' - identified {:,} species spanning {:,} genomes designated as type strain of species by GTDB'.format(
            len(gtdb_type_sp),
            sum([len(gids) for gids in gtdb_type_sp.values()])))
        self.log.info(' - identified {:,} species spanning {:,} genomes designated as type strain of subspecies or synonym by GTDB'.format(
            len(gtdb_type_subsp),
            sum([len(gids) for gids in gtdb_type_subsp.values()])))
        self.log.info(' - identified {:,} species spanning {:,} genomes designated as assembled from type material by NCBI'.format(
            len(ncbi_type_sp),
            sum([len(gids) for gids in ncbi_type_sp.values()])))
        self.log.info(' - identified {:,} species spanning {:,} genomes designated as assembled from proxytype material by NCBI'.format(
            len(ncbi_proxy),
            sum([len(gids) for gids in ncbi_proxy.values()])))
        self.log.info(' - identified {:,} species spanning {:,} genomes designated as assembly from synonym type material by NCBI'.format(
            len(ncbi_type_subsp),
            sum([len(gids) for gids in ncbi_type_subsp.values()])))
        self.log.info(' - identified {:,} species spanning {:,} genomes designated as a reference or representative by NCBI'.format(
            len(ncbi_reps),
            sum([len(gids) for gids in ncbi_reps.values()])))

        fout = open(os.path.join(self.output_dir,
                                 'type_material_stats.tsv'), 'w')
        fout.write('Type\tNo. taxa\tNo. genomes\n')
        fout.write('%s\t%d\t%d\n' % ('GTDB type species of genus', len(
            gtdb_type_genus), sum([len(gids) for gids in gtdb_type_genus.values()])))
        fout.write('%s\t%d\t%d\n' % ('GTDB type strain of species', len(
            gtdb_type_sp), sum([len(gids) for gids in gtdb_type_sp.values()])))
        fout.write('%s\t%d\t%d\n' % ('GTDB type strain of subspecies', len(
            gtdb_type_subsp), sum([len(gids) for gids in gtdb_type_subsp.values()])))
        fout.write('%s\t%d\t%d\n' % ('NCBI assembled from type material', len(
            ncbi_type_sp), sum([len(gids) for gids in ncbi_type_sp.values()])))
        fout.write('%s\t%d\t%d\n' % ('NCBI assembled from proxytype material', len(
            ncbi_proxy), sum([len(gids) for gids in ncbi_proxy.values()])))
        fout.write('%s\t%d\t%d\n' % ('NCBI assembly from synonym type material', len(
            ncbi_type_subsp), sum([len(gids) for gids in ncbi_type_subsp.values()])))
        fout.write('%s\t%d\t%d\n' % ('NCBI reference or representative', len(
            ncbi_reps), sum([len(gids) for gids in ncbi_reps.values()])))
        fout.close()

        self.validate_type_designations(cur_genomes,
                                        gtdb_type_sp,
                                        gtdb_type_subsp,
                                        ncbi_type_sp,
                                        ncbi_proxy,
                                        ncbi_type_subsp,
                                        ncbi_reps)

        # select representative for new named NCBI species
        self.log.info('Selecting representative for each of the {:,} unrepresentative named NCBI species:'.format(
            len(unrepresented_ncbi_sp)))

        if True:  # ***
            new_rep_genomes = self.select_rep_genomes(unrepresented_ncbi_sp,
                                                      cur_genomes,
                                                      gtdb_type_sp,
                                                      gtdb_type_subsp,
                                                      ncbi_type_sp,
                                                      ncbi_proxy,
                                                      ncbi_type_subsp,
                                                      ncbi_reps)

            pickle.dump(new_rep_genomes, open(os.path.join(
                self.output_dir, 'new_rep_genomes.pkl'), 'wb'))
        else:
            new_rep_genomes = pickle.load(
                open(os.path.join(self.output_dir, 'new_rep_genomes.pkl'), 'rb'))

        self.log.info(
            f' - identified representatives for {len(new_rep_genomes):,} new NCBI named species')

        # report species names present twice
        dup_sp = set(updated_sp_clusters.species_names.values()
                     ).intersection(new_rep_genomes.values())
        for cur_sp in dup_sp:
            # This does happen and indicates that an existing GTDB species name will need to be updated.
            # e.g., GTDB R89 has a Idiomarina aestuarii species cluster, but this is not represented by a
            # type strain genome. In R95, there is a type strain genome for I. aestuarii, but this does
            # not cluster in the existing species cluster. The current I. aestuarii cluster needs a new
            # name since it is incorrect based on the type strain genome.
            self.log.warning('{} is represented by both {} and {}.'.format(
                cur_sp,
                [(gid, cur_genomes[gid].is_effective_type_strain(
                )) for gid, sp in updated_sp_clusters.species_names.items() if sp == cur_sp],
                [(gid, cur_genomes[gid].is_effective_type_strain()) for gid, sp in new_rep_genomes.items() if sp == cur_sp]))

        # calculate ANI between representatives and resolve cases where type genomes have close ANI neighbours
        all_rep_genomes = {}
        for gid, sp in updated_sp_clusters.species_names.items():
            all_rep_genomes[gid] = sp

        for gid, sp in new_rep_genomes.items():
            all_rep_genomes[gid] = sp

        self.log.info(
            f'Calculate ANI between {len(all_rep_genomes):,} representatives of GTDB species clusters.')

        ani_af = self.ani_reps(cur_genomes, all_rep_genomes)

        self.log.info('Establishing ANI neighbours.')
        ani_neighbours = self.ani_neighbours(updated_sp_clusters,
                                             cur_genomes,
                                             ani_af)

        self.log.info('Resolving ANI neighbours.')
        excluded_gids = self.resolve_close_ani_neighbours(cur_genomes,
                                                          ani_neighbours,
                                                          gtdb_type_genus,
                                                          gtdb_type_sp,
                                                          gtdb_type_subsp,
                                                          ncbi_type_sp,
                                                          ncbi_proxy,
                                                          ncbi_type_subsp,
                                                          ncbi_reps,
                                                          sp_priority_ledger,
                                                          genus_priority_ledger,
                                                          lpsn_gss_file)

        self.write_final_reps(cur_genomes,
                              all_rep_genomes,
                              excluded_gids)

        self.write_ani_neighbour_table(updated_sp_clusters,
                                       cur_genomes,
                                       ani_af,
                                       ani_neighbours,
                                       excluded_gids)
