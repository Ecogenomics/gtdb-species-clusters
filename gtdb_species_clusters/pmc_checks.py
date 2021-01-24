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

import dendropy

from biolib.taxonomy import Taxonomy

from taxon_utils import specific_epithet


class PMC_Checks():
    """Perform post-manual curation checks"""

    def __init__(self, output_dir):
        """Initialization."""

        self.output_dir = output_dir

        self.logger = logging.getLogger('timestamp')

    def manual_species(self, init_taxonomy, manually_curated_tree):
        """Identify species names manually set by curators."""

        # read initial and manually curated taxonomy
        self.logger.info('Reading initial species names.')
        init_taxonomy = Taxonomy().read(init_taxonomy, use_canonical_gid=True)
        init_num_gids = sum(
            [1 for gid in init_taxonomy if not gid.startswith('D-')])
        self.logger.info(
            ' - read taxonomy for {:,} genomes.'.format(init_num_gids))

        self.logger.info('Reading manually-curated species names from tree.')
        mc_tree = dendropy.Tree.get_from_path(manually_curated_tree,
                                              schema='newick',
                                              rooting='force-rooted',
                                              preserve_underscores=True)
        mc_taxonomy = Taxonomy().read_from_tree(mc_tree)

        mc_specific = {}
        for gid, taxa in mc_taxonomy.items():
            if gid.startswith('D-'):
                continue

            mc_sp = taxa[-1]
            if not mc_sp.startswith('s__') or mc_sp == 's__':
                self.logger.error(
                    'Most specific classification for {} is {}.'.format(gid, taxa))
                continue

            mc_specific[gid] = specific_epithet(mc_sp)

        self.logger.info(
            ' - read taxonomy for {:,} genomes.'.format(len(mc_specific)))

        # report genomes with modified specific name assignment
        self.logger.info(
            'Identifying genomes with manually-curated species names.')
        fout = open(os.path.join(self.output_dir,
                                 'manual_species_names.tsv'), 'w')
        fout.write('Genome ID\tInitial species\tManually-curated species\n')
        num_mc = 0
        for gid, mc_sp in mc_specific.items():
            init_species = init_taxonomy[gid][Taxonomy.SPECIES_INDEX]
            init_specific = specific_epithet(init_species)

            if init_specific != mc_sp:
                mc_generic = mc_taxonomy[gid][Taxonomy.GENUS_INDEX].replace(
                    'g__', '')
                mc_species = 's__{} {}'.format(mc_generic, mc_sp)
                num_mc += 1
                fout.write('{}\t{}\t{}\n'.format(
                    gid,
                    init_species,
                    mc_species))

        fout.close()

        self.logger.info(
            ' - identified {:,} manually-curated species names.'.format(num_mc))

    def replace_generic(self, manual_species_names, manual_taxonomy):
        """Replace generic names with genus assignment."""

        # read manually-curated species names
        self.logger.info('Reading manually-curated species names.')
        mc_species = {}
        with open(manual_species_names) as f:
            f.readline()

            for line in f:
                tokens = line.strip().split('\t')
                mc_species[tokens[0]] = tokens[2]
        self.logger.info(
            ' - read manually-curated species for {:,} genomes.'.format(len(mc_species)))

        # read manual taxonomy file
        self.logger.info('Reading manually-curated taxonomy.')
        mc_taxonomy = Taxonomy().read(manual_taxonomy, use_canonical_gid=True)
        mc_num_gids = sum(
            [1 for gid in mc_taxonomy if not gid.startswith('D-')])
        self.logger.info(
            ' - read taxonomy for {:,} genomes.'.format(mc_num_gids))

        # replace generic names with genus names
        self.logger.info('Creating taxonomy file with updated species names.')
        fout = open(os.path.join(self.output_dir,
                                 'taxonomy_updated_sp.tsv'), 'w')
        num_genomes = 0
        for gid, taxa in mc_taxonomy.items():
            if gid.startswith('D-'):
                continue

            genus = taxa[Taxonomy.GENUS_INDEX]
            generic = genus.replace('g__', '')
            if not generic:
                self.logger.error(
                    'Genome is missing genus assignment: {}'.format(gid))
                sys.exit(-1)

            species = taxa[Taxonomy.SPECIES_INDEX]

            if gid in mc_species:
                if generic not in species and species != 's__':
                    self.logger.error('Genus assignment does not agree with manually-curated species assignment: {} {} {}'.format(
                        gid,
                        mc_species[gid],
                        '; '.join(mc_taxonomy[gid])))

            sp_tokens = species.split()
            if len(sp_tokens) < 2:
                self.logger.error(
                    'Species name appear to be erroneous: {} {}'.format(gid, species))
                specific = '<unassigned>'
            else:
                specific = species.split()[-1]

            final_sp = 's__{} {}'.format(generic, specific)
            taxa[Taxonomy.SPECIES_INDEX] = final_sp
            fout.write('{}\t{}\n'.format(
                gid,
                ';'.join(taxa)))

            num_genomes += 1

        fout.close()
        self.logger.info(' - processed {:,} genomes.'.format(num_genomes))
