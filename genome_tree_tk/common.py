###############################################################################
#
# common.py - utility functions used in many places in CheckM
#
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


def readTreeModel(reportFile):
    for line in open(reportFile):
        if 'Model of evolution:' in line:
            modelStr = line[line.find(':') + 1:].strip()

    return modelStr


def read_genome_ids(genomeFile):
    """Read genome ids from file."""

    genome_ids = set()
    img_genome_ids = set()
    user_genome_ids = set()
    for line in open(genomeFile):
        if '\t' in line:
            genome_id = line.split('\t')[0].strip()
        else:
            genome_id = line.split()[0].strip()
        genome_ids.add(genome_id)

        if genome_id.startswith('IMG_'):
            img_genome_ids.add(genome_id.replace('IMG_', ''))
        else:
            user_genome_ids.add(genome_id)

    return genome_ids, img_genome_ids, user_genome_ids
