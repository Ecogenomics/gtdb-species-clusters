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

import logging

class Update_ANI_Cache():
    """Update ANI cache to remove entries for updates genomes."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

    def run(self, genomes_new_updated_file, prev_ani_cache, out_ani_cache):
        """Update ANI cache to remove entries for updates genomes."""

        # get list of new or updated genomes. In theory, the ANI cache
        # shouldn't contain any genomes marked as NEW so it seems wies
        # to filtered them here in case this situation does somehow
        # occur.
        self.logger.info('Identifying new or updated genomes:')
        new_updated_gids = set()
        with open(genomes_new_updated_file, encoding='utf-8') as f:
            header = f.readline().strip().split('\t')

            status_idx = header.index('Status')

            for line in f:
                tokens = line.strip().split('\t')

                gid = tokens[0]
                if tokens[status_idx].upper() in ['NEW', 'UPDATED']:
                    new_updated_gids.add(gid)

        self.logger.info(
            f' - identified {len(new_updated_gids):,} new or updated genomes to remove from ANI cache')

        # remove updated genomes from ANI cache
        self.logger.info('Updating ANI cache:')
        fout = open(out_ani_cache, 'w')
        filtered_rows = 0
        retained_rows = 0
        with open(prev_ani_cache) as f:
            for line in f:
                tokens = line.strip().split('\t')

                gid1 = tokens[0]
                if gid1 in new_updated_gids:
                    filtered_rows += 1
                    continue

                gid2 = tokens[1]
                if gid2 in new_updated_gids:
                    filtered_rows += 1
                    continue

                retained_rows += 1
                fout.write(line)

        fout.close()

        self.logger.info(f' - filtered {filtered_rows:,} rows and retained {retained_rows:,} rows')
                