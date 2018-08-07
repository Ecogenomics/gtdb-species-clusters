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
import csv
import logging

import biolib.seq_io as seq_io
from biolib.taxonomy import Taxonomy

from genometreetk.common import (read_gtdb_metadata,
                                    read_genome_dir_file,
                                    read_gtdb_taxonomy)


class Arb(object):
    """Methods for producing files for ARB."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger()
        
    def _record(self, fout, 
                    genome_id,
                    metadata_fields, 
                    metadata_values,
                    aligned_seq):
        """Write out ARB record for genome."""

        fout.write("BEGIN\n")
        fout.write("db_name=%s\n" % genome_id)
        for col_header, value in zip(metadata_fields, metadata_values):
            # replace equal signs as these are incompatible with the ARB parser
            if value:
                value = value.replace('=', '/')

            fout.write("%s=%s\n" % (col_header, value))
        
        fout.write("warning=\n")
        fout.write("aligned_seq=%s\n" % (aligned_seq))
        fout.write("END\n\n")

    def create_records(self, metadata_file, msa_file, taxonomy_file, genome_list, output_file):
        """Create ARB records from GTDB metadata."""
        
        seqs = {}
        if msa_file:
            seqs = seq_io.read(msa_file)
        
        taxonomy = {}
        if taxonomy_file:
            taxonomy = Taxonomy().read(taxonomy_file)
            
        genomes_to_keep = set()
        if genome_list:
            for line in open(genome_list):
                genomes_to_keep.add(line.strip())
                
        fout = open(output_file, 'w')
        
        delimiter = ','
        if metadata_file.endswith('.tsv'):
            delimiter = '\t'
        
        header = True
        for row in csv.reader(open(metadata_file, 'rb'), delimiter=delimiter):
            if header:
                fields = [f.lower().replace(' ', '_').replace('-', '_') for f in row[1:]]
                if taxonomy:
                    fields.append('gtdb_taxonomy')
                header = False
            else:
                genome_id = row[0]
                values = row[1:]
                if taxonomy:
                    values.append('; '.join(taxonomy[genome_id]))
                aligned_seq = seqs.get(genome_id, '')
                
                if not genomes_to_keep or genome_id in genomes_to_keep:
                    self._record(fout, genome_id, fields, values, aligned_seq)

        fout.close()
        
    def write(self, hashes, io):
        """Write data to a GreenGenes formatted files.

        Parameters
        ----------
        hashes : list of dict[feature] -> value
            GreenGenes style metadata.
        io : output stream
            Handle to output stream.
        """
        is_first = True
        for dahash in hashes:
            if is_first:
                is_first = False
            else:
                io.write('\n')

            io.write('BEGIN\n')
            for key in sorted(dahash.keys()):
                if key != 'warning' and key != 'aligned_seq':
                    io.write('='.join([key, dahash[key]]) + '\n')

            # the warning field must be the second to last as it is used
            # to indicates that the aligned sequence is to follow
            io.write('='.join(['warning', dahash.get('warning', '')]) + '\n')

            # the aligned sequence must be the last field
            io.write('='.join(['aligned_seq', dahash.get('aligned_seq', '')]) + '\n')

            io.write('END\n')
            