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
from collections import namedtuple


def same_assembly_version(ncbi_accn1, ncbi_accn2):
    """Check if NCBi accessions have same version number."""

    return int(ncbi_accn1.split('.')[1]) == int(ncbi_accn2.split('.')[1])


def select_highest_quality(gids, cur_genomes):
    """Select highest quality genome."""

    # sort by decreasing type strain score, followed by
    # genome ID in order to ensure ties are broken in a
    # deterministic fashion between runs
    q = {k: cur_genomes[k].score_type_strain() for k in gids}
    q_sorted = sorted(q.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)

    return q_sorted[0][0]


def read_genome_path(genome_path_file):
    """Determine path to genomic FASTA file for each genome."""

    genome_files = {}
    for line in open(genome_path_file):
        line_split = line.strip().split('\t')

        gid = line_split[0]
        gid = canonical_gid(gid)

        genome_path = line_split[1]
        accession = os.path.basename(os.path.normpath(genome_path))

        genome_files[gid] = os.path.join(
            genome_path, accession + '_genomic.fna')

    return genome_files


def read_gtdb_metadata(metadata_file, fields):
    """Parse genome quality from GTDB metadata.

    Parameters
    ----------
    metadata_file : str
        Metadata for all genomes in CSV file.
    fields : iterable
        Fields  to read.

    Return
    ------
    dict : d[genome_id] -> namedtuple
        Value for fields indicted by genome IDs.
    """

    gtdb_metadata = namedtuple('gtdb_metadata', ' '.join(fields))
    m = {}

    with open(metadata_file, encoding='utf-8') as f:
        headers = f.readline().strip().split('\t')

        genome_index = headers.index('accession')

        indices = []
        for field in fields:
            indices.append(headers.index(field))

        for line in f:
            line_split = line.strip().split('\t')
            genome_id = canonical_gid(line_split[genome_index])

            values = []
            for i in indices:
                # save values as floats or strings
                v = line_split[i]
                try:
                    values.append(float(v))
                except ValueError:
                    if v is None or v == '' or v == 'none':
                        values.append(None)
                    elif v == 'f' or v.lower() == 'false':
                        values.append(False)
                    elif v == 't' or v.lower() == 'true':
                        values.append(True)
                    else:
                        values.append(v)
            m[genome_id] = gtdb_metadata._make(values)

    return m


def exclude_from_refseq(genbank_assembly_file):
    """Parse exclude from RefSeq field from NCBI assembly files."""

    excluded_from_refseq_note = {}
    with open(genbank_assembly_file, encoding='utf-8') as f:
        for line in f:
            if line[0] == '#':
                if line.startswith('# assembly_accession'):
                    header = line.strip().split('\t')
                    exclude_index = header.index('excluded_from_refseq')
            else:
                line_split = line.strip('\n\r').split('\t')
                gid = canonical_gid(line_split[0])
                excluded_from_refseq_note[gid] = line_split[exclude_index]

    return excluded_from_refseq_note


def parse_ncbi_bioproject(genbank_assembly_file):
    """Parse BioProject field from NCBI assembly files."""

    bioproject = {}
    with open(genbank_assembly_file, encoding='utf-8') as f:
        for line in f:
            if line[0] == '#':
                if line.startswith('# assembly_accession'):
                    header = line.strip().split('\t')
                    bioproject_index = header.index('bioproject')
            else:
                line_split = line.strip('\n\r').split('\t')
                gid = canonical_gid(line_split[0])
                bioproject[gid] = line_split[bioproject_index]

    return bioproject


def read_qc_file(qc_file):
    """Read genomes passing QC from file."""

    passed_qc = set()
    with open(qc_file, encoding='utf-8') as f:
        f.readline()

        for line in f:
            line_split = line.strip().split('\t')
            gid = canonical_gid(line_split[0])
            passed_qc.add(gid)

    return passed_qc


def read_cur_new_updated(genomes_new_updated_file):
    """Determine new and updated genomes in current GTDB release."""

    cur_new = set()
    cur_updated = set()
    with open(genomes_new_updated_file, encoding='utf-8') as f:
        header = f.readline().strip().split('\t')

        status_index = header.index('Status')

        for line in f:
            line_split = line.strip().split('\t')

            gid = line_split[0]
            status = line_split[status_index]
            if status == 'NEW':
                cur_new.add(gid)
            elif status == 'UPDATED':
                cur_updated.add(gid)

    return cur_new, cur_updated


def read_gtdbtk_classifications(gtdbtk_classify_file):
    """Determine classification of genomes according to GTDB-Tk."""

    classification = {}
    with open(gtdbtk_classify_file) as f:
        header = f.readline().strip().split('\t')

        classification_index = header.index('classification')

        for line in f:
            line_split = line.strip().split('\t')

            gid = canonical_gid(line_split[0])
            taxa = [t.strip()
                    for t in line_split[classification_index].split(';')]

            # covert degenerate cases into 7 rank classifications
            if taxa[0] == 'Unclassified Bacteria':
                taxa = ['d__Bacteria', 'p__', 'c__',
                        'o__', 'f__', 'g__', 's__']
            elif taxa[0] == 'Unclassified Archaea':
                taxa = ['d__Archaea', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

            classification[gid] = taxa

    return classification
