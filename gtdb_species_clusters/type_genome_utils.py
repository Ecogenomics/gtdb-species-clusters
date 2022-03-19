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

import sys
import logging
from collections import namedtuple

from numpy import (mean as np_mean)

from gtdb_species_clusters.genome_utils import canonical_gid, read_gtdb_metadata


NCBI_TYPE_SPECIES = set(['assembly from type material',
                         'assembly from neotype material',
                         'assembly designated as neotype'])
NCBI_PROXYTYPE = set(['assembly from proxytype material'])
NCBI_TYPE_SUBSP = set(['assembly from synonym type material'])

GTDB_TYPE_SPECIES = set(['type strain of species', 'type strain of neotype'])
GTDB_TYPE_SUBSPECIES = set(['type strain of subspecies',
                            'type strain of heterotypic synonym'])
GTDB_NOT_TYPE_MATERIAL = set(['not type material'])


ClusteredGenome = namedtuple('ClusteredGenome', 'ani af gid')
GenomeRadius = namedtuple('GenomeRadius', 'ani af neighbour_gid')


def parse_disbanded_cluster_ledger(disbanded_cluster_ledger):
    """Parse file indicating GTDB species clusters to be disbanded."""

    disbanded = set()
    with open(disbanded_cluster_ledger) as f:
        f.readline()
        for line in f:
            tokens = line.strip().split('\t')
            disbanded.add(canonical_gid(tokens[0]))

    return disbanded


def parse_manual_sp_curation_files(manual_sp_names, pmc_custom_species):
    """Parse files indicating manually curated species names."""

    # read species names explicitly set via manual curation
    logging.getLogger('timestamp').info('Parsing manually-curated species.')
    mc_species = {}
    if not (manual_sp_names is None or manual_sp_names.lower() == 'none'):
        with open(manual_sp_names) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                mc_species[tokens[0]] = tokens[2]
    logging.getLogger('timestamp').info(' - identified manually-curated species names for {:,} genomes.'.format(
        len(mc_species)))

    # read post-curation, manually defined species
    logging.getLogger('timestamp').info(
        'Parsing post-curation, manually-curated species.')
    pmc_species = 0
    if not (pmc_custom_species is None or pmc_custom_species.lower() == 'none'):
        with open(pmc_custom_species) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                gid = tokens[0]
                species = tokens[1]
                if gid in mc_species:
                    logging.getLogger('timestamp').warning('Manually-curated genome {} reassigned from {} to {}.'.format(
                        gid, mc_species[gid], species))
                pmc_species += 1
                mc_species[gid] = species
    logging.getLogger('timestamp').info(' - identified post-curation, manually-curated species names for {:,} genomes.'.format(
        pmc_species))

    return mc_species


def parse_updated_species_reps(updated_species_reps):
    """Get map indicating the updating of GTDB representatives."""

    new_to_prev_rid = {}
    with open(updated_species_reps) as f:
        f.readline()
        for line in f:
            tokens = line.strip().split('\t')

            prev_rid = tokens[0]
            new_rid = tokens[1]

            if new_rid.lower() != 'none':
                new_to_prev_rid[new_rid] = prev_rid

    return new_to_prev_rid


def resolve_merged_prev_representatives(cur_rid, prev_genomes, updated_gtdb_rids, prev_rids_in_cluster):
    """Determine best representative to associated with new species cluster containing multiple previous representative.

    This association between previous and new representatives is used 
    to establish the most suitable name for a species cluster.
    """

    sel_prev_rid = None
    if cur_rid in prev_genomes:
        # if current representative was a previous representative,
        # it is the most natural choice
        sel_prev_rid = cur_rid
    elif cur_rid in updated_gtdb_rids:
        prev_rid = updated_gtdb_rids[cur_rid]
        if prev_rid in prev_rids_in_cluster:
            # makes sense to use the previous representative that was explicitly
            # updated to the new, current representative
            sel_prev_rid = prev_rid
        else:
            logging.getLogger('timestamp').error('Updated representative of cluster {} no longer contains the previous representative of cluster: {}'.format(
                cur_rid,
                prev_rid))
            sys.exit(-1)
    else:
        # the current representative of the cluster is a new genome which doesn't
        # provide any guidance on which of the previous representatives to select
        # so just pick the genome with the highest assembly score
        highest_prev_gid = None
        highest_score = 0
        for prev_rid in prev_rids_in_cluster:
            if prev_genomes[prev_rid].score_assembly() > highest_score:
                highest_score = prev_genomes[prev_rid].score_assembly()
                highest_prev_gid = prev_rid

        sel_prev_rid = highest_prev_gid

    return sel_prev_rid


def infer_prev_gtdb_reps(prev_genomes, cur_clusters, updated_gtdb_rids):
    """Infer previous GTDB representative for each current GTDB representative.

    This can't just be taken from the `updated_species_reps` file since the
    final results of de novo cluster can, in some rare cases, cause small movements
    in the genomes associated with each species clusters and the generation or loss of
    a representative.
    """

    prev_rids = set(prev_genomes.sp_clusters)
    new_to_prev_rid = {}
    for cur_rid, cur_cids in cur_clusters.items():
        prev_rids_in_cluster = prev_rids.intersection(cur_cids)
        if len(prev_rids_in_cluster) == 1:
            new_to_prev_rid[cur_rid] = prev_rids_in_cluster.pop()
        elif len(prev_rids_in_cluster) > 1:
            resolved_rid = resolve_merged_prev_representatives(cur_rid,
                                                               prev_genomes,
                                                               updated_gtdb_rids,
                                                               prev_rids_in_cluster)

            if resolved_rid:
                logging.getLogger('timestamp').info(
                    ' - cluster {} contains multiple previous representatives and has been associated with representative {}.'.format(
                        cur_rid, resolved_rid))

                new_to_prev_rid[cur_rid] = resolved_rid

    return new_to_prev_rid


def write_clusters(clusters, rep_radius, genomes, out_file):
    """Write out clustering information."""

    fout = open(out_file, 'w')
    fout.write('Representative\tGTDB species\tNCBI species')
    fout.write(
        '\tClosest GTDB species\tClosest representative\tANI radius\tAF closest')
    fout.write(
        '\tNo. clustered genomes\tMean ANI\tMin ANI\tMean AF\tMin AF\tClustered genomes\n')
    for gid in sorted(clusters, key=lambda x: len(clusters[x]), reverse=True):
        if clusters[gid]:
            mean_ani = '%.2f' % np_mean([d.ani for d in clusters[gid]])
            min_ani = '%.2f' % min([d.ani for d in clusters[gid]])
            mean_af = '%.2f' % np_mean([d.af for d in clusters[gid]])
            min_af = '%.2f' % min([d.af for d in clusters[gid]])
        else:
            mean_ani = min_ani = mean_af = min_af = 'N/A'
        fout.write('%s\t%s\t%s' % (gid,
                                   genomes[gid].gtdb_taxa.species,
                                   genomes[gid].ncbi_taxa.species))

        ani, af, closest_gid = rep_radius[gid]
        if not af:
            af = 0

        if not closest_gid or closest_gid == 'N/A':
            closest_gid = 'N/A'
            closest_sp = 'N/A'
        else:
            closest_sp = genomes[closest_gid].gtdb_taxa.species

        fout.write('\t%s\t%s\t%f\t%f' % (closest_sp,
                                         closest_gid,
                                         ani,
                                         af))

        fout.write('\t%d\t%s\t%s\t%s\t%s\t%s\n' % (
            len(clusters[gid]),
            mean_ani, min_ani,
            mean_af, min_af,
            ','.join([d.gid for d in clusters[gid]])))
    fout.close()


def read_clusters(cluster_file):
    """Read cluster file."""

    clusters = {}
    rep_radius = {}
    with open(cluster_file) as f:
        headers = f.readline().strip().split('\t')

        rep_index = headers.index('Representative')

        num_clustered_index = headers.index('No. clustered genomes')
        clustered_genomes_index = headers.index('Clustered genomes')
        closest_gid_index = headers.index('Closest representative')
        closest_ani_index = headers.index('ANI radius')
        closest_af_index = headers.index('AF closest')

        for line in f:
            line_split = line.strip().split('\t')

            rid = line_split[rep_index]

            num_clustered = int(line_split[num_clustered_index])
            clusters[rid] = set([rid])
            if num_clustered > 0:
                clusters[rid].update(
                    [g.strip() for g in line_split[clustered_genomes_index].split(',')])

            closest_ani = float(line_split[closest_ani_index])
            closest_af = float(line_split[closest_af_index])
            closest_gid = line_split[closest_gid_index]
            rep_radius[rid] = GenomeRadius(ani=closest_ani,
                                           af=closest_af,
                                           neighbour_gid=closest_gid)

    return clusters, rep_radius


def write_rep_radius(rep_radius, genomes, out_file):
    """Write out ANI radius for each representative genomes."""

    fout = open(out_file, 'w')
    fout.write(
        'Representative\tGTDB species\tNCBI species\tANI\tAF\tClosest species\tClosest representative\n')

    for gid in rep_radius:
        ani, af, neighbour_gid = rep_radius[gid]
        if not af:
            af = 0

        if not neighbour_gid or neighbour_gid == 'N/A':
            neighbour_gid = 'N/A'
            neighbour_sp = 'N/A'
        else:
            neighbour_sp = genomes[neighbour_gid].gtdb_taxa.species

        fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gid,
                                                         genomes[gid].gtdb_taxa.species,
                                                         genomes[gid].ncbi_taxa.species,
                                                         ani,
                                                         af,
                                                         neighbour_sp,
                                                         neighbour_gid))
    fout.close()


def read_quality_metadata(metadata_file):
    """Read statistics needed to determine genome quality."""

    return read_gtdb_metadata(metadata_file, ['gtdb_taxonomy',
                                              'checkm_completeness',
                                              'checkm_contamination',
                                              'checkm_strain_heterogeneity_100',
                                              'genome_size',
                                              'contig_count',
                                              'n50_contigs',
                                              'scaffold_count',
                                              'ambiguous_bases',
                                              'total_gap_length',
                                              'ssu_count',
                                              'ssu_length',
                                              'mimag_high_quality',
                                              'gtdb_type_designation',
                                              'ncbi_assembly_level',
                                              'ncbi_genome_representation',
                                              'ncbi_refseq_category',
                                              'ncbi_type_material_designation',
                                              'ncbi_molecule_count',
                                              'ncbi_unspanned_gaps',
                                              'ncbi_spanned_gaps',
                                              'ncbi_genome_category'])
