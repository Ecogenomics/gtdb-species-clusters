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
import re
import random
import shutil
import tempfile
import ntpath
import pickle
from itertools import combinations
from collections import defaultdict, namedtuple, Counter

from biolib.external.execute import check_dependencies

from numpy import (mean as np_mean,
                    std as np_std)

from gtdb_species_clusters.common import read_gtdb_metadata
                                    
from gtdb_species_clusters.genome_utils import (read_genome_path,
                                                exclude_from_refseq,
                                                read_qc_file)
                                                
from gtdb_species_clusters.taxon_utils import (read_gtdb_taxonomy,
                                            read_gtdb_ncbi_taxonomy)
                                    
from gtdb_species_clusters.type_genome_utils import (GenomeRadius,
                                            quality_score,
                                            read_quality_metadata,
                                            symmetric_ani,
                                            write_clusters,
                                            write_rep_radius)
                                    
from gtdb_species_clusters.fastani import FastANI
from gtdb_species_clusters.mash import Mash

class ClusterDeNovo(object):
    """Infer de novo species clusters and type genomes for remaining genomes."""

    def __init__(self, ani_sp, af_sp, ani_cache_file, cpus, output_dir):
        """Initialization."""
        
        check_dependencies(['fastANI', 'mash'])
        
        self.cpus = cpus
        self.output_dir = output_dir

        self.logger = logging.getLogger('timestamp')
        
        self.true_str = ['t', 'T', 'true', 'True']
        
        self.ani_sp = ani_sp
        self.af_sp = af_sp

        self.min_mash_ani = 90.0
        
        self.ClusteredGenome = namedtuple('ClusteredGenome', 'ani af gid')
        
        self.fastani = FastANI(ani_cache_file, cpus)
        
    def _parse_type_clusters(self, type_genome_cluster_file):
        """Parse type genomes clustering information."""
        
        type_species = set()
        species_type_gid = {}
        type_gids = set()
        type_clustered_gids = set()
        type_radius = {}
        with open(type_genome_cluster_file) as f:
            headers = f.readline().strip().split('\t')
            
            type_sp_index = headers.index('NCBI species')
            type_genome_index = headers.index('Type genome')
            num_clustered_index = headers.index('No. clustered genomes')
            clustered_genomes_index = headers.index('Clustered genomes')
            closest_type_index = headers.index('Closest type genome')
            ani_radius_index = headers.index('ANI radius')
            af_index = headers.index('AF closest')

            for line in f:
                line_split = line.strip().split('\t')
                
                type_sp = line_split[type_sp_index]
                type_species.add(type_sp)
                
                type_gid = line_split[type_genome_index]
                type_gids.add(type_gid)
                
                species_type_gid[type_gid] = type_sp
                
                num_clustered = int(line_split[num_clustered_index])
                if num_clustered > 0:
                    for gid in [g.strip() for g in line_split[clustered_genomes_index].split(',')]:
                        type_clustered_gids.add(gid)
                        
                type_radius[type_gid] = GenomeRadius(ani = float(line_split[ani_radius_index]), 
                                                     af = float(line_split[af_index]),
                                                     neighbour_gid = line_split[closest_type_index])
                        
        return type_species, species_type_gid, type_gids, type_clustered_gids, type_radius
        
    def _parse_synonyms(self, type_genome_synonym_file):
        """Parse synonyms."""
        
        synonyms = set()
        with open(type_genome_synonym_file) as f:
            headers = f.readline().strip().split('\t')
            
            synonym_index = headers.index('Synonym')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                synonym = line_split[synonym_index]
                synonyms.add(synonym)
                
        return synonyms
        
    def _nontype_radius(self, unclustered_gids, type_gids, ani_af_nontype_vs_type):
        """Calculate circumscription radius for unclustered, nontype genomes."""
        
        # set type radius for all type genomes to default values
        nontype_radius = {}
        for gid in unclustered_gids:
            nontype_radius[gid] = GenomeRadius(ani = self.ani_sp, 
                                                     af = None,
                                                     neighbour_gid = None)

        # determine closest type ANI neighbour and restrict ANI radius as necessary
        ani_af = pickle.load(open(ani_af_nontype_vs_type, 'rb'))
        for nontype_gid in unclustered_gids:
            if nontype_gid not in ani_af:
                continue
                    
            for type_gid in type_gids:
                if type_gid not in ani_af[nontype_gid]:
                    continue
                    
                ani, af = symmetric_ani(ani_af, nontype_gid, type_gid)

                if ani > nontype_radius[nontype_gid].ani and af >= self.af_sp:
                    nontype_radius[nontype_gid] = GenomeRadius(ani = ani, 
                                                                 af = af,
                                                                 neighbour_gid = type_gid)
                    
        self.logger.info('ANI circumscription radius: min=%.2f, mean=%.2f, max=%.2f' % (
                                min([d.ani for d in nontype_radius.values()]), 
                                np_mean([d.ani for d in nontype_radius.values()]), 
                                max([d.ani for d in nontype_radius.values()])))
                        
        return nontype_radius
        
    def _mash_ani_unclustered(self, genome_files, gids):
        """Calculate pairwise Mash ANI estimates between genomes."""
        
        mash = Mash(self.cpus)
        
        # create Mash sketch for potential representative genomes
        mash_nontype_sketch_file = os.path.join(self.output_dir, 'gtdb_unclustered_genomes.msh')
        genome_list_file = os.path.join(self.output_dir, 'gtdb_unclustered_genomes.lst')
        mash.sketch(gids, genome_files, genome_list_file, mash_nontype_sketch_file)

        # get Mash distances
        mash_dist_file = os.path.join(self.output_dir, 'gtdb_unclustered_genomes.dst')
        mash.dist_pairwise( float(100 - self.min_mash_ani)/100, mash_nontype_sketch_file, mash_dist_file)

        # read Mash distances
        mash_ani = mash.read_ani(mash_dist_file)
        
        # report pairs above Mash threshold
        mash_ani_pairs = []
        for qid in mash_ani:
            for rid in mash_ani[qid]:
                if mash_ani[qid][rid] >= self.min_mash_ani:
                    if qid != rid:
                        mash_ani_pairs.append((qid, rid))
                        mash_ani_pairs.append((rid, qid))
                
        self.logger.info('Identified %d genome pairs with a Mash ANI >= %.1f%%.' % (len(mash_ani_pairs), self.min_mash_ani))

        return mash_ani
        
    def _selected_rep_genomes(self,
                                genome_files,
                                nontype_radius, 
                                unclustered_qc_gids, 
                                mash_ani,
                                quality_metadata,
                                rnd_type_genome):
        """Select representative genomes for species clusters in a  greedy fashion using species-specific ANI thresholds."""

        # sort genomes by quality score
        if rnd_type_genome:
            self.logger.info('Selecting random de novo type genomes.')
            sorted_gids = []
            for gid in random.sample(unclustered_qc_gids, len(unclustered_qc_gids)):
                sorted_gids.append((gid, 0))
        else:
            self.logger.info('Selecting de novo type genomes in a greedy manner based on quality.')
            qscore = quality_score(unclustered_qc_gids, quality_metadata)
            sorted_gids = sorted(qscore.items(), key=operator.itemgetter(1), reverse=True)

        # greedily determine representatives for new species clusters
        cluster_rep_file = os.path.join(self.output_dir, 'cluster_reps.tsv')
        clusters = set()
        if not os.path.exists(cluster_rep_file):
            self.logger.info('Clustering genomes to identify representatives.')
            clustered_genomes = 0
            max_ani_pairs = 0
            for idx, (cur_gid, _score) in enumerate(sorted_gids):

                # determine reference genomes to calculate ANI between
                ani_pairs = []
                if cur_gid in mash_ani:
                    for rep_gid in clusters:
                        if mash_ani[cur_gid].get(rep_gid, 0) >= self.min_mash_ani:
                            ani_pairs.append((cur_gid, rep_gid))
                            ani_pairs.append((rep_gid, cur_gid))

                # determine if genome clusters with representative
                clustered = False
                if ani_pairs:
                    if len(ani_pairs) > max_ani_pairs:
                        max_ani_pairs = len(ani_pairs)
                    
                    ani_af = self.fastani.pairs(ani_pairs, genome_files, report_progress=False)

                    closest_rep_gid = None
                    closest_rep_ani = 0
                    closest_rep_af = 0
                    for rep_gid in clusters:
                        ani, af = symmetric_ani(ani_af, cur_gid, rep_gid)

                        if af >= self.af_sp:
                            if ani > closest_rep_ani or (ani == closest_rep_ani and af > closest_rep_af):
                                closest_rep_gid = rep_gid
                                closest_rep_ani = ani
                                closest_rep_af = af

                        if ani > nontype_radius[cur_gid].ani and af >= self.af_sp:
                            nontype_radius[cur_gid] = GenomeRadius(ani = ani, 
                                                                         af = af,
                                                                         neighbour_gid = rep_gid)
                                                                         
                    if closest_rep_gid and closest_rep_ani > nontype_radius[closest_rep_gid].ani:
                        clustered = True
                    
                if not clustered:
                    # genome is a new species cluster representative
                    clusters.add(cur_gid)
                else:
                    clustered_genomes += 1
                
                if (idx+1) % 10 == 0 or idx+1 == len(sorted_gids):
                    statusStr = '-> Clustered %d of %d (%.2f%%) genomes [ANI pairs: %d; clustered genomes: %d; clusters: %d].'.ljust(96) % (
                                    idx+1, 
                                    len(sorted_gids), 
                                    float(idx+1)*100/len(sorted_gids),
                                    max_ani_pairs,
                                    clustered_genomes,
                                    len(clusters))
                    sys.stdout.write('%s\r' % statusStr)
                    sys.stdout.flush()
                    max_ani_pairs = 0
            sys.stdout.write('\n')
            
            # write out selected cluster representative
            fout = open(cluster_rep_file, 'w')
            for gid in clusters:
                fout.write('%s\n' % gid)
            fout.close()
        else:
            # read cluster reps from file
            self.logger.warning('Using previously determined cluster representatives.')
            for line in open(cluster_rep_file):
                gid = line.strip()
                clusters.add(gid)
                
        self.logger.info('Selected %d representative genomes for de novo species clusters.' % len(clusters))
        
        return clusters
        
    def _cluster_genomes(self, 
                            genome_files,
                            rep_genomes,
                            type_gids, 
                            passed_qc,
                            final_cluster_radius):
        """Cluster all non-type/representative genomes to selected type/representatives genomes."""

        all_reps = rep_genomes.union(type_gids)
        
        # calculate MASH distance between non-type/representative genomes and selected type/representatives genomes
        mash = Mash(self.cpus)
        
        mash_type_rep_sketch_file = os.path.join(self.output_dir, 'gtdb_rep_genomes.msh')
        type_rep_genome_list_file = os.path.join(self.output_dir, 'gtdb_rep_genomes.lst')
        mash.sketch(all_reps, genome_files, type_rep_genome_list_file, mash_type_rep_sketch_file)
        
        mash_none_rep_sketch_file = os.path.join(self.output_dir, 'gtdb_nonrep_genomes.msh')
        type_none_rep_file = os.path.join(self.output_dir, 'gtdb_nonrep_genomes.lst')
        mash.sketch(passed_qc - all_reps, genome_files, type_none_rep_file, mash_none_rep_sketch_file)

        # get Mash distances
        mash_dist_file = os.path.join(self.output_dir, 'gtdb_rep_vs_nonrep_genomes.dst')
        mash.dist(float(100 - self.min_mash_ani)/100, mash_type_rep_sketch_file, mash_none_rep_sketch_file, mash_dist_file)

        # read Mash distances
        mash_ani = mash.read_ani(mash_dist_file)
        
        # calculate ANI between non-type/representative genomes and selected type/representatives genomes
        clusters = {}
        for gid in all_reps:
            clusters[gid] = []
        
        genomes_to_cluster = passed_qc - set(clusters)
        ani_pairs = []
        for gid in genomes_to_cluster:
            if gid in mash_ani:
                for rep_gid in clusters:
                    if mash_ani[gid].get(rep_gid, 0) >= self.min_mash_ani:
                        ani_pairs.append((gid, rep_gid))
                        ani_pairs.append((rep_gid, gid))
                        
        self.logger.info('Calculating ANI between %d species clusters and %d unclustered genomes (%d pairs):' % (
                            len(clusters), 
                            len(genomes_to_cluster),
                            len(ani_pairs)))
        ani_af = self.fastani.pairs(ani_pairs, genome_files)

        # assign genomes to closest representatives 
        # that is within the representatives ANI radius
        self.logger.info('Assigning genomes to closest representative.')
        for idx, cur_gid in enumerate(genomes_to_cluster):
            closest_rep_gid = None
            closest_rep_ani = 0
            closest_rep_af = 0
            for rep_gid in clusters:
                ani, af = symmetric_ani(ani_af, cur_gid, rep_gid)
                
                if ani >= final_cluster_radius[rep_gid].ani and af >= self.af_sp:
                    if ani > closest_rep_ani or (ani == closest_rep_ani and af > closest_rep_af):
                        closest_rep_gid = rep_gid
                        closest_rep_ani = ani
                        closest_rep_af = af
                
            if closest_rep_gid:
                clusters[closest_rep_gid].append(self.ClusteredGenome(gid=cur_gid, 
                                                                        ani=closest_rep_ani, 
                                                                        af=closest_rep_af))
            else:
                self.logger.warning('Failed to assign genome %s to representative.' % cur_gid)
                if closest_rep_gid:
                    self.logger.warning(' ...closest_rep_gid = %s' % closest_rep_gid)
                    self.logger.warning(' ...closest_rep_ani = %.2f' % closest_rep_ani)
                    self.logger.warning(' ...closest_rep_af = %.2f' % closest_rep_af)
                    self.logger.warning(' ...closest rep radius = %.2f' % final_cluster_radius[closest_rep_gid].ani)
                else:
                    self.logger.warning(' ...no representative with an AF >%.2f identified.' % self.af_sp)
             
            statusStr = '-> Assigned %d of %d (%.2f%%) genomes.'.ljust(86) % (idx+1, 
                                                                                len(genomes_to_cluster), 
                                                                                float(idx+1)*100/len(genomes_to_cluster))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
        sys.stdout.write('\n')

        return clusters, ani_af
        
    def _assign_species_names(self, clusters, names_in_use, gtdb_taxonomy, gtdb_user_to_genbank):
        """Assign a species name to each species cluster."""
        
        orig_names_in_use = set(names_in_use)

        fout = open(os.path.join(self.output_dir, 'gtdb_assigned_sp.tsv'), 'w')
        fout.write('Representative genome\tAssigned species\tGTDB taxonomy\tNo. clustered genomes\tClustered GTDB genera\tClustered GTDB species\tSpecies name in use\tMost common name in use\tClustered genomes\n')
        cluster_sp_names = {}
        for rid in sorted(clusters, key=lambda x: len(clusters[x]), reverse=True):
            clustered_gids = [c.gid for c in clusters[rid]]
            
            # find most common genus name in cluster
            gtdb_genera = [gtdb_taxonomy[gid][5] for gid in clustered_gids] + [gtdb_taxonomy[rid][5]]
            gtdb_genus_counter = Counter(gtdb_genera)
            gtdb_common_genus = None 
            gtdb_common_genus_count = 0
            for genus, count in gtdb_genus_counter.most_common(): 
                if genus != 'g__':
                    gtdb_common_genus = genus
                    gtdb_common_genus_count = count
                    break
                    
            # in case of ties involving genus of representative genome, 
            # defer to classification of representative
            rep_genus = gtdb_taxonomy[rid][5]
            if gtdb_genus_counter[rep_genus] == gtdb_common_genus_count and rep_genus != 'g__':
                gtdb_common_genus = rep_genus
            
            # get most common GTDB species name 
            gtdb_sp = [gtdb_taxonomy[gid][6] for gid in clustered_gids] + [gtdb_taxonomy[rid][6]]
            gtdb_sp_counter = Counter(gtdb_sp)
            gtdb_common_sp = None
            gtdb_common_sp_count = 0
            for sp, count in gtdb_sp_counter.most_common(): 
                if sp != 's__':
                    gtdb_common_sp = sp
                    gtdb_common_sp_count = count
                    break
                    
            most_common_in_use = gtdb_common_sp in names_in_use

            min_req_genomes = 0.5*(sum(gtdb_sp_counter.values()) - gtdb_sp_counter.get('s__', 0))
            if gtdb_common_sp_count >= min_req_genomes and not most_common_in_use:
                # assign common species if it occurs in >=50% of the clustered genomes,
                # excluding genomes with no species assignment
                names_in_use.add(gtdb_common_sp)
                cluster_sp_names[rid] = gtdb_common_sp
            else:
                # derive new species name from genus, if possible, 
                # and accession number of representative genome
                genus = '{unresolved}'
                if gtdb_common_genus and gtdb_common_genus != 'g__':
                    genus = gtdb_common_genus[3:]
                
                acc = rid
                if rid.startswith('U_'):
                    if rid in gtdb_user_to_genbank:
                        acc = gtdb_user_to_genbank[rid]
                    else:
                        # create accession from GTDB User ID of the form:
                        # U_<number>u.0 which will give 'sp<number>u'
                        acc = 'U_' + rid.replace('U_', '') + 'u.0'

                derived_sp = 's__' + '%s sp%s' % (genus, acc[acc.rfind('_')+1:acc.rfind('.')])
                if derived_sp in names_in_use:
                    self.logger.error('Derived species name already in use: %s, %s' % (derived_sp, acc))
                    sys.exit(-1)

                names_in_use.add(derived_sp)
                cluster_sp_names[rid] = derived_sp
                
            fout.write('%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n' % (
                        rid, 
                        cluster_sp_names[rid],
                        '; '.join(gtdb_taxonomy[rid]),
                        len(clustered_gids),
                        ', '.join("%s=%r" % (genus, count) for (genus, count) in gtdb_genus_counter.most_common()),
                        ', '.join("%s=%r" % (sp, count) for (sp, count) in gtdb_sp_counter.most_common()),
                        ', '.join("%s=%s" % (sp, sp in names_in_use) for sp, _count in gtdb_sp_counter.most_common()),
                        '%s=%d' % (gtdb_common_sp, gtdb_common_sp_count) if most_common_in_use else 'n/a',
                        ', '.join(clustered_gids)))
                
        fout.close()
        
        return cluster_sp_names
        
    def _write_rep_info(self, 
                        clusters, 
                        cluster_sp_names, 
                        quality_metadata, 
                        genome_quality,
                        excluded_from_refseq_note,
                        ani_af,
                        output_file):
        """Write out information about selected representative genomes."""
                                            
        fout = open(output_file, 'w')
        fout.write('Species\tType genome\tNCBI assembly level\tNCBI genome category')
        fout.write('\tGenome size (bp)\tQuality score\tCompleteness (%)\tContamination (%)\tNo. scaffolds\tNo. contigs\tN50 contigs\tAmbiguous bases\tSSU count\tSSU length (bp)')
        fout.write('\tNo. genomes in cluster\tMean ANI\tMean AF\tMin ANI\tMin AF\tNCBI exclude from RefSeq\n')
        
        for gid in clusters:
            fout.write('%s\t%s\t%s\t%s' % (
                        cluster_sp_names[gid], 
                        gid, 
                        quality_metadata[gid].ncbi_assembly_level,
                        quality_metadata[gid].ncbi_genome_category))

            fout.write('\t%d\t%.2f\t%.2f\t%.2f\t%d\t%d\t%.1f\t%d\t%d\t%d' % (
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
                            
            anis = []
            afs = []
            for cluster_id in clusters[gid]:
                ani, af = symmetric_ani(ani_af, gid, cluster_id)
                anis.append(ani)
                afs.append(af)
            
            if anis:
                fout.write('\t%d\t%.1f\t%.2f\t%.1f\t%.2f\t%s\n' % (len(clusters[gid]),
                                                                    np_mean(anis), np_mean(afs),
                                                                    min(anis), min(afs),
                                                                    excluded_from_refseq_note.get(gid, '')))
            else:
                fout.write('\t%d\t%s\t%s\t%s\t%s\t%s\n' % (len(clusters[gid]),
                                                            'n/a', 'n/a', 'n/a', 'n/a',
                                                            excluded_from_refseq_note.get(gid, '')))
        fout.close()
        
    def _gtdb_user_genomes(self, gtdb_user_genomes_file, metadata_file):
        """Get map between GTDB User genomes and GenBank accessions."""
        
        uba_to_genbank = {}
        for line in open(gtdb_user_genomes_file):
            line_split = line.strip().split('\t')
            gb_acc = line_split[0]
            uba_id = line_split[4]
            uba_to_genbank[uba_id] = gb_acc
        
        user_to_genbank = {}
        m = read_gtdb_metadata(metadata_file, ['organism_name'])
        for gid, metadata in m.items():
            if '(UBA' in str(metadata.organism_name):
                uba_id = metadata.organism_name[metadata.organism_name.find('(')+1:-1]
                if uba_id in uba_to_genbank:
                    user_to_genbank[gid] = uba_to_genbank[uba_id]

        return user_to_genbank

    def run(self, qc_file,
                metadata_file,
                gtdb_user_genomes_file,
                genome_path_file,
                type_genome_cluster_file,
                type_genome_synonym_file,
                ncbi_refseq_assembly_file,
                ncbi_genbank_assembly_file,
                ani_af_nontype_vs_type,
                species_exception_file,
                rnd_type_genome):
        """Infer de novo species clusters and type genomes for remaining genomes."""
        
        # identify genomes failing quality criteria
        self.logger.info('Reading QC file.')
        passed_qc = read_qc_file(qc_file)
        self.logger.info('Identified %d genomes passing QC.' % len(passed_qc))
        
        # get NCBI taxonomy strings for each genome
        self.logger.info('Reading NCBI taxonomy from GTDB metadata file.')
        ncbi_taxonomy, ncbi_update_count = read_gtdb_ncbi_taxonomy(metadata_file, species_exception_file)
        gtdb_taxonomy = read_gtdb_taxonomy(metadata_file)
        self.logger.info('Read NCBI taxonomy for %d genomes with %d manually defined updates.' % (len(ncbi_taxonomy), ncbi_update_count))
        self.logger.info('Read GTDB taxonomy for %d genomes.' % len(gtdb_taxonomy))
        
        # parse NCBI assembly files
        self.logger.info('Parsing NCBI assembly files.')
        excluded_from_refseq_note = exclude_from_refseq(ncbi_refseq_assembly_file, ncbi_genbank_assembly_file)

        # get path to genome FASTA files
        self.logger.info('Reading path to genome FASTA files.')
        genome_files = read_genome_path(genome_path_file)
        self.logger.info('Read path for %d genomes.' % len(genome_files))
        for gid in set(genome_files):
            if gid not in passed_qc:
                genome_files.pop(gid)
        self.logger.info('Considering %d genomes as potential representatives after removing unwanted User genomes.' % len(genome_files))
        assert(len(genome_files) == len(passed_qc))
        
        # determine type genomes and genomes clustered to type genomes
        type_species, species_type_gid, type_gids, type_clustered_gids, type_radius = self._parse_type_clusters(type_genome_cluster_file)
        assert(len(type_species) == len(type_gids))
        self.logger.info('Identified %d type genomes.' % len(type_gids))
        self.logger.info('Identified %d clustered genomes.' % len(type_clustered_gids))
        
        # calculate quality score for genomes
        self.logger.info('Parse quality statistics for all genomes.')
        quality_metadata = read_quality_metadata(metadata_file)
        
        # calculate genome quality score
        self.logger.info('Calculating genome quality score.')
        genome_quality = quality_score(quality_metadata.keys(), quality_metadata)

        # determine genomes left to be clustered
        unclustered_gids = passed_qc - type_gids - type_clustered_gids
        self.logger.info('Identified %d unclustered genomes passing QC.' % len(unclustered_gids))

        # establish closest type genome for each unclustered genome
        self.logger.info('Determining ANI circumscription for %d unclustered genomes.' % len(unclustered_gids))
        nontype_radius = self._nontype_radius(unclustered_gids, type_gids, ani_af_nontype_vs_type)
        
        # calculate Mash ANI estimates between unclustered genomes
        self.logger.info('Calculating Mash ANI estimates between unclustered genomes.')
        mash_anis = self._mash_ani_unclustered(genome_files, unclustered_gids)

        # select species representatives genomes in a greedy fashion based on genome quality
        rep_genomes = self._selected_rep_genomes(genome_files,
                                                    nontype_radius, 
                                                    unclustered_gids, 
                                                    mash_anis,
                                                    quality_metadata,
                                                    rnd_type_genome)
        
        # cluster all non-type/non-rep genomes to species type/rep genomes
        final_cluster_radius = type_radius.copy()
        final_cluster_radius.update(nontype_radius)
        
        final_clusters, ani_af = self._cluster_genomes(genome_files,
                                                        rep_genomes,
                                                        type_gids, 
                                                        passed_qc,
                                                        final_cluster_radius)
        rep_clusters = {}
        for gid in rep_genomes:
            rep_clusters[gid] = final_clusters[gid]

        # get list of synonyms in order to restrict usage of species names
        synonyms = self._parse_synonyms(type_genome_synonym_file)
        self.logger.info('Identified %d synonyms.' % len(synonyms))
        
        # determine User genomes with NCBI accession number that may form species names
        gtdb_user_to_genbank = self._gtdb_user_genomes(gtdb_user_genomes_file, metadata_file)
        self.logger.info('Identified %d GTDB User genomes with NCBI accessions.' % len(gtdb_user_to_genbank))
        
        # assign species names to de novo species clusters
        names_in_use = synonyms.union(type_species)
        self.logger.info('Identified %d species names already in use.' % len(names_in_use))
        self.logger.info('Assigning species name to each de novo species cluster.')
        cluster_sp_names = self._assign_species_names(rep_clusters, 
                                                        names_in_use, 
                                                        gtdb_taxonomy,
                                                        gtdb_user_to_genbank)
        
         # write out file with details about selected representative genomes
        self._write_rep_info(rep_clusters, 
                                cluster_sp_names,
                                quality_metadata,
                                genome_quality,
                                excluded_from_refseq_note,
                                ani_af,
                                os.path.join(self.output_dir, 'gtdb_rep_genome_info.tsv'))
                                             
        # remove genomes that are not representatives of a species cluster and then write out representative ANI radius
        for gid in set(final_cluster_radius) - set(final_clusters):
            del final_cluster_radius[gid]
            
        all_species = cluster_sp_names
        all_species.update(species_type_gid)

        self.logger.info('Writing %d species clusters to file.' % len(all_species))
        self.logger.info('Writing %d cluster radius information to file.' % len(final_cluster_radius))
        
        write_clusters(final_clusters, 
                        final_cluster_radius, 
                        all_species, 
                        os.path.join(self.output_dir, 'gtdb_clusters_final.tsv'))

        write_rep_radius(final_cluster_radius, 
                            all_species, 
                            os.path.join(self.output_dir, 'gtdb_ani_radius_final.tsv'))
        