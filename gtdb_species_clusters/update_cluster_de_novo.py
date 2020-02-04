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

from gtdb_species_clusters.mash import Mash
from gtdb_species_clusters.fastani import FastANI

from gtdb_species_clusters.genome import Genome
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_clusters import SpeciesClusters
                                    
from gtdb_species_clusters.genome_utils import (canonical_gid,
                                                read_qc_file,
                                                exclude_from_refseq)

from gtdb_species_clusters.type_genome_utils import (GenomeRadius,
                                                        symmetric_ani,
                                                        write_rep_radius,
                                                        write_clusters)

class UpdateClusterDeNovo(object):
    """Infer de novo species clusters and representatives for remaining genomes."""

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
        
    def _parse_named_clusters(self, named_cluster_file):
        """Parse named GTDB species clusters."""
        
        rep_species = set()
        species_rep_gid = {}
        rep_gids = set()
        rep_clustered_gids = set()
        rep_radius = {}
        with open(named_cluster_file) as f:
            headers = f.readline().strip().split('\t')
            
            rep_index = headers.index('Representative')
            sp_index = headers.index('Proposed species')
            num_clustered_index = headers.index('No. clustered genomes')
            clustered_genomes_index = headers.index('Clustered genomes')
            closest_type_index = headers.index('Closest representative')
            ani_radius_index = headers.index('ANI radius')
            af_index = headers.index('AF closest')

            for line in f:
                line_split = line.strip().split('\t')
                
                sp = line_split[sp_index]
                if sp in rep_species:
                    self.logger.error('Proposed species represented by multiple genomes: {}'.format(sp))
                    sys.exit(-1)
                rep_species.add(sp)
                
                rep_gid = line_split[rep_index]
                rep_gids.add(rep_gid)
                
                species_rep_gid[rep_gid] = sp
                
                num_clustered = int(line_split[num_clustered_index])
                if num_clustered > 0:
                    for gid in [g.strip() for g in line_split[clustered_genomes_index].split(',')]:
                        rep_clustered_gids.add(gid)
                        
                rep_radius[rep_gid] = GenomeRadius(ani = float(line_split[ani_radius_index]), 
                                                     af = float(line_split[af_index]),
                                                     neighbour_gid = line_split[closest_type_index])
                        
        return rep_species, species_rep_gid, rep_gids, rep_clustered_gids, rep_radius
        
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
        
    def _nonrep_radius(self, unclustered_gids, rep_gids, ani_af_rep_vs_nonrep):
        """Calculate circumscription radius for unclustered, nontype genomes."""
        
        # set radius for genomes to default values
        nonrep_radius = {}
        for gid in unclustered_gids:
            nonrep_radius[gid] = GenomeRadius(ani = self.ani_sp, 
                                                     af = None,
                                                     neighbour_gid = None)

        # determine closest type ANI neighbour and restrict ANI radius as necessary
        ani_af = pickle.load(open(ani_af_rep_vs_nonrep, 'rb'))
        for nonrep_gid in unclustered_gids:
            if nonrep_gid not in ani_af:
                continue
                    
            for rep_gid in rep_gids:
                if rep_gid not in ani_af[nonrep_gid]:
                    continue
                    
                ani, af = symmetric_ani(ani_af, nonrep_gid, rep_gid)

                if ani > nonrep_radius[nonrep_gid].ani and af >= self.af_sp:
                    nonrep_radius[nonrep_gid] = GenomeRadius(ani = ani, 
                                                             af = af,
                                                             neighbour_gid = rep_gid)
                    
        self.logger.info('ANI circumscription radius: min={:.2f}, mean={:.2f}, max={:.2f}'.format(
                                min([d.ani for d in nonrep_radius.values()]), 
                                np_mean([d.ani for d in nonrep_radius.values()]), 
                                max([d.ani for d in nonrep_radius.values()])))
                        
        return nonrep_radius
        
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
                                cur_genomes,
                                genome_files,
                                nonrep_radius, 
                                unclustered_qc_gids, 
                                mash_ani):
        """Select representative genomes for species clusters in a  greedy fashion using species-specific ANI thresholds."""

        # sort genomes by quality score
        self.logger.info('Selecting de novo representatives in a greedy manner based on quality.')
        q = {gid:cur_genomes[gid].score_type_strain() for gid in unclustered_qc_gids}
        q_sorted = sorted(q.items(), key=operator.itemgetter(1), reverse=True)

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

                        if ani > nonrep_radius[cur_gid].ani and af >= self.af_sp:
                            nonrep_radius[cur_gid] = GenomeRadius(ani = ani, 
                                                                         af = af,
                                                                         neighbour_gid = rep_gid)
                                                                         
                    if closest_rep_gid and closest_rep_ani > nonrep_radius[closest_rep_gid].ani:
                        clustered = True
                    
                if not clustered:
                    # genome is a new species cluster representative
                    clusters.add(cur_gid)
                else:
                    clustered_genomes += 1
                
                if (idx+1) % 10 == 0 or idx+1 == len(sorted_gids):
                    statusStr = '-> Clustered {:,} of {:,} ({:.2f}%) genomes [ANI pairs: {:,}; clustered genomes: {:,}; clusters: {:,}].'.format(
                                    idx+1, 
                                    len(sorted_gids), 
                                    float(idx+1)*100/len(sorted_gids),
                                    max_ani_pairs,
                                    clustered_genomes,
                                    len(clusters)).ljust(96)
                    sys.stdout.write('{}\r'.format(statusStr))
                    sys.stdout.flush()
                    max_ani_pairs = 0
            sys.stdout.write('\n')
            
            # write out selected cluster representative
            fout = open(cluster_rep_file, 'w')
            for gid in clusters:
                fout.write('{}\n'.format(gid))
            fout.close()
        else:
            # read cluster reps from file
            self.logger.warning('Using previously determined cluster representatives.')
            for line in open(cluster_rep_file):
                gid = line.strip()
                clusters.add(gid)
                
        self.logger.info('Selected {:,} representative genomes for de novo species clusters.'.format(len(clusters)))
        
        return clusters
        
    def _cluster_genomes(self,
                            cur_genomes,
                            genome_files,
                            rep_genomes,
                            named_rep_gids, 
                            passed_qc,
                            final_cluster_radius):
        """Cluster new representatives to representatives of named GTDB species clusters."""

        all_reps = rep_genomes.union(named_rep_gids)
        
        # calculate MASH distance between non-type/representative genomes and selected type/representatives genomes
        mash = Mash(self.cpus)
        
        mash_rep_sketch_file = os.path.join(self.output_dir, 'gtdb_rep_genomes.msh')
        rep_genome_list_file = os.path.join(self.output_dir, 'gtdb_rep_genomes.lst')
        mash.sketch(all_reps, genome_files, rep_genome_list_file, mash_rep_sketch_file)
        
        mash_none_rep_sketch_file = os.path.join(self.output_dir, 'gtdb_nonrep_genomes.msh')
        non_rep_file = os.path.join(self.output_dir, 'gtdb_nonrep_genomes.lst')
        mash.sketch(cur_genomes.genomes.keys() - all_reps, genome_files, non_rep_file, mash_none_rep_sketch_file)

        # get Mash distances
        mash_dist_file = os.path.join(self.output_dir, 'gtdb_rep_vs_nonrep_genomes.dst')
        mash.dist(float(100 - self.min_mash_ani)/100, mash_rep_sketch_file, mash_none_rep_sketch_file, mash_dist_file)

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
                        
        self.logger.info('Calculating ANI between {:,} species clusters and {:,} unclustered genomes ({:,} pairs):'.format(
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
                self.logger.warning('Failed to assign genome {} to representative.'.format(cur_gid))
                if closest_rep_gid:
                    self.logger.warning(' ...closest_rep_gid = {}'.format(closest_rep_gid))
                    self.logger.warning(' ...closest_rep_ani = {:.2f}'.format(closest_rep_ani))
                    self.logger.warning(' ...closest_rep_af = {:.2f}'.format(closest_rep_af))
                    self.logger.warning(' ...closest rep radius = {:.2f}'.format(final_cluster_radius[closest_rep_gid].ani))
                else:
                    self.logger.warning(' ...no representative with an AF >{:.2f} identified.'.format(self.af_sp))
             
            statusStr = '-> Assigned {:,} of {:,} ({:.2f%}) genomes.'.format(idx+1, 
                                                                                len(genomes_to_cluster), 
                                                                                float(idx+1)*100/len(genomes_to_cluster)).ljust(86)
            sys.stdout.write('{}\r'.format(statusStr))
            sys.stdout.flush()
        sys.stdout.write('\n')

        return clusters, ani_af
        
    def _assign_species_names(self, cur_genomes, clusters, names_in_use):
        """Assign a species name to each species cluster."""
        
        orig_names_in_use = set(names_in_use)

        fout = open(os.path.join(self.output_dir, 'gtdb_assigned_sp.tsv'), 'w')
        fout.write('Representative genome\tAssigned species\tGTDB taxonomy\tNo. clustered genomes\tClustered GTDB genera\tClustered GTDB species\tSpecies name in use\tMost common name in use\tClustered genomes\n')
        cluster_sp_names = {}
        for rid in sorted(clusters, key=lambda x: len(clusters[x]), reverse=True):
            clustered_gids = [c.gid for c in clusters[rid]]
            
            # find most common genus name in cluster
            gtdb_genera = [cur_genomes[gid].gtdb_taxa.genus for gid in clustered_gids] + [cur_genomes[rid].gtdb_taxa.genus]
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
            rep_genus = cur_genomes[rid].gtdb_taxa.genus
            if gtdb_genus_counter[rep_genus] == gtdb_common_genus_count and rep_genus != 'g__':
                gtdb_common_genus = rep_genus
            
            # get most common GTDB species name 
            gtdb_sp = [cur_genomes[gid].gtdb_taxa.species for gid in clustered_gids] + [cur_genomes[rid].gtdb_taxa.species]
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
                    self.logger.error('User genome was selected as the representative for a new species cluster.')
                    sys.exit(-1)

                derived_sp = 's__' + '{} sp{}'.format(genus, acc[acc.rfind('_')+1:acc.rfind('.')])
                if derived_sp in names_in_use:
                    self.logger.error('Derived species name already in use: {}, {}'.format(derived_sp, acc))
                    sys.exit(-1)

                names_in_use.add(derived_sp)
                cluster_sp_names[rid] = derived_sp
                
            fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        rid, 
                        cluster_sp_names[rid],
                        cur_genomes[rid].gtdb_taxa,
                        len(clustered_gids),
                        ', '.join("%s=%r" % (genus, count) for (genus, count) in gtdb_genus_counter.most_common()),
                        ', '.join("%s=%r" % (sp, count) for (sp, count) in gtdb_sp_counter.most_common()),
                        ', '.join("%s=%s" % (sp, sp in names_in_use) for sp, _count in gtdb_sp_counter.most_common()),
                        '%s=%d' % (gtdb_common_sp, gtdb_common_sp_count) if most_common_in_use else 'n/a',
                        ', '.join(clustered_gids)))
                
        fout.close()
        
        return cluster_sp_names
        
    def _write_rep_info(self, 
                        cur_genomes,
                        clusters, 
                        cluster_sp_names, 
                        excluded_from_refseq_note,
                        ani_af,
                        output_file):
        """Write out information about selected representative genomes."""
                                            
        fout = open(output_file, 'w')
        fout.write('Species\tRepresentative\tNCBI assembly level\tNCBI genome category')
        fout.write('\tGenome size (bp)\tQuality score\tCompleteness (%)\tContamination (%)\tNo. scaffolds\tNo. contigs\tN50 contigs\tAmbiguous bases\tSSU count\tSSU length (bp)')
        fout.write('\tNo. genomes in cluster\tMean ANI\tMean AF\tMin ANI\tMin AF\tNCBI exclude from RefSeq\n')
        
        for gid in clusters:
            fout.write('{}\t{}\t{}\t{}'.format(
                        cluster_sp_names[gid], 
                        gid, 
                        cur_genomes[gid].ncbi_assembly_level,
                        cur_genomes[gid].ncbi_genome_category))

            fout.write('\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{}\t{:.1f}\t{}\t{}\t{}'.format(
                            cur_genomes[gid].length,
                            cur_genomes[gid].score_assembly(), 
                            cur_genomes[gid].comp,
                            cur_genomes[gid].cont,
                            cur_genomes[gid].scaffold_count,
                            cur_genomes[gid].contig_count,
                            cur_genomes[gid].contig_n50,
                            cur_genomes[gid].ambiguous_bases,
                            cur_genomes[gid].ssu_count,
                            cur_genomes[gid].ssu_length))
                            
            anis = []
            afs = []
            for cluster_id in clusters[gid]:
                ani, af = symmetric_ani(ani_af, gid, cluster_id)
                anis.append(ani)
                afs.append(af)
            
            if anis:
                fout.write('\t{}\t{:.1f}\t{:.2f}\t{:.1f}\t{:.2f}\t{}\n'.format(len(clusters[gid]),
                                                                    np_mean(anis), np_mean(afs),
                                                                    min(anis), min(afs),
                                                                    excluded_from_refseq_note.get(gid, '')))
            else:
                fout.write('\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(len(clusters[gid]),
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

    def run(self, qc_passed_file,
                    cur_gtdb_metadata_file,
                    cur_genomic_path_file,
                    uba_genome_paths,
                    named_cluster_file,
                    synonym_file,
                    ncbi_genbank_assembly_file,
                    ani_af_rep_vs_nonrep,
                    species_exception_file,
                    genus_exception_file,
                    gtdb_type_strains_ledger):
        """Infer de novo species clusters and representatives for remaining genomes."""
        
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

        # determine type genomes and genomes clustered to type genomes
        rep_species, species_rep_gid, named_rep_gids, rep_clustered_gids, rep_radius = self._parse_named_clusters(named_cluster_file)
        assert(len(rep_species) == len(named_rep_gids))
        self.logger.info('Identified {:,} type genomes.'.format(len(named_rep_gids)))
        self.logger.info('Identified {:,} clustered genomes.'.format(len(rep_clustered_gids)))
        
        # determine genomes left to be clustered
        unclustered_gids = cur_genomes.genomes.keys() - named_rep_gids - rep_clustered_gids
        self.logger.info('Identified {:,} unclustered genomes passing QC.'.format(len(unclustered_gids)))

        # establish closest representative for each unclustered genome
        self.logger.info('Determining ANI circumscription for {:,} unclustered genomes.' % len(unclustered_gids))
        nonrep_radius = self._nonrep_radius(unclustered_gids, named_rep_gids, ani_af_rep_vs_nonrep)

        # calculate Mash ANI estimates between unclustered genomes
        self.logger.info('Calculating Mash ANI estimates between unclustered genomes.')
        mash_anis = self._mash_ani_unclustered(cur_genomes.genomic_files, unclustered_gids)

        # select species representatives genomes in a greedy fashion based on genome quality
        rep_genomes = self._selected_rep_genomes(cur_genomes,
                                                    cur_genomes.genomic_files,
                                                    nonrep_radius, 
                                                    unclustered_gids, 
                                                    mash_anis)

        # cluster all non-type/non-rep genomes to species type/rep genomes
        final_cluster_radius = rep_radius.copy()
        final_cluster_radius.update(nonrep_radius)
        
        final_clusters, ani_af = self._cluster_genomes(cur_genomes,
                                                        cur_genomes.genomic_files,
                                                        rep_genomes,
                                                        named_rep_gids, 
                                                        passed_qc,
                                                        final_cluster_radius)
        rep_clusters = {}
        for gid in rep_genomes:
            rep_clusters[gid] = final_clusters[gid]

        # get list of synonyms in order to restrict usage of species names
        synonyms = self._parse_synonyms(synonym_file)
        self.logger.info('Identified {:,} synonyms.'.format(len(synonyms)))
        
        # assign species names to de novo species clusters
        names_in_use = synonyms.union(rep_species)
        self.logger.info('Identified {:,} species names already in use.'.format(len(names_in_use)))
        self.logger.info('Assigning species name to each de novo species cluster.')
        cluster_sp_names = self._assign_species_names(cur_genomes,
                                                        rep_clusters, 
                                                        names_in_use)
                                                        
        # parse NCBI assembly files
        self.logger.info('Parsing NCBI assembly files.')
        excluded_from_refseq_note = exclude_from_refseq(ncbi_genbank_assembly_file)
        
         # write out file with details about selected representative genomes
        self._write_rep_info(cur_genomes,
                                rep_clusters, 
                                cluster_sp_names,
                                excluded_from_refseq_note,
                                ani_af,
                                os.path.join(self.output_dir, 'gtdb_rep_genome_info.tsv'))
                                             
        # remove genomes that are not representatives of a species cluster and then write out representative ANI radius
        for gid in set(final_cluster_radius) - set(final_clusters):
            del final_cluster_radius[gid]
            
        all_species = cluster_sp_names
        all_species.update(species_rep_gid)

        self.logger.info('Writing {:,} species clusters to file.'.format(len(all_species)))
        self.logger.info('Writing {:,} cluster radius information to file.'.format(len(final_cluster_radius)))
        
        write_clusters(final_clusters, 
                        final_cluster_radius, 
                        all_species,
                        cur_genomes,
                        os.path.join(self.output_dir, 'gtdb_clusters_final.tsv'))

        write_rep_radius(final_cluster_radius, 
                            all_species, 
                            cur_genomes,
                            os.path.join(self.output_dir, 'gtdb_ani_radius_final.tsv'))

        