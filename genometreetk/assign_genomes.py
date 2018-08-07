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
import uuid
import time
import tempfile
import argparse
import operator
import ntpath
import logging
import multiprocessing as mp
from collections import defaultdict

from biolib.taxonomy import Taxonomy
from biolib.external.execute import check_dependencies

csv.field_size_limit(sys.maxsize)


class AssignGenomes(object):
    """Assign genomes to canonical genomes comprising GTDB reference tree."""

    def __init__(self, cpus, output_dir):
        """Initialization."""
        
        check_dependencies(['mash', 'fastANI'])
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        self.output_dir = output_dir
        self.cpus = cpus

        self.logger = logging.getLogger('timestamp')
        
        # genome assignment parameters
        self.ani_threshold = 95             # assign genomes above this ANI value
        self.af_threshold = 0.65            # assign genomes above this alignment fraction AF

        self.mash_ani_threshold = 96.5      # assign genomes with Mash above this ANI threshold
        
    def _genomes_to_process(self, full_gtdb_taxonomy, metadata_file, user_genomes):
        """Read GTDB metadata to determine genomes to process."""
        
        genomes_to_process = set()
        num_uba_genomes = 0
        with open(metadata_file) as f:
            header = f.readline().strip().split('\t')
            org_name_index = header.index('organism_name')

            for line in f:
                line_split = line.strip().split('\t')
                
                gid = line_split[0]
                if not gid in full_gtdb_taxonomy:
                    continue
                    
                org_name = line_split[org_name_index]
                
                if gid.startswith('RS_') or gid.startswith('GB_'):
                    genomes_to_process.add(gid)
                elif user_genomes or '(UBA' in org_name:
                    genomes_to_process.add(gid)
                    num_uba_genomes += 1
                                                
        return genomes_to_process, num_uba_genomes
        
    def _get_genome_id(self, genome_path):
        """Extract genome ID from path to genomic file."""
        
        genome_id = ntpath.basename(genome_path)
        if genome_id.startswith('GCA_') or genome_id.startswith('GCF_'):
            genome_id = '_'.join(genome_id.split('_')[0:2])
            if genome_id.startswith('GCA_'):
                genome_id = 'GB_' + genome_id
            else:
                genome_id = 'RS_' + genome_id
        else:
            genome_id = '_'.join(genome_id.split('_')[0:2])
            
        return genome_id

    def _mash_ani(self, gids, rep_ids, genomic_files):
        """Calculate ANI between genomes and representatives using Mash."""
        
        # create Mash sketches for representatives
        tmp_ref_file = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
        fout = open(tmp_ref_file, 'w')
        for gid in rep_ids:
            fout.write('%s\n' % genomic_files[gid])
        fout.close()
        
        # calculate pairwise Mash distances
        tmp_ref_sketch_file = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()) + '.msh')
        cmd = 'mash sketch -l -p %d -k 16 -s 5000 -o %s %s 2> /dev/null' % (self.cpus,
                                                                            tmp_ref_sketch_file, 
                                                                            tmp_ref_file)
        os.system(cmd)
        
        # create Mash sketches for genomes
        tmp_genome_file = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
        fout = open(tmp_genome_file, 'w')
        for gid in gids:
            fout.write('%s\n' % genomic_files[gid])
        fout.close()
        
        # calculate pairwise Mash distances
        tmp_genome_sketch_file = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()) + '.msh')
        cmd = 'mash sketch -l -p %d -k 16 -s 5000 -o %s %s 2> /dev/null' % (self.cpus,
                                                                            tmp_genome_sketch_file, 
                                                                            tmp_genome_file)
        os.system(cmd)
        
        # calculate distances between references and genomes
        tmp_mash_file = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
        cmd = 'mash dist -p %d -d %f -v %f %s %s > %s 2> /dev/null' % (self.cpus,
                                                                        (1.0 - self.ani_threshold/100.0),
                                                                        1e-3,
                                                                        tmp_ref_sketch_file, 
                                                                        tmp_genome_sketch_file, 
                                                                        tmp_mash_file)
        os.system(cmd)
        
        # read results
        ani_values = defaultdict(lambda: {})
        for line in open(tmp_mash_file):
            line_split = line.strip().split('\t')
            
            ref_genome = self._get_genome_id(line_split[0])
            query_genome = self._get_genome_id(line_split[1])
            dist = float(line_split[2])
            ani_values[query_genome][ref_genome] = (100.0 - 100*dist)

        os.remove(tmp_ref_file)
        os.remove(tmp_genome_file)
        os.remove(tmp_ref_sketch_file)
        os.remove(tmp_genome_sketch_file)
        os.remove(tmp_mash_file)

        return ani_values
        
    def _mash_assignments(self, gids, rep_ids, genome_files):
        """Assign genomes to representatives using Mash distances."""

        ani_values = self._mash_ani(gids, rep_ids, genome_files)

        rep_assignments = defaultdict(set)
        for gid in gids:
            if gid in rep_ids:
                self.logger.error('Representative %d also found in genome set: %s' % gid)
                sys.exit(-1)
                
            if not gid in ani_values:
                continue
                
            # find highest Mash ANI value
            max_ani = 0
            max_rep_id = None
            for rep_id in rep_ids:
                ani = ani_values[gid].get(rep_id, 0.0)
                if ani >= self.mash_ani_threshold and ani > max_ani:
                    max_ani = ani
                    max_rep_id = rep_id
                    
            if max_rep_id:
                rep_assignments[max_rep_id].add((gid, 'Mash', max_ani, 0))
                    
        return rep_assignments
        
    def _fastani_assignments(self, gids, rep_ids, genome_files):
        """Assign genomes to representatives using FastANI distances."""

        ani_af = self._fastani_parallel(gids, rep_ids, genome_files)

        rep_assignments = defaultdict(set)
        for gid in gids:
            if gid in rep_ids:
                self.logger.error('Representative %d also found in genome set: %s' % gid)
                sys.exit(-1)
                
            # find highest Mash ANI value
            max_ani = 0
            max_rep_id = None
            for rep_id in rep_ids:
                try:
                    ani, af = ani_af[gid][rep_id]
                except:
                    ani, af = 0, 0
                    self.logger.warning('No ANI value calculated between %s and representative %s.' % (gid, rep_id))
                
                if ani >= self.ani_threshold and af >= self.af_threshold and ani > max_ani:
                    max_ani = ani
                    max_af = af
                    max_rep_id = rep_id
                    
            if max_rep_id:
                rep_assignments[max_rep_id].add((gid, 'FastANI', max_ani, max_af))
                    
        return rep_assignments

    def _fastani_parallel(self, gids, rep_ids, genome_files):
        """Calculate FastANI between genomes and representatives in parallel."""
        
        ani_af = mp.Manager().dict()
        
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()
        
        rep_ids = list(rep_ids)
        for query_gid in gids:
            # process representatives in batches of 100 to keep
            # memory requirements in check
            for start_pos in range(0, len(rep_ids), 100):
                end_pos = min(start_pos + 100, len(rep_ids))
                worker_queue.put((query_gid, rep_ids[start_pos:end_pos]))

        for _ in range(self.cpus):
            worker_queue.put((None, None))

        try:
            workerProc = [mp.Process(target = self.__fastani_worker, args = (genome_files,
                                                                                worker_queue, 
                                                                                writer_queue)) for _ in range(self.cpus)]
            writeProc = mp.Process(target = self.__fastani_writer, args = (ani_af, writer_queue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writer_queue.put(None)
            writeProc.join()
        except:
            for p in workerProc:
                p.terminate()
            writeProc.terminate()
        
        return dict(ani_af)

    def __fastani_worker(self, genomic_files, queue_in, queue_out):
        """Process each data item in parallel."""

        while True:
            gid, rep_ids = queue_in.get(block=True, timeout=None)
            if gid == None:
                break

            ani_af = self._fastani(gid, rep_ids, genomic_files)

            queue_out.put(ani_af)

    def __fastani_writer(self, all_ani_af, queue_writer):
        """Store or write results of worker threads in a single thread."""
        
        full_results = {}
        while True:
            ani_af = queue_writer.get(block=True, timeout=None)
            if ani_af == None:
                for qid in full_results:
                    all_ani_af[qid] = full_results[qid]
                    
                break

            for qid in ani_af:
                if qid not in full_results:
                    full_results[qid] = ani_af[qid]
                else:
                    full_results[qid].update(ani_af[qid])
                
    def _fastani(self, query_gid, rep_ids, genomic_files):
        """Calculate ANI between genomes and representatives genomes using FastANI."""
        
        # create file pointing to representative genome files
        tmp_rep_file = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
        fout = open(tmp_rep_file, 'w')
        for gid in rep_ids:
            fout.write('%s\n' % genomic_files[gid])
        fout.close()
        
        # ANI must be calculated between each genome and each representative
        tmp_fastani_file = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
        cmd = 'fastANI --minFrag -1 -q %s --rl %s -o %s 2> /dev/null' % (
                    genomic_files[query_gid], 
                    tmp_rep_file, 
                    tmp_fastani_file)
        os.system(cmd)

        ani_af = {}
        if os.path.exists(tmp_fastani_file) and os.stat(tmp_fastani_file).st_size > 0:
            for line in open(tmp_fastani_file):
                line_split = line.strip().split()

                query_genome = self._get_genome_id(line_split[0])
                ref_genome = self._get_genome_id(line_split[1])
                if query_genome not in ani_af:
                    ani_af[query_genome] = {}
                    
                ani = float(line_split[2])
                af = float(line_split[3])/int(line_split[4])
                ani_af[query_genome][ref_genome] = (ani, af)
        else:
            self.logger.warning('FastANI results file not created: %s' % query_gid)
            os.remove(tmp_rep_file)
            return ani_af
        
        os.remove(tmp_fastani_file)
        os.remove(tmp_rep_file)
        
        return ani_af

    def _genome_genus_clusters(self, genomes_to_process, full_gtdb_taxonomy):
        """Create clusters of genomes in the same genus."""
        
        genus_clusters = defaultdict(set)
        for gid in genomes_to_process:
            taxa = full_gtdb_taxonomy[gid]
            
            # *** [DEBUG] Useful case for debugging.
            #if taxa[5] != 'g__Pyrobaculum': #'g__Lactobacillus':
            #    continue

            if taxa[5] == 'g__':
                self.logger.error('Genome does not have a defined genus: %s' % gid)
                sys.exit(-1)
            else:
                genus_clusters[taxa[5]].add(gid)

        return genus_clusters
        
    def _assign_genomes(self, genus_clusters, canonical_gtdb_taxonomy, genomic_files):
        """Assign genomes in genus to best canonical genome."""
        
        # perform initial assignments using Mash
        self.logger.info('Performing initial assignments using Mash.')
        mash_genus_clusters = {}
        mash_genus_reps = {}
        mash_rep_assignments = {}
        for i, genus in enumerate(genus_clusters):
            # get representatives for this genus
            genus_rep_ids = set()
            for gid, taxa in canonical_gtdb_taxonomy.items():
                if taxa[5] == genus:
                    genus_rep_ids.add(gid)
            
            rep_assignments = self._mash_assignments(genus_clusters[genus] - genus_rep_ids, 
                                                            genus_rep_ids, 
                                                            genomic_files)
                                                            
            mash_rep_assignments.update(rep_assignments)
            
            assigned_genomes = set([d[0] for gid_data in rep_assignments.values() for d in gid_data])
            mash_genus_clusters[genus] = genus_clusters[genus] - assigned_genomes - genus_rep_ids
            
            mash_genus_reps[genus] = genus_rep_ids

            statusStr = 'Assigned %d of %d %s genomes to %d canonical genomes using Mash [%d of %d (%.2f%%)].'.ljust(104) % (
                            len(assigned_genomes),
                            len(genus_clusters[genus]) - len(genus_rep_ids),
                            genus,
                            len(genus_rep_ids),
                            i+1, 
                            len(genus_clusters), 
                            float(i+1)*100/len(genus_clusters))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

        # finalizes assignment using FastANI
        self.logger.info('Performing additional assignments using FastANI.')
        fout = open(os.path.join(self.output_dir, 'report_fastani_genus_sizes.tsv'), 'w')
        fout.write('Genus\tNo. genomes\tNo. representatives\n')
        fastani_rep_assignments = {}
        for i, genus in enumerate(mash_genus_clusters):
            rep_assignments = self._fastani_assignments(mash_genus_clusters[genus], 
                                                            mash_genus_reps[genus], 
                                                            genomic_files)
            fout.write('%s\t%d\t%d\n' % (genus, len(mash_genus_clusters[genus]), len(mash_genus_reps[genus])))
            
            fastani_rep_assignments.update(rep_assignments)

            statusStr = 'Assigned %d of %d %s genomes to %d canonical genomes using FastANI [%d of %d (%.2f%%)].'.ljust(104) % (
                            sum([len(gid_data) for gid_data in rep_assignments.values()]),
                            len(mash_genus_clusters[genus]),
                            genus,
                            len(mash_genus_reps[genus]),
                            i+1, 
                            len(mash_genus_clusters), 
                            float(i+1)*100/len(mash_genus_clusters))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')
        
        fout.close()
        
        return mash_rep_assignments, fastani_rep_assignments
        
    def _write_clusters(self, clusters, output_file):
        """Write out clustering information."""
        
        fout = open(output_file, 'w')
        for rep_id in sorted(clusters, key=lambda x: len(clusters[x]), reverse=True):  
            fout.write('%s\t%d\t%s\n' % (rep_id, len(clusters[rep_id]), ','.join([d[0] for d in clusters[rep_id]])))
        fout.close()

    def run(self, 
                canonical_taxonomy_file,
                full_taxonomy_file,
                metadata_file,
                genome_path_file,
                user_genomes):
        """Determine mapped reads between pairs of genomes at different ANI values."""
        
        t1 = time.time()
        
        self.logger.info('Assigning genomes to canonical genome set.')

        # read taxonomy of all genomes
        self.logger.info('Reading taxonomy of all genomes.')
        full_gtdb_taxonomy = Taxonomy().read(full_taxonomy_file)
        self.logger.info('Identified taxonomy for %d genomes.' % len(full_gtdb_taxonomy))

        # determining genomes to process
        self.logger.info('Determining genomes to process.')
        genomes_to_process, num_uba_genomes = self._genomes_to_process(full_gtdb_taxonomy, metadata_file, user_genomes)
        self.logger.info('Processing %d genomes, include %d UBA or User genomes.' % (len(genomes_to_process), num_uba_genomes))
        
        # get path to genomic files for ANI calculation
        self.logger.info('Reading path to genomic FASTA files.')
        genomic_files = {}
        for line in open(genome_path_file):
            line_split = line.strip().split('\t')
            gid = line_split[0]
            if gid in genomes_to_process:
                genome_path = line_split[1]
                assembly_id = os.path.basename(os.path.normpath(genome_path))
                genome_file = os.path.join(genome_path, assembly_id + '_genomic.fna')
                if user_genomes and gid.startswith('U_'): # sanity check genome file
                    if not os.path.exists(genome_file):
                        genomes_to_process.remove(gid)
                        self.logger.info('Skipping %s as genome file is missing.' % gid)
                        continue
                genomic_files[gid] = genome_file

        # read taxonomy of canonical genomes
        self.logger.info('Reading taxonomy of canonical genomes.')
        canonical_gtdb_taxonomy = Taxonomy().read(canonical_taxonomy_file)
        for gid in canonical_gtdb_taxonomy:
            if gid not in genomic_files or gid not in genomes_to_process:
                self.logger.error('Canonical genome is not in genome path file or list of genomes to process: %s' % gid)
                sys.exit(-1)
        self.logger.info('Identified taxonomy for %d canonical genomes.' % len(canonical_gtdb_taxonomy))
                
        # cluster genomes by genus
        self.logger.info('Assigning genomes to genus clusters.')
        genus_clusters = self._genome_genus_clusters(genomes_to_process, full_gtdb_taxonomy)
        self.logger.info('Identified %d genus clusters.' % len(genus_clusters))

        # create species clusters ordered by quality
        self.logger.info('Assigning genomes in each genera to genome in canonical genome set.')
        mash_rep_assignments, fastani_rep_assignments = self._assign_genomes(genus_clusters, canonical_gtdb_taxonomy, genomic_files)
        
        self._write_clusters(mash_rep_assignments, os.path.join(self.output_dir, 'assignments_mash.tsv'))
        self._write_clusters(fastani_rep_assignments, os.path.join(self.output_dir, 'assignments_fastani.tsv'))
        
        # get assigned genomes
        assigned_genomes = {}
        for assignments in [mash_rep_assignments, fastani_rep_assignments]:
            for rep_id, gid_data in assignments.items():
                for d in gid_data:
                    gid, source, ani, af = d
                    if gid in assigned_genomes:
                        print(gid, rep_id)
                    assigned_genomes[gid] = (rep_id, source, ani, af)
                    
        # add in canonical genomes
        for rep_id in canonical_gtdb_taxonomy:
            assigned_genomes[rep_id] = (rep_id, 'Representative', 100.0, 1.0)

        # write out taxonomy file
        fout = open(os.path.join(self.output_dir, 'gtdb_taxonomy.tsv'), 'w')
        num_assigned = 0
        num_unassigned = 0
        for gid in genomes_to_process:
            if gid in assigned_genomes:
                rep_id, source, ani, af = assigned_genomes[gid]
                if source == 'Mash':
                    ani_str = 'Mash: %.2f' % ani
                else:
                    ani_str = '%s: %.2f, %.3f' % (source, ani, af)
                    
                fout.write('%s\t%s\t%s\t%s\n' % (gid, '; '.join(canonical_gtdb_taxonomy[rep_id]), rep_id, ani_str))
                num_assigned += 1
            else:
                taxa = full_gtdb_taxonomy[gid][0:6] + ['s__']
                fout.write('%s\t%s\t%s\n' % (gid, '; '.join(taxa), 'unassigned'))
                num_unassigned += 1
        fout.close()
        
        total_genomes = num_assigned + num_unassigned
        self.logger.info('Assigned %d (%.2f%%) genomes and was unable to assign %d (%.2f%%) genomes.' % (num_assigned,
                                                                                                            num_assigned*100.0/total_genomes,
                                                                                                            num_unassigned,
                                                                                                            num_unassigned*100.0/total_genomes))
        
        t2 = time.time()
        self.logger.info('Elapsed time to assign genomes (seconds): %.2f' % (t2 - t1))
        self.logger.info('Done.')
