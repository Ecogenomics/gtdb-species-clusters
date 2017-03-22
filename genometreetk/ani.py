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
import shutil
import logging
import operator
import itertools
import tempfile
from collections import defaultdict
import multiprocessing as mp

from biolib.taxonomy import Taxonomy
from biolib.external.execute import check_on_path

from genometreetk.common import (read_gtdb_taxonomy,
                                    read_gtdb_metadata,
                                    filter_genomes)
                                    
from numpy import (mean as np_mean,
                    percentile as np_percentile)


class ANI(object):
    """Calculate ANI between named species."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Maximum number of cpus/threads to use.
        """

        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus
        
        check_on_path('ani_calculator')

    def __worker(self,
                    metadata_file,
                    nt_files,
                    max_genomes,
                    queue_in,
                    queue_out):
        """Process each species in parallel."""
        
        metadata = read_gtdb_metadata(metadata_file, ['checkm_completeness', 
                                                        'checkm_contamination'])
            
        genome_quality = {}
        for genome_id, m in metadata.iteritems():
            genome_quality[genome_id] = m.checkm_completeness - 5*m.checkm_contamination
                
        while True:
            species, genome_ids = queue_in.get(block=True, timeout=None)
            if species == None:
                break
                
            # select highest quality genomes
            if len(genome_ids) > max_genomes:
                t = [(genome_id, q) 
                        for genome_id, q in genome_quality.iteritems() 
                        if genome_id in genome_ids]
                hq_genomes = sorted(t, key=lambda x: x[1], reverse=True)[0:max_genomes]
                genome_ids = [x[0] for x in hq_genomes]

            ani = []
            af = []
            results = ''
            tmp_dir = tempfile.mkdtemp()
            for gi, gj in itertools.combinations(genome_ids, 2):
                tmp_file = tempfile.NamedTemporaryFile(delete=False, dir=tmp_dir)
                tmp_file.close()
                cmd = ('ani_calculator ' + 
                        '-genome1fna %s ' + 
                        '-genome2fna %s ' +
                        '-outfile %s -outdir %s ' +
                        '> /dev/null') % (nt_files[gi], 
                                            nt_files[gj],
                                            tmp_file.name,
                                            tmp_dir)
      
                os.system(cmd)
                with open(tmp_file.name) as f:
                    f.readline()
                    for line in f:
                        results += '%s\t%s' % (species, line)
                        line_split = line.strip().split('\t')
                        ani.append(0.5*(float(line_split[2]) + float(line_split[3])))
                        af.append(0.5*(float(line_split[4]) + float(line_split[5])))

            shutil.rmtree(tmp_dir)

            queue_out.put((species, ani, af, genome_ids, results))

    def __writer(self, num_species, output_dir, writer_queue):
        """Write results for each species."""
        
        # gather results for each genome
        output_file = os.path.join(output_dir, 'ani_species.tsv')
        fout = open(output_file, 'w')
        fout.write('Species\tNo. Sampled Genomes\tMean ANI\tMedian ANI\t5th Percentile\t95th Percentile')
        fout.write('\tMean AF\tMedian AF\t5th Percentile\t95th Percentile')
        fout.write('\tSampled Genomes\n')
        
        output_file = os.path.join(output_dir, 'ani.tsv')
        fout_pw = open(output_file, 'w')
        fout_pw.write('Species\tGenome 1\tGenome 2\tANI(1->2)\tANI(2->1)\tAF(1->2)\tAF(2->1)\n')
        processed = 0
        while True:
            species, ani, af, genome_ids, results = writer_queue.get(block=True, timeout=None)
            if species == None:
              break

            processed += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) species.' % (processed,
                                                                            num_species,
                                                                            float(processed) * 100 / num_species)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            fout_pw.write(results)
            
            row = '%s\t%d' % (species, len(genome_ids))
            mean_ani = np_mean(ani)
            p5, median, p95 = np_percentile(ani, [5, 50, 95])
            row += '\t%.2f\t%.2f\t%.2f\t%.2f' % (mean_ani,
                                                    median,
                                                    p5, p95)
            mean_af = np_mean(af)
            p5, median, p95 = np_percentile(af, [5, 50, 95])
            row += '\t%.2f\t%.2f\t%.2f\t%.2f' % (mean_af*100,
                                                    median*100,
                                                    p5*100, p95*100)
            fout.write('%s\t%s\n' % (row, ','.join(genome_ids)))

        sys.stdout.write('\n')

        fout.close()
        fout_pw.close()

    def run(self,
                input_taxonomy,
                genome_path_file,
                metadata_file, 
                max_genomes,
                min_comp,
                max_cont,
                min_quality, 
                max_contigs, 
                min_N50, 
                max_ambiguous, 
                max_gap_length, 
                output_dir):
        """Calculate ANI for named species."""
        
        # get genomes passing filtering criteria
        filtered_genome_ids = filter_genomes(metadata_file,
                                                min_comp,
                                                max_cont,
                                                min_quality, 
                                                max_contigs, 
                                                min_N50, 
                                                max_ambiguous, 
                                                max_gap_length)
                                                
        # get species in each named species
        taxonomy = Taxonomy().read(input_taxonomy)
        genome_ids_to_remove = set(taxonomy.keys()) - filtered_genome_ids
        for genome_id in genome_ids_to_remove:
            del taxonomy[genome_id]
            
        named_species = Taxonomy().extant_taxa_for_rank('species', taxonomy)
        
        # get path to nucleotide files
        nt_files = {}
        for line in open(genome_path_file):
            line_split = line.strip().split('\t')

            gtdb_id = line_split[0]
            genome_id = gtdb_id.replace('GB_', '').replace('RS_', '')

            genome_dir = line_split[1]

            nt_file = os.path.join(genome_dir, 'prodigal', genome_id + '_protein.fna')
            nt_files[gtdb_id] = nt_file

        # populate worker queue with data to process
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()

        num_species = 0
        for species, genome_ids in named_species.iteritems():
            if len(genome_ids) > 1:
                worker_queue.put((species, genome_ids))
                num_species += 1

        for _ in range(self.cpus):
          worker_queue.put((None, None))

        try:
          worker_proc = [mp.Process(target=self.__worker, args=(metadata_file,
                                                                    nt_files,
                                                                    max_genomes,
                                                                    worker_queue,
                                                                    writer_queue)) for _ in range(self.cpus)]
          write_proc = mp.Process(target=self.__writer, args=(num_species,
                                                                  output_dir,
                                                                  writer_queue))

          write_proc.start()

          for p in worker_proc:
              p.start()

          for p in worker_proc:
              p.join()

          writer_queue.put((None, None, None, None, None))
          write_proc.join()
        except:
          for p in worker_proc:
            p.terminate()

          write_proc.terminate()
