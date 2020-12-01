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
import argparse
import logging
import re
import csv
import datetime
from collections import defaultdict

from numpy import (mean as np_mean, std as np_std)

from gtdb_species_clusters.genome import Genome
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.species_priority_manager import SpeciesPriorityManager
from gtdb_species_clusters.taxon_utils import is_placeholder_taxon, generic_name, specific_epithet
from gtdb_species_clusters.genome_utils import parse_ncbi_bioproject


class Sandbox(object):
    """Play to explore new ideas or calculate one off statistics relted to species clusters."""

    def __init__(self, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
    def parse_lpsn_scrape_sp_priorities(self):
        """Get priority of species from scraping LPSN website."""
        
        test_years = {}
        with open('../../strains_table/test/year_table.tsv') as f:
            for line in f:
                tokens = line.strip().split('\t')
                
                test_years['s__' + tokens[0]] = int(tokens[1])
        
        lpsn_sp_priorities = {}
        with open('../../lpsn/20201118/parse_html/lpsn_species.tsv') as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                
                sp = tokens[0]
                if sp == 's__':
                    # *** hack to skip bad case in file
                    # Pierre to fix
                    print('bad', line)
                    continue 
                    
                lpsn_species_authority = tokens[2]
                ref_str = lpsn_species_authority.split(', ')[0]
                references = ref_str.replace('(', '').replace(')', '')
                years = re.sub(r'emend\.[^\d]*\d{4}', '', references)
                years = re.sub(r'ex [^\d]*\d{4}', ' ', years)
                years = re.findall('[1-3][0-9]{3}', years, re.DOTALL)
                years = [int(y) for y in years if int(y) <= datetime.datetime.now().year]
                
                if len(years) == 0:
                    # assume this name is validated through ICN and just take the first 
                    # date given as the year of priority
                    years = re.findall('[1-3][0-9]{3}', references, re.DOTALL)
                    years = [int(y) for y in years if int(y) <= datetime.datetime.now().year]
                    
                lpsn_sp_priorities[sp] = years[0]
                
                if years[0] != test_years[sp]:
                    print('error', sp, years[0], test_years[sp])
                
        return lpsn_sp_priorities
        
    def parse_lpsn_sp_priorities(self, lpsn_gss_file):
        """Get priority of species from LPSN GSS file."""
        
        self.logger.info('Reading priority references from LPSN.')
        lpsn_sp_priorities = {}
        illegitimate_names = set()
        with open(lpsn_gss_file, encoding='utf-8', errors='ignore') as f:
            csv_reader = csv.reader(f)

            for line_num, tokens in enumerate(csv_reader):
                if line_num == 0:
                    genus_idx = tokens.index('genus_name')
                    specific_idx = tokens.index('sp_epithet')
                    subsp_idx = tokens.index('subsp_epithet')
                    status_idx = tokens.index('status')
                    author_idx = tokens.index('authors')
                else:
                    generic = tokens[genus_idx].strip().replace('"', '')
                    specific = tokens[specific_idx].strip().replace('"', '')
                    subsp = tokens[subsp_idx].strip().replace('"', '')
                    if not specific or subsp:
                        # skip genus names and subspecies
                        continue
                        
                    species = 's__{} {}'.format(generic, specific)
                    
                    status = tokens[status_idx].strip().replace('"', '')
                    status_tokens = [t.strip() for t in status.split(';')]
                    status_tokens = [tt.strip() for t in status_tokens for tt in t.split(',') ]
                    
                    if 'illegitimate name' in status_tokens:
                        illegitimate_names.add(species)
                        if species in lpsn_sp_priorities:
                            continue

                    # get priority references, ignoring references if they are
                    # marked as being a revied name as indicated by a 'ex' or 'emend'
                    # (e.g. Holospora (ex Hafkine 1890) Gromov and Ossipov 1981)
                    ref_str = tokens[author_idx]
                    references = ref_str.replace('(', '').replace(')', '')
                    years = re.sub(r'emend\.[^\d]*\d{4}', '', references)
                    years = re.sub(r'ex [^\d]*\d{4}', ' ', years)
                    years = re.findall('[1-3][0-9]{3}', years, re.DOTALL)
                    years = [int(y) for y in years if int(y) <= datetime.datetime.now().year]

                    if (species not in illegitimate_names
                        and species in lpsn_sp_priorities 
                        and years[0] != lpsn_sp_priorities[species]):
                            # conflict that can't be attributed to one of the entries being
                            # considered an illegitimate name
                            self.logger.error('Conflicting priority references for {}: {} {}'.format(
                                                species, years, lpsn_sp_priorities[species]))

                    lpsn_sp_priorities[species] = years[0]
                        
        self.logger.info(f' - establish priority for {len(lpsn_sp_priorities):,} species using LPSN GSS metadata.')
        
        return lpsn_sp_priorities
        
    def eval_bacdive_priority_info(self, cur_genomes, lpsn_gss_file):
        """Determine situations where BacDive priority conflicts with LPSN, 
            and if BacDive should be retained as source of priority."""
            
        # get priority references for species at LPSN
        lpsn_scrape_sp_priorities = self.parse_lpsn_scrape_sp_priorities()
        lpsn_sp_priorities = self.parse_lpsn_sp_priorities(lpsn_gss_file)

        # determine cases where GTDB type strain genomes has different priority
        # at BacDive then in LPSN GSS file
        diff_bacdive_priority = 0
        diff_lpsn_priority = 0
        identical_priority = 0
        both_no_priority = 0
        for gid in cur_genomes:
            if not cur_genomes[gid].is_gtdb_type_strain():
                continue
                
            bacdive_priority = cur_genomes[gid].dsmz_priority_year
            
            lpsn_priority = Genome.NO_PRIORITY_YEAR
            ncbi_sp = cur_genomes[gid].ncbi_taxa.species.replace('[', '').replace(']', '')
            if ncbi_sp in lpsn_sp_priorities:
                lpsn_priority = lpsn_sp_priorities[ncbi_sp]
            elif ncbi_sp in lpsn_scrape_sp_priorities:
                lpsn_priority = lpsn_scrape_sp_priorities[ncbi_sp]
            
            if bacdive_priority != Genome.NO_PRIORITY_YEAR and bacdive_priority != lpsn_priority:
                if bacdive_priority < lpsn_priority:
                    diff_bacdive_priority += 1
                    print(ncbi_sp, bacdive_priority, lpsn_priority)
                else:
                    diff_lpsn_priority += 1
            elif bacdive_priority == lpsn_priority:
                if bacdive_priority == Genome.NO_PRIORITY_YEAR:
                    both_no_priority += 1
                else:
                    identical_priority += 1
        
        print('diff_bacdive_priority', diff_bacdive_priority)
        print('diff_lpsn_priority', diff_lpsn_priority)
        print('identical_priority', identical_priority)
        print('both_no_priority', both_no_priority)
        
    def eval_priority_gtdb_clusters(self, prev_genomes, lpsn_gss_file):
        """Evaluate correct naming priority of GTDB species clusters.
        
        Brought about due to change in rules use to establish the 
        priority of names.
        """
        
        # get priority references for species at LPSN
        lpsn_sp_priorities = self.parse_lpsn_sp_priorities(lpsn_gss_file)
        
        # identify GTDB species clusters with type strain genomes from multiple species
        bad_rids = {}
        good_rids = set()
        for rid, cids in prev_genomes.sp_clusters.clusters():
            if prev_genomes[rid].is_gtdb_type_strain():
                rep_gtdb_sp = prev_genomes[rid].gtdb_taxa.species
                rep_ncbi_sp = prev_genomes[rid].ncbi_taxa.species
                if rep_ncbi_sp in lpsn_sp_priorities:
                    for cid in cids:
                        if cid == rid:
                            continue
                            
                        if prev_genomes[cid].is_gtdb_type_strain():
                            ncbi_sp = prev_genomes[cid].ncbi_taxa.species
                            
                            if ncbi_sp in lpsn_sp_priorities:
                                # if only one of the genomes has an NCBI species assignment
                                # that matches the proposed GTDB species this is the name
                                # the should be used
                                if generic_name(rep_gtdb_sp) == generic_name(rep_ncbi_sp) and generic_name(rep_gtdb_sp) != generic_name(ncbi_sp):
                                    # species name follows rules of priority
                                    good_rids.add(rid)
                                else:
                                    # priority is based on the species name
                                    if ((lpsn_sp_priorities[ncbi_sp] < lpsn_sp_priorities[rep_ncbi_sp] and specific_epithet(rep_gtdb_sp) != specific_epithet(ncbi_sp))
                                            or (lpsn_sp_priorities[rep_ncbi_sp] < lpsn_sp_priorities[ncbi_sp] and specific_epithet(rep_gtdb_sp) != specific_epithet(rep_ncbi_sp))):
                                            print('Incorrect naming priority', rid, rep_ncbi_sp, lpsn_sp_priorities[rep_ncbi_sp], cid, ncbi_sp, lpsn_sp_priorities[ncbi_sp])
                                            
                                            if rid not in bad_rids or lpsn_sp_priorities[ncbi_sp] < bad_rids[rid][4]:
                                                bad_rids[rid] = (rep_gtdb_sp, 
                                                                    rep_ncbi_sp, 
                                                                    lpsn_sp_priorities[rep_ncbi_sp], 
                                                                    ncbi_sp, 
                                                                    lpsn_sp_priorities[ncbi_sp])
                                    else:
                                        good_rids.add(rid)
                                        
        print('good_rids', len(good_rids))
                                    
        # write out bad cases
        fout = open(os.path.join(self.output_dir, 'gtdb_r95_incorrect_sp_priority.tsv'), 'w')
        fout.write('GTDB representative\tGTDB species\tNCBI species\tYear of priority\tCorrected NCBI species with priority\tYear of priority\n')
        for rid in bad_rids:
            rep_gtdb_sp, rep_ncbi_sp, rep_priority, cor_ncbi_sp, cor_priority = bad_rids[rid]
            fout.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        rid,
                        rep_gtdb_sp,
                        rep_ncbi_sp,
                        rep_priority,
                        cor_ncbi_sp,
                        cor_priority))
                        
        fout.close()
        
    def test_sp_priority_rules(self, cur_genomes, sp_priority_mngr, lpsn_gss_file):
        """Test implications of different species priority rules."""
        
        # Test cases:
        #references = "(Aizawa et al. 2010) Patel and Gupta 2020"
        #references = "Pukall et al. 2006 emend. Nouioui et al. 2018"
        #references = "Robinson and Ritchie 1981"
        #references = "Patel et al. 1980 emend. Robinson and Ritchie 1981 emend. Khan et al. 1984 emend. Murray 1986 emend. Tindall 2019"
        #references = "(Pasteur 1864) De Ley and Frateur 1974 (Approved Lists 1980)"
        #references = "(ex Frankland and Frankland 1887) Vandamme et al. 2016"
        #references = "(van Niel 1928) Scholz and Kilian 2016 emend. Nouioui et al. 2018"
        #references = "(ex Ahrens 1968) Uchino et al. 1999 emend. Liu et al. 2016"

        # get priority references for species at LPSN
        self.logger.info('Reading priority references from LPSN.')
        sp_approved_lists_1980 = set()
        lpsn_sp_priorities = defaultdict(list)
        with open(lpsn_gss_file, encoding='utf-8', errors='ignore') as f:
            csv_reader = csv.reader(f)

            for line_num, tokens in enumerate(csv_reader):
                if line_num == 0:
                    genus_idx = tokens.index('genus_name')
                    specific_idx = tokens.index('sp_epithet')
                    subsp_idx = tokens.index('subsp_epithet')
                    status_idx = tokens.index('status')
                    author_idx = tokens.index('authors')
                else:
                    generic = tokens[genus_idx].strip().replace('"', '')
                    specific = tokens[specific_idx].strip().replace('"', '')
                    subsp = tokens[subsp_idx].strip().replace('"', '')
                    if not specific or subsp:
                        # skip genus names and subspecies
                        continue
                        
                    species = 's__{} {}'.format(generic, specific)
                    
                    status = tokens[status_idx].strip().replace('"', '')
                    status_tokens = [t.strip() for t in status.split(';')]
                    status_tokens = [tt.strip() for t in status_tokens for tt in t.split(',') ]

                    if 'illegitimate name' in status_tokens:
                        continue

                    # get priority references, ignoring references if they are
                    # marked as being a revied name as indicated by a 'ex' or 'emend'
                    # (e.g. Holospora (ex Hafkine 1890) Gromov and Ossipov 1981)
                    ref_str = tokens[author_idx]
                    references = ref_str.replace('(', '').replace(')', '')
                    years = re.sub(r'emend\.[^\d]*\d{4}', '', references)
                    years = re.sub(r'ex [^\d]*\d{4}', ' ', years)
                    years = re.findall('[1-3][0-9]{3}', years, re.DOTALL)
                    years = [int(y) for y in years if int(y) <= datetime.datetime.now().year]

                    if species in lpsn_sp_priorities and years != lpsn_sp_priorities[species]:
                        self.logger.error('Conflicting priority references for {}: {} {}'.format(
                                            species, years, lpsn_sp_priorities[species]))
                    
                    lpsn_sp_priorities[species] = years 
                    
                    if 'Approved Lists 1980' in ref_str:
                        sp_approved_lists_1980.add(species)
                        
        self.logger.info(f' - establish priority for {len(lpsn_sp_priorities):,} species using LPSN GSS metadata.')
        
        # get current date of priority for type strain of species
        self.logger.info('Determing NCBI species with year of priority in GTDB metadata.')
        fout = open(os.path.join(self.output_dir, 'gtdb_sp_priority.tsv'), 'w')
        fout.write('Representative ID\tNCBI species\tGTDB species\tGTDB type strain\tGTDB type subspecies\tPriority\n')
        cur_priority_year = {}
        for rid in cur_genomes:
            ncbi_sp = cur_genomes[rid].ncbi_taxa.species
            priority = cur_genomes[rid].year_of_priority()

            fout.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        rid,
                        ncbi_sp, 
                        cur_genomes[rid].gtdb_taxa.species,
                        cur_genomes[rid].is_gtdb_type_strain(),
                        cur_genomes[rid].is_gtdb_type_subspecies(),
                        priority))
            
            if priority == Genome.NO_PRIORITY_YEAR:
                continue
                
            if not cur_genomes[rid].is_gtdb_type_strain():
                continue
            
            if ncbi_sp in cur_priority_year and priority != cur_priority_year[ncbi_sp][1]:
                prev_rid, prev_priority = cur_priority_year[ncbi_sp]
                self.logger.error('NCBI species has conflicting priority dates: {}/{}/{}/{}, {}/{}/{}/{}'.format(
                                    rid,
                                    ncbi_sp, 
                                    cur_genomes[rid].gtdb_taxa.species,
                                    priority,
                                    prev_rid, 
                                    cur_genomes[rid].ncbi_taxa.species,
                                    cur_genomes[rid].gtdb_taxa.species,
                                    prev_priority))
                                    
            cur_priority_year[ncbi_sp] = (rid, priority)
            
        fout.close()
            
        self.logger.info(' - identified year of priority for {:,} NCBI species.'.format(len(cur_priority_year)))
        
        # get number of species that follow "first reference" and "second reference" priority rule
        self.logger.info('Evaluating priorities used in GTDB compared to LPSN:')
        single_year_agree = 0
        single_year_disagree = 0
        first_year_agree = 0
        first_year_disagree = 0
        last_year_agree = 0
        last_year_disagree = 0
        last_year_from_dsmz = 0
        for sp in lpsn_sp_priorities:
            if sp not in cur_priority_year:
                continue
                
            rid, gtdb_priority_year = cur_priority_year[sp]
                
            if len(lpsn_sp_priorities[sp]) == 1:
                if lpsn_sp_priorities[sp][0] == gtdb_priority_year:
                    single_year_agree += 1
                else:
                    single_year_disagree += 1
                    #***print('single priority reference', sp, lpsn_sp_priorities[sp][0], gtdb_priority_year)
            elif sp not in sp_approved_lists_1980:
                # check agreement with first date for species not on AL 1980
                if lpsn_sp_priorities[sp][0] == gtdb_priority_year:
                    first_year_agree += 1
                else:
                    first_year_disagree += 1
                
                # check agreement with last date for species not on AL 1980
                if lpsn_sp_priorities[sp][-1] == gtdb_priority_year:
                    last_year_agree += 1
                else:
                    last_year_disagree += 1
                    
                    if cur_priority_year[sp][1] == cur_genomes[rid].dsmz_priority_year:
                        last_year_from_dsmz += 1
                        
                    print('multiple priority references', sp, lpsn_sp_priorities[sp], gtdb_priority_year, gtdb_priority_year == cur_genomes[rid].dsmz_priority_year)
                    
        self.logger.info(' - single priority reference, same priority year: {:,} ({:.1f}%)'.format(
                            single_year_agree, single_year_agree*100.0 / (single_year_agree+single_year_disagree)))
        self.logger.info(' - single priority reference, different priority year: {:,} ({:.1f}%)'.format(
                            single_year_disagree, single_year_disagree*100.0 / (single_year_agree+single_year_disagree)))
                            
        self.logger.info(' - multiple priority references, same first LPSN priority year: {:,} ({:.1f}%)'.format(
                            first_year_agree, first_year_agree*100.0 / (first_year_agree+first_year_disagree)))
        self.logger.info(' - multiple priority references, different first LPSN priority year: {:,} ({:.1f}%)'.format(
                            first_year_disagree, first_year_disagree*100.0 / (first_year_agree+first_year_disagree)))
                            
        self.logger.info(' - multiple priority references, same last LPSN priority year: {:,} ({:.1f}%)'.format(
                            last_year_agree, last_year_agree*100.0 / (last_year_agree+last_year_disagree)))
        self.logger.info(' - multiple priority references, different last LPSN priority year: {:,} ({:.1f}%)'.format(
                            last_year_disagree, last_year_disagree*100.0 / (last_year_agree+last_year_disagree)))
        self.logger.info(' - multiple priority references, different last LPSN priority year result of using DSMZ priority: {:,} ({:.1f}%)'.format(
                            last_year_from_dsmz, last_year_from_dsmz*100.0 / last_year_disagree))

    def run(self, 
            prev_gtdb_metadata_file,
            cur_gtdb_metadata_file,
            qc_passed_file,
            ncbi_genbank_assembly_file,
            untrustworthy_type_file,
            gtdb_type_strains_ledger,
            sp_priority_ledger,
            genus_priority_ledger,
            ncbi_env_bioproject_ledger,
            lpsn_gss_file):
        """Play to explore new ideas or calculate one off statistics relted to species clusters."""
        
        lpsn_scrape_sp_priorities = self.parse_lpsn_scrape_sp_priorities()
        lpsn_sp_priorities = self.parse_lpsn_sp_priorities(lpsn_gss_file)
        
        diff_p = 0
        same_p = 0
        for sp in set(lpsn_sp_priorities).intersection(lpsn_scrape_sp_priorities):
            if lpsn_sp_priorities[sp] != lpsn_scrape_sp_priorities[sp]:
                print('diff', sp, lpsn_sp_priorities[sp], lpsn_scrape_sp_priorities[sp])
                diff_p += 1
            else:
                same_p += 1
                
        print('diff_p', diff_p)
        print('same_p', same_p)
        
        print('exclusive to scraping file', len(set(lpsn_scrape_sp_priorities) - set(lpsn_sp_priorities)))
        print('exclusive to gss file', len(set(lpsn_sp_priorities) - set(lpsn_scrape_sp_priorities)))
        
        return #***

        # create previous and current GTDB genome sets
        if False:
            self.logger.info('Creating previous GTDB genome set.')
            prev_genomes = Genomes()
            prev_genomes.load_from_metadata_file(prev_gtdb_metadata_file,
                                                    gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                    ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                                    untrustworthy_type_ledger=untrustworthy_type_file,
                                                    ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)
            self.logger.info(' - previous genome set has {:,} species clusters spanning {:,} genomes.'.format(
                                len(prev_genomes.sp_clusters),
                                prev_genomes.sp_clusters.total_num_genomes()))

            # evalute priority issues with GTDB R95 clusters
            self.eval_priority_gtdb_clusters(prev_genomes, lpsn_gss_file)

        # create current GTDB genome sets
        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                create_sp_clusters=False,
                                                qc_passed_file=qc_passed_file,
                                                ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                                untrustworthy_type_ledger=untrustworthy_type_file,
                                                ncbi_env_bioproject_ledger=ncbi_env_bioproject_ledger)

        # evalute priority information from BacDive
        self.eval_bacdive_priority_info(cur_genomes, lpsn_gss_file)
        
        return 
                            
        # initialize species priority manager
        sp_priority_mngr = SpeciesPriorityManager(sp_priority_ledger,
                                                        genus_priority_ledger,
                                                        lpsn_gss_file,
                                                        self.output_dir)
                                                        
        # test implications of different species priority rules
        self.test_sp_priority_rules(cur_genomes, 
                                    sp_priority_mngr,
                                    lpsn_gss_file)
