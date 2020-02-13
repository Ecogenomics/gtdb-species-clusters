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
import pickle
from collections import defaultdict, namedtuple
from itertools import combinations

from numpy import (mean as np_mean, std as np_std)

from biolib.external.execute import check_dependencies

from gtdb_species_clusters.fastani import FastANI
from gtdb_species_clusters.genomes import Genomes
from gtdb_species_clusters.genome_utils import canonical_gid, exclude_from_refseq
from gtdb_species_clusters.taxon_utils import generic_name, specific_epithet, canonical_taxon
from gtdb_species_clusters.type_genome_utils import symmetric_ani


class ResolveTypes(object):
    """Resolve cases where a species has multiple genomes assembled from the type strain."""

    def __init__(self, ani_cache_file, cpus, output_dir):
        """Initialization."""
        
        self.ltp_dir = 'rna_ltp_132'
        self.ltp_results_file = 'ssu.taxonomy.tsv'
        self.LTP_METADATA = namedtuple('LTP_METADATA', 'taxonomy taxa species ssu_len evalue bitscore aln_len perc_iden perc_aln')
        
        self.ltp_pi_threshold = 99.0
        self.ltp_pa_threshold = 90.0
        self.ltp_ssu_len_threshold = 900
        self.ltp_evalue_threshold = 1e-10
        
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        self.cpus = cpus
        
        self.fastani = FastANI(ani_cache_file, cpus)
        
        self.ani_pickle_dir = os.path.join(self.output_dir, 'ani_pickles')
        if not os.path.exists(self.ani_pickle_dir):
            os.makedirs(self.ani_pickle_dir)
            
    def _parse_ltp_taxonomy_str(self, ltp_taxonomy_str):
        """Parse taxa and species from LTP taxonomy string."""
        
        if ';type sp.|' in ltp_taxonomy_str:
            taxa = ltp_taxonomy_str.split(';type sp.|')[0].split(';')
        elif ';|' in ltp_taxonomy_str:
            taxa = ltp_taxonomy_str.split(';|')[0].split(';')
        elif '|' in ltp_taxonomy_str:
            taxa = ltp_taxonomy_str.split('|')[0].split(';')
        elif ltp_taxonomy_str[-1] == ';':
            taxa = ltp_taxonomy_str[0:-1].split(';')
        else:
            taxa = ltp_taxonomy_str.split(';')
            
        sp = taxa[-1]
        if ' subsp. ' in sp:
            sp = ' '.join(sp.split()[0:2])
            
        # validate that terminal taxon appears to be a 
        # valid binomial species name
        if (sp[0].islower() or
            any(c.isdigit() for c in sp)
            or any(c.isupper() for c in sp[1:])):
            print(ltp_taxonomy_str, taxa)
            assert False
            
        return taxa, 's__' + sp
            
    def parse_ltp_metadata(self, type_gids, cur_genomes):
        """Parse Living Tree Project 16S rRNA metadata."""

        metadata = defaultdict(list)
        for gid in type_gids:
            genome_path = os.path.dirname(os.path.abspath(cur_genomes[gid].genomic_file))
            ltp_file = os.path.join(genome_path, self.ltp_dir, self.ltp_results_file)
            if os.path.exists(ltp_file):
                with open(ltp_file) as f:
                    header = f.readline().strip().split('\t')
                    
                    query_id_index = header.index('query_id')
                    taxonomy_index = header.index('taxonomy')
                    ssu_len_index = header.index('length')
                    evalue_index = header.index('blast_evalue')
                    bitscore_index = header.index('blast_bitscore')
                    aln_len_index = header.index('blast_align_len')
                    pi_index = header.index('blast_perc_identity')
                    
                    for line in f:
                        tokens = line.strip().split('\t')
                        
                        query_id = tokens[query_id_index]
                        taxonomy = tokens[taxonomy_index]
                        ssu_len = int(tokens[ssu_len_index])
                        evalue = float(tokens[evalue_index])
                        bitscore = float(tokens[bitscore_index])
                        aln_len = int(tokens[aln_len_index])
                        pi = float(tokens[pi_index])
                        
                        taxa, sp = self._parse_ltp_taxonomy_str(taxonomy)
                        
                        metadata[gid].append(self.LTP_METADATA(
                                                taxonomy=taxonomy,
                                                taxa=taxa,
                                                species=sp,
                                                ssu_len=ssu_len,
                                                evalue=evalue,
                                                bitscore=bitscore,
                                                aln_len=aln_len,
                                                perc_iden=pi,
                                                perc_aln=aln_len*100.0/ssu_len))
                                                
        return metadata
        
    def ltp_defined_species(self, ltp_taxonomy_file):
        """Get all species present in the LTP database."""
        
        ltp_species = set()
        with open(ltp_taxonomy_file, encoding = 'utf-8') as f:
            for line in f:
                tokens = line.strip().split('\t')
                
                taxonomy = tokens[1]
                taxa, sp = self._parse_ltp_taxonomy_str(taxonomy)
                ltp_species.add(sp)
                
        return ltp_species
        
    def ltp_species(self, gid, ltp_metadata):
        """Get high confident species assignments."""
        
        sp = set()
        for hit in ltp_metadata[gid]:
            # check if hit should be trusted
            if (hit.perc_iden >= self.ltp_pi_threshold 
                and hit.perc_aln >= self.ltp_pa_threshold 
                and hit.ssu_len >= self.ltp_ssu_len_threshold 
                and hit.evalue < self.ltp_evalue_threshold):
                sp.add(hit.species)
                
        return sp

    def check_strain_ani(self, gid_anis, untrustworthy_gids):
        """Check if genomes meet strain ANI criteria."""
        
        all_similar = True
        for gid1, gid2 in combinations(gid_anis, 2):
            if gid1 in untrustworthy_gids or gid2 in untrustworthy_gids:
                continue

            if gid_anis[gid1][gid2] < 99:
                return False
                
        return True
            
    def resolve_by_intra_specific_ani(self, gid_anis):
        """Resolve by removing intra-specific genomes with divergent ANI values."""
        
        if len(gid_anis) <= 2:
            return False, {}
        
        # consider most divergent genome as untrustworthy
        untrustworthy_gids = {}
        while True:
            # find most divergent genome
            min_ani = 100
            untrustworthy_gid = None
            for gid in gid_anis:
                if gid in untrustworthy_gids:
                    continue
                    
                anis = [ani for cur_gid, ani in gid_anis[gid].items() if cur_gid not in untrustworthy_gids]
                if np_mean(anis) < min_ani:
                    min_ani = np_mean(anis)
                    untrustworthy_gid = gid
            
            untrustworthy_gids[untrustworthy_gid] = f'{min_ani:.2f}% ANI to other type strain genomes'
            
            all_similar = self.check_strain_ani(gid_anis, untrustworthy_gids)
    
            if all_similar:
                return True, untrustworthy_gids
            
            remaining_genomes = len(gid_anis) - len(untrustworthy_gids)
            if remaining_genomes <= 2 or len(untrustworthy_gids) >= len(gid_anis):
                return False, {}

    def resolve_by_ncbi_types(self, gid_anis, type_gids, cur_genomes):
        """Resolve by consulting NCBI type material metadata."""

        untrustworthy_gids = {}
        ncbi_type_count = 0
        for gid in type_gids:
            if not cur_genomes[gid].is_ncbi_type_strain():
                untrustworthy_gids[gid] = 'Not classified as assembled from type material at NCBI'
            else:
                ncbi_type_count += 1
            
        all_similar = self.check_strain_ani(gid_anis, untrustworthy_gids)
        
        if all_similar and len(untrustworthy_gids) > 0 and ncbi_type_count > 0:
            return True, untrustworthy_gids
            
        return False, {}
        
    def resolve_gtdb_family(self, gid_anis, ncbi_sp, type_gids, cur_genomes):
        """Resolve by identifying genomes with a conflicting GTDB family assignment."""
        
        genus = 'g__' + generic_name(ncbi_sp)
        gtdb_genus_rep = cur_genomes.gtdb_type_species_of_genus(genus)
        if not gtdb_genus_rep:
            return False, {}
        
        expected_gtdb_family = cur_genomes[gtdb_genus_rep].gtdb_taxa.family
        
        untrustworthy_gids = {}
        matched_family = 0
        for gid in type_gids:
            if cur_genomes[gid].gtdb_taxa.family == expected_gtdb_family:
                matched_family += 1
            else:
                # genome is classified to a different GTDB family than
                # expected for this species
                 untrustworthy_gids[gid] = f'Conflicting GTDB family assignment of {cur_genomes[gid].gtdb_taxa.family}, expected {expected_gtdb_family}'

        all_similar = self.check_strain_ani(gid_anis, untrustworthy_gids)
        
        # conflict is resolved if remaining genomes pass ANI similarity test, 
        if all_similar and len(untrustworthy_gids) > 0 and matched_family > 0:
            return True, untrustworthy_gids

        return False, {}
        
    def resolve_gtdb_genus(self, gid_anis, ncbi_sp, type_gids, cur_genomes):
        """Resolve by identifying genomes with a conflicting GTDB genus assignments."""
        
        ncbi_genus = 'g__' + generic_name(ncbi_sp)
        
        untrustworthy_gids = {}
        matched_genus = 0
        for gid in type_gids:
            canonical_gtdb_genus = canonical_taxon(cur_genomes[gid].gtdb_taxa.genus)

            if ncbi_genus == canonical_gtdb_genus:
                matched_genus += 1
            else:
                untrustworthy_gids[gid] = f'Conflicting GTDB genus assignment of {cur_genomes[gid].gtdb_taxa.genus}, expected {ncbi_genus}'

        all_similar = self.check_strain_ani(gid_anis, untrustworthy_gids)
        
        if all_similar and len(untrustworthy_gids) > 0 and matched_genus > 0:
            return True, untrustworthy_gids

        return False, {}
        
    def resolve_gtdb_species(self, gid_anis, ncbi_sp, type_gids, cur_genomes):
        """Resolve by identifying genomes with a conflicting GTDB species assignments to different type material."""
        
        ncbi_sp_epithet = specific_epithet(ncbi_sp)
        
        untrustworthy_gids = {}
        matched_sp_epithet = 0
        for gid in type_gids:
            if ncbi_sp_epithet == cur_genomes[gid].gtdb_taxa.specific_epithet:
                matched_sp_epithet += 1
            else:
                # check if genome is classified to a GTDB species cluster supported
                # by a type strain genome in which case we should consider this 
                # genome untrustworthy
                gtdb_sp = cur_genomes[gid].gtdb_taxa.species
                if gtdb_sp != 's__':
                    gtdb_sp_rid = cur_genomes.gtdb_sp_rep(gtdb_sp)
                    if cur_genomes[gtdb_sp_rid].is_effective_type_strain():
                        # genome has been assigned to another species
                        # defined by a type strain genome
                        ani, af = self.fastani.symmetric_ani_cached(gid, 
                                                        gtdb_sp_rid, 
                                                        cur_genomes[gid].genomic_file, 
                                                        cur_genomes[gtdb_sp_rid].genomic_file)
                        untrustworthy_gids[gid] = f'Conflicting GTDB species assignment of {cur_genomes[gid].gtdb_taxa.species} [ANI={ani:.2f}%; AF={af:.2f}%]'

        all_similar = self.check_strain_ani(gid_anis, untrustworthy_gids)
        
        # conflict is resolved if remaining genomes pass ANI similarity test, 
        if all_similar and len(untrustworthy_gids) > 0 and matched_sp_epithet > 0:
            return True, untrustworthy_gids

        return False, {}
        
    def resolve_validated_untrustworthy_ncbi_genomes(self, gid_anis, ncbi_sp, type_gids, ltp_metadata, ltp_defined_species, cur_genomes):
        """Resolve by identifying genomes marked as `untrustworthy as type` at NCBI and with conflicting LTP assignments."""

        if ncbi_sp not in ltp_defined_species:
            return False, {}
        
        untrustworthy_gids = {}
        for gid in type_gids:
            if 'untrustworthy as type' in cur_genomes[gid].excluded_from_refseq_note.lower():
                ltp_species = self.ltp_species(gid, ltp_metadata)

                if ncbi_sp not in ltp_species and len(ltp_species) > 0:
                    untrustworthy_gids[gid] = f"Conflicting 16S rRNA hits to LTP database of {' / '.join(set(ltp_species))}"
                    
        all_similar = self.check_strain_ani(gid_anis, untrustworthy_gids)

        # conflict is resolved if remaining genomes pass ANI similarity test, 
        if all_similar and len(untrustworthy_gids) > 0:
            return True, untrustworthy_gids
                
        return False, {}
        
    def resolve_ltp_conflict(self, gid_anis, ncbi_sp, type_gids, ltp_metadata, require_conflict_sp):
        """Resolve by considering BLAST hits of 16S rRNA genes to LTP database."""
        
        untrustworthy_gids = {}
        genomes_matching_expected_sp = 0
        for gid in type_gids:
            expected_sp_count = 0
            match_unexpected_sp = []
            for hit in ltp_metadata[gid]:
                # check if hit should be trusted
                if (hit.perc_iden >= self.ltp_pi_threshold 
                    and hit.perc_aln >= self.ltp_pa_threshold 
                    and hit.ssu_len >= self.ltp_ssu_len_threshold 
                    and hit.evalue < self.ltp_evalue_threshold):
                    ltp_sp = hit.species
                    if ltp_sp == ncbi_sp:
                        expected_sp_count += 1
                    else:
                        match_unexpected_sp.append(ltp_sp)
 
            if expected_sp_count == 0 and len(match_unexpected_sp) >= require_conflict_sp:
                if len(match_unexpected_sp) > 0:
                    untrustworthy_gids[gid] = f"Conflicting 16S rRNA hits to LTP database of {' / '.join(set(match_unexpected_sp))}"
                else:
                    untrustworthy_gids[gid] = "Lack of 16S rRNA hits to LTP database"
            elif expected_sp_count > len(match_unexpected_sp):
                genomes_matching_expected_sp += 1

        all_similar = self.check_strain_ani(gid_anis, untrustworthy_gids)
        
        if all_similar and len(untrustworthy_gids) > 0 and genomes_matching_expected_sp > 0:
            return True, untrustworthy_gids

        return False, {}
                
    def run(self, 
                cur_gtdb_metadata_file,
                cur_genomic_path_file,
                qc_passed_file,
                gtdbtk_classify_file,
                ncbi_genbank_assembly_file,
                ltp_taxonomy_file,
                gtdb_type_strains_ledger,
                untrustworthy_type_ledger):
        """Resolve cases where a species has multiple genomes assembled from the type strain."""
        
        # get species in LTP reference database
        self.logger.info('Determining species defined in LTP reference database.')
        ltp_defined_species = self.ltp_defined_species(ltp_taxonomy_file)
        self.logger.info(f' ... identified {len(ltp_defined_species):,} species.')
        
        # create current GTDB genome sets
        self.logger.info('Creating current GTDB genome set.')
        cur_genomes = Genomes()
        cur_genomes.load_from_metadata_file(cur_gtdb_metadata_file,
                                                gtdb_type_strains_ledger=gtdb_type_strains_ledger,
                                                create_sp_clusters=False,
                                                uba_genome_file=None,
                                                qc_passed_file=qc_passed_file,
                                                ncbi_genbank_assembly_file=ncbi_genbank_assembly_file,
                                                untrustworthy_type_ledger=untrustworthy_type_ledger,
                                                gtdbtk_classify_file=gtdbtk_classify_file)
        cur_genomes.load_genomic_file_paths(cur_genomic_path_file)
        self.logger.info(f' ... current genome set contains {len(cur_genomes):,} genomes.')
        
        # parsing genomes manually established to be untrustworthy as type
        self.logger.info('Determining genomes manually annotated as untrustworthy as type.')
        manual_untrustworthy_types = {}
        with open(untrustworthy_type_ledger) as f:
            header = f.readline().strip().split('\t')
            
            ncbi_sp_index = header.index('NCBI species')
            reason_index = header.index('Reason for declaring untrustworthy')
            
            for line in f:
                tokens = line.strip().split('\t')
                
                gid = canonical_gid(tokens[0])
                manual_untrustworthy_types[gid] = (tokens[ncbi_sp_index], tokens[reason_index])
        self.logger.info(f' ... identified {len(manual_untrustworthy_types):,} genomes manually annotated as untrustworthy as type.')

        # identify NCBI species with multiple genomes assembled from type strain of species
        self.logger.info('Determining number of type strain genomes in each NCBI species.')
        sp_type_strain_genomes = defaultdict(set)
        for gid in cur_genomes:
            if cur_genomes[gid].is_effective_type_strain():
                ncbi_sp = cur_genomes[gid].ncbi_taxa.species
                if ncbi_sp != 's__':
                    # yes, NCBI has genomes marked as assembled from type material
                    # that do not actually have a binomial species name
                    sp_type_strain_genomes[ncbi_sp].add(gid)

        multi_type_strains_sp = [ncbi_sp for ncbi_sp, gids in sp_type_strain_genomes.items() if len(gids) > 1]
        self.logger.info(f' ... identified {len(multi_type_strains_sp):,} NCBI species with multiple assemblies indicated as being type strain genomes.')
        
        # sort by number of genome assemblies
        self.logger.info('Calculating ANI between type strain genomes in each species.')
        
        fout = open(os.path.join(self.output_dir, 'multi_type_strain_species.tsv'), 'w')
        fout.write('NCBI species\tNo. type strain genomes\t>=99% ANI\tMean ANI\tStd ANI\tMean AF\tStd AF\tResolution\tGenome IDs\n')
        
        fout_genomes = open(os.path.join(self.output_dir, 'type_strain_genomes.tsv'), 'w')
        fout_genomes.write('Genome ID\tUntrustworthy\tNCBI species\tGTDB genus\tGTDB species\tLTP species\tConflict with prior GTDB assignment')
        fout_genomes.write('\tMean ANI\tStd ANI\tMean AF\tStd AF\tExclude from RefSeq\tNCBI taxonomy\tGTDB taxonomy\n')
        
        fout_unresolved = open(os.path.join(self.output_dir, 'unresolved_type_strain_genomes.tsv'), 'w')
        fout_unresolved.write('Genome ID\tNCBI species\tGTDB genus\tGTDB species\tLTP species')
        fout_unresolved.write('\tMean ANI\tStd ANI\tMean AF\tStd AF\tExclude from RefSeq\tNCBI taxonomy\tGTDB taxonomy\n')
        
        fout_high_divergence = open(os.path.join(self.output_dir, 'highly_divergent_type_strain_genomes.tsv'), 'w')
        fout_high_divergence.write('Genome ID\tNCBI species\tGTDB genus\tGTDB species\tLTP species\tMean ANI\tStd ANI\tMean AF\tStd AF\tExclude from RefSeq\tNCBI taxonomy\tGTDB taxonomy\n')
        
        fout_untrustworthy = open(os.path.join(self.output_dir, 'untrustworthy_type_material.tsv'), 'w')
        fout_untrustworthy.write('Genome ID\tNCBI species\tGTDB species\tLTP species\tReason for declaring untrustworthy\n')
        for gid in manual_untrustworthy_types:
            ncbi_sp, reason = manual_untrustworthy_types[gid]
            fout_untrustworthy.write('{}\t{}\t{}\t{}\t{}\n'.format(
                                        gid, 
                                        ncbi_sp, 
                                        cur_genomes[gid].gtdb_taxa.species,
                                        '<not tested>',
                                        'n/a',
                                        'Manual curation: ' + reason))
        
        processed = 0
        num_divergent = 0
        unresolved_sp_count = 0
        
        ncbi_ltp_resolved = 0
        intra_ani_resolved = 0
        ncbi_type_resolved = 0
        gtdb_family_resolved = 0
        gtdb_genus_resolved = 0
        gtdb_sp_resolved = 0
        ltp_resolved = 0
        
        use_pickled_results = False #***
        if use_pickled_results:
            self.logger.warning('Using previously calculated ANI results in: {}'.format(self.ani_pickle_dir))
        
        prev_gtdb_sp_conflicts = 0
        for ncbi_sp, type_gids in sorted(sp_type_strain_genomes.items(), key=lambda kv: len(kv[1])):
            if len(type_gids) == 1:
                continue
                
            status_str = '-> Processing {} with {:,} type strain genomes [{:,} of {:,} ({:.2f}%)].'.format(
                                ncbi_sp, 
                                len(type_gids),
                                processed+1, 
                                len(multi_type_strains_sp),
                                (processed+1)*100.0/len(multi_type_strains_sp)).ljust(128)
            sys.stdout.write('{}\r'.format(status_str))
            sys.stdout.flush()
            processed += 1

            # calculate ANI between type strain genomes
            ncbi_sp_str = ncbi_sp[3:].lower().replace(' ', '_')
            if not use_pickled_results: #***
                ani_af = self.fastani.pairwise(type_gids, cur_genomes.genomic_files)
                pickle.dump(ani_af, open(os.path.join(self.ani_pickle_dir, f'{ncbi_sp_str}.pkl'), 'wb'))
            else:
                ani_af = pickle.load(open(os.path.join(self.ani_pickle_dir, f'{ncbi_sp_str}.pkl'), 'rb'))
            
            anis = []
            afs = []
            gid_anis = defaultdict(lambda: {})
            gid_afs = defaultdict(lambda: {})
            all_similar = True
            for gid1, gid2 in combinations(type_gids, 2):
                ani, af = symmetric_ani(ani_af, gid1, gid2)
                if ani < 99 or af < 0.65:
                    all_similar = False
                    
                anis.append(ani)
                afs.append(af)
                
                gid_anis[gid1][gid2] = ani
                gid_anis[gid2][gid1] = ani
                
                gid_afs[gid1][gid2] = af
                gid_afs[gid2][gid1] = af
                
            note = 'All type strain genomes have ANI >99% and AF >65%.'
            unresolved_species = False
            
            # read LTP metadata for genomes
            ltp_metadata = self.parse_ltp_metadata(type_gids, cur_genomes)

            untrustworthy_gids = {}
            gtdb_resolved_sp_conflict = False
            if not all_similar:
                # need to establish which genomes are untrustworthy as type
                num_divergent += 1
                unresolved_species = True
                
                # write out highly divergent cases for manual inspection; 
                # these should be compared to the automated selection
                if np_mean(anis) < 95:
                    for gid in type_gids:
                        ltp_species = self.ltp_species(gid, ltp_metadata)
                            
                        fout_high_divergence.write('{}\t{}\t{}\t{}\t{}\t{:.2f}\t{:.3f}\t{:.3f}\t{:.4f}\t{}\t{}\t{}\n'.format(
                                                        gid,
                                                        ncbi_sp,
                                                        cur_genomes[gid].gtdb_taxa.genus,
                                                        cur_genomes[gid].gtdb_taxa.species,
                                                        ' / '.join(ltp_species),
                                                        np_mean(list(gid_anis[gid].values())),
                                                        np_std(list(gid_anis[gid].values())),
                                                        np_mean(list(gid_afs[gid].values())),
                                                        np_std(list(gid_afs[gid].values())),
                                                        cur_genomes[gid].excluded_from_refseq_note,
                                                        cur_genomes[gid].ncbi_taxa,
                                                        cur_genomes[gid].gtdb_taxa))
                
                # filter genomes marked as `untrustworthy as type` at NCBI and where the LTP
                # assignment also suggest the asserted type material is incorrect
                resolved, untrustworthy_gids = self.resolve_validated_untrustworthy_ncbi_genomes(gid_anis, 
                                                                                                    ncbi_sp, 
                                                                                                    type_gids, 
                                                                                                    ltp_metadata, 
                                                                                                    ltp_defined_species,
                                                                                                    cur_genomes)
                if resolved:
                    note = "Species resolved by removing genomes considered `untrustworthy as type` and with a LTP BLAST hit confirming the assembly is likely untrustworthy"
                    ncbi_ltp_resolved += 1

                # try to resolve by LTP 16S BLAST results
                if not resolved:
                    resolved, untrustworthy_gids = self.resolve_ltp_conflict(gid_anis, ncbi_sp, type_gids, ltp_metadata, 0)
                    if resolved:
                        note = 'Species resolved by identifying conflicting or lack of LTP BLAST results'
                        ltp_resolved += 1

                # try to resolve species using intra-specific ANI test
                if not resolved:
                    resolved, untrustworthy_gids = self.resolve_by_intra_specific_ani(gid_anis)
                    if resolved:
                        note = 'Species resolved by intra-specific ANI test'
                        intra_ani_resolved += 1

                # try to resolve by GTDB family assignment
                if not resolved:
                    resolved, untrustworthy_gids = self.resolve_gtdb_family(gid_anis, ncbi_sp, type_gids, cur_genomes)
                    if resolved:
                        note = 'Species resolved by consulting GTDB family classifications'
                        gtdb_family_resolved += 1
                
                # try to resolve by GTDB genus assignment
                if not resolved:
                    resolved, untrustworthy_gids = self.resolve_gtdb_genus(gid_anis, ncbi_sp, type_gids, cur_genomes)
                    if resolved:
                        note = 'Species resolved by consulting GTDB genus classifications'
                        gtdb_genus_resolved += 1
                           
                # try to resolve by GTDB species assignment
                if not resolved:
                    resolved, untrustworthy_gids = self.resolve_gtdb_species(gid_anis, ncbi_sp, type_gids, cur_genomes)
                    if resolved:
                        note = 'Species resolved by consulting GTDB species classifications'
                        gtdb_sp_resolved += 1
                        
                # try to resolve by considering genomes annotated as type material at NCBI,
                # which includes considering if genomes are marked as untrustworthy as type
                if not resolved:
                    resolved, untrustworthy_gids = self.resolve_by_ncbi_types(gid_anis, type_gids, cur_genomes)
                    if resolved:
                        note = 'Species resolved by consulting NCBI assembled from type metadata'
                        ncbi_type_resolved += 1

                if resolved:
                    unresolved_species = False
                    
                    # check if type strain genomes marked as trusted or untrusted conflict
                    # with current GTDB species assignment
                    untrustworthy_gtdb_sp_match = False
                    trusted_gtdb_sp_match = False
                    for gid in type_gids:
                        gtdb_canonical_epithet = canonical_taxon(specific_epithet(cur_genomes[gid].gtdb_taxa.species))
                        if gtdb_canonical_epithet == specific_epithet(ncbi_sp):
                            if gid in untrustworthy_gids:
                                untrustworthy_gtdb_sp_match = True
                            else:
                                trusted_gtdb_sp_match = True

                    if untrustworthy_gtdb_sp_match and not trusted_gtdb_sp_match:
                        prev_gtdb_sp_conflicts += 1
                        gtdb_resolved_sp_conflict = True

                    # write results to file
                    for gid, reason in untrustworthy_gids.items():
                        ltp_species = self.ltp_species(gid, ltp_metadata)
                        
                        if 'untrustworthy as type' in cur_genomes[gid].excluded_from_refseq_note:
                            reason += "; considered `untrustworthy as type` at NCBI"
                        fout_untrustworthy.write('{}\t{}\t{}\t{}\t{}\n'.format(gid,
                                                                                ncbi_sp,
                                                                                cur_genomes[gid].gtdb_taxa.species,
                                                                                ' / '.join(ltp_species),
                                                                                reason))
                                                                                
                        # Sanity check that if the untrustworthy genome has an LTP to only the
                        # expected species, that all other genomes also have a hit to the 
                        # expected species (or potentially no hit). Otherwise, more consideration
                        # should be given to the genome with the conflicting LTP hit.
                        if len(ltp_species) == 1 and ncbi_sp in ltp_species:
                            other_sp = set()
                            for test_gid in type_gids:
                                ltp_species = self.ltp_species(test_gid, ltp_metadata)
                                if ltp_species and ncbi_sp not in ltp_species:
                                    other_sp.update(ltp_species)
                                
                            if other_sp:
                                self.logger.warning(f'Genome {gid} marked as untrustworthy, but this conflicts with high confidence LTP 16S rRNA assignment.')
                                
                    num_ncbi_untrustworthy = sum([1 for gid in type_gids if 'untrustworthy as type' in cur_genomes[gid].excluded_from_refseq_note])
                    if num_ncbi_untrustworthy != len(type_gids):
                        for gid in type_gids:
                            if (gid not in untrustworthy_gids 
                                and 'untrustworthy as type' in cur_genomes[gid].excluded_from_refseq_note):
                                self.logger.warning("Retaining genome {} from {} despite being marked as `untrustworthy as type` at NCBI [{:,} of {:,} considered untrustworthy].".format(
                                                        gid, 
                                                        ncbi_sp,
                                                        num_ncbi_untrustworthy,
                                                        len(type_gids)))
                else:
                    note = 'Species is unresolved; manual curation is required!'
                    unresolved_sp_count += 1
                    
                if unresolved_species:
                    for gid in type_gids:
                        ltp_species = self.ltp_species(gid, ltp_metadata)
                            
                        fout_unresolved.write('{}\t{}\t{}\t{}\t{}\t{:.2f}\t{:.3f}\t{:.3f}\t{:.4f}\t{}\t{}\t{}\n'.format(
                                    gid,
                                    ncbi_sp,
                                    cur_genomes[gid].gtdb_taxa.genus,
                                    cur_genomes[gid].gtdb_taxa.species,
                                    ' / '.join(ltp_species),
                                    np_mean(list(gid_anis[gid].values())),
                                    np_std(list(gid_anis[gid].values())),
                                    np_mean(list(gid_afs[gid].values())),
                                    np_std(list(gid_afs[gid].values())),
                                    cur_genomes[gid].excluded_from_refseq_note,
                                    cur_genomes[gid].ncbi_taxa,
                                    cur_genomes[gid].gtdb_taxa))

            for gid in type_gids:
                ltp_species = self.ltp_species(gid, ltp_metadata)
                    
                fout_genomes.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{:.3f}\t{:.3f}\t{:.4f}\t{}\t{}\t{}\n'.format(
                            gid,
                            gid in untrustworthy_gids,
                            ncbi_sp,
                            cur_genomes[gid].gtdb_taxa.genus,
                            cur_genomes[gid].gtdb_taxa.species,
                            ' / '.join(ltp_species),
                            gtdb_resolved_sp_conflict,
                            np_mean(list(gid_anis[gid].values())),
                            np_std(list(gid_anis[gid].values())),
                            np_mean(list(gid_afs[gid].values())),
                            np_std(list(gid_afs[gid].values())),
                            cur_genomes[gid].excluded_from_refseq_note,
                            cur_genomes[gid].ncbi_taxa,
                            cur_genomes[gid].gtdb_taxa))

            fout.write('{}\t{}\t{}\t{:.2f}\t{:.3f}\t{:.3f}\t{:.4f}\t{}\t{}\n'.format(
                        ncbi_sp,
                        len(type_gids),
                        all_similar,
                        np_mean(anis),
                        np_std(anis),
                        np_mean(afs),
                        np_std(afs),
                        note,
                        ', '.join(type_gids)))

        sys.stdout.write('\n')
        fout.close()
        fout_unresolved.close()
        fout_high_divergence.close()
        fout_genomes.close()
        fout_untrustworthy.close()
        
        self.logger.info(f'Identified {num_divergent:,} species with 1 or more divergent type strain genomes.')
        self.logger.info(f' ... resolved {ncbi_ltp_resolved:,} species by removing NCBI `untrustworthy as type` genomes with a conflicting LTP 16S rRNA classifications.')
        self.logger.info(f' ... resolved {ltp_resolved:,} species by considering conflicting LTP 16S rRNA classifications.')
        self.logger.info(f' ... resolved {intra_ani_resolved:,} species by considering intra-specific ANI values.')
        self.logger.info(f' ... resolved {gtdb_family_resolved:,} species by considering conflicting GTDB family classifications.')
        self.logger.info(f' ... resolved {gtdb_genus_resolved:,} species by considering conflicting GTDB genus classifications.')
        self.logger.info(f' ... resolved {gtdb_sp_resolved:,} species by considering conflicting GTDB species classifications.')
        self.logger.info(f' ... resolved {ncbi_type_resolved:,} species by considering type material designations at NCBI.')

        if unresolved_sp_count > 0:
            self.logger.warning(f'There are {unresolved_sp_count:,} unresolved species with multiple type strain genomes.')
            self.logger.warning('These should be handled before proceeding with the next step of GTDB species updating.')
            self.logger.warning("This can be done by manual curation and adding genomes to 'untrustworthy_type_ledger'.")
        
        self.logger.info(f'Identified {prev_gtdb_sp_conflicts:,} cases where resolved type strain conflicts with prior GTDB assignment.')
