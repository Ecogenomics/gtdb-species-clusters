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
import subprocess
import shutil

import tomli 

from gtdblib.util.shell.execute import check_dependencies


def compress(input_file: str) -> str:
    """Gzip compress file."""

    logger = logging.getLogger('timestamp')

    cmd = ['pigz', '-f', input_file]

    logger.info(f"Executing: {' '.join(cmd)}")

    try:
        proc = subprocess.Popen(cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                encoding='utf-8')

        while True:
            out = proc.stdout.readline()
            if not out and proc.poll() is not None:
                break

        if proc.returncode != 0:
            logger.error(f'Return code: {proc.returncode}')
            sys.exit(1)
    except OSError as e:
        print(e)
        logger.error('Failed to compress file.')
        sys.exit(1)

    return input_file + '.gz'


def compress_dir(input_dir: str, remove_dir: bool = False) -> str:
    """Gzip compress directory."""

    logger = logging.getLogger('timestamp')

    compressed_dir_file = f'{input_dir}.tar.gz'
    cmd = ['tar', '--use-compress-program=pigz',
           '-cf', compressed_dir_file, input_dir]

    logger.info(f"Executing: {' '.join(cmd)}")

    try:
        proc = subprocess.Popen(cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                encoding='utf-8')

        while True:
            out = proc.stdout.readline()
            if not out and proc.poll() is not None:
                break

        if proc.returncode != 0:
            logger.error(f'Return code: {proc.returncode}')
            sys.exit(1)
    except OSError as e:
        print(e)
        logger.error('Failed to compress directory.')
        sys.exit(1)

    if remove_dir:
        shutil.rmtree(input_dir)

    return compressed_dir_file


class PMC_Cleanup(object):
    """Remove temporary files and compress large files."""

    def __init__(self):
        """Initialization."""

        check_dependencies(['pigz'])
        self.log = logging.getLogger('rich')

    def check_output_dir_exists(self, output_dir: str) -> bool:
        """Check if output directory exists."""

        if not os.path.exists(output_dir):
            self.log.warning(f'Missing output directory: {output_dir}')
            return False

        return True

    def clean_root(self, output_dir: str) -> None:
        """Compress files in root directory."""

        if not self.check_output_dir_exists(output_dir):
            self.log.warning('Invalid root directory.')
            sys.exit(1)

        for f in os.listdir(output_dir):
            fp = os.path.join(output_dir, f)
            if os.path.isdir(fp):
                continue

            if f.endswith('.gz'):
                continue
            
            if os.path.getsize(fp) > 1_000_000:
                # compress files > 1MB
                compress(fp)

    def clean_u_new_genomes(self, output_dir: str) -> None:
        """Clean up files for u_new_genomes step."""

        if not self.check_output_dir_exists(output_dir):
            return

        compress_dir(output_dir, remove_dir=True)

    def clean_u_qc_genomes(self, output_dir: str) -> None:
        """Clean up files for u_qc_genomes step."""

        if not self.check_output_dir_exists(output_dir):
            return

        compress_dir(output_dir, remove_dir=True)

    def clean_u_gtdbtk(self, output_dir: str) -> None:
        """Clean up files for u_gtdbtk step."""

        if not self.check_output_dir_exists(output_dir):
            return

        # remove temporary directories
        for gtdb_tk_dir in os.listdir(output_dir):
            if gtdb_tk_dir == 'gtdbtk_batch0':
                # keep the first directory as an example of what was executed
                continue

            if gtdb_tk_dir.startswith('gtdbtk_batch'):
                shutil.rmtree(os.path.join(output_dir, gtdb_tk_dir))

        compress_dir(output_dir, remove_dir=True)

    def clean_u_lpsn_rna_types(self, output_dir: str) -> None:
        """Clean up files for u_lpsn_rna_types step."""

        if not self.check_output_dir_exists(output_dir):
            return

        compress_dir(output_dir, remove_dir=True)

    def clean_u_resolve_types(self, output_dir: str) -> None:
        """Clean up files for u_resolve_types step."""

        if not self.check_output_dir_exists(output_dir):
            return

        # remove temporary ANI results directory
        tmp_ani_dir = os.path.join(output_dir, 'ani_pickles')
        if os.path.exists(tmp_ani_dir):
            shutil.rmtree(tmp_ani_dir)

        compress_dir(output_dir, remove_dir=True)

    def clean_u_rep_changes(self, output_dir: str) -> None:
        """Clean up files for u_rep_changes step."""

        if not self.check_output_dir_exists(output_dir):
            return

        compress_dir(output_dir, remove_dir=True)

    def clean_u_rep_actions(self, output_dir: str) -> None:
        """Clean up files for u_rep_actions step."""

        if not self.check_output_dir_exists(output_dir):
            return

        compress_dir(output_dir, remove_dir=True)

    def clean_u_sel_reps(self, output_dir: str) -> None:
        """Clean up files for u_sel_reps step."""

        if not self.check_output_dir_exists(output_dir):
            return

        # remove temporary ANI and Mash files
        tmp_ani_file = os.path.join(output_dir, 'reps_ani_af.pkl')
        if os.path.exists(tmp_ani_file):
            os.remove(tmp_ani_file)

        tmp_mash_file = os.path.join(output_dir, 'gtdb_reps.msh')
        if os.path.exists(tmp_mash_file):
            os.remove(tmp_mash_file)

        compress_dir(output_dir, remove_dir=True)

    def clean_u_cluster_named_reps(self, output_dir: str) -> None:
        """Clean up files for u_cluster_named_reps step."""

        if not self.check_output_dir_exists(output_dir):
            return

        # remove temporary ANI and Mash files
        tmp_ani_file = os.path.join(output_dir, 'ani_af_nonrep_vs_rep.pkl')
        if os.path.exists(tmp_ani_file):
            os.remove(tmp_ani_file)

        tmp_mash_nonrep_file = os.path.join(output_dir, 'gtdb_nonreps.msh')
        if os.path.exists(tmp_mash_nonrep_file):
            os.remove(tmp_mash_nonrep_file)

        tmp_mash_rep_file = os.path.join(output_dir, 'gtdb_reps.msh')
        if os.path.exists(tmp_mash_rep_file):
            os.remove(tmp_mash_rep_file)

        compress_dir(output_dir, remove_dir=True)

    def clean_u_cluster_de_novo(self, output_dir: str) -> None:
        """Clean up files for u_cluster_de_novo step."""

        if not self.check_output_dir_exists(output_dir):
            return

        # remove temporary ANI and Mash files
        tmp_ani_file = os.path.join(output_dir, 'ani_af_rep_vs_nonrep.de_novo.pkl')
        if os.path.exists(tmp_ani_file):
            os.remove(tmp_ani_file)

        tmp_mash_nonrep_file = os.path.join(output_dir, 'gtdb_nonrep_genomes.msh')
        if os.path.exists(tmp_mash_nonrep_file):
            os.remove(tmp_mash_nonrep_file)

        tmp_mash_rep_file = os.path.join(output_dir, 'gtdb_rep_genomes.msh')
        if os.path.exists(tmp_mash_rep_file):
            os.remove(tmp_mash_rep_file)

        tmp_mash_unclustered_file = os.path.join(output_dir, 'gtdb_unclustered_genomes.msh')
        if os.path.exists(tmp_mash_unclustered_file):
            os.remove(tmp_mash_unclustered_file)

        tmp_skani_sketch_dir = os.path.join(output_dir, 'skani_ref_sketches')
        if os.path.exists(tmp_skani_sketch_dir):
            shutil.rmtree(tmp_skani_sketch_dir)

        compress_dir(output_dir, remove_dir=True)

    def clean_u_cluster_stats(self, output_dir: str) -> None:
        """Clean up files for u_cluster_stats step."""

        if not self.check_output_dir_exists(output_dir):
            return

        compress_dir(output_dir, remove_dir=True)

    def clean_u_ncbi_erroneous(self, output_dir: str) -> None:
        """Clean up files for u_ncbi_erroneous step."""

        if not self.check_output_dir_exists(output_dir):
            return

        compress_dir(output_dir, remove_dir=True)

    def clean_u_synonyms(self, output_dir: str) -> None:
        """Clean up files for u_synonyms step."""

        if not self.check_output_dir_exists(output_dir):
            return

        compress_dir(output_dir, remove_dir=True)

    def clean_u_curation_trees(self, output_dir: str) -> None:
        """Clean up files for u_curation_trees step."""

        if not self.check_output_dir_exists(output_dir):
            return

        compress_dir(output_dir, remove_dir=True)

    def clean_u_species_init(self, output_dir: str) -> None:
        """Clean up files for u_species_init step."""

        if not self.check_output_dir_exists(output_dir):
            return

        compress_dir(output_dir, remove_dir=True)

    def clean_u_pmc_species_names_ar(self, output_dir: str) -> None:
        """Clean up files for u_pmc_species_names_ar step."""

        if not self.check_output_dir_exists(output_dir):
            return

        compress_dir(output_dir, remove_dir=True)

    def clean_u_pmc_species_names_bac(self, output_dir: str) -> None:
        """Clean up files for u_pmc_species_names_bac step."""

        if not self.check_output_dir_exists(output_dir):
            return

        compress_dir(output_dir, remove_dir=True)

    def run(self, output_toml_file: str):
        """Remove temporary files and compress large files."""

        with open(output_toml_file, mode="rb") as fp:
            output_dirs = tomli.load(fp)

        self.clean_root(output_dirs["root_dir"])
        self.clean_u_new_genomes(os.path.join(output_dirs["root_dir"], output_dirs["u_new_genomes"]))
        self.clean_u_qc_genomes(os.path.join(output_dirs["root_dir"], output_dirs["u_qc_genomes"]))
        self.clean_u_gtdbtk(os.path.join(output_dirs["root_dir"], output_dirs["u_gtdbtk"]))
        self.clean_u_lpsn_rna_types(os.path.join(output_dirs["root_dir"], output_dirs["u_lpsn_rna_types"]))
        self.clean_u_resolve_types(os.path.join(output_dirs["root_dir"], output_dirs["u_resolve_types"]))
        self.clean_u_rep_changes(os.path.join(output_dirs["root_dir"], output_dirs["u_rep_changes"]))
        self.clean_u_rep_actions(os.path.join(output_dirs["root_dir"], output_dirs["u_rep_actions"]))
        self.clean_u_sel_reps(os.path.join(output_dirs["root_dir"], output_dirs["u_sel_reps"]))
        self.clean_u_cluster_named_reps(os.path.join(output_dirs["root_dir"], output_dirs["u_cluster_named_reps"]))
        self.clean_u_cluster_de_novo(os.path.join(output_dirs["root_dir"], output_dirs["u_cluster_de_novo"]))
        self.clean_u_cluster_stats(os.path.join(output_dirs["root_dir"], output_dirs["u_cluster_stats"]))
        self.clean_u_ncbi_erroneous(os.path.join(output_dirs["root_dir"], output_dirs["u_ncbi_erroneous"]))
        self.clean_u_synonyms(os.path.join(output_dirs["root_dir"], output_dirs["u_synonyms"]))
        self.clean_u_curation_trees(os.path.join(output_dirs["root_dir"], output_dirs["u_curation_trees"]))
        self.clean_u_species_init(os.path.join(output_dirs["root_dir"], output_dirs["u_species_init"]))
        self.clean_u_pmc_species_names_ar(os.path.join(output_dirs["root_dir"], output_dirs["u_pmc_species_names_ar"]))
        self.clean_u_pmc_species_names_bac(os.path.join(output_dirs["root_dir"], output_dirs["u_pmc_species_names_bac"]))
