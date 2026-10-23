# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Internal GTDB toolkit for establishing, updating, and validating the species clusters of each Genome Taxonomy Database release (e.g., R232). It is a single CLI, `bin/gtdb_species_clusters`, whose subcommands are run sequentially by a curator during a release cycle. Distributed via the `ace-internal` conda channel.

## Working conventions

- **Always track the current task in `CURRENT_TASK.md`** at the repo root, so work can be resumed after a machine crash or a lost connection. Create it when starting a task. Record the goal, the branch, steps done, steps remaining, and any non-obvious notes (gotchas found, decisions the user made). Update it as each step completes, not only at the end. When resuming a session, read it first. The file is git-ignored. Clear it or replace its contents when a new task starts.
- **Every PR bumps the version in `gtdb_species_clusters/VERSION`.** The first line is the version string read by `__init__.py` and `setup.py`. Add the new version at the top of the file, followed by brief `- ` bullet lines summarizing the changes and then a blank line, matching the existing entries. Put the detailed description in the PR body.

## Commands

```bash
# Run the CLI from a checkout (no install needed)
python bin/gtdb_species_clusters                  # list all subcommands, grouped by release stage
python bin/gtdb_species_clusters <command> -h     # command-specific options

# Most commands take a TOML params file as their first positional argument
python bin/gtdb_species_clusters u_qc_genomes params.toml

# Tests (pytest; covers the Genome and Genomes classes and LTP taxonomy parsing)
python -m pytest gtdb_species_clusters/tests
python -m pytest gtdb_species_clusters/tests/test_genome.py::TestGenome::test_pass_qc
```

There is no lint/format config in the repo. The system `python` is 3.6 with no dependencies, so activate the appropriate conda env before running anything. For R237 that env is `gtdb_species_clustering-r237`, which has Python 3.12, skani 0.4.0 and an editable install of this checkout (`pip install -e .`). Runtime dependencies are declared in `setup.py`. They include `gtdblib`, the shared GTDB utility library (`check_file_exists`, `canonical_gid`, `Taxonomy`, `check_dependencies`, etc.), which is available only from PyPI. They also include `numpy`, `dendropy`, `Levenshtein`, `tqdm`, and `toml` on Python < 3.11. External binaries are `skani`, `mash`, and `gtdbtk`. Use GTDB-Tk 2.8.0 (env `gtdbtk-2.8.0`, which sets `GTDBTK_DATA_PATH` to the R232 reference data) for `u_gtdbtk`, stacked on the toolkit env: `mamba activate gtdb_species_clustering-r237 && mamba activate --stack gtdbtk-2.8.0`. `gtdb_species_clusters` keeps using the toolkit env's Python via its shebang, while `gtdbtk` and its dependencies (pplacer, skani, etc.) come from the GTDB-Tk env. Deactivate the stacked env afterwards.

## Architecture

**CLI → OptionsParser → step class.** `bin/gtdb_species_clusters` defines every argparse subparser and the help menu. `gtdb_species_clusters/main.py:OptionsParser` has one method per subcommand. Each method parses the TOML file, calls `check_file_exists` on the inputs, and then creates a step class (e.g., `UpdateSelectRepresentatives`) and calls its `.run(...)`. `OptionsParser.run()` dispatches through a long `if/elif` chain. Adding a command means touching all three places: the subparser, the `OptionsParser` method, and the dispatch branch. Each `u_*` method also calls `log_command_inputs` with the TOML keys and earlier-step outputs it uses. This logs them as written in the TOML, using the raw copy that `parse_toml_file` keeps in `params['toml_raw']`. When you add or remove an input of a `u_*` command, update that list too.

**Release pipeline, in module-name order:**
- `update_*.py` (`u_*` commands): identify new/updated genomes → QC → GTDB-Tk → resolve type material → representative changes/actions → select reps → cluster to named reps → de novo clusters → synonyms/erroneous NCBI → initial species names.
- Manual curation happens outside this tool.
- `pmc_*.py` ("post-manual curation"): final names, type species/strain checks, `pmc_validate`, stats, cleanup.
- Inspection/one-off commands (`merge_test`, `ani_sp_pair`, `intra_sp_derep`, etc.) take explicit file arguments instead of TOML.

**Steps communicate through files in per-step output dirs.** The TOML has the sections `data_files`, `ledgers`, `output_dirs` (keyed by subcommand name), `ltp`, `release`, and `root_dir`. A later step finds an earlier step's output by joining `params['output_dirs'][<earlier step>]` with a hard-coded filename. Examples: `qc_passed.tsv` from `u_qc_genomes`, `untrustworthy_type_material.tsv` from `u_resolve_types`, `gtdb_clusters_de_novo.tsv` from `u_cluster_de_novo`. Renaming an output file breaks every downstream consumer, so grep `main.py` for the filename first. `main.parse_toml_file` resolves every `data_files`, `ledgers` and `output_dirs` path, plus `ltp.ltp_taxonomy_file`, against `root_dir`. Absolute paths are left unchanged. `ltp.ltp_dir` is deliberately not resolved, because it is a subdirectory name inside each genome directory. The remaining `os.path.join(params['root_dir'], ...)` calls in the `pmc_*` methods are therefore no-ops. The bin script also writes a `gtdb_species_clusters.log` into the command's output dir.

**Per-domain and gzipped inputs.** The previous release's metadata and the NCBI GenBank assembly summaries are split by domain in the TOML: `prev_gtdb_{ar,bac}_metadata_file` and `ncbi_assembly_summary_{ar,bac}_genbank_file`. `main.py` combines each pair into a list with `prev_gtdb_metadata_files(params)` or `ncbi_assembly_summary_genbank_files(params)`. That list is passed through unchanged, so parameters named `*_file` in the step classes may hold a list. Read these files with `genome_utils.MultiFileReader`, which takes one file or a list, concatenates them, emits the header once and checks the headers match. Read other single inputs with `genome_utils.open_file`. Both readers handle `.gz` transparently.

**Ledgers** are curator-maintained TSVs that override NCBI/LPSN data. Examples: `gtdb_type_strains_ledger`, `untrustworthy_type_ledger`, `sp_priority_ledger`, `genus_priority_ledger`, `ncbi_env_bioproject_ledger`, `qc_exception_ledger`, `species_classification_ledger`, and `specific_epithet_ledger_{ar,bac}`. They are applied when genomes are loaded. `ncbi_env_bioproject_ledger` has an `Assembly type` column (`Environmental` or `Isolates`) read by `genome_utils.parse_ncbi_env_bioproject_ledger`: only Environmental BioProjects force genomes to be treated as MAGs/SAGs, while Isolates BioProjects only suppress the `u_rep_actions` warning about BioProjects replacing many representatives with presumed isolates. In `u_resolve_types` and `u_rep_changes`, the species in `gtdb_type_strains_ledger` takes the place of a genome's NCBI species, and ledger genomes are never untrustworthy as type material. `Genome.is_gtdb_type_strain_ledger()` identifies these genomes; loading marks them with the type designation source `GTDB curator`.

**Core data model:**
- `genomes.py:Genomes` is a dict-like container of `Genome` objects. It is built with `load_from_metadata_file(...)`, which takes the GTDB metadata TSV plus optional QC-passed file and ledgers, and also builds species clusters. Most steps load both a `prev_genomes` (last release) and a `cur_genomes` (this release).
- `genome.py:Genome` is a dataclass holding the type-material logic (`is_gtdb_type_strain`, `is_effective_type_strain`, ...), QC (`pass_qc` with `QcCriteria`), and representative scoring (`score_assembly`, `score_type_strain`). Rep selection is driven by these scores.
- `species_clusters.py:SpeciesClusters` maps a representative ID to its member genome IDs.
- `taxa.py` and `taxon_utils.py` handle rank-prefixed taxonomy strings (`d__;...;s__`). The `*_manager.py` modules handle nomenclature: name priority, specific epithets, taxon suffixes, and NCBI species handling.
- `defaults.py` holds the clustering/QC thresholds: ANI 95% / AF 50% species, 97% synonyms, 99% same strain, skani prefilter 85%, and the QC cutoffs used as CLI defaults.
- Shigella species are synonyms of E. coli under the GTDB (`Defaults.SHIGELLA_NCBI_SPECIES`). Their genomes are never selected as representatives by `u_rep_actions` or `u_sel_reps`, and `u_cluster_named_reps` and `u_cluster_de_novo` log an error unless all E. coli and Shigella type strain genomes are in a single cluster not represented by a Shigella genome (`type_genome_utils.check_ecoli_shigella_type_strains`).

**ANI/AF calculation** goes through `skani.py:Skani`. It wraps the `skani` binary and runs it in a multiprocessing pool. Results are cached in an SQLite table `ani_table` that is loaded into memory, and `ani_cache_file` is passed to the inspection commands. ANI/AF dicts are `ani_af[qid][rid] = (ani, af_r, af_q)`; use `Skani.symmetric_ani_af` for order-independent lookups. `Skani.pairs` runs each pair in the pool and waits for all of them, so calling it once per species leaves most CPUs idle when species have few genomes. Calculate all pairs in one call first, which fills the in-memory cache, and then process species (see `RepActions.calculate_rep_anis`). The `mash.py` wrapper is used only by `intra_sp_derep`. `u_gtdbtk` shells out to `gtdbtk classify_wf` in batches. Several instances can run at once on different servers that share the output dir. An instance claims a batch by atomically creating `canary.running` (host, PID, start time) in `gtdbtk_batch<N>/`. It writes `done` when GTDB-Tk succeeds and then removes the canary. Other instances skip batches that have either file. A canary left on the same host by a dead PID is reclaimed. Results are combined into `gtdbtk_classify.tsv` only once every batch is done. Batch dirs from before 0.7.4 count as done if they have a summary file. GTDB-Tk writes `gtdbtk.{ar53,bac120}.summary.tsv` only for domains present in the batch, so either file may be missing.

**Debug toggles for resuming expensive steps.** Several steps wrap a costly computation in `if True:  # ***DEBUGGING` and save the result to a `.pkl` in the output dir. The `else:` branch reloads that pickle. See `update_select_reps.py`, `update_cluster_named_reps.py`, `update_cluster_de_novo.py`, `update_rep_actions.py`, and `pmc_cluster_stats.py`. Developers flip these to `False` by hand to rerun from cached results. Leave them as `True` in committed code.

`sandbox.py` is a scratch command for one-off exploration, not part of the release pipeline.
