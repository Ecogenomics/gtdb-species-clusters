# SeqCode taxonomy support: implementation plan

Status: plan agreed 2026-11-03; PRs 1-9 merged (0.9.10), steps 9-12 verified on the R237 run; PR 10 (0.9.11) in review. Tick the PR checklist at the end as work lands.
Baselines for verification: the pre-SeqCode R237 step outputs are in `sp_cluster_update/deprecated/<N>_<step>/`.

## 1. Problem

The toolkit treats the NCBI taxonomy (`Genome.ncbi_taxa`) as the only external classification of a genome.
GTDB now registers genomes as type material under SeqCode, and the R237 metadata flags them as such
(`gtdb_type_designation_ncbi_taxa = type strain of species`, source `SeqCode`;
`gtdb_type_species_of_genus = t` from `seqcode_type_species_of_genus`), but their *names* still come from
NCBI. 1,193 of the 1,783 SeqCode genomes have no NCBI species (`s__`) and 62 have an NCBI species that
differs from the SeqCode name. Consequences seen in R237:

- `u_species_init` exits on the first type species genome with no NCBI species (G963970155, *Acidiparvus
  lacustris*); 614 representatives hit this.
- `u_resolve_types` never checks the 1,193 SeqCode-only type strains (skipped as `s__`).
- `u_synonyms` writes synonym rows with `s__` as the name, which `u_species_init` then reads back as a
  synonym of every unnamed cluster.
- SeqCode type genomes lose every naming-priority contest (no priority year, `s__` name) in
  `u_rep_actions`, `u_sel_reps` and `u_synonyms`.
- `u_species_init` leaves ~580 SeqCode-named clusters with placeholder names until `pmc_species_names`.

## 2. Decisions (agreed 2026-11-03)

| Topic | Decision |
|---|---|
| Source of SeqCode data | `data_files.seqcode_file` (`seqcode_table.tsv`) for **both** `prev_genomes` and `cur_genomes`. The cur metadata `seqcode_*` columns are cross-checked and any disagreement is logged as a warning. |
| Previous release | Only `ref_taxa` of previous genomes is populated from the table. `prev_genomes[*].gtdb_taxa` and the R232 type-designation flags stay exactly as in the R232 metadata: comparisons of GTDB classifications between releases must keep seeing the R232 classification. |
| Ranks | Full lineage: `seqcode_classification` (d__ to g__) plus the species from `seqcode_name`. |
| Resolved accessor | `Genome.ref_taxa` returns `seqcode_taxa` when present, else `ncbi_taxa`. `Genome.ref_taxa_source` is `'SeqCode'` or `'NCBI'`. `ncbi_taxa` stays for code that must mean NCBI specifically. |
| Type strains ledger | Ledger wins over SeqCode. `ref_taxa` stays SeqCode > NCBI; a new `Genome.type_species()` returns the ledger species for ledger genomes, else `ref_taxa.species`, and replaces the per-step special cases. Loading warns on each ledger/SeqCode conflict and on each ledger row that is redundant (changes nothing). |
| `ncbi_untrustworthy_sp` ledger | Unchanged: a genome in it is never an effective type strain, whatever its SeqCode name. |
| Priority year | `year_of_priority()` returns the earliest of the LPSN and SeqCode years when the name is validly published under both codes (SeqCode Rule 23c: priority dates from the earliest valid publication under either code), else whichever year exists. Decided 2026-11-03, replacing the earlier "SeqCode year first" rule. (SeqCode dates are validation dates: later than the LPSN year for 181 of the 380 genomes that have both, so for those the LPSN year is used.) |
| Output tables | Keep existing `NCBI species`/`NCBI genus` headers. Append `SeqCode species`, `SeqCode genus` and `Reference source used` (`NCBI` or `SeqCode`) columns at the end of the row. |
| Curation trees | `u_curation_trees` uses `ref_taxa` for prev and cur, so taxa new to SeqCode since the last release are highlighted. |
| Packaging | PR 1 is the core: data model, loading, accessors, warnings, tests; no step reads the new accessors, so outputs are unchanged. Then one PR per pipeline step, in pipeline order, so each step can be re-run and verified before the next is switched. |

## 3. Data facts (R237, for verification)

- `seqcode_table.tsv`: 1,783 rows, all `seqcode_rank = species`, all `seqcode_species_status = Valid (SeqCode)`,
  one row per accession, no species name listed twice. Accessions are versioned and mix GCA/GCF (203 rows are
  GCF where the metadata has GCA); `canonical_gid` reconciles them. `seqcode_classification` always ends in an
  empty `s__`; its genus always equals the generic name of `seqcode_name`; domain always equals the GTDB domain;
  `seqcode_name` never carries `Candidatus`. Genus status is `Valid (SeqCode)` for 1,549 and `Valid (ICNP)` for 234.
- Same 1,783 genomes carry `seqcode_*` columns in `gtdb_r237_metadata.tsv.gz`; names agree 100%. R232 metadata
  has no `seqcode_*` columns.
- Type flags already set by the metadata pipeline: 1,401 `type strain of species | SeqCode`, 381
  `type strain of species | LPSN; SeqCode`, 1 `not used as type | none`.
- NCBI names: 1,193 have no NCBI species; of the 590 that do, 528 agree with SeqCode, 10 differ in epithet only,
  52 differ in genus (27 same epithet, 25 both). 28 have an NCBI genus but no species.
- 1,717 are R232 representatives; 66 are not. 8 fail R237 QC. 1,539 MAGs, 9 SAGs, 235 isolates.
- Ledgers: 99 `gtdb_type_strains` rows are SeqCode genomes (97 same name; 2 differ: G001534665
  `Arsenophonus_A lipoptenae` vs `Arsenophonus lipoptenae`, G002072215 `Methanobrevibacter arboriphilicus` vs
  `Methanobinarius arboriphilus`); 41 `species_classification` rows (all same name); 1 in
  `ncbi_untrustworthy_sp_assignments`; 1 in `gtdb_untrustworthy_type_genomes`.
- Priority years: 1,403 have no LPSN year. SeqCode priority years: 2022 (18), 2023 (80), 2024 (315), 2025 (1,295),
  2026 (75).
- Of the 890 `seqcode_type_species_of_genus` genomes, 271 clusters already carry the SeqCode name, 576 have a
  placeholder, 26 a different Latin name, 12 are new unnamed clusters, 5 are unclustered (details in the
  2026-11-02 analysis; 542 of the placeholders are Sereika et al. 2025 names).

## 4. Core design (PR 1, version 0.9.0)

### 4.1 `taxa.py`

No change to `Taxa`. Add a helper `seqcode_taxa(classification_str, species_name)` (or do this in the loader)
that parses the classification string and sets `species = 's__' + species_name`.

### 4.2 `genome.py`

New fields (dataclass, after `ncbi_unfiltered_taxa`; all callers of `Genome(...)` are in `genomes.py` and tests):

```python
seqcode_taxa: Optional[Taxa]        # None when the genome is not SeqCode type material
seqcode_priority_year: int          # Genome.NO_PRIORITY_YEAR when not SeqCode
```

New methods/properties:

```python
@property
def ref_taxa(self) -> Taxa:            # seqcode_taxa if not None else ncbi_taxa
@property
def ref_taxa_source(self) -> str:      # 'SeqCode' or 'NCBI'
def is_seqcode_type_strain(self):      # seqcode_taxa is not None and not gtdb_untrustworthy_as_type
def type_species(self) -> str:         # ledger species if is_gtdb_type_strain_ledger() else ref_taxa.species
def year_of_priority(self):            # min(lpsn_priority_year, seqcode_priority_year if is_seqcode_type_strain())
def ref_taxa_columns(self) -> str:     # '\t'.join([seqcode species or '', seqcode genus or '', ref_taxa_source])
```

`year_of_priority()` keeps its current body in PR 1 (see §5, step 7) so PR 1 changes no output.

`is_gtdb_type_strain()` / `is_effective_type_strain()` stay metadata-driven; `is_seqcode_type_strain()` is an
additional source check used by `u_resolve_types` and reports. For current genomes the two agree
(1,782 of 1,783).

`__str__` gains `ref_sp` and the source.

### 4.3 `genomes.py`

- `load_from_metadata_file(..., seqcode_file=None)`: parse `seqcode_table.tsv` once
  (`genome_utils.parse_seqcode_table(path) -> {gid: SeqCodeRecord}` with `lru_cache`, like the NCBI assembly
  parser); keys are `canonical_gid`. For each genome in the table build `seqcode_taxa`, set
  `seqcode_priority_year` from `seqcode_priority_date[:4]`.
- If the metadata has `seqcode_name` (current release), compare with the table: warn per genome when the
  metadata has a SeqCode name and the table does not (or vice versa), or the names differ. Log counts:
  genomes with SeqCode taxa, with both NCBI and SeqCode species, with differing names.
- Ledger warnings (only when `gtdb_type_strains_ledger` given): for each ledger genome with `seqcode_taxa`,
  warn `conflict` when the ledger species differs from the SeqCode species, and warn `redundant` when the
  ledger species equals the SeqCode species (ledger row changes nothing). Also warn `redundant` when the ledger
  species equals the NCBI species *and* the metadata already marks the genome as a type strain.
- Replace the three ad-hoc "ledger species" look-ups (`u_resolve_types`, `u_ncbi_erroneous`, and the loader's
  own grouping in `ncbi_sp_effective_type_genomes()`) with `type_species()` **in the step PRs**, not here.
- New helpers, used by step PRs: `named_ref_species()` (same as `named_ncbi_species()` on `ref_taxa`) and
  `type_strain_genomes_by_species()` (groups effective type strains by `type_species()`). The old helpers stay
  until their last caller is switched, then are removed.

### 4.4 `main.py`

- Every `u_*` command that loads genomes gets `seqcode_file` from `params['data_files']` (already in the TOML,
  used today only by the pmc commands), passes it to both `prev_genomes` and `cur_genomes` loads, and adds
  `'seqcode_file'` to its `log_command_inputs` list. Step classes gain a `seqcode_file` parameter.
- No TOML change needed (`seqcode_file` exists under `[data_files]`).

### 4.5 Tests

- `test_genome.py`: `ref_taxa`/`ref_taxa_source` with and without SeqCode; `type_species()` for ledger, SeqCode
  and NCBI genomes; `year_of_priority()` for SeqCode/LPSN/none; `is_seqcode_type_strain()` respects
  `gtdb_untrustworthy_as_type`.
- `test_genomes.py`: load a small metadata file plus a SeqCode table (GCF accession in the table, GCA in the
  metadata); check `seqcode_taxa` lineage and species; conflict/redundant ledger warnings (`caplog`).
- Existing tests must pass unchanged.

### 4.6 Verification of PR 1

Re-run `u_qc_genomes` on R237 into a scratch directory with the new code: every output file must be
byte-identical to `2_u_qc_genomes/`. The log gains the SeqCode loading lines and ledger warnings only.

## 5. Per-step changes and behaviour

Conventions for every step PR:

- Replace `ncbi_taxa.species/genus/...` with `ref_taxa.*` unless the code must mean NCBI specifically
  (`is_ncbi_subspecies`, `is_ncbi_type_strain`, `ncbi_unfiltered_taxa`, and places that explicitly report the
  NCBI name next to the resolved one).
- Where type strains are grouped or compared by species, use `type_species()`.
- Output tables: keep the existing header text; append `SeqCode species`, `SeqCode genus`,
  `Reference source used` as the last columns (`Genome.ref_taxa_columns()`). Readers that index by header name
  (`read_clusters`, `parse_synonyms_table`, `load_from_metadata_file` for `qc_passed.tsv`) are unaffected; grep
  the readers of any file whose columns change before changing it.
- Bump `VERSION`; PR body states the step, the behaviour change and the verification done.
- Verification = the step is run on R237 and its output diffed against the `deprecated/` baseline; every difference
  must be explained by the behaviour change listed below. From PR 8 on the user runs the step after the PR is
  merged and the comparison is done on that output (decided 2026-11-03).

### Step 1: `u_new_genomes` (no PR)

No taxonomy is used. No change.

### Step 2: `u_qc_genomes` (PR 2)

Changes: load with `seqcode_file`; `qc_passed.tsv`/`qc_failed.tsv` and the type-strain QC tables
(`update_qc_genomes.py:246, 268, 432, 476`) report `ref_taxa` and gain the three columns. The "species lost
all type genomes" tally groups by `type_species()`.

Behaviour: QC decisions are unchanged (the type-strain exception already uses `is_gtdb_type_strain()`, which is
true for SeqCode genomes). SeqCode genomes that fail QC (8 in R237) are now listed under their SeqCode species
instead of `s__`. `qc_passed.tsv` is read downstream by its first column only.

### Step 3: `u_gtdbtk` (no PR)

No taxonomy is used. No change.

### Step 4: `u_lpsn_rna_types` (PR 3)

Changes: genomes are grouped under both their reference species and their NCBI species (when different) when
matching LPSN type 16S sequences, so a SeqCode genome is still found under the LPSN basonym (e.g.
*Methanobrevibacter smithii* for the SeqCode name *Methanocatella smithii*); `gtdb_type_strains[...]` is
grouped the same way using `type_species()`; output table gains the columns.

Behaviour: the 62 genomes whose SeqCode name differs from NCBI, and the 1,193 with no NCBI species, are now
also searched under their SeqCode species. They only match if LPSN lists that name with a type 16S sequence, so
expect only additional rows, never removed ones. The "add N genomes to the type strains ledger" warning is
unaffected (SeqCode genomes are already GTDB type strains).

### Step 5: `u_resolve_types` (PR 4)

Changes:

- `sp_with_mult_type_strains` (`update_resolve_types.py:700-712`) groups by `type_species()`; drop the local
  ledger special case.
- `resolve_by_ncbi_types` (line 377) becomes "resolve by source type designation": a genome counts as type
  material at its source when `is_ncbi_type_strain()` **or** `is_seqcode_type_strain()`; the reason string
  becomes "Not classified as type material at NCBI or SeqCode".
- LTP conflict checks and report tables use `ref_taxa` and gain the columns.
- `check_ledger_conflicts` keeps checking the ledgers; the SeqCode conflict/redundancy warnings come from the
  loader.

Behaviour: SeqCode-only type strains (1,193) are no longer skipped. Each SeqCode species has one genome in the
table, so new multi-type-strain species arise only where a SeqCode genome shares `type_species()` with an
NCBI/LPSN type strain of the same name (possible for the 380 names LPSN also lists). Those species get ANI
checks and may add rows to `untrustworthy_type_material.tsv`; without the `is_seqcode_type_strain()` change,
every SeqCode genome without NCBI type annotation would be declared untrustworthy, which must not happen.
Review every new/removed row of `untrustworthy_type_material.tsv` and `type_strain_ani.tsv`.

### Step 6: `u_rep_changes` (PR 5)

Changes: species-change detection (`update_rep_changes.py:165-175`) compares `prev.ref_taxa` with
`cur.ref_taxa` (prev gets SeqCode names from the same table). Detailed message says which source changed
("Reference species reassigned from X (NCBI) to Y (SeqCode)"). Keep the ledger exemption. Type-strain
gain/loss logic unchanged (it reads the metadata type flags, which for prev are R232's).

Behaviour: representatives whose NCBI name changed but whose SeqCode name is stable are now UNCHANGED instead
of REASSIGNED; representatives whose NCBI species is `s__` in both releases but now have a SeqCode name are
UNCHANGED (both sides resolve to SeqCode). Expect `rep_change_summary.tsv` to differ only in the
`NCBI_SPECIES_CHANGE` column for SeqCode genomes.

### Step 7: `u_rep_actions` (PR 6)

Changes:

- `year_of_priority()` switches to the earliest of the LPSN and SeqCode years (§2, §4.2). This is the first step that uses it.
- `action_naming_priority` (lines 834-975): group type strains by `type_species()`, compare with
  `cur_genomes[updated_rid].type_species()`; `violates_naming_priority` compares `ref_taxa.specific_epithet`.
- `rep_candidates` Shigella check uses `ref_taxa` (identical result).
- `species_priority_manager.species_priority()` and `species_priority_year()` use `ref_taxa`; the
  `lpsn_valid_names` test is extended so a SeqCode species counts as a valid name (SeqCode names are validly
  published under SeqCode).
- `action_log.tsv` parameters keep `prev_ncbi_species`/`new_ncbi_species` keys and add
  `prev_species_source`/`new_species_source`; `updated_naming_priority.tsv` and the untrustworthy-at-NCBI
  table gain the columns.

Behaviour: SeqCode type genomes take part in naming priority under their SeqCode name and year. Where a
cluster's representative is not type material and a SeqCode type genome is a member, the representative is
replaced (`NOMENCLATURE_PRIORITY:REPLACED`). 1,717 SeqCode genomes are already representatives, so this mostly
affects the 66 non-representatives. Contests involving the 1,403 SeqCode genomes with no LPSN year now use
the SeqCode year instead of no year; the 380 with both keep the earlier (usually LPSN) year. Review every
`NOMENCLATURE_PRIORITY` row that differs.

### Step 8: `u_sel_reps` (PR 7)

Changes: `unrepresented_ncbi_sp`, `named_ncbi_species`, the `gtdb_type_sp`/`ncbi_type_sp` tabulation (lines
1081-1098), `select_rep`, ANI-neighbour resolution and all report tables use `ref_taxa`/`type_species()`;
tables gain the columns. `validate_type_designations` reports SeqCode as a type source where applicable.

Behaviour: SeqCode species whose genome has no NCBI species now count as named species with a type strain.
Those not yet represented get a representative selected under the SeqCode name (12 new clusters and 5
unclustered genomes in R237 instead of waiting for de novo clustering). ANI-neighbour exclusion
(`gtdb_excluded_ani_neighbours.tsv`) uses the SeqCode year and name. `validate_type_strain_of_species.tsv`
lists SeqCode species.

### Steps 9-11: `u_cluster_named_reps`, `u_cluster_de_novo`, `u_cluster_stats` (PR 8, one PR)

Changes: only reporting via `type_genome_utils.write_clusters`/`write_rep_radius` (lines 236, 372) and
`check_ecoli_shigella_type_strains` (`ref_taxa`, identical result). Cluster tables gain the columns.

Behaviour: clustering unchanged. `gtdb_clusters_de_novo.tsv` and `gtdb_named_rep_clusters.tsv` are read by
`read_clusters` (header-indexed) so the extra columns are safe.

### Step 12: `u_ncbi_erroneous` (PR 9)

Changes: species of every genome from `ref_taxa`; anchoring type strains by `type_species()` (replaces the
local ledger species look-up); table gains the columns. The fallback that lets a ledger genome also anchor its
NCBI species when no other type strain does (added in 4a22c6d for ledger rows that correct spelling or add a
genus suffix) is kept and, from 0.9.10, applies to every type strain genome: its `ref_taxa` and NCBI species
also anchor to its cluster when they differ from `type_species()` and have no other type strain genome. The `NCBI species` column keeps the NCBI name;
the appended columns give the SeqCode name the genome was grouped by.

Behaviour: SeqCode type genomes anchor their SeqCode species. Any genome whose reference species equals that
name but sits in another cluster is now reported as misclassified (previously these species had no anchor).
Genomes whose NCBI species was wrong but whose SeqCode name is right stop being flagged for that name.
In 0.9.9 the R237 run dropped 76 rows because SeqCode type genomes with a different SeqCode name (e.g.
*Methanobrevibacter arboriphilus* → *Methanobinarius arboriphilus*, *Pelagibacter communis* → *P. ubiqueversans*)
no longer anchored their NCBI name; the 0.9.10 fallback restores them, since genomes keeping that NCBI name in
other clusters must still be excluded from NCBI-name votes in `u_synonyms` and `u_species_init`.

### Step 13: `u_synonyms` (PR 10)

Changes in `ncbi_species_manager.py`: type strain genomes are grouped by `type_species()`
(`Genomes.type_strain_genomes_by_species()`, replacing `ncbi_sp_effective_type_genomes()`, which is removed);
`identify_type_strain_synonyms` names type strain genomes by `type_species()`, and the representative's
`ref_taxa` species is never a synonym of its `type_species()`; `identify_consensus_synonyms` groups genomes by
`ref_taxa` and skips species that are the `type_species()` or `ref_taxa` species of a type strain genome.
`write_synonym_table` keeps the NCBI species columns and appends `Reference species`, `Reference source used`,
`Synonym reference species` and `Synonym reference source used`, where the reference species is the
`type_species()` of type strain genomes (source `GTDB curator` when it is the ledger species) and the `ref_taxa`
species otherwise; the `Warnings` field is now always written. `parse_synonyms_table` reads the reference
species columns (NCBI columns for older tables) and skips, with a warning, rows where either name is `s__`.
`classify_ncbi_species` and `ncbi_sp_gtdb_cluster_table` are only used by the `pmc_*` commands and are switched
with them.

Behaviour: the four R237 rows with `s__` disappear or get real names (the Palauibacter row vanishes because
representative and synonym are type material of the same species; Gammaporibacter and Caproiciproducens rows
show their SeqCode names; the Kariarchaeum row names the ledger species *Kariarchaeum pelagium* of the
representative). SeqCode type genomes in another type genome's cluster become `TYPE_STRAIN_SYNONYM`s under
their SeqCode name, with the SeqCode year in the priority check. The NCBI species of a type genome renamed
under the SeqCode (e.g. *Pelagibacter communis*) is no longer a type strain species, so it can become a
`MAJORITY_VOTE_SYNONYM` of the SeqCode name when all its correctly classified genomes are in that cluster.
The two "incorrect priority" warnings may change.

### Step 14: `u_curation_trees` (PR 11)

Changes: `new_ncbi_taxa` uses `ref_taxa` for prev and cur at every rank; output files keep their names, the
table gains a `Reference source used` column.

Behaviour: SeqCode taxa (genus and above) that were not in the previous release's reference taxa are listed
and placed in the curation trees. Since prev genomes also get SeqCode taxa from the same table, only taxa
whose genomes are new to GTDB, or newly registered in SeqCode, appear.

### Step 15: `u_species_init` (PR 12)

Changes:

- `_sort_by_naming_priority`, `name_type_species_clusters`, `name_type_strain_clusters`,
  `name_binomial_clusters`, `_determine_sp_epithet` and `update_species_names` use `type_species()` for type
  material and `ref_taxa` elsewhere. The `sys.exit` at line 189 becomes a warning for genomes that are type
  species of a genus but still resolve to `s__` (none expected in R237).
- `species_priority_manager.genus_priority_year` falls back to the SeqCode year of the genus's type species
  when LPSN has no entry for the genus.
- Curation and update logs record the source of the proposed name.

Behaviour: this is the largest change. Clusters represented by SeqCode type species/type strains are named
from the SeqCode name at initialisation instead of after manual curation: in R237, ~576 placeholder clusters
get their SeqCode names, 4 suffixed names are checked against the SeqCode name, and the 26 clusters whose
current Latin name differs from SeqCode are written to `conflicting_type_species_of_genus.tsv` /
the curation log for a decision. Genus names for the 563 placeholder-genus clusters come from
`seqcode_classification`. The `species_classification_ledger` is not read here, so no conflict with its 41
SeqCode rows.

### Post-curation commands (PR 13+, after the `u_*` pipeline is verified)

- `pmc_species_names`: keep the existing `parse_seqcode_classification`; `set_type_species`/`set_type_strain`
  and the NCBI classification use `ref_taxa`/`type_species()`. The "final taxonomy does not reflect SeqCode
  specific name" warning stays.
- `pmc_validation` (56 `ncbi_taxa` uses): genus/family/domain consistency checks and type-material checks
  use `ref_taxa`; the three `seqcode_species.get(rid, '')` comparisons use `seqcode_taxa.species`.
- `pmc_check_type_species`: compare GTDB genus with `ref_taxa.genus`; genomes whose reference genus is empty
  are skipped (removes the noise rows for SeqCode genomes with no NCBI genus).
- `pmc_check_type_strains`: compare with `ref_taxa.species`; drop the `s__` skip for SeqCode genomes.
- `pmc_cluster_stats`, `pmc_checks`, `pmc_cleanup`: no taxonomy use beyond reporting; add columns where
  `NCBI species` is written.

## 6. Open points to settle when the step is reached

- Step 5: decided 2026-11-03: a genome that is type material via SeqCode is treated exactly like one that is type
  material via LPSN; neither source gets a ledger-like override, since neither is fully reliable. Conflicts between
  SeqCode and LPSN/NCBI type genomes (e.g. *Aalborgiella aggregata*, *Planktophila lacus*) go to curation via the
  untrustworthy type ledger.
- Step 7: decided 2026-11-03: SeqCode names count as valid names in `species_priority`.
- Step 15: whether the 26 Latin-name conflicts should be auto-renamed to SeqCode or always left for curation
  (plan: curation log only, no automatic rename).

## 7. PR checklist

- [x] PR 1: core (0.9.0): fields, loader, `ref_taxa`, `type_species()`, warnings, tests, this document
- [x] PR 2: `u_qc_genomes` (0.9.1)
- [x] PR 3: `u_lpsn_rna_types` (0.9.2)
- [x] PR 4: `u_resolve_types` (0.9.3)
- [x] PR 5: `u_rep_changes` (0.9.4)
- [x] PR 6: `u_rep_actions` (0.9.5; `year_of_priority`, `species_priority_manager`)
- [x] PR 7: `u_sel_reps` (0.9.7)
- [x] PR 8: `u_cluster_named_reps`, `u_cluster_de_novo`, `u_cluster_stats` (0.9.8)
- [x] PR 9: `u_ncbi_erroneous` (0.9.9; NCBI-name fallback for all type strain genomes in 0.9.10)
- [x] PR 10: `u_synonyms` (0.9.11)
- [ ] PR 11: `u_curation_trees`
- [ ] PR 12: `u_species_init`
- [ ] PR 13+: `pmc_*`
