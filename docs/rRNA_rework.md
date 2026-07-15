# rRNA table rework — change log & debugging guide

This documents the rework of the rRNA PostgreSQL tables from a Human-only,
source-specific layout into a **generic, multi-species schema** where new
versions / species / rRNAs are just more rows (no new columns, no model changes).

Scope of the change: the Snakemake pipeline + PostgreSQL tables only.
The cross-species `rRNA_alignment_*` tables are **not** part of this rework.

---

## 1. Why

Before, the rRNA side was rigid:

- `rRNA_modifications` had three parallel position columns
  (`pos_snoRNABAse`, `pos_incarnato`, `pos_snOPY`) — one per source. Adding a
  4th source meant adding a column and changing every consumer.
- `conversion_18S` and `conversion_28S` were two near-identical wide tables
  (one column per source). Adding 5.8S or mouse meant new tables.
- `rRNAs` had no species, so it was implicitly Human-only.

After, everything is **long-format and keyed by a universal `pos_id`**:

- `rRNA_modifications` is version-agnostic (one `pos_id` per modification).
- `conversion_18S` + `conversion_28S` are merged into one `rRNA_conversion`.
- `rRNAs` carries `species` and `rRNA_type` (ref/alt).
- Adding mouse / 5.8S / a new source version = adding rows, nothing else.

---

## 2. Schema: before → after

### `rRNAs`

```diff
 CREATE TABLE "rRNAs" (
     rRNA_id varchar(50) NOT NULL PRIMARY KEY,
     rRNA_name varchar(50) NOT NULL,
+    species varchar(50) NOT NULL,          -- NEW: 'Human' for existing rows
+    rRNA_type varchar(30) NOT NULL,        -- NEW: 'ref' | 'alt'
     sequence varchar(7500) NOT NULL,
     sources varchar(200) NOT NULL,
-    modif_version varchar(30) NOT NULL
+    modif_version varchar(30) NOT NULL,
+    UNIQUE (species, rRNA_name, modif_version)
 );
```

### `rRNA_modifications`

```diff
 CREATE TABLE "rRNA_modifications" (
     modification_id SMALLSERIAL NOT NULL PRIMARY KEY,
-    rRNA_name varchar(20) NOT NULL,
-    pos_snoRNABAse INTEGER NOT NULL,
-    pos_incarnato INTEGER NOT NULL,
-    pos_snOPY INTEGER NOT NULL,
+    species varchar(50) NOT NULL,
+    rRNA_name varchar(20) NOT NULL,
+    pos_id INTEGER NOT NULL,               -- universal alignment position
     nucleotide varchar(1) NOT NULL,
     modification_type varchar(30) NOT NULL,
     modification_name varchar(30) NOT NULL
 );
```

The three `pos_*` columns are gone. A modification is now located by
`(species, rRNA_name, pos_id)`. Native per-version positions are recovered via
`rRNA_conversion` (see §4). This drops per-source attestation ("who reported
this mod") — accepted by design.

### `conversion_18S` + `conversion_28S`  →  `rRNA_conversion`

The two wide tables are replaced by one long table:

```sql
CREATE TABLE "rRNA_conversion" (
    id bigserial PRIMARY KEY,
    species varchar(50) NOT NULL,
    rRNA_name varchar(50) NOT NULL,         -- scopes pos_id (18S/28S/5.8S, per species)
    pos_id INTEGER NOT NULL,                -- universal alignment column
    pos INTEGER,                            -- native position in THIS version
    base varchar(1),                        -- base in THIS version
    modif_version varchar(30) NOT NULL      -- -> rRNAs (species + modif_version)
);
```

One row per **(species, rRNA_name, pos_id, modif_version)**. The distinct set of
`(species, rRNA_name, pos_id)` is the alignment grid.

---

## 3. What `pos_id` is (and is not)

`pos_id` is the **column of a multiple-version alignment**, scoped per
`(species, rRNA_name)`. It is NOT:

- **not** the snoRNABase position — alt-only insertions shift it. E.g. the
  Human 18S modification `Um1668` (snoRNABase position 1668) has `pos_id = 1670`
  because 2 snOPY-specific insertions precede it.
- **not** global — `pos_id` is only meaningful together with `(species, rRNA_name)`.

The reference (`rRNA_type = 'ref'`, i.e. `modif_version = 'snoRNABase'`) is the
backbone: its positions run 1..N, and alt-only insertions are slotted between them.

### The threading algorithm (in `scripts/build_rRNA_conversion.py`)

The source conversion CSVs are wide tables where **each row aligns corresponding
positions across versions** (`NA` where a version is absent), but the rows are
**not** in biological order (28S is sorted by Incarnato; both files append orphan
rows at the end). The algorithm reorders them:

1. **Backbone**: take the reference (`snoRNABase`) positions in ascending order.
2. **For each alt-only insertion** (a row with no reference position but a
   position in some alt version): find, **for that alt version**, the largest
   alt position `< this one` that *does* map to a reference position; the
   insertion is placed immediately after that reference position.
3. Tiebreak insertions sharing the same anchor by their alt position (deterministic).
4. Assign `pos_id = 1..M` (1-based) in that order.
5. Emit one long row per **(column, version present at that column)**:
   `(species, rRNA_name, pos_id, modif_version, pos, base)`.

Reference columns always sort before insertions sharing the same anchor, so an
insertion lands strictly between `ref_pos` and `ref_pos + 1`.

This is implemented as a pure, testable function `compute_alignment(df, versions,
species, rRNA_name)`; see §8 for how to exercise it.

---

## 4. Data flow (compute once, load twice)

`rRNA_modifications.pos_id` must agree with `rRNA_conversion.pos_id`, so the
alignment is computed **once** into an intermediate, then both loaders consume it:

```
conversion_18S.csv ┐
                   ├─ build_rRNA_conversion ─→ data/rRNA_processed/rRNA_conversion_long.tsv
conversion_28S.csv ┘            │                       (species, rRNA_name, pos_id,
                                │                        modif_version, pos, base — WITH header)
              ┌─────────────────┴──────────┐
              ▼                            ▼
   psql_rRNA_conversion          psql_rRNA_modifications
   (load verbatim)               (join ref rows: modif_version='snoRNABase',
                                  pos == pos_snoRNABase → pos_id)
```

- The intermediate TSV is written **with a header** so it's inspectable.
- `rRNA_modifications` left-joins the ref rows; rRNAs with **no** conversion
  source (5.8S is ref-only and has no conversion CSV) fall back to
  `pos_id = pos_snoRNABase` (their own position — correct, since there are no
  alt insertions to shift anything).

---

## 5. File-by-file change log

| File | Action | Notes |
|---|---|---|
| `scripts/build_rRNA_conversion.py` | **NEW** | Threading algorithm; writes the intermediate long TSV. Pure pandas, **no docker**. Core fn `compute_alignment()`. |
| `scripts/psql_rRNA_conversion.py` | **NEW** | Loads the intermediate → `rRNA_conversion` (bigserial `id` skipped in `\copy`). |
| `scripts/psql_rRNAs.py` | MODIFIED | Adds `species='Human'`, `rRNA_type` (`ref` iff `modif_version=='snoRNABase'` else `alt`), renames `references`→`sources`, explicit column order before `to_csv`. Keeps T→U + modif_version logic. |
| `scripts/psql_rRNA_modifications.py` | MODIFIED | Drops the 3 `pos_*` columns; adds `species`; renames `rRNA`→`rRNA_name`, `modif_name`→`modification_name`; derives `pos_id` by joining the intermediate ref rows (with fallback). |
| `rules/psql.smk` | MODIFIED | Removed `psql_conversion_18S` + `psql_conversion_28S`; added `build_rRNA_conversion` + `psql_rRNA_conversion`; `psql_rRNA_modifications` now takes the intermediate as input; grant list swaps conversion_18S/28S → rRNA_conversion. |
| `config.yaml` | MODIFIED | `path.rRNA_processed: data/rRNA_processed`; `rRNA_modifications.rRNA_conversion_long: rRNA_conversion_long.tsv`; `psql:` drops `conversion_18S`/`conversion_28S`, adds `rRNA_conversion`. The input CSV keys `conversion_18S`/`conversion_28S` under `rRNA_modifications:` are **kept** (the build rule still reads them). |
| `Snakefile` | MODIFIED | `rule all`: `conversion_18S_psql`/`conversion_28S_psql` → `rRNA_conversion_psql`. |
| `scripts/psql_conversion_18S.py` | **DELETED** | Replaced. |
| `scripts/psql_conversion_28S.py` | **DELETED** | Replaced. |

---

## 6. Source data quirks (things that bite)

The inputs under `data/Hermes_data/processed/` are hand-curated and have quirks
the code now accounts for:

- **`conversion_18S.csv`** (1,871 rows): `snoRNABase_pos` dense 1..1869, then
  **2 snOPY-only rows appended at the end** (`NA,NA,1648` and `NA,NA,1640`).
  Missing = literal `NA`. snOPY numbering is snoRNABase + 2 past the insertions.
- **`conversion_28S.csv`** (5,075 rows): **sorted by Incarnato_pos**; contains
  **40 Incarnato-only insertions** (`snoRNABase_pos=NA`, interspersed) and a
  **5-row tail of snoRNABase-only positions** (`3136,3188,3303,3317,3336` with
  `Incarnato_pos=NA`) appended after Incarnato ends. **No snOPY columns at all.**
- **`rRNA.csv`** (6 rows, Human): 5.8S (1 ref), 18S (ref snoRNABase + alt
  incarnato + alt snOPY), 28S (ref snoRNABase + alt incarnato). No mouse.
- **`all_rRNA_modifications.csv`** (218 rows, tab-delimited despite `.csv`):
  `pos_snoRNABase` is **fully populated (0 NA)**, no duplicate
  `(rRNA, pos_snoRNABase)`. modification_type ∈ {`2'-O-methylation`,
  `pseudouridylation`, `both`}.

### Bugs fixed as a side effect

- The old `conversion_18S.base_zavolan` column was actually loaded from
  `snOPY_seq` (mislabeled; there is no "Zavolan" source). Gone — snOPY's base
  now correctly populates the `modif_version='snOPY'` rows of `rRNA_conversion`.
- The old wide `pos_snoRNABAse` typo (capital `BA`) is gone.

### NULL sentinels

Inconsistent across the codebase (not normalized here): conversions used `'-1'`,
rRNAs/modifications use `'.'`, specie uses `''`. The new `rRNA_conversion` loader
uses `'.'` for consistency with rRNAs/modifications. The build step drops
`NA`/missing rather than writing `-1`.

---

## 7. How to debug

### "The positions look wrong / off by N"
`pos_id` is an alignment column, not a source position — it is **shifted by
preceding alt insertions**. Do not compare `pos_id` directly to a snoRNABase
position. Always recover a version's native position via the round-trip join
(§9). To see the shift, inspect the intermediate:

```bash
# Um1668: snoRNABase pos 1668 -> pos_id 1670 (2 snOPY insertions precede it)
awk -F'\t' 'NR==1 || ($2=="18S" && $3>=1668 && $3<=1672)' \
    data/rRNA_processed/rRNA_conversion_long.tsv
```

### "A modification has the wrong/missing pos_id"
Check the build output and the join in `psql_rRNA_modifications.py`:

- The modification's `pos_snoRNABase` must exist in the intermediate's
  `modif_version='snoRNABase'` rows for that `(species, rRNA_name)`.
- For 5.8S (no conversion CSV) the join intentionally misses and `pos_id`
  falls back to `pos_snoRNABase`. That is correct.

### "Columns shifted after a load" (silent data corruption)
`\copy` is **positional** (the TSV has no header). If a `df.to_csv` column order
ever drifts from the `\copy (...)` list, data loads into the wrong columns with
no error. Every modified script builds an explicit column list right before
`to_csv`, kept beside the `\copy` list — keep them in sync on any future edit.

### "DAG won't resolve / MissingInputError"
Check `config.yaml`: `path.rRNA_processed` and
`rRNA_modifications.rRNA_conversion_long` must exist; `psql.rRNA_conversion`
must exist. The build rule inputs `conversion_18S`/`conversion_28S` CSV keys
must still be present.

### Inspect the generated SQL before loading
Each `data/psql/<table>/data_script.sql` is plain text — open it to confirm the
`CREATE TABLE` and `\copy` column list match.

### Exercise the algorithm in isolation
`compute_alignment` is a pure function; run it directly without snakemake:

```python
import sys; sys.path.insert(0, 'scripts')
import pandas as pd
from build_rRNA_conversion import compute_alignment, VERSIONS
df = pd.read_csv('data/Hermes_data/processed/conversion_18S.csv', na_values=['NA'])
long = compute_alignment(df, VERSIONS['18S'], 'Human', '18S')
print(long[long.pos_id.between(1668, 1672)])
```

---

## 8. Verification performed (all passed, DB load excluded)

These checks were run against the real data **without loading the DB** (the
load is deferred until the Django cutover — see `rRNA_django_migration.md`):

- **DAG**: `snakemake -n` resolves; `build_rRNA_conversion`, `psql_rRNA_conversion`,
  `psql_rRNA_modifications`, `psql_rRNAs` present; no `conversion_18S`/`28S`.
- **Alignment**: 18S = 1,871 columns (1,869 ref + 2 snOPY insertions); 28S =
  5,075 (5,035 ref + 40 Incarnato insertions); `pos_id` contiguous 1..M, no gaps.
- **Round-trip** (`Um1668`): ref pos 1668 → `pos_id 1670`; at 1670 the native
  positions recover snoRNABase=1668, incarnato=1668, **snOPY=1670** (matches the
  source row; also confirms the zavolan mislabel is fixed).
- **Threading**: 18S snOPY insertions land at ~ref 1639 (not at the end); 28S
  ref-only tail rows (3136–3336) back in biological order (pos_id 3153–3360).
- **rRNAs**: 6 rows; `rRNA_type` rule correct; `UNIQUE(species, rRNA_name,
  modif_version)` satisfied.
- **rRNA_modifications**: 218 rows (28S=128, 18S=86, 5.8S=4); no duplicate
  `(species, rRNA_name, pos_id)`; 5.8S `pos_id` ∈ {14, 55, 69, 75}.
- **rRNA_conversion**: 15,714 rows, counts per (rRNA_name, modif_version) as expected.

### Live-DB verification (run after the coordinated load)

```bash
source /home/danx/.zshrc && conda activate snakemake7.24
snakemake -c1 --use-conda \
    data/psql/rRNA_conversion/data_table.tsv \
    data/psql/rRNAs/data_table.tsv \
    data/psql/rRNA_modifications/data_table.tsv data/psql/permission.tok
```

```sql
-- counts
SELECT species, rRNA_name, rRNA_type, COUNT(*) FROM raw_data."rRNAs" GROUP BY 1,2,3;
SELECT rRNA_name, COUNT(*) FROM raw_data."rRNA_modifications" GROUP BY 1;
SELECT rRNA_name, modif_version, COUNT(*) FROM raw_data."rRNA_conversion" GROUP BY 1,2;

-- round-trip: Um1668 -> each version's native position
SELECT c.modif_version, c.pos AS native_pos, c.base
FROM raw_data."rRNA_modifications" m
JOIN raw_data."rRNA_conversion" c
  ON c.species=m.species AND c.rRNA_name=m.rRNA_name AND c.pos_id=m.pos_id
WHERE m.species='Human' AND m.rRNA_name='18S' AND m.modification_name='Um1668'
ORDER BY c.modif_version;
-- expect: incarnato=1668, snOPY=1670, snoRNABase=1668
```

---

## 9. Follow-ups / known limitations

- **No mouse / 5.8S-alt input data yet.** The schema supports them; they'll be
  rows once a mouse `rRNA.csv` / conversion source is supplied. Ref-only rRNAs
  (5.8S today) correctly get no `rRNA_conversion` rows.
- **`rRNA_percentage_modification`** still uses the old `pos_snoRNABAse` column
  (typo, lowercase `a`). Untouched — it is now the only table left on the old
  per-source-pos model. Consider aligning it to `pos_id` later.
- **NULL sentinels** remain inconsistent across scripts (see §6).
- **`database_organisation/snoDB3_schema.svg`** is now stale and should be
  updated to show `rRNA_conversion` and the new `rRNAs`/`rRNA_modifications` columns.
- **Django cutover** — the downstream web project must migrate to the new schema
  in lockstep with the DB load; see `rRNA_django_migration.md`.
