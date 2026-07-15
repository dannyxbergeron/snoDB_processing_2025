# Django migration guide — rRNA schema rework

> **Audience:** the AI (future Claude) that will modify the **snoDB3 Django web
> project** (a separate repo) after the pipeline reworked the rRNA tables.
> This file lives in the *pipeline* repo (`snodb3`) and describes exactly what
> changed in the database so the Django side can be updated to match.
>
> **Companion doc:** `docs/rRNA_rework.md` — full pipeline change log and the
> `pos_id` threading algorithm. Read §3 there for *why* `pos_id` is what it is.

---

## 0. TL;DR — the mental model shift

The rRNA tables moved from **wide + source-specific** to **long + generic**:

- A modification is no longer stored with three per-source position columns.
  It now has a single **universal `pos_id`** (an alignment column).
- To get a modification's position **in a specific version's coordinates**, you
  join `rRNA_modifications` → `rRNA_conversion` on `(species, rRNA_name, pos_id)`
  and filter by `modif_version`.
- The two old `conversion_18S` / `conversion_28S` tables are gone, replaced by
  one `rRNA_conversion` table that works for any rRNA and any species.

**Golden rule:** `pos_id` is **not** a source position and **not** global. It is
shifted by alt-version insertions (e.g. Human-18S `Um1668` is at `pos_id = 1670`,
not 1668). Never assume `pos_id == snoRNABase position`. Always resolve a real
coordinate through `rRNA_conversion`.

---

## 1. Exact schema: before → after (schema `raw_data`)

### `rRNAs`  (PK `rRNA_id`)

| before | after |
|---|---|
| `rRNA_id, rRNA_name, sequence, sources, modif_version` | `rRNA_id, rRNA_name, **species**, **rRNA_type**, sequence, sources, modif_version` + `UNIQUE(species, rRNA_name, modif_version)` |

- `species` — `varchar(50)`, e.g. `'Human'` (matches the snoRNA `specie` table).
- `rRNA_type` — `varchar(30)`, `'ref'` or `'alt'`. The **reference** version
  (`rRNA_type='ref'`) is the one whose numbering anchored `pos_id`; currently
  that is always `modif_version='snoRNABase'`.

### `rRNA_modifications`  (PK `modification_id`)

| before | after |
|---|---|
| `modification_id, rRNA_name, pos_snoRNABAse, pos_incarnato, pos_snOPY, nucleotide, modification_type, modification_name` | `modification_id, **species**, rRNA_name, **pos_id**, nucleotide, modification_type, modification_name` |

- The three `pos_*` columns are **removed**. Replaced by `pos_id`.
- Note the old typo `pos_snoRNABAse` (capital `BA`) is gone.

### `conversion_18S`, `conversion_28S`  →  **DELETED**

Replaced by a single long table:

```sql
CREATE TABLE "rRNA_conversion" (
    id bigserial PRIMARY KEY,
    species varchar(50) NOT NULL,
    rRNA_name varchar(50) NOT NULL,
    pos_id INTEGER NOT NULL,
    pos INTEGER,            -- native position in THIS version
    base varchar(1),        -- base in THIS version
    modif_version varchar(30) NOT NULL
);
```

### Reference data (current contents)

- `rRNAs`: **6 rows**, all `species='Human'`:
  - 5.8S / ref / snoRNABase
  - 18S / ref / snoRNABase; 18S / alt / incarnato; 18S / alt / snOPY
  - 28S / ref / snoRNABase; 28S / alt / incarnato
- `rRNA_modifications`: **218 rows** (18S=86, 28S=128, 5.8S=4).
- `rRNA_conversion`: **15,714 rows** (18S: snoRNABase 1869, incarnato 1869, snOPY 1871; 28S: snoRNABase 5035, incarnato 5070). **5.8S has NO rows** (ref-only, no alt versions to convert).
- `modif_version` values in use: `'snoRNABase'`, `'incarnato'`, `'snOPY'`.
- `rRNA_name` values: `'5.8S'`, `'18S'`, `'28S'`.

---

## 2. Model changes (Django)

The pipeline loads the DB; Django reads it. Models almost certainly use
`managed = False` with explicit `db_table`. Update field definitions to mirror §1.

```python
class rRNA(models.Model):
    rRNA_id       = models.CharField(max_length=50, primary_key=True)
    rRNA_name     = models.CharField(max_length=50)
    species       = models.CharField(max_length=50)          # NEW
    rRNA_type     = models.CharField(max_length=30)          # NEW ('ref'|'alt')
    sequence      = models.CharField(max_length=7500)
    sources       = models.CharField(max_length=200)
    modif_version = models.CharField(max_length=30)
    class Meta:
        managed = False
        db_table = 'rRNAs'
        constraints = [models.UniqueConstraint(
            fields=['species', 'rRNA_name', 'modif_version'], name='unique_rRNA_version')]

class rRNAModification(models.Model):
    modification_id   = models.AutoField(primary_key=True)
    species           = models.CharField(max_length=50)      # NEW
    rRNA_name         = models.CharField(max_length=20)
    pos_id            = models.IntegerField()                # NEW (replaces the 3 pos_* cols)
    nucleotide        = models.CharField(max_length=1)
    modification_type = models.CharField(max_length=30)
    modification_name = models.CharField(max_length=30)
    class Meta:
        managed = False
        db_table = 'rRNA_modifications'

class rRNAConversion(models.Model):                          # NEW (replaces conversion_18S/28S)
    id            = models.AutoField(primary_key=True)
    species       = models.CharField(max_length=50)
    rRNA_name     = models.CharField(max_length=50)
    pos_id        = models.IntegerField()
    pos           = models.IntegerField(null=True)
    base          = models.CharField(max_length=1, null=True)
    modif_version = models.CharField(max_length=30)
    class Meta:
        managed = False
        db_table = 'rRNA_conversion'
        indexes = [models.Index(fields=['species', 'rRNA_name', 'pos_id']),
                   models.Index(fields=['species', 'rRNA_name', 'modif_version'])]
```

**Delete** the `conversion_18S` and `conversion_28S` models entirely.

> There are **no SQL `FOREIGN KEY` constraints** in this DB (the pipeline never
> declares `REFERENCES`). So there are no Django `ForeignKey` fields to define —
> all joins are explicit composite joins on `(species, rRNA_name, pos_id)` and
> `(species, rRNA_name, modif_version)`. Add the indexes above yourself; they
> are essential for the join-heavy queries in §3 to be fast.

---

## 3. Query rewrite patterns (the core of the migration)

Below, each pattern shows the **old** way (now broken) and the **new** way, as
SQL plus a Django ORM sketch. Adapt model/field names to the real codebase.

### Q1 — A modification's position in a given version's coordinates

Old: read `rRNA_modifications.pos_snOPY` directly.
New: join through `rRNA_conversion` and filter `modif_version`.

```sql
SELECT c.pos AS native_pos, c.base
FROM raw_data."rRNA_modifications" m
JOIN raw_data."rRNA_conversion" c
  ON c.species = m.species AND c.rRNA_name = m.rRNA_name AND c.pos_id = m.pos_id
WHERE m.species='Human' AND m.rRNA_name='18S'
  AND m.modification_name='Um1668'
  AND c.modif_version='snOPY';
-- -> native_pos = 1670
```

```python
pos = (rRNAConversion.objects
       .filter(species='Human', rRNA_name='18S', modif_version='snOPY',
               pos_id__in=rRNAModification.objects
                   .filter(species='Human', rRNA_name='18S', modification_name='Um1668')
                   .values('pos_id'))
       .values_list('pos', 'base'))
```

### Q2 — All modifications for an rRNA, in a chosen version's coordinates

```sql
SELECT m.modification_name, m.nucleotide, m.modification_type, c.pos AS pos_in_version
FROM raw_data."rRNA_modifications" m
JOIN raw_data."rRNA_conversion" c
  ON c.species=m.species AND c.rRNA_name=m.rRNA_name AND c.pos_id=m.pos_id
WHERE m.species='Human' AND m.rRNA_name='18S' AND c.modif_version='snoRNABase'
ORDER BY c.pos;
```

### Q3 — Reconstruct the OLD wide format (the "derived view")

This is how you reproduce the old `pos_snoRNABAse / pos_incarnato / pos_snOPY`
columns per modification — useful to keep existing views/templates working
with minimal change. Pivot `rRNA_conversion` across versions:

```sql
SELECT m.modification_id, m.modification_name, m.nucleotide, m.modification_type,
       MAX(CASE WHEN c.modif_version='snoRNABase' THEN c.pos END) AS pos_snoRNABase,
       MAX(CASE WHEN c.modif_version='incarnato'  THEN c.pos END) AS pos_incarnato,
       MAX(CASE WHEN c.modif_version='snOPY'      THEN c.pos END) AS pos_snOPY
FROM raw_data."rRNA_modifications" m
JOIN raw_data."rRNA_conversion" c
  ON c.species=m.species AND c.rRNA_name=m.rRNA_name AND c.pos_id=m.pos_id
WHERE m.species='Human' AND m.rRNA_name='18S'
GROUP BY m.modification_id, m.modification_name, m.nucleotide, m.modification_type
ORDER BY pos_snoRNABase;
```

**Recommended:** wrap this as a Django manager method / raw SQL helper so views
keep calling `Modification.objects.wide(species='Human', rRNA_name='18S')` and
get back rows that *look* like the old model. This is the cheapest path to
"generate the same data" the user asked for.

```python
# sketch: a manager method returning the pivoted wide rows
def wide(self, species, rRNA_name):
    qs = (self.get_queryset()
            .filter(species=species, rRNA_name=rRNA_name)
            .values('modification_id', 'modification_name', 'nucleotide',
                    'modification_type', 'pos_id'))
    # fetch conversion rows for these pos_ids and pivot in Python, or use raw SQL above.
```

### Q4 — Convert a position between two versions (replaces conversion_18S/28S)

Old: look up `conversion_18S` by `pos_snoRNABase` to get `pos_incarnato`/`pos_snOPY`.
New: go version A → `pos_id` → version B.

```sql
-- given (species, rRNA_name, native pos in version A), find the pos in version B
SELECT cb.pos AS pos_in_B
FROM raw_data."rRNA_conversion" ca
JOIN raw_data."rRNA_conversion" cb
  ON cb.species=ca.species AND cb.rRNA_name=ca.rRNA_name AND cb.pos_id=ca.pos_id
WHERE ca.species=? AND ca.rRNA_name=? AND ca.modif_version=?  -- version A
  AND ca.pos=?                                               -- native pos in A
  AND cb.modif_version=?;                                    -- version B
```

### Q5 — The reference / all versions of an rRNA

```python
# reference sequence for a (species, rRNA_name)
rRNA.objects.get(species='Human', rRNA_name='18S', rRNA_type='ref')

# all versions of a (species, rRNA_name)
rRNA.objects.filter(species='Human', rRNA_name='18S').values('modif_version', 'rRNA_type', 'rRNA_id')
```

### Q6 — Which rRNAs exist for a species

```python
rRNA.objects.filter(species='Human').values_list('rRNA_name', flat=True).distinct()
# -> ['5.8S', '18S', '28S']
```

---

## 4. Critical gotchas

1. **`pos_id != source position`.** Alt insertions shift it (`Um1668` → 1670).
   Always resolve coordinates via `rRNA_conversion`. Never display `pos_id` as if
   it were a sequence position without converting.
2. **5.8S (and any ref-only rRNA) has NO `rRNA_conversion` rows.** Joining
   `rRNA_modifications` → `rRNA_conversion` for 5.8S yields **nothing**. For such
   rRNAs, `pos_id` already equals the (single) native position, so use a **LEFT
   join** and fall back to `m.pos_id` when the join misses. In practice:
   ```sql
   COALESCE(c.pos, m.pos_id) AS native_pos
   ```
   (For 18S/28S, the ref row always exists, so `c.pos` is present.)
3. **`species = 'Human'`** (not `'Homo_sapiens'`), to match the snoRNA `specie`
   table. The `rRNA_alignment_*` tables use `'Homo_sapiens'` — do not assume a
   single species vocabulary; keep them separate unless you unify deliberately.
4. **`modif_version` casing**: `'snoRNABase'` (capital B), `'incarnato'`,
   `'snOPY'` (caps). The reference is `modif_version='snoRNABase'`
   (equivalently `rRNA_type='ref'`).
5. **Composite joins, no FKs.** Every join keys on the full composite:
   - modifications ↔ conversion: `(species, rRNA_name, pos_id)`
   - conversion ↔ rRNAs: `(species, rRNA_name, modif_version)` (all three —
     `modif_version` repeats across `rRNA_name`, e.g. snoRNABase exists for 5.8S,
     18S, and 28S).
6. **Per-source attestation is gone.** The old row implicitly meant "all sources
   report this mod." The new row just says "this universal position is modified."
   If you ever need "reported by source X but not Y," that information is **no
   longer in the DB** — flag to the user if a view relied on it.
7. **`base_zavolan` never existed correctly.** The old column was mislabeled
   (it held snOPY data). Don't try to resurrect a "Zavolan" base; it was snOPY
   all along.

---

## 5. What to grep for in the Django code

Search the Django repo for these tokens — every hit is likely a change site:

```
pos_snoRNABAse      # old typo'd column on rRNA_modifications (note capital BA)
pos_snoRNABase      # the correctly-spelled variants if any
pos_incarnato
pos_snOPY
conversion_18S
conversion_28S
base_zavolan        # mislabel; should map to snOPY
base_snoRNABase
base_incarnato
rRNA_modifications
rRNAs               # the model/table — gain species, rRNA_type
modif_version
```

Typical places: `models.py`, `views.py`, any `queries.sql` / raw SQL, serializers,
templates that render modification positions, and any "rRNA browser" / target
view.

---

## 6. Suggested implementation order

1. **Update models** (§2): add `species`/`rRNA_type` to the rRNA model; rewrite
   the modifications model; add the `rRNAConversion` model; delete the two
   conversion models. Make migrations (`managed=False`, so migrations just record
   the schema state — the *real* schema is owned by the pipeline).
2. **Add a `wide()` helper** (§3 Q3) and a `position_in_version()` helper (§3 Q1)
   so view code can keep asking the same questions.
3. **Rewrite views/templates** that referenced the dropped columns or the deleted
   tables, routing them through the helpers. Mind 5.8S (§4.2).
4. **Add the composite indexes** (§2) so the joins are fast at ~15k conversion rows.
5. **Validate** against the counts in §1 and the `Um1668` round-trip (§3 Q1 → 1670).

---

## 7. Cutover coordination (important)

The DB load and the Django deploy must happen **together**:

- If the **DB loads first**, the old Django queries non-existent columns/tables
  (`pos_snoRNABAse`, `conversion_18S`) → 500s on any rRNA page.
- If **Django deploys first**, the new queries hit the old schema → 500s.

Recommended sequence:
1. Finish + test the Django changes against a **copy/staging** of the DB loaded
   with the new schema.
2. In one maintenance window: re-run the pipeline load (see below) **and** deploy
   Django.

Pipeline load command (from the `snodb3` repo):
```bash
source /home/danx/.zshrc && conda activate snakemake7.24
snakemake -c1 --use-conda \
    data/psql/rRNA_conversion/data_table.tsv \
    data/psql/rRNAs/data_table.tsv \
    data/psql/rRNA_modifications/data_table.tsv data/psql/permission.tok
```
This drops + recreates those three tables and re-grants `scottweb_surfer`.

---

## 8. Post-migration sanity checks (run against the DB)

```sql
-- counts match §1
SELECT COUNT(*) FROM raw_data."rRNAs";                 -- 6
SELECT rRNA_name, COUNT(*) FROM raw_data."rRNA_modifications" GROUP BY 1;  -- 18S=86,28S=128,5.8S=4
SELECT rRNA_name, modif_version, COUNT(*) FROM raw_data."rRNA_conversion" GROUP BY 1,2;

-- 5.8S has no conversion rows
SELECT COUNT(*) FROM raw_data."rRNA_conversion" WHERE rRNA_name='5.8S';   -- 0

-- Um1668 round-trip
SELECT c.modif_version, c.pos
FROM raw_data."rRNA_modifications" m
JOIN raw_data."rRNA_conversion" c
  ON c.species=m.species AND c.rRNA_name=m.rRNA_name AND c.pos_id=m.pos_id
WHERE m.species='Human' AND m.rRNA_name='18S' AND m.modification_name='Um1668'
ORDER BY c.modif_version;
-- incarnato=1668, snOPY=1670, snoRNABase=1668
```

---

## 9. Open questions to confirm with the user before/while migrating

- Does any view rely on **per-source reporting** of a modification (Q: "which
  sources annotated this mod")? That is no longer stored (§4.6).
- Should the site expose **multiple species** in the rRNA browser now (mouse is
  coming), or stay Human-only for now? The schema is ready either way.
- Is `species='Human'` acceptable long-term, or should the codebase standardize
  on `'Homo_sapiens'` (would require migrating the snoRNA `specie` table too)?
