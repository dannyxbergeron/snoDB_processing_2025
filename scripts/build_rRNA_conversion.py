"""Build the generic, long-format rRNA conversion table.

Reads the wide conversion CSVs (``conversion_18S.csv`` / ``conversion_28S.csv``),
where each row aligns corresponding positions across the rRNA versions of a single
(species, rRNA_name), and produces a *long* intermediate TSV keyed by a universal
``pos_id`` alignment column.

``pos_id`` is NOT a copy of any source's numbering -- it is the column of an
alignment built from all versions:

  * the reference (``snoRNABase``) positions form the backbone, in ascending order;
  * every alt-only insertion (a position present in an alt version but absent from
    the reference) is placed immediately after the nearest preceding reference
    position it is bounded by, so it gets its own ``pos_id``;
  * ``pos_id`` is then 1..M in that biological order.

This intermediate is consumed by both ``psql_rRNA_conversion`` (load it verbatim)
and ``psql_rRNA_modifications`` (map each modification's reference position onto
``pos_id``), guaranteeing both tables share the same coordinate grid.

The conversion is *within* a (species, rRNA_name), across its versions. A
ref-only rRNA (e.g. 5.8S, which has no conversion CSV) simply contributes no rows
here -- that is correct, not missing data.
"""

import bisect
import os

import pandas as pd

SPECIES = 'Human'

# Per-rRNA version layout: (modif_version, position_column, sequence_column).
# The FIRST entry of each list is the reference (backbone) version.
VERSIONS = {
    '18S': [
        ('snoRNABase', 'snoRNABase_pos', 'snoRNABase_seq'),
        ('incarnato',  'Incarnato_pos',  'Incarnato_seq'),
        ('snOPY',      'snOPY_pos',      'snOPY_seq'),
    ],
    '28S': [
        ('snoRNABase', 'snoRNABase_pos', 'snoRNABase_seq'),
        ('incarnato',  'Incarnato_pos',  'Incarnato_seq'),
    ],
}


def _is_na(value):
    return pd.isna(value) or value == 'NA'


def compute_alignment(df, versions, species, rRNA_name):
    """Return the long-format conversion rows for one (species, rRNA_name).

    Parameters
    ----------
    df : DataFrame
        Wide conversion table (one row = one alignment column across versions).
    versions : list[tuple[str, str, str]]
        ``(modif_version, position_column, sequence_column)``; ``versions[0]``
        is the reference backbone.
    """
    ref_pos_col = versions[0][1]

    # --- 1. Parse each wide row into an alignment column. ------------------
    # Each column holds the reference position (if any) and, per version present,
    # that version's native position + base.
    raw_cols = []
    for _, row in df.iterrows():
        rp = row[ref_pos_col]
        ref_pos = None if _is_na(rp) else int(rp)
        col = {'ref_pos': ref_pos, 'versions': {}}
        for modif_version, pos_col, seq_col in versions:
            pos = row[pos_col]
            if _is_na(pos):
                continue
            base = row[seq_col]
            base = '' if _is_na(base) else str(base).replace('T', 'U')
            col['versions'][modif_version] = (int(pos), base)
        raw_cols.append(col)

    # --- 2. Per alt version, map alt_pos -> ref_pos for the mapped columns. -
    # Used to anchor alt-only insertions after their nearest preceding ref pos.
    alt_maps = {}
    for modif_version, _, _ in versions[1:]:
        mapping = {col['versions'][modif_version][0]: col['ref_pos']
                   for col in raw_cols
                   if col['ref_pos'] is not None and modif_version in col['versions']}
        alt_maps[modif_version] = (sorted(mapping.keys()), mapping)

    def lower_bound_ref(col):
        """Reference position this column should sort after."""
        if col['ref_pos'] is not None:
            return col['ref_pos']
        # alt-only insertion: anchor on the first alt version present here
        for modif_version, _, _ in versions[1:]:
            if modif_version in col['versions']:
                alt_pos = col['versions'][modif_version][0]
                positions, mapping = alt_maps[modif_version]
                i = bisect.bisect_right(positions, alt_pos) - 1
                return 0 if i < 0 else mapping[positions[i]]
        return 0

    def sort_key(idx, col):
        # Reference columns sort before alt insertions sharing the same anchor,
        # so an insertion lands strictly between ref_pos and ref_pos + 1.
        if col['ref_pos'] is not None:
            return (col['ref_pos'], 0, 0, idx)
        anchor = lower_bound_ref(col)
        alt_pos = next(iter(col['versions'].values()))[0]
        return (anchor, 1, alt_pos, idx)

    ordered = [col for _, col in sorted(enumerate(raw_cols),
                                        key=lambda t: sort_key(t[0], t[1]))]

    # --- 3. Assign pos_id (1-based) and emit one row per (column, version). -
    rows = []
    for pos_id, col in enumerate(ordered, start=1):
        for modif_version, _, _ in versions:
            if modif_version in col['versions']:
                pos, base = col['versions'][modif_version]
                rows.append((species, rRNA_name, pos_id, modif_version, pos, base))

    return pd.DataFrame(rows, columns=['species', 'rRNA_name', 'pos_id',
                                       'modif_version', 'pos', 'base'])


def main():

    sources = {
        '18S': snakemake.input.conversion_18S,
        '28S': snakemake.input.conversion_28S,
    }

    frames = []
    for rRNA_name, path in sources.items():
        df = pd.read_csv(path, na_values=['NA'])
        frames.append(compute_alignment(df, VERSIONS[rRNA_name], SPECIES, rRNA_name))

    long = pd.concat(frames, ignore_index=True)
    long = (long.sort_values(['rRNA_name', 'pos_id', 'modif_version'])
                .reset_index(drop=True))

    out_path = snakemake.output.rRNA_conversion_long
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    long.to_csv(out_path, sep='\t', index=False)

    print(f'Wrote {len(long)} conversion rows to {out_path}')
    print(long.groupby(['rRNA_name', 'modif_version']).size())


if __name__ == '__main__':
    main()
