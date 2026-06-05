import os

import numpy as np
import pandas as pd

from snakemake.shell import shell

cd_file = snakemake.input.final_cd_box
cd_guide = snakemake.input.final_cd_guide
haca_file = snakemake.input.final_ha_box
haca_guide = snakemake.input.final_haca_guide
sca_file = snakemake.input.sca_atlas_boxes
ids_file = snakemake.input.final_snodb_ids

psql_script = snakemake.output.snoRNA_boxes_script
psql_cd_data = snakemake.output.cd_boxes_psql
psql_haca_data = snakemake.output.haca_boxes_psql
psql_sca_data = snakemake.output.sca_boxes_psql

host_script = snakemake.params.host_script


def validate(serie):
    unique_val = set(serie)
    # Make sure that there is not more than one value
    if len(unique_val) != 1:
        return ';'.join(sorted(list(unique_val)))
    return list(unique_val)[0]

def create_psql_script():

    # --------------- FOR target TABLE ---------------------------
    search_path = 'SET SEARCH_PATH TO raw_data;\n'
    del_tables = """
DROP TABLE IF EXISTS "cd_boxes_location";
DROP TABLE IF EXISTS "haca_boxes_location";
DROP TABLE IF EXISTS "sca_boxes_location";\n
"""
    create_tables = """
CREATE TABLE "cd_boxes_location" (
    unique_id varchar(50) NOT NULL PRIMARY KEY,
    c_seq varchar(10),
    c_start INTEGER NOT NULL,
    d_seq varchar(10),
    d_start INTEGER NOT NULL,
    c_prime_seq varchar(10),
    c_prime_start INTEGER NOT NULL,
    d_prime_seq varchar(10),
    d_prime_start INTEGER NOT NULL,
    guide1_seq varchar(20),
    guide1_start INTEGER NOT NULL,
    guide2_seq varchar(20),
    guide2_start INTEGER NOT NULL
);

CREATE TABLE "haca_boxes_location" (
    unique_id varchar(50) NOT NULL PRIMARY KEY,
    h_seq varchar(10),
    h_start INTEGER NOT NULL,
    aca_seq varchar(10),
    aca_start INTEGER NOT NULL,
    guide1A_seq varchar(15),
    guide1A_start INTEGER NOT NULL,
    guide1B_seq varchar(15),
    guide1B_start INTEGER NOT NULL,
    guide2A_seq varchar(15),
    guide2A_start INTEGER NOT NULL,
    guide2B_seq varchar(15),
    guide2B_start INTEGER NOT NULL
);

CREATE TABLE "sca_boxes_location" (
    unique_id varchar(50) NOT NULL PRIMARY KEY,
    box1_seq varchar(10),
    box1_start INTEGER NOT NULL,
    box2_seq varchar(10),
    box2_start INTEGER NOT NULL,
    box3_seq varchar(10),
    box3_start INTEGER NOT NULL,
    box4_seq varchar(10),
    box4_start INTEGER NOT NULL
);
    """
    import_data = """
\\copy cd_boxes_location \
FROM '/sql/cd_data_table.tsv' \
WITH (DELIMITER E'\\t', NULL '.');\n

\\copy haca_boxes_location \
FROM '/sql/haca_data_table.tsv' \
WITH (DELIMITER E'\\t', NULL '.');\n

\\copy sca_boxes_location \
FROM '/sql/sca_data_table.tsv' \
WITH (DELIMITER E'\\t', NULL '.');\n
"""

    with open(psql_script, 'w') as f:
        f.write(search_path)
        f.write(del_tables)
        f.write(create_tables)
        f.write(import_data)


def replace_ns(df_):

    df = df_.copy(deep=True)
    for i in df.index:
        for col in df.columns:
            if set(str(df.at[i, col])) == {'N',}:
                df.at[i, col] = np.nan
    return df


def add_missing_snoRNAs(box_df, ids_df, box_type_filter, id_col='gene_id'):
    """Add empty box rows for snoRNAs that don't have box annotation data."""
    matching_ids = ids_df[ids_df['box_type'].isin(box_type_filter)]
    existing_ids = set(box_df[id_col])
    missing = matching_ids[~matching_ids['unique_id'].isin(existing_ids)]

    if missing.empty:
        return box_df

    # Build empty rows with the same columns as box_df
    empty_rows = []
    for uid in missing['unique_id']:
        row = {col: np.nan for col in box_df.columns}
        row[id_col] = uid
        for col in box_df.columns:
            if 'start' in col:
                row[col] = -1
        empty_rows.append(row)

    return pd.concat([box_df, pd.DataFrame(empty_rows)], ignore_index=True)


def main():

    cd_df = pd.read_csv(cd_file, sep='\t')
    cd_guide_df = pd.read_csv(cd_guide, sep='\t')
    haca_df = pd.read_csv(haca_file, sep='\t')
    haca_guide_df = pd.read_csv(haca_guide, sep='\t')
    sca_df = pd.read_csv(sca_file, sep='\t')
    ids_df = pd.read_csv(ids_file, sep='\t')[['unique_id', 'box_type']]

    # Remove snoRNA atlas id for unprocessed scaRNA df
    sca_df['gene_id'] = sca_df.unique_id
    sca_df.drop(columns=['unique_id'], inplace=True)

    # Merge guide region with boxes for cd
    for col in cd_guide_df.columns[1:]:
        cd_df[col] = cd_df.gene_id.map(dict(zip(cd_guide_df.unique_id, cd_guide_df[col])))
        if 'start' in col:
            cd_df[col] = cd_df[col].fillna(-1).astype(int)

    # Merge guide region with boxes for haca
    for col in haca_guide_df.columns[1:]:
        haca_df[col] = haca_df.gene_id.map(dict(zip(haca_guide_df.unique_id, haca_guide_df[col])))
        if 'start' in col:
            haca_df[col] = haca_df[col].fillna(-1).astype(int)

    # Add empty box rows for snoRNAs without annotation data
    cd_df = add_missing_snoRNAs(cd_df, ids_df, ['C/D', 'SNORD-like'])
    haca_df = add_missing_snoRNAs(haca_df, ids_df, ['H/ACA', 'AluACA', 'telomerase RNA'])
    sca_df = add_missing_snoRNAs(sca_df, ids_df, ['scaRNA'])

    # Replace all Ns by np.nan
    cd_df = replace_ns(cd_df)
    haca_df = replace_ns(haca_df)
    sca_df = replace_ns(sca_df)

    create_psql_script()

    cd_df.to_csv(psql_cd_data, sep='\t', index=False, header=False)
    haca_df.to_csv(psql_haca_data, sep='\t', index=False, header=False)
    sca_df.to_csv(psql_sca_data, sep='\t', index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_cd_data)
    shell(f'cp scripts/psql_container.sh {data_dir}')

    # Run the docker host script
    shell(f'./{host_script} {data_dir}')


if __name__ == '__main__':
    main()
