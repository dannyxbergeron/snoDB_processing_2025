import os

import numpy as np
import pandas as pd

from snakemake.shell import shell

file = snakemake.input.targets_psql

psql_script = snakemake.output.target_grouped_by_script
psql_data = snakemake.output.target_grouped_by_psql

host_script = snakemake.params.host_script


def agg_biotypes(serie):
    res = []
    unique_val = set(serie)
    if 'rRNA' in unique_val:
        res.append('rRNA')
        unique_val.remove('rRNA')
    if 'snRNA' in unique_val:
        res.append('snRNA')
        unique_val.remove('snRNA')
    if len(unique_val) > 0:
        res.append('Others')
    return ';'.join(res)

def concat_targets(serie):
    return ';'.join(sorted(list(serie)))


def create_psql_script():

    # --------------- FOR target TABLE ---------------------------
    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "targets_grouped_by";\n\n'
    create_table = """
CREATE TABLE "targets_grouped_by" (
    unique_id varchar(50) NOT NULL PRIMARY KEY,
    target_count INTEGER,
    target_biotypes varchar(200),
    rRNA_targets varchar(100),
    snRNA_targets varchar(100),
    lncRNA_targets varchar(100),
    protein_coding_targets varchar(250),
    snoRNA_targets varchar(250),
    miRNA_targets varchar(100),
    tRNA_targets varchar(100),
    ncRNA_targets varchar(100),
    pseudogene_targets varchar(100),
    other_targets varchar(100),
    all_targets varchar(500)
);
    """
    import_data = """
\\copy targets_grouped_by (unique_id, target_count, target_biotypes, rRNA_targets, snRNA_targets, lncRNA_targets, protein_coding_targets, snoRNA_targets, miRNA_targets, tRNA_targets, ncRNA_targets, pseudogene_targets, other_targets, all_targets) \
FROM '/sql/data_table.tsv' \
WITH (DELIMITER E'\\t', NULL '.');\n
    """

    with open(psql_script, 'w') as f:
        f.write(search_path)
        f.write(del_table)
        f.write(create_table)
        f.write(import_data)




def main():

    header = ['target', 'target_biotype', 'target_gene_id', 'sources', 'RISE_name', 'unique_id']
    df = pd.read_csv(file, sep='\t', names=header)
    df = df[[
        'unique_id', 'target', 'target_biotype'
    ]]
    target_biotypes = [
        'rRNA',
        'snRNA',
        'lncRNA',
        'protein_coding',
        'snoRNA',
        'miRNA',
        'tRNA',
        'ncRNA',
        'pseudogene',
        'others',
    ]

    df_gb = df[['unique_id', 'target', 'target_biotype']].groupby('unique_id').agg({
        'target': 'count',
        'target_biotype': agg_biotypes
    }).reset_index()

    for biotype in target_biotypes:
        tmp = df.copy(deep=True)
        tmp = tmp.loc[tmp.target_biotype == biotype]
        tmp_gb = tmp[['unique_id', 'target']].groupby('unique_id').agg(concat_targets).reset_index()
        tmp_dict = dict(zip(tmp_gb.unique_id, tmp_gb.target))
        df_gb[biotype] = df_gb.unique_id.map(tmp_dict)

    # Deal with rRNA having 18S and 18S-something
    for subunit in ['18S', '28S']:
        df_gb['rRNA'] = [
            x.replace(f'{subunit};', '') if not pd.isnull(x) and f'{subunit}-' in x
            else x
            for x in df_gb['rRNA'].values
        ]

    # Create a all_targets columns for the search in snoDB
    all_targets_gb = df[['unique_id', 'target']].groupby('unique_id').agg(lambda x: ';'.join(x)).reset_index()
    all_targets_dict = dict(zip(all_targets_gb.unique_id, all_targets_gb.target))
    df_gb['all_targets']  = df_gb['unique_id'].map(all_targets_dict)

    # fill na with . for sql query
    df_gb.fillna('.', inplace=True)

    # create_psql_script()
    create_psql_script()

    print(df_gb)
    print(df_gb.columns)
    df_gb.to_csv(psql_data, sep='\t', index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_data)
    shell(f'cp scripts/psql_container.sh {data_dir}')

    # Run the docker host script
    shell(f'./{host_script} {data_dir}')


if __name__ == '__main__':
    main()
