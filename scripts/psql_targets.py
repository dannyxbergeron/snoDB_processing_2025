import os

import numpy as np
import pandas as pd

from snakemake.shell import shell

file = snakemake.input.all_targets
ext_file = snakemake.input.cd_extended

psql_script = snakemake.output.targets_script
psql_data = snakemake.output.targets_psql

host_script = snakemake.params.host_script


def validate(serie):
    unique_val = set(serie)
    # Make sure that there is not more than one value
    if len(unique_val) != 1:
        return ';'.join(sorted(list(unique_val)))
    return list(unique_val)[0]

def create_psql_script():

    # --------------- FOR target TABLE ---------------------------
    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "targets";\n\n'
    create_table = """
CREATE TABLE "targets" (
    target_id SMALLSERIAL NOT NULL PRIMARY KEY,
    target varchar(50) NOT NULL,
    target_biotype varchar(50) NOT NULL,
    target_gene_id varchar(50),
    sources varchar(100),
    RISE_name varchar(50),
    unique_id varchar(50) NOT NULL
);
    """
    import_data = """
\\copy targets (target, target_biotype, target_gene_id, sources, RISE_name, unique_id) \
FROM '/sql/data_table.tsv' \
WITH (DELIMITER E'\\t', NULL '.');\n
    """

    with open(psql_script, 'w') as f:
        f.write(search_path)
        f.write(del_table)
        f.write(create_table)
        f.write(import_data)




def main():

    df = pd.read_csv(file, sep='\t')
    ext_df = pd.read_csv(ext_file, sep='\t')
    ext_df.drop(columns=['gene_name'], inplace=True)
    df = pd.concat([df, ext_df])

    df = df[[
        'unique_id', 'source', 'target', 'target_biotype', 'target_gene_id',
        'RISE_name'
    ]]

    '''
        To deal with the RISE data without id and without coordinates.
        there is only one problematic interaction for target_gene_id:
        ens_ENSG00000262074_U3. For RISE_name, there was 5 problematic
        entries, I kept the more general name, ex: SNORD31 instead of SNORD31B,
        But it was only for target: 18S, 28S or U6...
        Not the best solution, but there was no best solution...
    '''
    df.sort_values(by=['target_gene_id', 'RISE_name'], inplace=True)
    df.drop_duplicates(subset=['unique_id', 'target', 'source'], inplace=True)

    df['u_target'] = df.unique_id + '_' + df.target

    df_gb = df.groupby('u_target').agg({
        'unique_id': validate,
        'source': list,
        'target': validate,
        'target_biotype': validate,
        'target_gene_id': validate,
        'RISE_name': validate,
    }).reset_index()

    # Validation that there is no duplicate valies
    for col in df_gb.columns:
        if col != 'source':
            df_gb[col] = df_gb[col].fillna('.')
            # print(df_gb.loc[df_gb[col].str.contains(';')])

    df_gb['source'] = [
        ';'.join(sources)
        for sources in df_gb.source.values
    ]

    df_gb = df_gb[[
        'target', 'target_biotype',
        'target_gene_id', 'source', 'RISE_name', 'unique_id'
    ]]
    # remove nonsens RRNA from RISE
    df_gb['target'] = df_gb['target'].str.replace('_RRNA', '')

    create_psql_script()

    df_gb.to_csv(psql_data, sep='\t', index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_data)
    shell(f'cp scripts/psql_container.sh {data_dir}')

    # Run the docker host script
    shell(f'./{host_script} {data_dir}')


if __name__ == '__main__':
    main()
