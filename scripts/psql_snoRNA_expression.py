import os

import numpy as np
import pandas as pd

from snakemake.shell import shell

file = snakemake.input.snoRNA_mapped_matrix

psql_script = snakemake.output.snoRNA_expression_script
psql_data = snakemake.output.snoRNA_expression_psql

host_script = snakemake.params.host_script

def get_cell_and_tissues(df):
    return '\n'.join([
        f'    {col} DOUBLE PRECISION NOT NULL,'
        for col in df.columns[1:-1]
    ])

def create_psql_script(df):

    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "snorna_expression";\n\n'
    cell_and_tissues = get_cell_and_tissues(df)
    create_table = f"""
CREATE TABLE "snorna_expression" (
    unique_id varchar(50) NOT NULL PRIMARY KEY,
{cell_and_tissues}
    is_expressed BOOLEAN NOT NULL
);
    """
    import_data = """
\\copy snorna_expression FROM '/sql/data_table.tsv' WITH (DELIMITER E'\\t', NULL '');
    """

    with open(psql_script, 'w') as f:
        f.write(search_path)
        f.write(del_table)
        f.write(create_table)
        f.write(import_data)



def main():

    df = pd.read_csv(file, sep="\t")
    print(df.columns)
    df = df[[
        'unique_id',
        'Breast_1', 'Breast_2', 'Breast_3',
        'Ovary_1', 'Ovary_2', 'Ovary_3',
        'Prostate_2', 'Prostate_1', 'Prostate_3',
        'Brain_1', 'Brain_2', 'Brain_3',
        'Liver_1', 'Liver_2', 'Liver_3',
        'Testis_1', 'Testis_2', 'Testis_3',
        'SkeletalMuscle_1', 'SkeletalMuscle_2', 'SkeletalMuscle_3',
        'HumanRef_1', 'HumanRef_2', 'HumanRef_3',
        'BrainLam_1', 'BrainLam_2', 'BrainLam_3',
        'HCT116_1', 'HCT116_2',
        'MCF7_1', 'MCF7_2',
        'PC3_1', 'PC3_2',
        'TOV112D_1', 'TOV112D_2',
        'SKOV_frg_1', 'SKOV_frg_2',
        'is_expressed'
    ]]
    # Rename lambowitz datasets
    df = df.rename(columns={
        'HumanRef_1': 'UHR_1',
        'HumanRef_2': 'UHR_2',
        'HumanRef_3': 'UHR_3',
        'BrainLam_1': 'HBR_1',
        'BrainLam_2': 'HBR_2',
        'BrainLam_3': 'HBR_3'
    })


    create_psql_script(df)

    print(df)
    df.to_csv(psql_data, sep='\t', index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_data)
    shell(f'cp scripts/psql_container.sh {data_dir}')

    # Run the docker host script
    shell(f'./{host_script} {data_dir}')


if __name__ == '__main__':
    main()
