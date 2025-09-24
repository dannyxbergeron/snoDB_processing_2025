import os

import numpy as np
import pandas as pd

from snakemake.shell import shell

file = snakemake.input.final_snodb_ids
cons_file = snakemake.input.final_ids_conservation
seq_file = snakemake.input.snoRNA_sequences

psql_script = snakemake.output.basic_features_script
psql_data = snakemake.output.basic_features_psql

host_script = snakemake.params.host_script

def create_psql_script():

    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "basic_features";\n\n'
    create_table = """
CREATE TABLE "basic_features" (
    unique_id varchar(50) NOT NULL PRIMARY KEY,
    gene_name varchar(50),
    synonyms varchar(100),
    box_type varchar(15) NOT NULL,
    conservation_snoRNA_atlas varchar(50),
    conservation_phastcons DOUBLE PRECISION NOT NULL,
    length INTEGER NOT NULL,
    sequence varchar(1000)
);
    """
    import_data = """
\\copy basic_features FROM '/sql/data_table.tsv' WITH (DELIMITER E'\\t', NULL '');
    """


    with open(psql_script, 'w') as f:
        f.write(search_path)
        f.write(del_table)
        f.write(create_table)
        f.write(import_data)



def main():

    df = pd.read_csv(file, sep='\t')
    cons_df = pd.read_csv(cons_file, sep='\t',)
    seq_df = pd.read_csv(seq_file, sep='\t',)

    print(seq_df.columns)
    length_list = df.end - df.start + 1
    df = df[[
        'unique_id', 'gene_name', 'synonyms', 'box_type',
        'conservation',
    ]]
    df['phastcons'] = df.unique_id.map(dict(zip(cons_df.unique_id, cons_df.mean_conservation)))
    df['phastcons'] = df['phastcons'].fillna(-1)
    df['length'] = length_list
    df['seq'] = df.unique_id.map(dict(zip(seq_df.unique_id, seq_df.seq)))

    create_psql_script()

    df.to_csv(psql_data, sep='\t', index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_data)
    shell(f'cp scripts/psql_container.sh {data_dir}')

    # Run the docker host script
    shell(f'./{host_script} {data_dir}')


if __name__ == '__main__':
    main()
