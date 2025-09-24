import os

import numpy as np
import pandas as pd

from snakemake.shell import shell

file = snakemake.input.final_sno_host_extended

psql_script = snakemake.output.host_features_script
psql_data = snakemake.output.host_features_psql

host_script = snakemake.params.host_script

def create_psql_script():

    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "host_features";\n\n'
    create_table = """
CREATE TABLE "host_features" (
    gene_id varchar(50) NOT NULL PRIMARY KEY,
    chr varchar(50) NOT NULL,
    start INTEGER NOT NULL,
    "end" INTEGER NOT NULL,
    strand varchar(1) NOT NULL,
    gene_name varchar(50),
    synonyms varchar(250),
    biotype varchar(50),
    function varchar(200)
);
    """
    import_data = """
\\copy host_features FROM '/sql/data_table.tsv' WITH (DELIMITER E'\\t', NULL '');
    """


    with open(psql_script, 'w') as f:
        f.write(search_path)
        f.write(del_table)
        f.write(create_table)
        f.write(import_data)



def main():

    df = pd.read_csv(file)
    print(df.columns)
    df = df[[
        'gene_id', 'seqname', 'host_start', 'host_end', 'strand',
        'gene_name', 'synonyms', 'gene_biotype', 'function'
    ]]
    df = df.drop_duplicates()

    create_psql_script()

    df.to_csv(psql_data, sep='\t', index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_data)
    shell(f'cp scripts/psql_container.sh {data_dir}')

    # Run the docker host script
    shell(f'./{host_script} {data_dir}')


if __name__ == '__main__':
    main()
