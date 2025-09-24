import os

import numpy as np
import pandas as pd

from snakemake.shell import shell

file = snakemake.input.all_rRNA_sample_percentage_modification

psql_data = snakemake.output.rRNA_sample_percentage_modification_psql
psql_script = snakemake.output.rRNA_sample_percentage_modification_script

host_script = snakemake.params.host_script


def create_psql_script():

    # --------------- FOR target TABLE ---------------------------
    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "rRNA_sample_percentage_modification";\n\n'
    create_table = """
CREATE TABLE "rRNA_sample_percentage_modification" (
    sample_id varchar(150) NOT NULL PRIMARY KEY,
    sample varchar(150) NOT NULL,
    sample_name varchar(50) NOT NULL,
    sample_rep varchar(15) NOT NULL,
    sample_type varchar(30) NOT NULL,
    modification_investigated varchar(30) NOT NULL,
    source varchar(100) NOT NULL
);
    """
    import_data = """
\\copy "rRNA_sample_percentage_modification" (sample_id, sample, sample_name, sample_rep, sample_type, modification_investigated, source) \
FROM '/sql/data_table.tsv' \
WITH (DELIMITER E'\\t', NULL '');\n
    """

    with open(psql_script, 'w') as f:
        f.write(search_path)
        f.write(del_table)
        f.write(create_table)
        f.write(import_data)



def main():

    df = pd.read_csv(file, dtype={'sample_rep': str})
    df['sorting_order'] = df.sample_rep.map(float)
    df.sort_values(by=['sample_id', 'sorting_order'], inplace=True)
    df.drop(columns=['sorting_order'], inplace=True)

    create_psql_script()

    df.to_csv(psql_data, sep='\t', index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_data)
    shell(f'cp scripts/psql_container.sh {data_dir}')

    # Run the docker host script
    shell(f'./{host_script} {data_dir}')


if __name__ == '__main__':
    main()
