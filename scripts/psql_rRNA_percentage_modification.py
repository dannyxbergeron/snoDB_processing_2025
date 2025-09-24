import os

import numpy as np
import pandas as pd

from snakemake.shell import shell

file = snakemake.input.all_rRNA_percentage_modification

psql_data = snakemake.output.rRNA_percentage_modification_psql
psql_script = snakemake.output.rRNA_percentage_modification_script

host_script = snakemake.params.host_script


def create_psql_script():

    # --------------- FOR target TABLE ---------------------------
    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "rRNA_percentage_modification";\n\n'
    create_table = """
CREATE TABLE "rRNA_percentage_modification" (
    modif_id SMALLSERIAL NOT NULL PRIMARY KEY,
    sample_id varchar(150) NOT NULL,
    sample varchar(150) NOT NULL,
    rRNA_name varchar(20) NOT NULL,
    pos_snoRNABAse INTEGER NOT NULL,
    percentage_modif DOUBLE PRECISION,
    modification_type varchar(30) NOT NULL
);
    """
    import_data = """
\\copy "rRNA_percentage_modification" (sample_id, sample, rRNA_name, pos_snoRNABAse, percentage_modif, modification_type) \
FROM '/sql/data_table.tsv' \
WITH (DELIMITER E'\\t', NULL '');\n
    """

    with open(psql_script, 'w') as f:
        f.write(search_path)
        f.write(del_table)
        f.write(create_table)
        f.write(import_data)




def main():

    df = pd.read_csv(file)
    df.sort_values(by=['rRNA_name', 'pos_snoRNABase'], inplace=True)

    create_psql_script()

    df.to_csv(psql_data, sep='\t', index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_data)
    shell(f'cp scripts/psql_container.sh {data_dir}')

    # Run the docker host script
    shell(f'./{host_script} {data_dir}')


if __name__ == '__main__':
    main()
