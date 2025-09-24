import os

import numpy as np
import pandas as pd

from snakemake.shell import shell

file = snakemake.input.all_rRNA_modifications

psql_data = snakemake.output.rRNA_modifications_psql
psql_script = snakemake.output.rRNA_modifications_script

host_script = snakemake.params.host_script


def create_psql_script():

    # --------------- FOR target TABLE ---------------------------
    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "rRNA_modifications";\n\n'
    create_table = """
CREATE TABLE "rRNA_modifications" (
    modification_id SMALLSERIAL NOT NULL PRIMARY KEY,
    rRNA_name varchar(20) NOT NULL,
    pos_snoRNABAse INTEGER NOT NULL,
    pos_incarnato INTEGER NOT NULL,
    pos_snOPY INTEGER NOT NULL,
    nucleotide varchar(1) NOT NULL,
    modification_type varchar(30) NOT NULL,
    modification_name varchar(30) NOT NULL
);
    """
    import_data = """
\\copy "rRNA_modifications" (rRNA_name, \
                             pos_snoRNABAse, \
                             pos_incarnato, \
                             pos_snOPY, \
                             nucleotide, \
                             modification_type, \
                             modification_name) \
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
    df.sort_values(by=['rRNA', 'pos_snoRNABase'], inplace=True)

    create_psql_script()

    df.to_csv(psql_data, sep='\t', index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_data)
    shell(f'cp scripts/psql_container.sh {data_dir}')

    # Run the docker host script
    shell(f'./{host_script} {data_dir}')


if __name__ == '__main__':
    main()
