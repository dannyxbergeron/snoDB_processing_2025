
import os

import numpy as np
import pandas as pd

from snakemake.shell import shell

file = snakemake.input.final_snodb_ids

psql_script = snakemake.output.specie_script
psql_data = snakemake.output.specie_psql

host_script = snakemake.params.host_script

def create_psql_script():

    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "specie";\n\n'
    create_table = """
CREATE TABLE "specie" (
    unique_id varchar(50) NOT NULL PRIMARY KEY,
    specie varchar(50)
);
    """
    import_data = """
\\copy specie FROM '/sql/data_table.tsv' WITH (DELIMITER E'\\t', NULL '');
    """


    with open(psql_script, 'w') as f:
        f.write(search_path)
        f.write(del_table)
        f.write(create_table)
        f.write(import_data)



def main():

    df = pd.read_csv(file, sep='\t')


    df = df[[
        'unique_id',
    ]]
    df['specie'] = 'Human'

    create_psql_script()

    df.to_csv(psql_data, sep='\t', index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_data)
    shell(f'cp scripts/psql_container.sh {data_dir}')

    # Run the docker host script
    shell(f'./{host_script} {data_dir}')


if __name__ == '__main__':
    main()
