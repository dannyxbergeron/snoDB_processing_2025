import os

import numpy as np
import pandas as pd

from snakemake.shell import shell

file = snakemake.input.mapped_ENCODE

psql_script = snakemake.output.encode_eclip_script
psql_data = snakemake.output.encode_eclip_psql

host_script = snakemake.params.host_script


def create_psql_script():

    # --------------- FOR target TABLE ---------------------------
    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "encode_eclip";\n\n'
    create_table = """
CREATE TABLE "encode_eclip" (
    eclip_id SMALLSERIAL NOT NULL PRIMARY KEY,
    rbp varchar(50) NOT NULL,
    cell_line varchar(50) NOT NULL,
    pvalue DOUBLE PRECISION NOT NULL,
    score DOUBLE PRECISION NOT NULL,
    unique_id varchar(50) NOT NULL
);
    """
    import_data = """
\\copy encode_eclip (rbp, cell_line, pvalue, score, unique_id) \
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
    print(df.columns)
    df = df[[
        'RBP', 'cell_line', 'pValue', 'score', 'unique_id'
    ]]

    create_psql_script()

    df.to_csv(psql_data, sep='\t', index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_data)
    shell(f'cp scripts/psql_container.sh {data_dir}')

    # Run the docker host script
    shell(f'./{host_script} {data_dir}')


if __name__ == '__main__':
    main()
