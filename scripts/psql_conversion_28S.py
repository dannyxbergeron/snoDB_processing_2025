import os

import numpy as np
import pandas as pd

from snakemake.shell import shell

file = snakemake.input.conversion_28S

psql_data = snakemake.output.conversion_28S_psql
psql_script = snakemake.output.conversion_28S_script

host_script = snakemake.params.host_script


def create_psql_script():

    # --------------- FOR target TABLE ---------------------------
    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "conversion_28S";\n\n'
    create_table = """
CREATE TABLE "conversion_28S" (
    pos_id SMALLSERIAL NOT NULL PRIMARY KEY,
    pos_snoRNABase INTEGER,
    pos_incarnato INTEGER,
    base_snoRNABase varchar(1),
    base_incarnato varchar(1)
);
    """
    import_data = """
\\copy "conversion_28S" (pos_snoRNABase, pos_incarnato, base_snoRNABase, base_incarnato) \
FROM '/sql/data_table.tsv' \
WITH (DELIMITER E'\\t', NULL '-1');\n
    """

    with open(psql_script, 'w') as f:
        f.write(search_path)
        f.write(del_table)
        f.write(create_table)
        f.write(import_data)




def main():

    df = pd.read_csv(file)
    df = df.fillna(-1)
    df['snoRNABase_pos'] = df.snoRNABase_pos.map(int)
    df['Incarnato_pos'] = df.Incarnato_pos.map(int)

    # Replace T by U everywhere
    df.snoRNABase_seq = df.snoRNABase_seq.str.replace('T', 'U')
    df.Incarnato_seq = df.Incarnato_seq.str.replace('T', 'U')

    create_psql_script()

    df.to_csv(psql_data, sep='\t', index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_data)
    shell(f'cp scripts/psql_container.sh {data_dir}')

    # Run the docker host script
    shell(f'./{host_script} {data_dir}')


if __name__ == '__main__':
    main()
