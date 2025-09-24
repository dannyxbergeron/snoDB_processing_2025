import os

import numpy as np
import pandas as pd

from snakemake.shell import shell

file = snakemake.input.rRNAs

psql_data = snakemake.output.rRNAs_psql
psql_script = snakemake.output.rRNAs_script

host_script = snakemake.params.host_script


def create_psql_script():

    # --------------- FOR target TABLE ---------------------------
    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "rRNAs";\n\n'
    create_table = """
CREATE TABLE "rRNAs" (
    rRNA_id varchar(50) NOT NULL PRIMARY KEY,
    rRNA_name varchar(50) NOT NULL,
    sequence varchar(7500) NOT NULL,
    sources varchar(200) NOT NULL,
    modif_version varchar(30) NOT NULL
);
    """
    import_data = """
\\copy "rRNAs" (rRNA_id, rRNA_name, sequence, sources, modif_version) \
FROM '/sql/data_table.tsv' \
WITH (DELIMITER E'\\t', NULL '.');\n
    """

    with open(psql_script, 'w') as f:
        f.write(search_path)
        f.write(del_table)
        f.write(create_table)
        f.write(import_data)




def main():

    df = pd.read_csv(file)

    modif_version = []
    for sources in df['references'].values:
        if 'snoRNABase' in sources:
            modif_version.append('snoRNABase')
        elif 'Incarnato' in sources:
            modif_version.append('incarnato')
        elif 'snOPY' in sources:
            modif_version.append('snOPY')
        else:
            print(sources)
            print('ERROR !!!')
            exit()

    df['modif_version'] = modif_version

    # Change DNA sequence to RNA sequence
    df.sequence = df.sequence.str.replace('T', 'U')

    create_psql_script()

    df.to_csv(psql_data, sep='\t', index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_data)
    shell(f'cp scripts/psql_container.sh {data_dir}')

    # Run the docker host script
    shell(f'./{host_script} {data_dir}')


if __name__ == '__main__':
    main()
