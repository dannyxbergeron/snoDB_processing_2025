import os

import pandas as pd

from snakemake.shell import shell

file = snakemake.input.rRNA_conversion_long

psql_data = snakemake.output.rRNA_conversion_psql
psql_script = snakemake.output.rRNA_conversion_script

host_script = snakemake.params.host_script


def create_psql_script():

    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "rRNA_conversion";\n\n'
    create_table = """
CREATE TABLE "rRNA_conversion" (
    id bigserial PRIMARY KEY,
    species varchar(50) NOT NULL,
    rRNA_name varchar(50) NOT NULL,
    pos_id INTEGER NOT NULL,
    pos INTEGER,
    base varchar(1),
    modif_version varchar(30) NOT NULL
);
    """
    import_data = """
\\copy "rRNA_conversion" (species, rRNA_name, pos_id, pos, base, modif_version) \
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

    # Column order is the contract for the positional \copy below.
    df = df[['species', 'rRNA_name', 'pos_id', 'pos', 'base', 'modif_version']]
    df['pos'] = df['pos'].astype(int)

    create_psql_script()

    df.to_csv(psql_data, sep='\t', index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_data)
    shell(f'cp scripts/psql_container.sh {data_dir}')

    # Run the docker host script
    shell(f'./{host_script} {data_dir}')


if __name__ == '__main__':
    main()
