import os

import pandas as pd

from snakemake.shell import shell

file = snakemake.input.lookup_info

psql_script = snakemake.output.lookup_script
psql_data = snakemake.output.lookup_psql

host_script = snakemake.params.host_script


def create_psql_script():

    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "lookup";\n\n'
    create_table = """
CREATE TABLE "lookup" (
    id bigserial PRIMARY KEY,
    unique_id CHARACTER VARYING (50),
    snoDB_id CHARACTER VARYING (50),
    gene_name CHARACTER VARYING (50),
    source_id CHARACTER VARYING (50),
    target_id CHARACTER VARYING (50),
    specie CHARACTER VARYING (50),
    perc_id DOUBLE PRECISION,
    sequence CHARACTER VARYING (1000)
);
    """
    import_data = """
\\copy lookup(unique_id, snoDB_id, gene_name, source_id, target_id, specie, perc_id, sequence) FROM '/sql/data_table.tsv' WITH (DELIMITER E'\\t', NULL '');
    """

    with open(psql_script, 'w') as f:
        f.write(search_path)
        f.write(del_table)
        f.write(create_table)
        f.write(import_data)


def main():

    df = pd.read_csv(file, sep='\t')
    df = df[[
        'unique_id', 'snoDB_id', 'gene_name', 'source_id',
        'target_id', 'specie', 'perc_id', 'sequence',
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
