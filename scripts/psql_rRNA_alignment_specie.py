import os

import pandas as pd

from snakemake.shell import shell

file = snakemake.input.specie_data

psql_script = snakemake.output.specie_script
psql_data = snakemake.output.specie_psql

host_script = snakemake.params.host_script


def create_psql_script():

    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "rRNA_alignment_specie";\n\n'
    create_table = """
CREATE TABLE "rRNA_alignment_specie" (
    id bigserial PRIMARY KEY,
    unique_id CHARACTER VARYING (50) NOT NULL,
    accession CHARACTER VARYING (50),
    type CHARACTER VARYING (20) NOT NULL,
    specie CHARACTER VARYING (60) NOT NULL,
    source CHARACTER VARYING (200),
    sequence TEXT
);
    """
    import_data = """
\\copy "rRNA_alignment_specie" (unique_id, accession, type, specie, source, sequence) FROM '/sql/data_table.tsv' WITH (DELIMITER E'\\t', NULL '');
    """

    with open(psql_script, 'w') as f:
        f.write(search_path)
        f.write(del_table)
        f.write(create_table)
        f.write(import_data)


def main():

    df = pd.read_csv(file, sep='\t')

    # Rename 'id' (GenBank accession) to match SQL column name
    df = df.rename(columns={'id': 'accession'})

    df = df[[
        'unique_id', 'accession', 'type', 'specie', 'source', 'sequence',
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
