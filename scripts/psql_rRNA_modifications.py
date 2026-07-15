import os

import numpy as np
import pandas as pd

from snakemake.shell import shell

file = snakemake.input.all_rRNA_modifications
conversion_long = snakemake.input.rRNA_conversion_long

psql_data = snakemake.output.rRNA_modifications_psql
psql_script = snakemake.output.rRNA_modifications_script

host_script = snakemake.params.host_script


def create_psql_script():

    # --------------- FOR rRNA_modifications TABLE ---------------------------
    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "rRNA_modifications";\n\n'
    create_table = """
CREATE TABLE "rRNA_modifications" (
    modification_id SMALLSERIAL NOT NULL PRIMARY KEY,
    species varchar(50) NOT NULL,
    rRNA_name varchar(20) NOT NULL,
    pos_id INTEGER NOT NULL,
    nucleotide varchar(1) NOT NULL,
    modification_type varchar(30) NOT NULL,
    modification_name varchar(30) NOT NULL
);
    """
    import_data = """
\\copy "rRNA_modifications" (species, \
                             rRNA_name, \
                             pos_id, \
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
    df['species'] = 'Human'
    df = df.rename(columns={'rRNA': 'rRNA_name',
                            'modif_name': 'modification_name'})

    # Map each modification's reference (snoRNABase) position onto the universal
    # pos_id computed in build_rRNA_conversion. rRNAs without a conversion source
    # (e.g. 5.8S, ref-only) are absent from the intermediate, so their pos_id is
    # simply their own snoRNABase position (no alt insertions to shift it).
    conv = pd.read_csv(conversion_long, sep='\t')
    ref_map = (conv[conv['modif_version'] == 'snoRNABase']
               [['species', 'rRNA_name', 'pos', 'pos_id']]
               .rename(columns={'pos': 'pos_snoRNABase'}))
    df = df.merge(ref_map, on=['species', 'rRNA_name', 'pos_snoRNABase'], how='left')
    df['pos_id'] = df['pos_id'].fillna(df['pos_snoRNABase']).astype(int)

    df.sort_values(by=['species', 'rRNA_name', 'pos_id'], inplace=True)

    create_psql_script()

    # Column order must match the \copy list in create_psql_script().
    df = df[['species', 'rRNA_name', 'pos_id', 'nucleotide',
             'modification_type', 'modification_name']]
    df.to_csv(psql_data, sep='\t', index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_data)
    shell(f'cp scripts/psql_container.sh {data_dir}')

    # Run the docker host script
    shell(f'./{host_script} {data_dir}')


if __name__ == '__main__':
    main()
