import os

import numpy as np
import pandas as pd

from snakemake.shell import shell

file = snakemake.input.final_snodb_ids
host_file = snakemake.input.final_sno_host
snoRNABase_file = snakemake.input.snoRNABase_hg38

psql_script = snakemake.output.ext_ids_psql_script
psql_data = snakemake.output.external_ids_psql

host_script = snakemake.params.host_script


def get_gtf_id(df):
    gtf_ids = []
    for snoDB_id, ensembl_id, unique_id in df[['snoDB_id', 'ensembl', 'unique_id']].values:
        if not pd.isnull(ensembl_id):
            gtf_ids.append(ensembl_id.split(';')[0])
        else:
            gtf_ids.append(snoDB_id)
    return gtf_ids


def get_snoRNABase_name(snoRNABase_df, df):

    snoRNAbase_dict = dict(zip(snoRNABase_df.gene_id, snoRNABase_df.gene_name))
    snoRNAbase_name = []
    for snornabase_ids in df['snoRNABase'].values:
        if pd.isnull(snornabase_ids):
            snoRNAbase_name.append(np.nan)
        else:
            tmp = []
            for snornabase_id in snornabase_ids.split(';'):
                tmp.append(snoRNAbase_dict[snornabase_id])
            snoRNAbase_name.append(';'.join(tmp))
    return snoRNAbase_name

def get_rnacentral_simplified(df):

    return [
        ';'.join([field.split('_')[0] for field in x.split(';')])
        if not pd.isnull(x)
        else x
        for x in df.rnacentral.values
    ]



def create_psql_script():

    search_path = 'SET SEARCH_PATH TO raw_data;\n\n'
    del_table = 'DROP TABLE IF EXISTS "external_ids";\n\n'
    create_table = """
CREATE TABLE "external_ids" (
    unique_id varchar(50) NOT NULL PRIMARY KEY,
    snoDB_id varchar(50) UNIQUE NOT NULL,
    Ensembl varchar(50) UNIQUE ,
    RefSeq varchar(50),
    HGNC varchar(50),
    NCBI varchar(50),
    snoRNABase varchar(50) UNIQUE,
    snoRNA_Atlas varchar(50),
    snoPY varchar(50),
    RNA_Central varchar(200),
    RFAM varchar(50),
    snoDB_legacy varchar(50) UNIQUE,
    host_gene_id varchar(50),
    snoRNABase_name varchar(50),
    snopy_orthologs_name varchar(50),
    all_ids varchar(500)
);
    """
    import_data = """
\\copy external_ids FROM '/sql/data_table.tsv' WITH (DELIMITER E'\\t', NULL '');
    """


    with open(psql_script, 'w') as f:
        f.write(search_path)
        f.write(del_table)
        f.write(create_table)
        f.write(import_data)


def create_external_ids(df):

    external_ids = []
    for row in df[df.columns[2:-2]].values:
        tmp = []
        for field in row:
            if not pd.isnull(field):
                tmp.append(field)
        if tmp == []:
            external_ids.append(np.nan)
        else:
            external_ids.append(';'.join(tmp))

    return external_ids



def main():

    initial_df = pd.read_csv(file, sep='\t', dtype={'ncbi': str})
    df = initial_df.copy(deep=True)
    df['gtf'] = get_gtf_id(df)
    host_df = pd.read_csv(host_file)
    df['host_gene_id'] = df.gtf.map(dict(zip(host_df.sno_id, host_df.gene_id)))

    df = df[[
        'unique_id', 'snoDB_id', 'ensembl', 'refseq', 'hgnc',
        'ncbi', 'snoRNABase', 'snoRNA_atlas', 'snopy',
        'rnacentral', 'rfam', 'snoDB', 'host_gene_id'
    ]]

    # CHANGE RNA_central to remove the coordinates
    df['rnacentral'] = get_rnacentral_simplified(df)

    # Get snoRNABAse name (use to link with snoRNABase)
    snoRNABase_df = pd.read_csv(snoRNABase_file, sep='\t')
    df['snoRNABase_name'] = get_snoRNABase_name(snoRNABase_df, df)

    # Re-introduce in the right order snopy_orthologs_name
    df['snopy_orthologs_name'] = initial_df['snopy_orthologs_name']

    df['all_external_ids'] = create_external_ids(df)

    print(df)
    create_psql_script()

    df.to_csv(psql_data, sep='\t', index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_data)
    shell(f'cp scripts/psql_container.sh {data_dir}')

    # Run the docker host script
    shell(f'./{host_script} {data_dir}')


if __name__ == '__main__':
    main()
