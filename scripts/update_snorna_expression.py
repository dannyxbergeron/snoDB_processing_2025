import os
import subprocess

import numpy as np
import pandas as pd

new_data = snakemake.input.new_data
matrix = snakemake.input.snoRNA_mapped_matrix

id_file = snakemake.input.new_final_snodb_ids

out = snakemake.output.new_snoRNA_mapped_matrix


def main():

    new_data_df = pd.read_csv(new_data, sep="\t")
    matrix_df = pd.read_csv(matrix, sep="\t")

    id_df = pd.read_csv(id_file, sep="\t")
    snodb_dict = dict(zip(id_df['unique_id'], id_df['snoDB_id']))
    snodb_unique_id_dict = dict(zip(id_df['snoDB_id'], id_df['unique_id']))
    id_name = dict(zip(id_df['unique_id'], id_df['gene_name']))

    new_data_cond = list(new_data_df.columns[1:])
    matrix_cond = list(matrix_df.columns[3:-1])

    old_more = sorted(list(set(matrix_cond) - set(new_data_cond)))
    new_more = sorted(list(set(new_data_cond) - set(matrix_cond)))

    print('-------------')
    print('Conditions that are only in the old matrix:')
    print(old_more)
    print('Conditions that are only in the new matrix:')
    print(new_more)
    print('-------------')


    ## Format the new matrix
    new_matrix = new_data_df.copy(deep=True)

    ## Keep only the wanted rows in the new matrix
    new_matrix = new_matrix.loc[
        (new_matrix.gene_id.str.startswith('snoID'))
        | (new_matrix.gene_id.str.startswith('snoDB'))
        | (new_matrix.gene_id.str.startswith('CD_'))
        | (new_matrix.gene_id.isin(['ENSG00000270141', 'ENSG00000284078'])) # TERC and SNORD138 ?
    ]
    print('Number of rows in the new matrix after filtering: ', new_matrix.shape[0])
    snoDB = new_matrix[new_matrix['gene_id'].str.startswith('snoDB')].shape[0]
    snoID = new_matrix[new_matrix['gene_id'].str.startswith('snoID')].shape[0]
    CD = new_matrix[new_matrix['gene_id'].str.startswith('CD_')].shape[0]
    missing = new_matrix[new_matrix['gene_id'].isin(['ENSG00000270141', 'ENSG00000284078'])].shape[0]
    print('Number of snoDB rows: ', snoDB)
    print('Number of snoID rows: ', snoID)
    print('Number of CD_ rows: ', CD)
    print('Number of TEC and SNORD138: ', missing)
    print('-------------------')

    conds = list(new_matrix.columns[1:])

    # Get the value of is_expressed
    new_matrix['is_expressed'] = (new_matrix.iloc[:, 1:] >= 1).any(axis=1)

    ## Change the value of the cell where the gene_id == 'ENSG00000270141' to 'ens_ENSG00000270141'
    new_matrix.loc[new_matrix['gene_id'] == 'ENSG00000270141', 'gene_id'] = 'ens_ENSG00000270141'
    new_matrix.loc[new_matrix['gene_id'] == 'ENSG00000284078', 'gene_id'] = 'ens_ENSG00000284078'

    new_matrix['snoDB_id'] = np.where(new_matrix.gene_id.str.startswith('snoDB'), new_matrix['gene_id'], new_matrix['gene_id'].map(snodb_dict))
    new_matrix['gene_id'] = np.where(new_matrix.gene_id.str.startswith('snoDB'), new_matrix['gene_id'].map(snodb_unique_id_dict), new_matrix['gene_id'])

    new_matrix['gene_name'] = new_matrix['gene_id'].map(id_name)
    new_matrix.rename(columns={'gene_id': 'unique_id'}, inplace=True)
    new_matrix['gtf'] = new_matrix['snoDB_id']


    new_cols = ['snoDB_id', 'unique_id', 'gtf', 'gene_name'] + conds + ['is_expressed']
    new_matrix = new_matrix[new_cols]

    new_matrix.sort_values(by='snoDB_id', inplace=True)

    print(new_matrix)

    new_matrix.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main()
