import os
import subprocess

import numpy as np
import pandas as pd

new_data = snakemake.input.new_data
full_info = snakemake.input.final_snodb_ids

out = snakemake.output.new_final_snodb_ids


def main():

    cols = [
        'chr', 'start', 'end', 'unique_id', 'score', 'strand', 'box_type'
    ]
    new_data_df = pd.read_csv(new_data, sep="\t", names=cols)
    full_info_df = pd.read_csv(full_info, sep="\t")

    # Format the new data
    new_data_df['chr'] = 'chr' + new_data_df['chr']

    ## Merge the two dataframes
    updated_df = pd.concat([full_info_df, new_data_df], ignore_index=True)
    updated_df.drop(columns=['score'], inplace=True)

    updated_df.sort_values(by=['snoDB_id', 'chr', 'start'], inplace=True)


    current = ''
    for i in updated_df.index:
        snoDB_id = updated_df.at[i, 'snoDB_id']
        if not pd.isna(snoDB_id):
            current = snoDB_id

        else:
            num = int(current.replace('snoDB', ''))
            new_val = 'snoDB' + str(num + 1)
            current = new_val
            updated_df.at[i, 'snoDB_id'] = new_val

    print(updated_df)

    updated_df.to_csv(out, sep="\t", index=False)



if __name__ == "__main__":
    main()
