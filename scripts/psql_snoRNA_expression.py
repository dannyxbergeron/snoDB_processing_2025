import os

import numpy as np
import pandas as pd

from snakemake.shell import shell

file = snakemake.input.snoRNA_mapped_matrix

psql_script = snakemake.output.snoRNA_expression_script
psql_data = snakemake.output.snoRNA_expression_psql

host_script = snakemake.params.host_script


def create_psql_script(df):
    search_path = "SET SEARCH_PATH TO raw_data;\n\n"
    del_table = 'DROP TABLE IF EXISTS "snorna_expression";\n\n'
    create_table = """
CREATE TABLE "snorna_expression" (
    id SERIAL PRIMARY KEY,
    unique_id varchar(50) NOT NULL,
    value DOUBLE PRECISION NOT NULL,
    sample_type VARCHAR(50) NOT NULL
);
    """
    import_data = """
\\copy snorna_expression FROM '/sql/data_table.tsv' WITH (DELIMITER E'\\t', NULL '');
    """

    with open(psql_script, "w") as f:
        f.write(search_path)
        f.write(del_table)
        f.write(create_table)
        f.write(import_data)


def main():
    df = pd.read_csv(file, sep="\t")
    df = df.drop(columns=["gtf", "gene_name", "snoDB_id"])

    # Rename lambowitz datasets
    df = df.rename(
        columns={
            "HumanRef_1": "UHR_1",
            "HumanRef_2": "UHR_2",
            "HumanRef_3": "UHR_3",
            "BrainLam_1": "HBR_1",
            "BrainLam_2": "HBR_2",
            "BrainLam_3": "HBR_3",
        }
    )

    df = df.melt(
        id_vars=["unique_id", "is_expressed"],
        var_name="sample_type",
        value_name="value",
    )
    df = df[["unique_id", "value", "sample_type"]]

    # Get a unique bigint as the index
    df = df.reset_index()
    print(df)

    create_psql_script(df)

    print(df)
    df.to_csv(psql_data, sep="\t", index=False, header=False)

    # Copy the docker container script in the folder
    data_dir = os.path.dirname(psql_data)
    shell(f"cp scripts/psql_container.sh {data_dir}")

    # Run the docker host script
    shell(f"./{host_script} {data_dir}")


if __name__ == "__main__":
    main()
