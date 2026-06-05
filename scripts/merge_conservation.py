import pandas as pd

# bigWigAverageOverBed output with 4-column BED (no header, 6 columns):
# name, size, covered, sum, mean0, mean
new_cons = pd.read_csv(
    snakemake.input.new_conservation,
    sep='\t',
    header=None,
    names=['name', 'size', 'covered', 'sum', 'mean0', 'mean']
)

# Format to match final_ids_conservation.tsv columns:
# unique_id, span, tot_cons, mean_conservation
new_formatted = pd.DataFrame({
    'unique_id': new_cons['name'],
    'span': new_cons['size'],
    'tot_cons': new_cons['sum'].round(3),
    'mean_conservation': (new_cons['sum'] / new_cons['size']).round(3)
})

# Read existing conservation data (has header)
existing_cons = pd.read_csv(
    snakemake.input.existing_conservation,
    sep='\t'
)

# Combine existing + new
combined = pd.concat([existing_cons, new_formatted], ignore_index=True)

# Write output with header, no index
combined.to_csv(
    snakemake.output[0],
    sep='\t',
    index=False
)
