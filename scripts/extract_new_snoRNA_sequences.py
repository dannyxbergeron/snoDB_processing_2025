import os
import subprocess
import tempfile

import pandas as pd

new_data = snakemake.input.new_data
new_final_snodb_ids = snakemake.input.new_final_snodb_ids
existing_sequences = snakemake.input.existing_sequences
reference_genome = snakemake.input.reference_genome

output_file = snakemake.output.new_snoRNA_sequences


def main():
    # Read the new final snoDB IDs to build unique_id -> snoDB_id mapping
    final_ids_df = pd.read_csv(new_final_snodb_ids, sep="\t")

    # Get only the new snoRNAs (those with snoDB_id >= snoDB2121)
    # Extract numeric part of snoDB_id for comparison
    final_ids_df['snoDB_num'] = final_ids_df['snoDB_id'].str.replace('snoDB', '').astype(int)
    new_ids_df = final_ids_df[final_ids_df['snoDB_num'] >= 2121]

    # Build mapping: unique_id -> snoDB_id
    id_mapping = dict(zip(new_ids_df['unique_id'], new_ids_df['snoDB_id']))

    # Read BED file with new snoRNA coordinates
    cols = ['chr', 'start', 'end', 'unique_id', 'score', 'strand', 'box_type']
    bed_df = pd.read_csv(new_data, sep="\t", names=cols)

    # Keep chromosome names as-is (reference genome uses "1", "2", etc. without "chr" prefix)

    # Create temporary directory and files
    tmp_dir = "data/tmp"
    os.makedirs(tmp_dir, exist_ok=True)

    tmp_bed = os.path.join(tmp_dir, "new_snorna.bed")
    tmp_fasta = os.path.join(tmp_dir, "new_snorna.fasta")

    # Write temporary BED file
    # bedtools getfasta expects: chrom, start, end, name, score, strand
    bed_output = bed_df[['chr', 'start', 'end', 'unique_id', 'score', 'strand']]
    bed_output.to_csv(tmp_bed, sep="\t", header=False, index=False)

    # Call bedtools getfasta with -s flag for strand-aware extraction
    # -name flag preserves the name column (unique_id) in FASTA header
    cmd = [
        "bedtools", "getfasta",
        "-fi", reference_genome,
        "-bed", tmp_bed,
        "-fo", tmp_fasta,
        "-s",   # Strand-aware: reverse complement for minus strand
        "-name"  # Use name column in FASTA header
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"bedtools getfasta failed: {result.stderr}")

    # Parse FASTA output
    new_sequences = []
    with open(tmp_fasta, 'r') as f:
        current_id = None
        current_seq = []

        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_id is not None:
                    seq = ''.join(current_seq)
                    # Get snoDB_id from mapping
                    if current_id in id_mapping:
                        new_sequences.append({
                            'snoDB_id': id_mapping[current_id],
                            'unique_id': current_id,
                            'seq': seq
                        })
                # Parse new header: >unique_id::chr:start-end(strand)
                # Extract just the unique_id (part before ::)
                header = line[1:].split()[0]  # Remove '>' and get first part
                current_id = header.split('::')[0]  # Get part before '::'
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget the last sequence
        if current_id is not None:
            seq = ''.join(current_seq)
            if current_id in id_mapping:
                new_sequences.append({
                    'snoDB_id': id_mapping[current_id],
                    'unique_id': current_id,
                    'seq': seq
                })

    # Create DataFrame with new sequences
    new_seq_df = pd.DataFrame(new_sequences)

    # Read existing sequences
    existing_df = pd.read_csv(existing_sequences, sep="\t")

    # Merge old and new sequences
    combined_df = pd.concat([existing_df, new_seq_df], ignore_index=True)

    # Write output
    combined_df.to_csv(output_file, sep="\t", index=False)

    print(f"Extracted {len(new_sequences)} new sequences")
    print(f"Total sequences in output: {len(combined_df)}")

    # Clean up temporary files
    os.remove(tmp_bed)
    os.remove(tmp_fasta)


if __name__ == "__main__":
    main()
