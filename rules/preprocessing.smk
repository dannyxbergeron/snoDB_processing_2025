##############################################
## Adding the new snoRNAs
##############################################
rule add_alphonse_snoRNA:
    input:
        new_data = "data/new_information/snoRNAs/new_snorna_abt_efc.bed",
        final_snodb_ids = join(config['path']['processed'],
                               config['processed']['final_snodb_ids']),
    output:
        new_final_snodb_ids = join(config['path']['processed_snodb3'],
                                config['processed']['final_snodb_ids']),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/add_alphonse_snoRNA.py"


##############################################
## Updating the snoRNA expression
##############################################
rule update_snorna_expression:
    input:
        new_data = "data/new_information/expression/snorna_abundance_tpm.tsv",
        snoRNA_mapped_matrix = join(config['path']['expression'],
                                    config['expression']['snoRNA_mapped_matrix']),
        new_final_snodb_ids = join(config['path']['processed_snodb3'],
                                config['processed']['final_snodb_ids']),
    output:
        new_snoRNA_mapped_matrix = join(config['path']['expression'],
                                        config['expression']['snoRNA_mapped_matrix_snodb3']),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/update_snorna_expression.py"


##############################################
## Extracting the new snoRNA sequences
##############################################
rule extract_new_snoRNA_sequences:
    input:
        new_data = "data/new_information/snoRNAs/new_snorna_abt_efc.bed",
        new_final_snodb_ids = join(config['path']['processed_snodb3'],
                                config['processed']['final_snodb_ids']),
        existing_sequences = join(config['path']['fasta'],
                                  config['snoRNA_sequences']),
        reference_genome = "data/references/Homo_sapiens.GRCh38.dna.primary_assembly_V101.fa",
    output:
        new_snoRNA_sequences = join(config['path']['fasta'],
                                    config['snoRNA_sequences_snodb3']),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_new_snoRNA_sequences.py"


##############################################
## Getting the snoRNA conservation
##############################################

##############################################
## Generating the final basic features file
##############################################

##############################################
## Creating final external ids
##############################################

##############################################
## Creating the final genomic location file
##############################################

##############################################
## Creating the final genomic location file
##############################################

##############################################
## Getting the boxes for all the new snoRNA
##############################################

##############################################
## Creating the new images for the new snoRNAs
##############################################
