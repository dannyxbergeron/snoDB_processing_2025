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
rule compute_new_conservation:
    input:
        bigwig = "data/phastcons/hg38.phastCons100way.bw",
        bed = "data/new_information/snoRNAs/new_snorna_abt_efc.bed",
    output:
        join(config['path']['conservation'],
             config['conservation']['new_raw_conservation']),
    conda:
        "../envs/phastcons.yaml"
    shell:
        """awk '{{print "chr"$1"\\t"$2"\\t"$3"\\t"$4}}' {input.bed} > {output}.tmp.bed && """
        "bigWigAverageOverBed {input.bigwig} {output}.tmp.bed {output} && """
        "rm {output}.tmp.bed"""

rule merge_conservation:
    input:
        new_conservation = join(config['path']['conservation'],
                                config['conservation']['new_raw_conservation']),
        existing_conservation = join(config['path']['conservation'],
                                     config['conservation']['final_ids_conservation']),
    output:
        join(config['path']['conservation'],
             config['conservation']['new_ids_conservation']),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_conservation.py"

##############################################
## Basic features, external ids, and genomic location
## are handled directly by the updated psql rules
## (psql_basic_features, psql_external_ids, psql_genomic_location)
## which now use the snodb3 data files directly.
##############################################

##############################################
## Getting the boxes for all the new snoRNA
## TODO: Add box annotations for new snoRNAs.
## Currently, new snoRNAs are inserted with empty
## box values in psql_snoRNA_boxes. Update when
## box annotation data becomes available.
##############################################

##############################################
## Creating the new images for the new snoRNAs
## TODO: Generate images for new snoRNAs
##############################################
