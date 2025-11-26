rule psql_external_ids:
    input:
        final_snodb_ids = join(config['path']['processed'],
                               config['processed']['final_snodb_ids']),
        final_sno_host = join(config['path']['sno_host'],
                              config['sno_host']['final_sno_host']),
        snoRNABase_hg38 = join(config['path']['processed'],
                               config['processed']['snoRNABase']),
    output:
        external_ids_psql = join(config['path']['psql'],
                                 config['psql']['external_ids'],
                                 'data_table.tsv'),
        ext_ids_psql_script = join(config['path']['psql'],
                                 config['psql']['external_ids'],
                                 'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_external_ids.py"

# --------------------------------------------------------------------------

rule psql_genomic_location:
    input:
        final_snodb_ids = join(config['path']['processed'],
                               config['processed']['final_snodb_ids'])
    output:
        genomic_location_psql = join(config['path']['psql'],
                                     config['psql']['genomic_location'],
                                     'data_table.tsv'),
        genomic_location_script = join(config['path']['psql'],
                                       config['psql']['genomic_location'],
                                       'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_genomic_location.py"

# --------------------------------------------------------------------------

rule psql_basic_features:
    input:
        final_snodb_ids = join(config['path']['processed'],
                               config['processed']['final_snodb_ids']),
        final_ids_conservation = join(config['path']['conservation'],
                                      config['conservation']['final_ids_conservation']),
        snoRNA_sequences = join(config['path']['fasta'],
                                config['snoRNA_sequences']),
        snoRNA_mapped_matrix = join(config['path']['expression'],
                                    config['expression']['snoRNA_mapped_matrix']),
    output:
        basic_features_psql = join(config['path']['psql'],
                                     config['psql']['basic_features'],
                                     'data_table.tsv'),
        basic_features_script = join(config['path']['psql'],
                                       config['psql']['basic_features'],
                                       'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_basic_features.py"

# --------------------------------------------------------------------------

rule psql_specie:
    input:
        final_snodb_ids = join(config['path']['processed'],
                               config['processed']['final_snodb_ids']),
    output:
        specie_psql = join(config['path']['psql'],
                                     config['psql']['specie'],
                                     'data_table.tsv'),
        specie_script = join(config['path']['psql'],
                                       config['psql']['specie'],
                                       'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_specie.py"


# --------------------------------------------------------------------------

rule psql_host_features:
    input:
        final_sno_host_extended = join(config['path']['sno_host'],
                                       config['sno_host']['final_sno_host_extended']),
    output:
        host_features_psql = join(config['path']['psql'],
                                  config['psql']['host_features'],
                                  'data_table.tsv'),
        host_features_script = join(config['path']['psql'],
                                    config['psql']['host_features'],
                                    'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_host_features.py"


# --------------------------------------------------------------------------

rule psql_snoRNA_expression:
    input:
        snoRNA_mapped_matrix = join(config['path']['expression'],
                                    config['expression']['snoRNA_mapped_matrix']),
    output:
        snoRNA_expression_psql = join(config['path']['psql'],
                                             config['psql']['snoRNA_expression'],
                                             'data_table.tsv'),
        snoRNA_expression_script = join(config['path']['psql'],
                                        config['psql']['snoRNA_expression'],
                                        'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_snoRNA_expression.py"


# --------------------------------------------------------------------------

rule psql_host_expression:
    input:
        host_mapped_matrix = join(config['path']['expression'],
                                    config['expression']['host_mapped_matrix']),
    output:
        host_expression_psql = join(config['path']['psql'],
                                    config['psql']['host_expression'],
                                    'data_table.tsv'),
        host_expression_script = join(config['path']['psql'],
                                      config['psql']['host_expression'],
                                      'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_host_expression.py"


# --------------------------------------------------------------------------

rule psql_targets:
    input:
        all_targets = join(config['path']['all_interactions'],
                           config['all_interactions']),
        cd_extended = join(config['path']['all_interactions_extended'],
                           config['all_interactions_extended']['cd_extended']),
    output:
        targets_psql = join(config['path']['psql'],
                                    config['psql']['targets'],
                                    'data_table.tsv'),
        targets_script = join(config['path']['psql'],
                                      config['psql']['targets'],
                                      'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_targets.py"


# --------------------------------------------------------------------------

rule psql_targets_grouped_by:
    input:
        targets_psql = join(config['path']['psql'],
                                    config['psql']['targets'],
                                    'data_table.tsv'),
    output:
        target_grouped_by_psql = join(config['path']['psql'],
                                    config['psql']['target_grouped_by'],
                                    'data_table.tsv'),
        target_grouped_by_script = join(config['path']['psql'],
                                      config['psql']['target_grouped_by'],
                                      'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_target_grouped_by.py"



# --------------------------------------------------------------------------

rule psql_encode_eclip:
    input:
        mapped_ENCODE = join(config['path']['ENCODE'],
                             config['ENCODE']['mapped'])
    output:
        encode_eclip_psql = join(config['path']['psql'],
                                    config['psql']['encode_eclip'],
                                    'data_table.tsv'),
        encode_eclip_script = join(config['path']['psql'],
                                      config['psql']['encode_eclip'],
                                      'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_encode_eclip.py"


# --------------------------------------------------------------------------

rule psql_snoRNA_boxes:
    input:
        final_cd_box = join(config['path']['boxes'],
                            config['boxes']['final_cd_boxes']),
        final_ha_box = join(config['path']['boxes'],
                            config['boxes']['final_ha_boxes']),
        sca_atlas_boxes = join(config['path']['boxes'],
                               config['boxes']['sca_atlas_boxes']),
        final_cd_guide = join(config['path']['all_interactions_extended'],
                               config['all_interactions_extended']['final_cd_guide_regions']),
        final_haca_guide = join(config['path']['all_interactions_extended'],
                                config['all_interactions_extended']['final_haca_guide_regions']),
    output:
        cd_boxes_psql = join(config['path']['psql'],
                             config['psql']['snoRNA_boxes'],
                             'cd_data_table.tsv'),
        haca_boxes_psql = join(config['path']['psql'],
                               config['psql']['snoRNA_boxes'],
                                    'haca_data_table.tsv'),
        sca_boxes_psql = join(config['path']['psql'],
                              config['psql']['snoRNA_boxes'],
                                    'sca_data_table.tsv'),
        snoRNA_boxes_script = join(config['path']['psql'],
                                   'snoRNA_boxes',
                                   'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_snoRNA_boxes.py"


# --------------------------------------------------------------------------

rule psql_rRNAs:
    input:
        rRNAs = join(config['path']['rRNA_modifications'],
                     config['rRNA_modifications']['rRNA']),
    output:
        rRNAs_psql = join(config['path']['psql'],
                          config['psql']['rRNAs'],
                          'data_table.tsv'),
        rRNAs_script = join(config['path']['psql'],
                           config['psql']['rRNAs'],
                           'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_rRNAs.py"


rule psql_rRNA_modifications:
    input:
        all_rRNA_modifications = join(config['path']['rRNA_modifications'],
                                      config['rRNA_modifications']['all_rRNA_modifications']),
    output:
        rRNA_modifications_psql = join(config['path']['psql'],
                                       config['psql']['rRNA_modifications'],
                                       'data_table.tsv'),
        rRNA_modifications_script = join(config['path']['psql'],
                                         config['psql']['rRNA_modifications'],
                                         'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_rRNA_modifications.py"


rule psql_conversion_18S:
    input:
        conversion_18S = join(config['path']['rRNA_modifications'],
                              config['rRNA_modifications']['conversion_18S']),
    output:
        conversion_18S_psql = join(config['path']['psql'],
                                       config['psql']['conversion_18S'],
                                       'data_table.tsv'),
        conversion_18S_script = join(config['path']['psql'],
                                         config['psql']['conversion_18S'],
                                         'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_conversion_18S.py"


rule psql_conversion_28S:
    input:
        conversion_28S = join(config['path']['rRNA_modifications'],
                              config['rRNA_modifications']['conversion_28S']),
    output:
        conversion_28S_psql = join(config['path']['psql'],
                                       config['psql']['conversion_28S'],
                                       'data_table.tsv'),
        conversion_28S_script = join(config['path']['psql'],
                                         config['psql']['conversion_28S'],
                                         'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_conversion_28S.py"
#____________________________________________________________________________
rule psql_rRNA_percentage_modification:
    input:
        all_rRNA_percentage_modification = join(config['path']['rRNA_modifications'],
                                             config['rRNA_modifications']['all_percentages_modifications']),
    output:
        rRNA_percentage_modification_psql = join(config['path']['psql'],
                                                 config['psql']['rRNA_percentage_modification'],
                                                 'data_table.tsv'),
        rRNA_percentage_modification_script = join(config['path']['psql'],
                                                   config['psql']['rRNA_percentage_modification'],
                                                   'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_rRNA_percentage_modification.py"


rule psql_rRNA_sample_percentage_modification:
    input:
        all_rRNA_sample_percentage_modification = join(config['path']['rRNA_modifications'],
                                                       config['rRNA_modifications']['sample_percentages_modifications']),
    output:
        rRNA_sample_percentage_modification_psql = join(config['path']['psql'],
                                                        config['psql']['rRNA_sample_percentage_modification'],
                                                        'data_table.tsv'),
        rRNA_sample_percentage_modification_script = join(config['path']['psql'],
                                                          config['psql']['rRNA_sample_percentage_modification'],
                                                          'data_script.sql'),
    params:
        host_script = 'scripts/psql_host.sh'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/psql_rRNA_sample_percentage_modification.py"




# --------- GRANT PERMISSION ON TABLES FOR scottweb_surfer ------------------

tables = [
     config['psql']['external_ids'],
     config['psql']['genomic_location'],
     config['psql']['basic_features'],
     config['psql']['specie'],
     config['psql']['host_features'],
     config['psql']['snoRNA_expression'],
     config['psql']['host_expression'],
     config['psql']['targets'],
     config['psql']['target_grouped_by'],
     config['psql']['encode_eclip'],
     config['psql']['rRNAs'],
     config['psql']['rRNA_modifications'],
     config['psql']['conversion_18S'],
     config['psql']['conversion_28S'],
     config['psql']['rRNA_percentage_modification'],
     config['psql']['rRNA_sample_percentage_modification'],
]
box_tables = [
    'cd_data_table',
    'haca_data_table',
    'sca_data_table'
]

rule grant_permission_on_tables:
    input:
        tables = expand('data/psql/{table}/data_table.tsv',
                        table=tables),
        box_tables = expand('data/psql/snoRNA_boxes/{box_tables}.tsv',
                            box_tables=box_tables)
    output:
        tok = 'data/psql/permission.tok'
    params:
        host_script = 'scripts/psql_host.sh',
        container_perm = 'psql_container_perm.sh',
    conda:
        "../envs/python.yaml"
    shell:
        "cd scripts && "
        "docker run "
        "-it "
        """--mount type=bind,source="$(pwd)",target=/sql """
        "psql13 "
        "'./sql/{params.container_perm}' && "
        "touch ../{output.tok}"
