source activate snakemake
set -e
# run all mmQTL jobs for maximum efficiency

# run in series due to massive numbers of files created.
#configs=$(input/GENCODE_transcript_config.yaml -n  input/isoseq_expression_config.yaml  input/isoseq_leafcutter_config.yaml  input/isoseq_SUPPA_RI_config.yaml  input/isoseq_transcript_config.yaml)

# make sure that the feature name in the config matches the dataCode in the config!
features=(
    GENCODE_expression 
    GENCODE_transcript 
    GENCODE_leafcutter
    GENCODE_SUPPA_AF
    GENCODE_SUPPA_RI
    GENCODE_SUPPA_AL
    GENCODE_SUPPA_SE
    GENCODE_SUPPA_A3
    GENCODE_SUPPA_A5
    #isoseq_expression 
    #isoseq_transcript 
    #isoseq_leafcutter
    #isoseq_SUPPA_RI 
    #isoseq_SUPPA_AF
    #isoseq_SUPPA_AL
    #isoseq_SUPPA_SE
    #isoseq_SUPPA_A3
    #isoseq_SUPPA_A5
    #union_expression
    #union_transcript
    #union_leafcutter
    #union_SUPPA_A3
    #union_SUPPA_A5
    #union_SUPPA_AL
    #union_SUPPA_AF
    #union_SUPPA_RI
    #union_SUPPA_SE
)

for i in ${features[@]};do
    echo $i
    snakejob -s Snakefile -c  input/${i}_config.yaml

    rm -rf results/$i/mmQTL/mmQTL_tmp
    rm -rf .snakemake/
    rm -rf cluster/*

done

