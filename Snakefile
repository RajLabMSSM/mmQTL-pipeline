# mmQTL pipeline 08/05/21
# Erica Brophy

#1. Make sample key
#2. Filtered plink genotype file
#3. List of gene positions - command out for now 
#4. GRM
#5. Normalised gene expression matrix
#6. Covariate matrix
#7. Combine Covariate matrix
#8. Harmonize phenotype input files across dataset
#9. Run mmQTL for features across datasets in dataKey
#10. Collate all mmQTL results   

import glob
import pandas as pd 
import os 

inFolder = config['inFolder']
outFolder = config['outFolder']
sampleKey = config['sampleKey']
VCF = config['VCF']
VCFstem = VCF.split(".vcf.gz")[0]
geneMatrix = config['geneMatrix']
geneAnno = config['geneAnno']
dataCode = config['dataCode']
PEER_values = config["PEER_values"]

matrix = pd.read_csv(geneMatrix, sep = ' ')

genename = matrix['Gene'] 

prefix =  dataCode

covariateFile = config['covariateFile']

####################################################

SNAKEDIR = os.path.dirname(workflow.snakefile) + "/"

chromosomes = [i for i in range(1,23)]

#Pipeline will run for each dataset in the dataKey created by the user
dataKey = config['dataKey']
meta = pd.read_csv(dataKey, sep = ';') 
dataset = meta['dataset']
metadata_dict = meta.set_index("dataset").T.to_dict()
#expand on PEER

rule all:
    output:
        #expand(outFolder + "/peer{PEER_N}/" + dataCode + "_peer{PEER_N}.combined_covariates.txt", PEER_N = PEER_values)
        expand(outFolder + "/peer{PEER_N}/" + dataCode + "_peer{PEER_N}.PEER_covariates.txt", PEER_N = PEER_values)
 
#1. Make sample key 

rule getParticipants:
    input:
        txt = sampleKey
    output:
        txt = prefix + "_participants.txt"
    run:
        sk = pd.read_csv(input.txt, sep = "\t")
        participants = sk[["participant_id"]]
        participants.to_csv(output.txt, index = False, header = False, sep = "\t")

#2. Filtered plink genotype file

rule VCFtoPLINK:
    input:
        vcf = VCFstem + ".vcf.gz",
        participants = prefix + "_participants.txt"
    output:
        prefix + "_genotypes.fam", 
        prefix + "_genotypes.afreq"
    params:
        stem = prefix + "_genotypes"
    shell:
        "ml plink2; "
        "plink2 --make-bed "
        "--output-chr chrM "
        "--max-alleles 2 "
        "--keep {input.participants} "
        "--maf 0.01 "
        "--allow-extra-chr "
        "--max-maf 0.9975 "
        "--vcf {input.vcf} "
        "--out {params.stem} "

#3. List of gene positions

rule VCF_chr_list:
    input:
        VCF = VCFstem + ".vcf.gz",
        index = VCFstem + ".vcf.gz.tbi"
    output:
        prefix + "_vcf_chr_list.txt"
    shell:
        "tabix -l {input.VCF} > {output}"

#4. GRM

rule generateGRM:
    input:
        genotypes = prefix + "_genotypes.afreq"
    output:
        prefix + "_genotypes.grm"
    params:
        stem = prefix + "_genotypes"
    conda:
        SNAKEDIR + "envs/plink2.yaml" 
    shell:
        "plink2 --bfile {params.stem} "
        "--read-freq {input.genotypes} "
        "--make-grm-list -out {params.stem}"

#5. Normalised gene expression matrix

rule geneNORM:
    input:
        vcf_chr_list = prefix + "_vcf_chr_list.txt",
        geneMatrix = geneMatrix,
	geneAnno = geneAnno,
        sampleKey = sampleKey
    output:
        prefix + ".expression.bed.gz"
    params:
        script = "scripts/eqtl_prepare_expression.py"
    shell:
        "module load python/3.7.3;"
        " python {params.script} {input.geneAnno} {input.geneMatrix}"
        " {input.sampleKey} {input.vcf_chr_list} {prefix} "
        " --tpm_threshold 0.1 "
        " --count_threshold 6 "
        " --sample_frac_threshold 0.2 "
        " --normalization_method tmm "


#6. Covariate matrix
#expand on PEER wildcard

rule runPEER:
    input:
        expand(outFolder + "/{dataset}/" + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.PEER_covariates.txt", dataset = dataset, PEER_N = PEER_values, allow_missing=True ),
        phenotype_matrix = prefix + ".expression.bed.gz" 
    params:
        script = "scripts/run_PEER.R",
        num_peer = "{PEER_N}"
    output:
        outFolder + "/{dataset}/" + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.PEER_covariates.txt"
    run:
        if int(wildcards.PEER_N) > 0:
            shell("ml R/3.6.0; ")
            shell("Rscript {params.script} {input} {outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer} {params.num_peer}")
        else:
            shell("touch {output}")

#7. Combine covariates 

rule combineCovariates:
    input:
        pheno = prefix + ".expression.bed.gz",
        peer =  outFolder + "/{dataset}/" + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.PEER_covariates.txt",
        covariates = covariateFile
    output:
        cov_df = outFolder + "/{dataset}/" + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.combined_covariates.txt"
    params:
        num_peer = "{PEER_N}",
        script = "scripts/combine_covariates.py",
        logNomFolder = outFolder + "peer{PEER_N}/logNomFolder",
        logPerFolder = outFolder + "peer{PEER_N}/logPerFolder"
    run:
        if int(wildcards.PEER_N) > 0:
            peerFile = "--add_covariates " + input.peer
        else:
            peerFile = ""
        shell("python {params.script} {peerFile} \
            --covariates {input.covariates} \
            {outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer}")
        # make sure combined covariate file has column names in same order as phenotype file
        phenotype_df = pd.read_csv(input.pheno, sep='\t', dtype={'#chr':str, '#Chr':str})
        covariate_df = pd.read_csv(output.cov_df, sep = "\t" )
        covariate_df = covariate_df[ ["ID"] + list(phenotype_df.columns[4:]) ]
        # write out
        covariate_df.to_csv(output.cov_df, header = True, index = False, sep = "\t")

#8. Harmonize phenoype files so that each file has the same features
#expand on dataset wildcard 
##this is the pheno file that will go into the path for pheno_file in the runMMQTL rule

rule phenoUnion: 
    input:
        expand(outFolder + "/{dataset}", dataset = dataset, allow_missing=True ),
        pheno = prefix + ".expression.bed.gz"
    output:
         prefix + "harmonized.expression.bed.gz" 
    params:
        script = "scripts/pheno_harmonize.R"
    shell: 
        "ml R 3.6;"
        "Rscript {params.script} {input.pheno}"

#9. Run mmQTL 
#expand on PEER and phenotype

rule runMMQTL: 
    input:
        #expand on phenotype and PEER
        expand(outFolder + "/{dataset}/" + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.PEER_covariates.txt", dataset = dataset, PEER_N = PEER_values, allow_missing=True ),
        pheno_file = lambda wildcards: os.path.join(inFolder,  metadata_dict[wildcards.dataset]["Phenotype"]),
        geno_file = lambda wildcards: os.path.join(inFolder, metadata_dict[wildcards.dataset]["Genotype"]),
        GRM_file = lambda wildcards: os.path.join(inFolder, metadata_dict[wildcards.dataset]["GRM"]),
        feature_annotation = lambda wildcards: os.path.join(inFolder, metadata_dict[wildcards.dataset]["Feature_annotation"]),
        covariate_file = lambda wildcards: os.path.join(inFolder, metadata_dict[wildcards.dataset]["Covariate_file"]),
    output:
        #mmQTL output files 
        directory( outFolder + "/{dataset}/" + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.PEER_covariates.txt" + "{CHROM}")
    params:
        #script = "scripts/"
        #run mmQTL through src and outputs   
    shell:
        "cd /sc/arion/projects/ad-omics/data/software/MMQTL/MMQTL"
        "MMQTL23 -b  -P  pheno_file.text   -Z  geno_file.txt   -R GRM_file.txt \
        -a feature_annotation.bed  -A random  -C covariate_matrix.csv  -gene  gene_name"

#10. Collate mmQTL results

rule mmQTLcollate: 
    input:
        expand(outFolder + "/{dataset}/" + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.PEER_covariates.txt" + "{CHROM}", dataset = dataset, PEER_N = PEER_values, CHROM = chromosomes, allow_missing=True )
    output:
        directory( outFolder + "/{dataset}/" + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.PEER_covariates.txt")
    params:
        script = "scripts/merge_results.R" 
    shell:
        "ml R 3.6;"
        "Rscript {params.script}" 

