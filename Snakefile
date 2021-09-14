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

#inFolder = config['inFolder']
outFolder = config['outFolder']
phenoMeta = config['phenoMeta']
#VCF = config['VCF']
#VCFstem = VCF.split(".vcf.gz")[0]
#geneMatrix = config['geneMatrix']
#geneAnno = config['geneAnno']
dataCode = config['dataCode']
#PEER_values = config["PEER_values"]


#genename = matrix['Gene'] 

prefix = outFolder + "/" + dataCode

#covariateFile = config['covariateFile']

####################################################

SNAKEDIR = os.path.dirname(workflow.snakefile) + "/"

#chromosomes = ["chr" + str(i) for i in range(1,23)]

N_CHUNKS = 1
CHUNKS = range(1, N_CHUNKS + 1)

#Pipeline will run for each dataset in the dataKey created by the user
dataKey = config['dataKey']

meta = pd.read_csv(dataKey, sep = '\t') 
print(meta)
datasets = meta['dataset']
metadata_dict = meta.set_index("dataset").T.to_dict()
#expand on PEER

prefix = outFolder + "{DATASET}/{DATASET}"

mmQTL_folder = outFolder + "mmQTL/"

def write_list_to_file(my_list, file_name):
    with open(file_name, 'w') as f:
        for item in my_list:
            full_path = os.path.join(SNAKEDIR, item)
            f.write("%s\n" % full_path)

rule all:
    input:
        #expand(outFolder + "/peer{PEER_N}/" + dataCode + "combined_covariates.txt", PEER_N = PEER_values)
        #expand(outFolder + "/peer{PEER_N}/" + dataCode + "PEER_covariates.txt", PEER_N = PEER_values)
        #expand(prefix + "_genotypes.grm", DATASET = datasets ) 
        expand(mmQTL_folder + "chunk_{CHUNK}_output.txt", CHUNK = CHUNKS )
        #mmQTL_folder + "all_results_collated.txt"

#1. Make sample key 

rule getParticipants:
    output:
        txt = prefix + "_participants.txt"
    run:
        sample_key = metadata_dict[wildcards.DATASET]["sample_key"]
        sk = pd.read_csv(sample_key, sep = "\t")
        participants = sk[["participant_id"]]
        participants.to_csv(output.txt, index = False, header = False, sep = "\t")

#2. Filtered plink genotype file

rule VCFtoPLINK:
    input:
        #vcf = VCFstem + ".vcf.gz",
        participants = prefix + "_participants.txt"
        #participants  = outFolder + "{DATASET}/{DATASET}_participants.txt"
    output:
        fam = prefix + "_genotypes.fam", 
        afreq = prefix + "_genotypes.afreq",
        chr_list = prefix + "_vcf_chr_list.txt"
    params:
        stem = prefix + "_genotypes"
    run:
        vcf = metadata_dict[wildcards.DATASET]["genotypes"]
        shell("ml tabix; \
            tabix -l {vcf} > {output.chr_list}"
        )
        shell("ml plink2; \
        plink2 --make-bed \
        --output-chr chrM \
        --max-alleles 2 \
        --keep {input.participants} \
        --maf 0.01 \
        --freq \
        --allow-extra-chr \
        --max-maf 0.9975 \
        --vcf {vcf} \
        --out {params.stem}"
        )

#4. GRM

rule generateGRM:
    input:
        genotypes = prefix + "_genotypes.afreq"
    output:
        prefix + "_genotypes.grm"
    params:
        stem = prefix + "_genotypes"
    shell:
        "ml plink2; "
        "plink2 --bfile {params.stem} "
        "--read-freq {input.genotypes} "
        "--make-grm-list -out {params.stem}"

#5. Normalise phenotype matrix

rule normalise_pheno:
    input:
        #vcf_chr_list = prefix + "_vcf_chr_list.txt",
        #geneMatrix = geneMatrix,
	    #geneAnno = geneAnno,
        #sampleKey = sampleKey
    output:
        prefix + "_pheno.tsv.gz"
    params:
        pheno_meta = phenoMeta,
        script = "scripts/prepare_phenotypes.R"
    run:
        pheno = metadata_dict[wildcards.DATASET]["phenotypes"]
        sk = metadata_dict[wildcards.DATASET]["sample_key"]
        #pheno_meta = metadata_dict[wildcards.DATASET]["phenotype_info"]
        threshold = 1,
        fraction = 0.5,
        group = False
        group_string = ""
        if( group == True):
            group_string = " --group "
        
        shell("ml R/4.0.3;\
        Rscript {params.script} \
        --key {sk} \
        --pheno_matrix {pheno} \
        --pheno_meta {params.pheno_meta} \
        {group_string} \
        --threshold {threshold} \
        --fraction {fraction} \
        --prefix {outFolder}/{wildcards.DATASET}/{wildcards.DATASET} ")
        #module load python/3.7.3;"
        #" python {params.script} {input.geneAnno} {input.geneMatrix}"
        #" {input.sampleKey} {input.vcf_chr_list} {prefix} "
        " --tpm_threshold 0.1 "
        " --count_threshold 6 "
        " --sample_frac_threshold 0.2 "
        " --normalization_method tmm "


#6. Covariate matrix
#expand on PEER wildcard

rule runPEER:
    input:
        phenotype_matrix = prefix + "_pheno.tsv.gz"
    params:
        script = "scripts/run_PEER.R",
    output:
        prefix  + ".PEER_covariates.txt"
    run:
        PEER_N = metadata_dict[wildcards.DATASET]["PEER"]
        if int(PEER_N) > 0:
            shell("ml R/3.6.0; ")
            shell("Rscript {params.script} {input} {outFolder}/{wildcards.DATASET}/{wildcards.DATASET} {PEER_N}")
        else:
            shell("touch {output}")

#7. Combine covariates 

rule combineCovariates:
    input:
        pheno = prefix + "_pheno.tsv.gz",
        peer = prefix  + ".PEER_covariates.txt"
        #peer =  outFolder + "/{dataset}/"  + dataCode + "PEER_covariates.txt",
        #covariates = ""
        #covariates = covariateFile
    output:
        cov_df = prefix + "_combined_covariates.txt"
    params:
        script = "scripts/combine_covariates.py",
    run:
        PEER_N = metadata_dict[wildcards.DATASET]["PEER"]
        if int(PEER_N) > 0:
            peerFile = "--add_covariates " + input.peer
        else:
            peerFile = ""
        #shell("python {params.script} {peerFile} \
        #    --covariates {input.covariates} \
        #    {outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer}")
        # make sure combined covariate file has column names in same order as phenotype file
        phenotype_df = pd.read_csv(input.pheno, sep='\t', dtype={'#chr':str, '#Chr':str})
        covariate_df = pd.read_csv(output.cov_df, sep = "\t" )
        covariate_df = covariate_df[ ["ID"] + list(phenotype_df.columns[4:]) ]
        # write out
        covariate_df.to_csv(output.cov_df, header = True, index = False, sep = "\t")

#8. Harmonize phenoype files so that each file has the same features
#expand on dataset wildcard 
##this is the pheno file that will go into the path for pheno_file in the runMMQTL rule

rule harmonise_phenotypes: 
    input:
        pheno_meta = phenoMeta,
        pheno = expand(prefix + "_pheno.tsv.gz", DATASET = datasets ),
    output:
         expand(mmQTL_folder +  "{DATASET}_pheno.harmonised.tsv", DATASET = datasets),
         mmQTL_folder + "phenotype_metadata.tsv"
    params:
        script = "scripts/pheno_harmonize.R",
        prefix = mmQTL_folder
    shell: 
        "ml R/4.0.3;"
        "Rscript {params.script} --prefix {params.prefix} --metadata {input.pheno_meta}  {input.pheno}"

## prepare inputs for mmQTL
rule prep_mmQTL:
        input:
           pheno = expand(mmQTL_folder +  "{DATASET}_pheno.harmonised.tsv", DATASET = datasets),
           geno = expand(prefix + "_genotypes.fam", DATASET = datasets),
           cov = expand(prefix + ".PEER_covariates.txt", DATASET = datasets),
           grm = expand(prefix + "_genotypes.grm", DATASET = datasets)
        output:
           pheno_txt = mmQTL_folder + "pheno_list.txt",
           geno_txt = mmQTL_folder + "geno_list.txt",
           grm_txt = mmQTL_folder + "grm_list.txt",
           cov_txt = mmQTL_folder + "cov_list.txt"
        run:
           # for genotypes, remove file extension
           plink_files = [ os.path.splitext(i)[0] for i in input.geno ]

           write_list_to_file(input.pheno, output.pheno_txt)
           write_list_to_file(plink_files, output.geno_txt)
           write_list_to_file(input.grm, output.grm_txt)
           write_list_to_file(input.cov, output.cov_txt) 

#9. Run mmQTL 
#expand on PEER and phenotype

rule runMMQTL: 
    input:
        pheno = mmQTL_folder + "pheno_list.txt",
        geno = mmQTL_folder + "geno_list.txt",
        grm = mmQTL_folder + "grm_list.txt",
        cov = mmQTL_folder + "cov_list.txt",
        pheno_meta = mmQTL_folder + "phenotype_metadata.tsv"
    params:
        script = "scripts/run_mmQTL.R",
        prefix = mmQTL_folder + "results/"
    output:
        mmQTL_folder + "chunk_{CHUNK}_output.txt"
        #directory( outFolder + "/{dataset}/"  + dataCode + "PEER_covariates.txt" + "{CHROM}")
    shell:
        "ml R/4.0.3;"
        "Rscript {params.script} "
        " --pheno_file {input.pheno} "
        " --geno_file {input.geno} "
        " --grm_file {input.grm} "
        " --cov_file {input.cov} "
        " --pheno_meta {input.pheno_meta} "
        " --prefix {mmQTL_folder} "
        " -i {wildcards.CHUNK} "
        " -n {N_CHUNKS} "

#10. Collate mmQTL results

rule mmQTLcollate: 
    input:
        expand(mmQTL_folder + "chunk_{CHUNK}_output.txt", CHUNK = CHUNKS )
    output:
        mmQTL_folder + "all_results_collated.txt"
    params:
        script = "scripts/merge_results.R" 
    shell:
        "ml R 3.6;"
        #"Rscript {params.script}" 

