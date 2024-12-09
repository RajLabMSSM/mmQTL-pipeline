# mmQTL pipeline 
# Jack Humphrey & Erica Brophy
# 2021 Raj lab
import glob
import pandas as pd 
import os 
import itertools

R_VERSION = "R/4.2.0"
mmQTL_bin = "/sc/arion/projects/als-omics/microglia_isoseq/mmQTL-pipeline/MMQTL_bin/MMQTL26a"


# set chunk number
#N_CHUNKS = 20
#CHUNKS = range(1, N_CHUNKS + 1)

# how many chunks?
chunk_factor = 15

## new cool way to do chunking proportional to chromosome size!
# multiply per-chr chunk number by chunk factor
# expand df so that each chromosome is repeated by the number of chunks
# and the chunk value goes from 1 to N
# these two lists will then be the input to the expand statement and zipped together
# chr21 gets 1 unit, chr1 gets 7 - 7 times as many genes
chunk_df = pd.read_csv("scripts/gencode_v38_chr_chunk_weights.tsv", sep = "\t")

chunk_df = chunk_df.assign(chunk = chunk_df['ceil'] * chunk_factor )
chunks = list(chunk_df[ "chunk" ])
chrom = list(chunk_df["chr"])
chr_zip = []
chunk_zip = []
# the magic part -
for i in range(len(chunks)):
    chr_zip.extend( [ chrom[i] ] * chunks[i] )
    chunk_zip.extend( [j for j in range(1, chunks[i] + 1) ] )

chunk_dict = chunk_df.set_index("chr").T.to_dict()    


chromosomes = ["chr" + str(i) for i in range(1,23) ]
# for testing
#chromosomes = ["chr21"]

#print(chromosomes)

# all feature-specific and QTL mapping data go to specific outFolder
outFolder = config['outFolder']

# genotype data is always shared between analyses so put in separate shared folder
genoFolder = os.path.join(outFolder,"genotypes/")

if "phenoMeta" not in config.keys():
    config['phenoMeta'] = "" 
phenoMeta = config['phenoMeta']
dataCode = config['dataCode']

# set defaults for filtering phenotypes
# default should be TPM >1 in >= 50% of samples
# genes should be >1 TPM, transcripts and SUPPA should be >0.1
# do not filter leafcutter
if "phenoThreshold" not in config.keys():
    config['phenoThreshold'] = 1
if "phenoFraction" not in config.keys():
    config['phenoFraction'] = 0.5

# minimum number of datasets a feature can appear in to be included in meta-analysis
if "minDatasets" not in config.keys():
    config["minDatasets"] = 2

min_datasets = config["minDatasets"]
pheno_threshold = config['phenoThreshold']
pheno_fraction = config['phenoFraction']

outFolder = os.path.join(outFolder, dataCode) + "/"
print( " output folder = " + outFolder)
####################################################

SNAKEDIR = os.path.dirname(workflow.snakefile) + "/"


#Pipeline will run for each dataset in the dataKey created by the user
dataKey = config['dataKey']

# if group = True; then divide each feature by group total
if "group" not in config.keys():
    config["group"] = False

group_features = bool(config['group'])

meta = pd.read_csv(dataKey, sep = '\t') 
print(meta)
datasets = meta['dataset']
metadata_dict = meta.set_index("dataset").T.to_dict()
#expand on PEER
GTF = config["GTF"]
prefix = outFolder + "{DATASET}/{DATASET}"
geno_prefix = genoFolder + "{DATASET}/{DATASET}"

pheno_matrix = prefix + "_pheno.tsv.gz"
leafcutter_string = ""
SUPPA_events = ""
mode_string = "normal"

## LEAFCUTTER SETTINGS
if "leafcutter" not in config.keys():
    config["leafcutter"] = False
leafcutter = config['leafcutter']
if leafcutter == True:
    print(" * Leafcutter mode!")
    # TESTING
    #pheno_matrix = prefix + ".leafcutter.phenotype_matrix.tsv"
    mode_string = "leafcutter"

## SUPPA SETTINGS
if "SUPPA" not in config.keys():
    config["SUPPA"] = False
SUPPA = config['SUPPA']
if SUPPA == True:
    print(" * SUPPA mode")
    SUPPA_events = config["SUPPA_events"]
    pheno_matrix = prefix + ".SUPPA.phenotype_matrix.tsv.gz"
    mode_string = "SUPPA"

if leafcutter == True and SUPPA == True:
    sys.exit(" * You can't run both SUPPA and Leafcutter simultaneously!" )

## RNA EDITING SETTINGS
if "edqtl" not in config.keys():
    config["edqtl"] = False
edqtl = config['edqtl']
if edqtl == True:
    print(" * edqtl mode")
    pheno_matrix = prefix + ".edqtl.phenotype_matrix.tsv.gz"

## SETTINGS FOR IF YOU'RE INCLUDING KNOWN COVARIATES AND PEER FACTORS
if "known_covars" not in config.keys():
    config["known_covars"] = False
known_covars = config["known_covars"]

mmQTL_folder = outFolder + "mmQTL/"
mmQTL_tmp_folder = outFolder + "mmQTL/mmQTL_tmp/"

def write_list_to_file(my_list, file_name):
    with open(file_name, 'w') as f:
        for item in my_list:
            full_path = os.path.join(SNAKEDIR, item)
            f.write("%s\n" % full_path)

rule all:
    input:
        #expand(mmQTL_folder + "output/{CHROM}_chunk_{CHUNK}_output.txt", CHUNK = CHUNKS, CHROM = chromosomes )
        #expand(mmQTL_folder + "output/{CHROM}_top_assoc.tsv", CHROM = chromosomes)
        mmQTL_folder + dataCode + "_full_assoc.tsv.gz"
        #mmQTL_folder + dataCode + "_top_assoc.tsv.gz"
        #expand(outFolder + "/peer{PEER_N}/" + dataCode + "combined_covariates.txt", PEER_N = PEER_values)
        #expand(outFolder + "/peer{PEER_N}/" + dataCode + "PEER_covariates.txt", PEER_N = PEER_values)
        #expand(prefix + "_genotypes_GRM.tsv", DATASET = datasets ) 
        #expand(mmQTL_folder + "chunk_{CHUNK}_output.txt", CHUNK = CHUNKS )
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

# 2. Filtered plink genotype file

# 2.1 Convert VCF to plink, remove multi-allelic SNPs, blacklisted regions of genome, individuals not in sample key, variants with nan allele frequency, and maintain allele-order to prevent allele-flipping
rule VCFtoPLINK:
    input:
        # vcf = VCFstem + ".vcf.gz",
        participants = prefix + "_participants.txt"
        # participants = outFolder + "{DATASET}/{DATASET}_participants.txt"
    output:
        bed = geno_prefix + "_genotypes.tmp2.bed",
        bim = geno_prefix + "_genotypes.tmp2.bim",
        fam = geno_prefix + "_genotypes.tmp2.fam", 
        afreq = geno_prefix + "_genotypes.tmp2.afreq",
        chr_list = geno_prefix + "_vcf_chr_list.txt"
    params:
        stem = geno_prefix + "_genotypes",
        blacklist = "scripts/Lifted_HighLDregion_hg38_RK_12_12_19.bed"
    run:
        vcf = metadata_dict[wildcards.DATASET]["genotypes"]
        shell("""
            ml tabix && tabix -l {vcf} > {output.chr_list}

            ml plink2/2.3
            plink2 --make-bed \
                   --output-chr chrM \
                   --max-alleles 2 \
                   --geno 0.1 \
                   --keep-allele-order \
                   --keep {input.participants} \
                   --exclude range {params.blacklist} \
                   --allow-extra-chr \
                   --vcf {vcf} \
                   --out {params.stem}.tmp

            plink2 --make-bed \
                   --output-chr chrM \
                   --max-alleles 2 \
                   --keep-allele-order \
                   --freq \
                   --bfile {params.stem}.tmp \
                   --out {params.stem}.tmp2
        
        rm {params.stem}.tmp.fam {params.stem}.tmp.bed {params.stem}.tmp.bim {params.stem}.tmp.log
        """)

# 2.2 Extract CHR POS positions from plink file
rule extract_chr_pos:
    input:
        bim = geno_prefix + "_genotypes.tmp2.bim"
    output:
        chr_pos = genoFolder + "{DATASET}/{DATASET}_genotypes_chr_pos.tsv",
    shell:
        """
        awk '{{print $1, $4}}' {input.bim} > {output.chr_pos}
        """

# 2.3 Combine all CHR POS files into one file. 
rule combine_chr_pos:
    input:
        chr_pos = expand(genoFolder + "{DATASET}/{DATASET}_genotypes_chr_pos.tsv", DATASET = datasets),
    output:
        combined_genotype_bed = genoFolder + dataCode + "_combined_chr_pos.bed"
    shell:
        """
        cat {input.chr_pos} >> {output.combined_genotype_bed}
        """

# 2.4 Count number of occurrences of each allele across DATASETS
rule count_chr_pos:
    input:
        combined_genotype_bed = genoFolder + dataCode + "_combined_chr_pos.bed"
    output: 
        combined_genotype_bed_with_counts = genoFolder + dataCode + "_combined_chr_pos_with_counts.bed",
        harmonized_genotype_bed = genoFolder + dataCode + "_combined_chr_pos_with_counts_without_singletons.bed"
    shell:
        """
        sort {input.combined_genotype_bed} | uniq -c | sort -nr > {output.combined_genotype_bed_with_counts}
        awk '{{if ($1 > 1) print $2 " " $3 " " $3}}' {output.combined_genotype_bed_with_counts} > {output.harmonized_genotype_bed} # separated CHR POS POS by using a space
        """

# 2.5 Filter out singletons using the harmonized genotype bed file, and remove rare variants
rule removeSingletons:
    input:
        harmonized_genotype_bed = genoFolder + dataCode + "_combined_chr_pos_with_counts_without_singletons.bed",
        bed = geno_prefix + "_genotypes.tmp2.bed",
        bim = geno_prefix + "_genotypes.tmp2.bim",
        fam = geno_prefix + "_genotypes.tmp2.fam", 
        afreq = geno_prefix + "_genotypes.tmp2.afreq"
    output:
        bed = geno_prefix + "_genotypes.bed",
        bim = geno_prefix + "_genotypes.bim",
        fam = geno_prefix + "_genotypes.fam", 
        afreq = geno_prefix + "_genotypes.afreq",
    params:
        stem = geno_prefix + "_genotypes"
    shell:
        """
        ml plink2/2.3

        plink2 --make-bed \
               --output-chr chrM \
               --max-alleles 2 \
               --maf 0.01 \
               --freq \
               --max-maf 0.9975 \
               --extract range {input.harmonized_genotype_bed} \
               --keep-allele-order \
               --allow-extra-chr \
               --bfile {params.stem}.tmp2 \
               --out {params.stem} 
        
        rm {params.stem}.tmp2.fam {params.stem}.tmp2.bed {params.stem}.tmp2.bim {params.stem}.tmp2.log {params.stem}.tmp2.afreq
        """

#3. Split plink file into chromosomes
rule splitPlinkChr:
    input: 
        fam = geno_prefix + "_genotypes.fam"
    params:
        stem = geno_prefix + "_genotypes"
    output: 
        expand( geno_prefix + "_genotypes_{CHROM}.fam", CHROM = chromosomes, allow_missing = True)
    run:
        for chrom in chromosomes:
            shell("ml plink2/2.3; \
            plink2 --bfile {params.stem} --make-bed --chr {chrom} --out {params.stem}_{chrom}" )

#4. GRM

rule generateGRM:
    input:
        genotypes = geno_prefix + "_genotypes.afreq"
    output:
        geno_prefix + "_genotypes_GRM.tsv"
    params:
        stem = geno_prefix + "_genotypes",
        script = "scripts/process_GRM.R"
    shell:
        "ml plink; ml {R_VERSION} ; ml gcta;"
        "plink --bfile {params.stem} --maf 0.05 --output-chr 26 --make-bed --out {params.stem}_GCTA;"
        "gcta64 --bfile {params.stem}_GCTA  --autosome --maf 0.05 --make-grm --out {params.stem}_GRM  --thread-num 20;"
        " Rscript {params.script} --prefix {params.stem}_GRM  "

#5. Normalise phenotype matrix
rule prepare_phenotypes:
    output:
        prefix + "_pheno.tsv.gz"
    params:
        pheno_meta = phenoMeta,
        script = "scripts/prepare_phenotypes.R"
    run:
        print(metadata_dict[wildcards.DATASET])
        pheno = metadata_dict[wildcards.DATASET]["phenotypes"]
        counts_string = ""
        if "counts" in meta.columns:
            counts_string = "--counts " + metadata_dict[wildcards.DATASET]["counts"]
        sk = metadata_dict[wildcards.DATASET]["sample_key"]
        #pheno_meta = metadata_dict[wildcards.DATASET]["phenotype_info"]
        threshold = pheno_threshold,
        fraction = pheno_fraction,
        group_string = ""
        if( group_features == True):
            group_string = " --group "
        
        shell("ml {R_VERSION};\
        Rscript {params.script} \
        --key {sk} \
        --pheno_matrix {pheno} \
        --pheno_meta {params.pheno_meta} \
        {counts_string} \
        {group_string} \
        --threshold {threshold} \
        --fraction {fraction} \
        --prefix {outFolder}/{wildcards.DATASET}/{wildcards.DATASET} ")

## Special case - Leafcutter junctions
rule prepare_leafcutter:
    input:
        # expecting gzipped junction files with extension {sample}.junc.gz
        exon_list = phenoMeta, # hg38 exons from gencode v30 with gene_id and gene_name
        genes_gtf = GTF # full GTF or just gene starts and ends?
    output:
        pheno_matrix = prefix + ".leafcutter.phenotype_matrix.tsv",
        pheno_meta= prefix + ".leafcutter.phenotype_meta.tsv"
    params: 
        leafcutter_dir = os.getcwd() + "/scripts/leafcutter/", # all leafcutter scripts hosted in a folder - some had to be converted py2 -> py3
        script = os.getcwd() + "/scripts/sqtl_prepare_splicing.py",
        min_clu_reads = 30,
        min_clu_ratio = 0.001,
        max_intron_len = 100000, # cut down to 100k to reduce SNP testing distance
        num_pcs = 10, # must be at least the number of samples!
        coord_mode = "cluster_middle"
        #coord_mode = "cluster_middle" # set coordinates to either "TSS" or "cluster_middle"
    run:
        junc_list = metadata_dict[wildcards.DATASET]["phenotypes"]
        sample_key = metadata_dict[wildcards.DATASET]["sample_key"] 
        shell(
        "ml {R_VERSION}; \
        ml tabix; \
        python {params.script}  \
                 {junc_list}  \
                 {input.exon_list}  \
                 {input.genes_gtf}  \
                 {wildcards.DATASET}  \
                 {sample_key}  \
         --min_clu_reads {params.min_clu_reads}  \
         --min_clu_ratio {params.min_clu_ratio}  \
         --max_intron_len {params.max_intron_len}  \
         --coord_mode {params.coord_mode}  \
         --num_pcs {params.num_pcs} \
         --leafcutter_dir {params.leafcutter_dir};  \
         mv -f {wildcards.DATASET}* {outFolder}/{wildcards.DATASET} " 
       )

rule prepare_SUPPA:
    input:
        SUPPA_events    
    output:
        pheno_matrix = prefix + ".SUPPA.phenotype_matrix.tsv.gz",
        pheno_meta = prefix + ".SUPPA.phenotype_meta.tsv.gz"
    params:
        script = "scripts/prepare_SUPPA.R"
    run:
        tx_matrix = metadata_dict[wildcards.DATASET]["phenotypes"]
        sample_key = metadata_dict[wildcards.DATASET]["sample_key"]
        shell(
        "ml {R_VERSION}; \
         Rscript {params.script} --pheno {tx_matrix} --key {sample_key} --events {input} --prefix {outFolder}{wildcards.DATASET}/{wildcards.DATASET}")

rule prepare_edqtl:
    output:
        pheno_matrix = prefix + ".edqtl.phenotype_matrix.tsv.gz"
    params:
        script = "scripts/prepare_edqtl.R"
    run:
        pheno = metadata_dict[wildcards.DATASET]["phenotypes"]
        sk = metadata_dict[wildcards.DATASET]["sample_key"]
        threshold = pheno_threshold,
        fraction = pheno_fraction,
        meta = phenoMeta,
        shell(
        "ml {R_VERSION}; \
        Rscript {params.script} \
        --key {sk} \
        --pheno_matrix {pheno} \
        --pheno_meta {meta} \
        --threshold {threshold} \
        --fraction {fraction} \
        --prefix {outFolder}/{wildcards.DATASET}/{wildcards.DATASET} "
        )

#6. Covariate matrix
# using PEER
rule run_PEER:
    input:
        phenotype_matrix = pheno_matrix
    params:
        script = "scripts/run_PEER_mmQTL.R",
    output:
        prefix  + "_PEER_covariates.txt"
    run:
        PEER_N = metadata_dict[wildcards.DATASET]["PEER"]
        if int(PEER_N) > 0:
            shell("ml R/4.0.3; Rscript {params.script} {input} {outFolder}/{wildcards.DATASET}/{wildcards.DATASET} {PEER_N}")
        else:
            shell("touch {output}")

# add known covariates to PEER factors
rule merge_covariates:
    input:
        prefix + "_PEER_covariates.txt"
    output:
        prefix + "_covariates.txt"
    params:
        script = "scripts/merge_covariates.R"
    run:
        covar_file = metadata_dict[wildcards.DATASET]["covariates"]
        if known_covars == True:
           shell(
            "ml {R_VERSION}; \
            Rscript {params.script} \
            --PEER_cov {input} \
            --known_cov {covar_file} \
            --prefix {outFolder}/{wildcards.DATASET}/{wildcards.DATASET} "
            )
        else:
            shell("cp {input} {output}")

# regress covariates from phenotype matrix
rule regress_covariates:
    input:
        pheno = pheno_matrix,
        cov = prefix + "_covariates.txt"
    output:
        prefix + "_pheno.regressed.tsv.gz"
    params:
        script = "scripts/regress_covariates_factor_sex_age.R"
    run:
        #known_covars = config["known_covars"]
        #if known_covars == True:
            shell("ml {R_VERSION}; Rscript {params.script} --pheno {input.pheno} --cov {input.cov} --out {output}")
        #else:
        #    shell("cp {input.pheno} {output}")

#8. Harmonize phenotype files so that each file has the same features
#expand on dataset wildcard 
##this is the pheno file that will go into the path for pheno_file in the runMMQTL rule

rule harmonise_phenotypes: 
    input:
        pheno_meta = phenoMeta,
        pheno = expand(prefix + "_pheno.regressed.tsv.gz", DATASET = datasets ),
    output:
         expand(mmQTL_folder +  "{DATASET}_pheno.regressed.harmonised.tsv", DATASET = datasets),
         mmQTL_folder + "phenotype_metadata.tsv"
    params:
        script = "scripts/pheno_harmonize.R",
        prefix = mmQTL_folder
    shell: 
        "ml {R_VERSION};"
        "Rscript {params.script} --prefix {params.prefix} --metadata {input.pheno_meta} --mode {mode_string} {input.pheno} --min_datasets {min_datasets}"

## prepare inputs for mmQTL
rule prep_mmQTL:
        input:
           pheno = expand(mmQTL_folder +  "{DATASET}_pheno.regressed.harmonised.tsv", DATASET = datasets),
           geno = expand(geno_prefix + "_genotypes_{CHROM}.fam", DATASET = datasets, CHROM = chromosomes),
           grm = expand(geno_prefix + "_genotypes_GRM.tsv", DATASET = datasets)
           #cov = expand(prefix + "_PEER_mmQTL.txt", DATASET = datasets)
        output:
           pheno_txt = mmQTL_folder + "pheno_list.txt",
           geno_txt = expand(mmQTL_folder + "{CHROM}_geno_list.txt", CHROM = chromosomes),
           grm_txt = mmQTL_folder + "grm_list.txt"
           #cov_txt = mmQTL_folder + "cov_list.txt"
        run:
           # for genotypes, remove file extension
           plink_files = [ os.path.splitext(i)[0] for i in input.geno ]
           # for each chromosome write a list of genotype files
           for CHROM in chromosomes:
               geno_chrom_files = [genoFolder + d + "/" + d + "_genotypes_" + CHROM for d in datasets ]     
               geno_chrom_out = mmQTL_folder + CHROM + "_geno_list.txt"
               write_list_to_file( geno_chrom_files, geno_chrom_out )

           write_list_to_file(input.pheno, output.pheno_txt)
           #write_list_to_file(plink_files, output.geno_txt)
           write_list_to_file(input.grm, output.grm_txt)
           #write_list_to_file(input.cov, output.cov_txt) 

#9. Run mmQTL 
rule runMMQTL: 
    input:
        pheno = mmQTL_folder + "pheno_list.txt",
        geno = expand(mmQTL_folder + "{CHROM}_geno_list.txt", allow_missing = True), #CHROM = chromosomes),
        grm = mmQTL_folder + "grm_list.txt",
        #cov = mmQTL_folder + "cov_list.txt",
        pheno_meta = mmQTL_folder + "phenotype_metadata.tsv"
    params:
        script = "scripts/run_mmQTL.R",
        prefix = mmQTL_tmp_folder
    output:
        mmQTL_tmp_folder + "{CHROM}_chunk_{CHUNK}_output.txt"
    run:
        max_chunk = chunk_dict[wildcards.CHROM]['chunk']
        shell(
        "ml {R_VERSION};\
        Rscript {params.script} \
         --chrom {wildcards.CHROM} \
         --pheno_file {input.pheno} \
         --geno_file {mmQTL_folder}/{wildcards.CHROM}_geno_list.txt \
         --grm_file {input.grm} \
         --pheno_meta {input.pheno_meta} \
         --prefix {params.prefix} \
         --mmQTL {mmQTL_bin} \
         -i {wildcards.CHUNK} \
         -n {max_chunk} " )
#10. Collate mmQTL results

rule mmQTLcollate: 
    input:
        # use zip to zip together different numbers of chunk per chromosome
        outputs = expand(mmQTL_tmp_folder + "{CHROM}_chunk_{CHUNK}_output.txt", zip, CHUNK = chunk_zip, CHROM = chr_zip ),
        meta = mmQTL_folder + "phenotype_metadata.tsv"
    output:
        mmQTL_tmp_folder + "{CHROM}_top_assoc.tsv"
    params:
        script = "scripts/collate_mmQTL.R",
        prefix = mmQTL_tmp_folder
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script} --prefix {params.prefix} --chrom {wildcards.CHROM} --metadata {input.meta} --geno {genoFolder}" 

rule topCollate:
    input:
        expand(mmQTL_tmp_folder + "{CHROM}_top_assoc.tsv", CHROM = chromosomes)
    output:
        mmQTL_folder + dataCode + "_top_assoc.tsv.gz"
    params:
        script = "scripts/collate_top_chrom.R"
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script} --output {output} {input}"

rule fullCollate:
    input:
        mmQTL_folder + dataCode + "_top_assoc.tsv.gz"
    output:
        #tsv = mmQTL_folder + dataCode + "_full_assoc.tsv",
        gz =  mmQTL_folder + dataCode + "_full_assoc.tsv.gz",
        tbi = mmQTL_folder + dataCode + "_full_assoc.tsv.gz.tbi"
    params:
        tsv = mmQTL_folder + dataCode + "_full_assoc.tsv",
        prefix = mmQTL_tmp_folder
    shell:
        "set +o pipefail;"
        "ml bcftools/1.9;"
        "echo time for sorting;"
        "for i in {chromosomes}; do echo $i;"
        "   cat {params.prefix}/${{i}}_*_all_nominal.tsv > {params.prefix}/${{i}}_full_assoc.tsv;"
        "   sort --parallel=4 -k 4,4n {params.prefix}/${{i}}_full_assoc.tsv > {params.prefix}/${{i}}_full_assoc.sorted.tsv;"
        "done;"
        "echo time for concatenating;"
        # qval column in top needs to be removed from header to make full assoc
        "zcat {input} | head -1 | awk 'BEGIN{{OFS=\"\t\"}}NF{{NF-=1}};1' > {params.prefix}/{dataCode}_full_assoc_header.txt;"
        "echo  {params.prefix}/chr{{1..22}}_full_assoc.sorted.tsv ;"
        "cat {params.prefix}/{dataCode}_full_assoc_header.txt {params.prefix}/chr{{1..22}}_full_assoc.sorted.tsv > {params.tsv};"
        "bgzip {params.tsv};"
        "echo time for tabixing;"
        "tabix -S 1 -s 3 -b 4 -e 4 {output.gz} "
