# mmQTL pipeline 
# Jack Humphrey, Erica Brophy, Kailash BP, Winston Cuddleston, Tatsuhiko Naito
# 2025 Raj lab
import glob
import pandas as pd 
import os 
import itertools

# MODULES AND SOFTWARE 
R_VERSION = config.get("R_version", "R/4.2.0")
PLINK_VERSION = config.get("PLINK_version", "plink2/2.3")
TABIX_VERSION = config.get("TABIX_version", "tabix/0.2.6")
BCFTOOLS_VERSION = config.get("BCFTOOLS_version", "bcftools/1.9")
GCTA_VERSION = config.get("GCTA_version", "gcta/1.94.1")
mmQTL_bin = config.get("mmQTL_bin", "/sc/arion/projects/bigbrain/MMQTL26a")

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

# all feature-specific and QTL mapping data go to specific outFolder
outFolder = config.get('outFolder', "results/")

# genotype data is always shared between analyses so put in separate shared folder
genoFolder = os.path.join(outFolder,"genotypes/")

# QTL-mapping settings
QTL_type = config.get("QTL_type", "cis")  # Default to "cis" if not defined, set "trans" to run trans-QTL pipeline
eQTL_number = config.get("eQTL_number", 1) # Default to primary QTL i.e. eQTL number = 1
variants_to_extract = config.get("variantsToExtract", "") # Default to all variants i.e. ""

phenoMeta = config.get('phenoMeta', "")
phenoMetaTrans = config.get('phenoMetaTrans', phenoMeta) # CHR START END FEEATURE of features to test for trans
dataCode = config.get('dataCode', "Default")

# set defaults for filtering phenotypes
# default should be TPM >1 in >= 50% of samples
# genes should be >1 TPM, transcripts and SUPPA should be >0.1
# do not filter leafcutter
pheno_threshold = config.get('phenoThreshold', 1)
pheno_fraction = config.get('phenoFraction', 0.5)

# minimum number of datasets a feature can appear in to be included in meta-analysis
min_datasets = config.get('minDatasets', 2)

outFolder = os.path.join(outFolder, dataCode) + "/"
print( " output folder = " + outFolder)

####################################################

SNAKEDIR = os.path.dirname(workflow.snakefile) + "/"

# Pipeline will run for each dataset in the dataKey created by the user
dataKey = config['dataKey']

# if group = True; then divide each feature by group total
group_features = bool(config.get('group', False))

meta = pd.read_csv(dataKey, sep = '\t') 
print(meta)
datasets = meta['dataset']
metadata_dict = meta.set_index("dataset").T.to_dict()

GTF = config.get("GTF", "")
prefix = outFolder + "{DATASET}/{DATASET}"
geno_prefix = genoFolder + "{DATASET}/{DATASET}"

pheno_matrix = prefix + "_pheno.tsv.gz"
leafcutter_string = ""
SUPPA_events = ""
mode_string = "normal"

## LEAFCUTTER SETTINGS
leafcutter = config.get('leafcutter', False)
if leafcutter == True:
    print(" * Leafcutter mode!")
    mode_string = "leafcutter"

## SUPPA SETTINGS
SUPPA = config.get('SUPPA', False)
if SUPPA == True:
    print(" * SUPPA mode")
    SUPPA_events = config["SUPPA_events"]
    pheno_matrix = prefix + ".SUPPA.phenotype_matrix.tsv.gz"
    mode_string = "SUPPA"

if leafcutter == True and SUPPA == True:
    sys.exit(" * You can't run both SUPPA and Leafcutter simultaneously!" )

## RNA EDITING SETTINGS
edqtl = config.get('edqtl', False)
if edqtl == True:
    print(" * edqtl mode")
    pheno_matrix = prefix + ".edqtl.phenotype_matrix.tsv.gz"

## SETTINGS FOR IF YOU'RE INCLUDING KNOWN COVARIATES AND PEER FACTORS
known_covars = config.get("known_covars", False) # If True, then provide a column in dataKey with path to covariates in long format

mmQTL_folder = outFolder + "mmQTL/"
mmQTL_tmp_folder = outFolder + "mmQTL/mmQTL_tmp/"

def write_list_to_file(my_list, file_name):
    with open(file_name, 'w') as f:
        for item in my_list:
            full_path = os.path.join(SNAKEDIR, item)
            f.write("%s\n" % full_path)

rule all:
    input:
        mmQTL_folder + dataCode + "_full_assoc.tsv.gz"

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
        participants = prefix + "_participants.txt"
    output:
        bed = temp(geno_prefix + "_genotypes.tmp2.bed"),
        bim = temp(geno_prefix + "_genotypes.tmp2.bim"),
        fam = temp(geno_prefix + "_genotypes.tmp2.fam"),
        afreq = temp(geno_prefix + "_genotypes.tmp2.afreq"),
        log = temp(geno_prefix + "_genotypes.tmp2.log"),
        chr_list = temp(geno_prefix + "_vcf_chr_list.txt")
    params:
        stem = geno_prefix + "_genotypes",
        blacklist = "scripts/Lifted_HighLDregion_hg38_RK_12_12_19.bed"
    run:
        vcf = metadata_dict[wildcards.DATASET]["genotypes"]
        shell("""
            ml {TABIX_VERSION} && tabix -l {vcf} > {output.chr_list}

            ml {PLINK_VERSION}
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

# 2.6 Extract select set of variants

rule extractVariants:
    input:
        geno_bed = geno_prefix + "_genotypes.tmp2.bed",
        geno_bim = geno_prefix + "_genotypes.tmp2.bim",
        geno_fam = geno_prefix + "_genotypes.tmp2.fam",
        variants_to_extract = variants_to_extract
    output:
        filtered_bed = temp(geno_prefix + "_genotypes.filtered.bed"),
        filtered_bim = temp(geno_prefix + "_genotypes.filtered.bim"),
        filtered_fam = temp(geno_prefix + "_genotypes.filtered.fam")
    params:
        stem = geno_prefix + "_genotypes"
    shell:
        """
        ml {PLINK_VERSION}

        echo "Running PLINK variant extraction..."

        if [ -s {input.variants_to_extract} ]; then
            plink2 --bfile {params.stem}.tmp2 \
                   --extract {input.variants_to_extract} \
                   --make-bed \
                   --out {params.stem}.filtered
        else
            echo "No variant list provided. Skipping extraction step."
            cp {params.stem}.tmp2.bed {params.stem}.filtered.bed
            cp {params.stem}.tmp2.bim {params.stem}.filtered.bim
            cp {params.stem}.tmp2.fam {params.stem}.filtered.fam
        fi
        """

# 2.6 Filter out singletons using the harmonized genotype bed file, and remove rare variants

rule removeSingletons:
    input:
        harmonized_genotype_bed = genoFolder + dataCode + "_combined_chr_pos_with_counts_without_singletons.bed",
        bed = geno_prefix + "_genotypes.filtered.bed",
        bim = geno_prefix + "_genotypes.filtered.bim",
        fam = geno_prefix + "_genotypes.filtered.fam"
    output:
        bed = geno_prefix + "_genotypes.bed",
        bim = geno_prefix + "_genotypes.bim",
        fam = geno_prefix + "_genotypes.fam",
        afreq = geno_prefix + "_genotypes.afreq"
    params:
        stem = geno_prefix + "_genotypes"
    shell:
        """
        ml {PLINK_VERSION}

        echo "Removing singletons and rare variants..."

        plink2 --make-bed \
            --output-chr chrM \
            --max-alleles 2 \
            --maf 0.01 \
            --freq \
            --max-maf 0.9975 \
            --extract range {input.harmonized_genotype_bed} \
            --keep-allele-order \
            --allow-extra-chr \
            --bfile {params.stem}.filtered \
            --out {params.stem}
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
            shell("ml {PLINK_VERSION}; \
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
        "ml {PLINK_VERSION}; ml {R_VERSION} ; ml {GCTA_VERSION};"
        "plink2 --bfile {params.stem} --maf 0.05 --make-bed --out {params.stem}_GCTA --human;"
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
        ml {TABIX_VERSION}; \
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
            shell("ml {R_VERSION}; Rscript {params.script} {input} {outFolder}/{wildcards.DATASET}/{wildcards.DATASET} {PEER_N}")
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
            shell("echo No Known covariates provided - check config or datakey if covariates were to be included; cp {input} {output}")

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
        shell("ml {R_VERSION}; Rscript {params.script} --pheno {input.pheno} --cov {input.cov} --out {output}")

#8. Harmonize phenotype files so that each file has the same features
# expand on dataset wildcard 
## this is the pheno file that will go into the path for pheno_file in the runMMQTL rule

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

# Using the harmonized phenotype metadata file as input, so that we only use harmonized features.
rule make_trans_pheno_meta:
    input:
        pheno_meta = phenoMeta, # original phenotype metadata file with all features
        pheno_meta_harmonized = mmQTL_folder + "phenotype_metadata.tsv", # harmonised phenotype metadata file with MIN_DATASET
        pheno_meta_trans_to_test = phenoMetaTrans # phenotype metadata file of features to test
    output:
        trans_pheno_meta = mmQTL_folder + "phenotype_metadata_trans.tsv"
    params:
        script = "scripts/make_trans_pheno_meta.R",
        out_folder = mmQTL_folder
    run:
        if QTL_type == "trans":  # Check if QTL_type is "trans"
            shell(
                  "ml {R_VERSION}; "
                  "Rscript {params.script} --outFolder {params.out_folder} "
                  " --pheno_meta {input.pheno_meta} "
                  " --pheno_meta_harmonized {input.pheno_meta_harmonized} "
                  " --pheno_meta_trans_to_test {input.pheno_meta_trans_to_test} ")
        else:
            print("Skipping 'make_trans_pheno_meta' because QTL_type is not 'trans'.")
            shell("touch {output.trans_pheno_meta}")  # Create an empty placeholder file

rule prepare_trans_phenotype:
    input:
        pheno = mmQTL_folder + "{DATASET}_pheno.regressed.harmonised.tsv",
        pheno_meta_trans = phenoMetaTrans
    output:
        trans_pheno_file = mmQTL_folder + "{DATASET}_pheno.regressed.harmonised.trans.tsv"
    params:
        script = "scripts/prepare_trans_phenotype.R",
        out_folder = mmQTL_folder
    run:
        if QTL_type == "trans":
            shell(
                "ml {R_VERSION}; "
                "Rscript {params.script} "
                "--pheno {input.pheno} "
                "--pheno_meta_trans {input.pheno_meta_trans} "
                "--outFolder {params.out_folder}"
            )
        else:
            print("Skipping 'prepare_trans_phenotype' because QTL_type is not 'trans'.")
            shell("touch {output.trans_pheno_file}")  # Create an empty file

rule cleanup_pheno_files:
    input:
        cis_pheno_files = expand(mmQTL_folder + "{DATASET}_pheno.regressed.harmonised.tsv", DATASET=datasets),
        trans_pheno_files = expand(mmQTL_folder + "{DATASET}_pheno.regressed.harmonised.trans.tsv", DATASET=datasets),
        trans_pheno_meta = mmQTL_folder + "phenotype_metadata_trans.tsv"
    output:
        cleanup_done = mmQTL_folder + "cleanup_complete.txt"  # Dummy file to signal cleanup completion
    run:
        if QTL_type == "cis":
            print("QTL_type is 'cis'. Removing trans-QTL files...")
            trans_files = " ".join(input.trans_pheno_files)  # Combine all trans files into one command
            shell(f"rm -f {trans_files} {input.trans_pheno_meta}")

        elif QTL_type == "trans":
            print("QTL_type is 'trans'. Removing cis-QTL files...")
            cis_files = " ".join(input.cis_pheno_files)  # Combine all cis files into one command
            shell(f"rm -f {cis_files}")

        # Create a dummy file to mark completion
        shell(f"echo 'Successfully cleaned up files' > {output.cleanup_done}")

## prepare inputs for mmQTL
rule prep_mmQTL:
        input:
           pheno = lambda wildcards: expand(mmQTL_folder + "{DATASET}_pheno.regressed.harmonised.trans.tsv", DATASET=datasets) if QTL_type == "trans" else expand(mmQTL_folder + "{DATASET}_pheno.regressed.harmonised.tsv", DATASET=datasets),
           geno = expand(geno_prefix + "_genotypes_{CHROM}.fam", DATASET = datasets, CHROM = chromosomes),
           grm = expand(geno_prefix + "_genotypes_GRM.tsv", DATASET = datasets),
           cleanup_done = mmQTL_folder + "cleanup_complete.txt"
        output:
           pheno_txt = mmQTL_folder + "pheno_list.txt",
           geno_txt = expand(mmQTL_folder + "{CHROM}_geno_list.txt", CHROM = chromosomes),
           grm_txt = mmQTL_folder + "grm_list.txt"
        run:
           # for genotypes, remove file extension
           plink_files = [ os.path.splitext(i)[0] for i in input.geno ]
           # for each chromosome write a list of genotype files
           for CHROM in chromosomes:
               geno_chrom_files = [genoFolder + d + "/" + d + "_genotypes_" + CHROM for d in datasets ]     
               geno_chrom_out = mmQTL_folder + CHROM + "_geno_list.txt"
               write_list_to_file( geno_chrom_files, geno_chrom_out )

           write_list_to_file(input.pheno, output.pheno_txt)
           write_list_to_file(input.grm, output.grm_txt)

#9. Run mmQTL 
rule runMMQTL: 
    input:
        pheno = mmQTL_folder + "pheno_list.txt",
        geno = expand(mmQTL_folder + "{CHROM}_geno_list.txt", allow_missing = True), #CHROM = chromosomes),
        grm = mmQTL_folder + "grm_list.txt",
        pheno_meta = lambda wildcards: mmQTL_folder + ("phenotype_metadata_trans.tsv" if QTL_type == "trans" else "phenotype_metadata.tsv")
    params:
        script = "scripts/run_mmQTL.R",
        prefix = mmQTL_tmp_folder,
        eQTL_number = eQTL_number
    output:
        mmQTL_tmp_folder + "{CHROM}_chunk_{CHUNK}_output.txt",
        mmQTL_tmp_folder + "{CHROM}_chunk_{CHUNK}_meta.tsv"
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
         --eQTL_number {params.eQTL_number} \
         --prefix {params.prefix} \
         --mmQTL {mmQTL_bin} \
         -i {wildcards.CHUNK} \
         -n {max_chunk} " )
         
rule gzip_results:
    input:
        chunk_success_check = mmQTL_tmp_folder + "{CHROM}_chunk_{CHUNK}_output.txt",
        chunk_meta = mmQTL_tmp_folder + "{CHROM}_chunk_{CHUNK}_meta.tsv"
    output:
        temp(mmQTL_tmp_folder + "{CHROM}_chunk_{CHUNK}_output_gzip.txt")
    params:
        script = "scripts/gzip_runMMQTL_output.R",
        eQTL_number = eQTL_number
    shell:
        "ml {R_VERSION}; Rscript {params.script} --chunk_meta {input.chunk_meta} --eQTL_number {params.eQTL_number} --output_file {output} "

#10. Collate mmQTL results

rule mmQTLcollate: 
    input:
        # use zip to zip together different numbers of chunk per chromosome
        outputs = expand(mmQTL_tmp_folder + "{CHROM}_chunk_{CHUNK}_output.txt", zip, CHUNK = chunk_zip, CHROM = chr_zip),
        outputs_gzip = expand(mmQTL_tmp_folder + "{CHROM}_chunk_{CHUNK}_output_gzip.txt", zip, CHUNK = chunk_zip, CHROM = chr_zip),
        meta = mmQTL_folder + "phenotype_metadata.tsv"
    output:
        temp(mmQTL_tmp_folder + "{CHROM}_top_assoc.tsv")
    params:
        script = "scripts/collate_mmQTL.R",
        prefix = mmQTL_tmp_folder
    run:
        if QTL_type == "cis":
            shell("""
                ml {R_VERSION};
                Rscript {params.script} --prefix {params.prefix} --chrom {wildcards.CHROM} --metadata {input.meta} --geno {genoFolder}
            """)
        elif QTL_type == "trans":
            print("Skipping 'mmQTLcollate' because QTL_type is not 'cis'.")
            shell("touch {output[0]}")

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
        gz =  mmQTL_folder + dataCode + "_full_assoc.tsv.gz",
        tbi = mmQTL_folder + dataCode + "_full_assoc.tsv.gz.tbi"
    params:
        tsv = mmQTL_folder + dataCode + "_full_assoc.tsv",
        prefix = mmQTL_tmp_folder
    shell:
        """
        set +o pipefail
        ml {BCFTOOLS_VERSION}
        ml {TABIX_VERSION}
        echo time for sorting
        chr=$(zcat {input} | tail -n +2 | cut -f3 | sort | uniq)
        for i in $chr; do echo $i
           zcat {params.prefix}/${{i}}_*_all_nominal.tsv.gz > {params.prefix}/${{i}}_full_assoc.tsv
           rm {params.prefix}/${{i}}_*_all_nominal.tsv.gz
           sort --parallel=4 -k 4,4n {params.prefix}/${{i}}_full_assoc.tsv > {params.prefix}/${{i}}_full_assoc.sorted.tsv
        done
        echo time for concatenating
        # qval column in top needs to be removed from header to make full assoc
        zcat {input} | head -1 | awk 'BEGIN{{OFS=\"\t\"}}NF{{NF-=1}};1' > {params.prefix}/{dataCode}_full_assoc_header.txt
        cat {params.prefix}/{dataCode}_full_assoc_header.txt {params.prefix}/chr*_full_assoc.sorted.tsv > {params.tsv}
        rm {params.prefix}/{dataCode}_full_assoc_header.txt
        rm {params.prefix}/chr*_full_assoc.tsv
        rm {params.prefix}/chr*_full_assoc.sorted.tsv
        bgzip {params.tsv}
        echo time for tabixing
        tabix -S 1 -s 3 -b 4 -e 4 {output.gz}
        find {params.prefix} -type d -empty -delete
        """
        
        

        
