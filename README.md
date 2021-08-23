# mmQTL-pipeline

Pipeline will be wrapped using a **dataset key** created by the user in the format of a dataframe. The dataset key contains 5 rows: **1. Dataset, 2. sampleKey, 3. VCF, 4. Phenotype file, 5. Covariate file**; where the column provides the path to each of these files. 

### User input data: 
The inputs required by the user are a *sampleKey* matching RNA samples to DNA samples, a *VCF file* containing the genotype information, a *gene expression matrix* (TPM or similar) and a *gene annotation file* (GENCODE or similar). 

### The snakemake pipeline will create the following files needed to run MMQTL:
1. Filtered plink genotype file
2. List of gene positions - need this for PLINK 
3. GRM
4. Normalized gene expression matrix
5. Covariate matrix
