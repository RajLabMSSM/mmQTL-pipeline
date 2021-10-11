# mmQTL-pipeline

Pipeline will be wrapped using a **dataset key** created by the user in the format of a dataframe. The dataset key contains 5 rows: **1. Dataset, 2. sampleKey, 3. VCF, 4. Phenotype file, 5. Covariate file**; where the column provides the path to each of these files. 

### User input data: 
The inputs required by the user are a *sampleKey* matching RNA samples to DNA samples, a *VCF file* containing the genotype information, a *gene expression matrix* (TPM or similar) and a *gene annotation file* (GENCODE or similar). 

The *sample key* must have two columns named participant_id (listing the names of the genotype IDs) and sample_id (listing the name of the phenotype sample IDs)

The *genotypes* must be in VCF format in a single file

The *phenotype matrix* must be a tab-separated table with one feature per row, the first column should be the unique identifier of the feature, followed by each sample ID in the dataset. The matrix should contain normalised values for each sample. See below for types of phenotypes that mmQTL can work with.

The *PEER factor* should be a number of PEER factors to be regressed out of the phenotype matrix in each dataset. Ideally this has already been determined using the `QTL-mapping-pipeline`.



## Types of phenotypes used for mmQTL:

* Gene expression

This must be TPM normalised. A phenotype metadata file should be created. This can be created from a GTF file using `scripts/get_gene_meta_from_gtf.R` in "gene" mode.

* Transcripts

The user must prepare a transcript expression matrix with TPM normalised transcript expression, from Kallisto, RSEM or similar.

If the user wishes to treat each transcript as it's own feature - (transcript QTLs), in the config.yaml set `group: False`

If the user wishes to divide each transcript by the total gene TPM, (transcript usage QTLs), set `group: True`

Transcript metadata must be created from the same GTF as used to create the transcript expression matrix using `scripts/get_gene_meta_from_gtf.R` in "transcript" mode.

* Other features (SUPPA, txrevise etc)

To be described.

* Junctions

If the user wants to perform junction usage QTLs (Leafcutter splicing QTLs) then provide a list of junction files, created by regtools.

The user must provide an exon metadata created from a GTF of the user's choice using `scripts/get_gene_meta_from_gtf.R` in "exon" mode. Phenotype metadata for the junctions will be created by the pipeline.

