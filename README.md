# mmQTL-pipeline

Jack Humphrey & Erica Brophy

Towfique Raj Lab

Mount Sinai, New York

2021-2023


Snakemake pipeline for running multi-cohort QTL meta-analysis using the linear mixed model [mmQTL](https://pubmed.ncbi.nlm.nih.gov/35058635/) to account for diverse ancestries and repeated donors.

Pipeline allows for multiple cohorts to be uniformly processed for joint QTL mapping and meta-analysis of multiple phenotypes, including expression, transcript usage, splice junction usage (Leafcutter) and splicing event usage (SUPPA).

## Input data for the meta-analysis:

User must assemble a  **dataset key** (see example) created by the user in the format of a dataframe. The dataset key contains 5 rows: **1. Dataset, 2. sampleKey, 3. VCF, 4. Phenotype file, 5. Covariate file**; where the column provides the path to each of these files. 

In this example, 3 datasets from 2 human brain cohorts are combined:

| dataset | sample_key |  genotypes |  phenotypes |  PEER | counts |
| ----------- | ----------- | ------| ------| ---- | ---- |
| cohort1_frontal_cortex  | cohort1_temporal_sample_key.tsv | cohort1.vcf | cohort1_frontal_gene_tpm.tsv | 10 | cohort1_frontal_gene_counts.tsv |
| cohort1_temporal_cortex  | cohort1_frontal_sample_key.tsv | cohort1.vcf | cohort1_temporal_gene_tpm.tsv | 15 | cohort1_temporal_gene_counts.tsv |
| cohort2_frontal_cortex  | cohort2_frontal_sample_key.tsv | cohort2.vcf | cohort2_frontal_gene_tpm.tsv | 15 | cohort2_frontal_gene_tpm.tsv |

### Input data per cohort:

The inputs required by the user are a *sampleKey* matching RNA samples to DNA samples, a *VCF file* containing the genotype information, a *phenotype matrix* and a *gene annotation file* (GENCODE or similar). 

The *sample key* must have two columns named **participant_id** listing the names of the genotype IDs and **sample_id** listing the name of the phenotype sample IDs.

The *genotypes* must be in VCF format with all chromsomes in a single file.

The *phenotype matrix* must be a tab-separated table with one feature per row, the first column should be the unique identifier of the feature, followed by each sample ID in the dataset. The matrix should contain normalised values for each sample. See below for types of phenotypes that mmQTL can work with.

The *PEER factor* should be a number of PEER factors to be regressed out of the phenotype matrix in each dataset. Ideally this has already been determined using the [QTL-mapping-pipeline](https://github.com/RajLabMSSM/QTL-mapping-pipeline)

*counts* For mapping gene expression QTLs, TPM values are used for removing lowly expressed genes, and read counts are used as the phenotype. For all other QTL types a counts column is not necessary.

## Types of phenotypes used for mmQTL:

* Gene expression

This must be TPM normalised. A phenotype metadata file should be created. This can be created from a GTF file using `scripts/get_gene_meta_from_gtf.R` in "gene" mode.

* Transcripts

The user must prepare a transcript expression matrix with TPM normalised transcript expression, from Kallisto, RSEM or similar.

If the user wishes to treat each transcript as it's own feature - (transcript QTLs), in the config.yaml set `group: False`

If the user wishes to divide each transcript by the total gene TPM, (transcript usage QTLs), set `group: True`

Transcript metadata must be created from the same GTF as used to create the transcript expression matrix using `scripts/get_gene_meta_from_gtf.R` in "transcript" mode.

* SUPPA local splicing events

```
python3.4 suppa.py generateEvents -i <input-file.gtf> -o <output-file> -f ioe -e <list-of-events>
```

* Junction usage QTLs (Leafcutter)

If the user wants to perform junction usage QTLs (Leafcutter splicing QTLs) then provide a list of junction files, created by regtools.

The user must provide an exon metadata created from a GTF of the user's choice using `scripts/get_gene_meta_from_gtf.R` in "exon" mode. Phenotype metadata for the junctions will be created by the pipeline.

Then prepare the full leafcutter junction matrix for the entire cohort:

```
Rscript scripts/prepare_leafcutter_matrix.R --prefix <prefix> --junctions <junction list> --exons <exon metadata> 
```

Where prefix is a directory followed by the stem of a filename, ie "junctions/test"


## NEW - Wrapper script to create all metadata!

`prepare_qtl_metadata.lsf` creates all QTL metadata files based on your GTF and list of junctions


## Config options

* minDatasets 

This sets the minimum number of datasets that a feature must be present in to be used for QTL meta-analysis (default = 2)

* phenoThreshold

Minimum transcripts per million (TPM) threshold for keeping a gene or transcript for QTL mapping. Default = 1

* phenoFraction

Minimum fraction of the samples in each cohort that must meet the `phenoThreshold`. Default = 0.5

