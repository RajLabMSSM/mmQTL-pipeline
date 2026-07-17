# mmQTL-pipeline

Jack Humphrey, Erica Brophy, Kailash BP, Winston Cuddleston, Beomjin Jang, Tatsuhiko Naito

Raj Lab, Icahn School of Medicine at Mount Sinai, New York, NY

2021-2026

The **mmQTL pipeline** is a Snakemake workflow for performing multi-cohort QTL meta-analysis using the linear mixed model [**mmQTL**](https://pubmed.ncbi.nlm.nih.gov/35058635/), which accounts for both diverse ancestries and repeated donors.

The pipeline uniformly processes multiple cohorts before jointly mapping QTLs across studies. It supports both *cis*- and *trans*-QTL mapping and multiple molecular phenotypes, including:

- Gene expression
- Transcript expression
- Transcript usage
- Junction usage (Leafcutter)
- Splicing event usage (SUPPA)
- RNA editing

The workflow performs phenotype preprocessing, normalization, covariate regression, genotype harmonization, mixed-model association testing, and meta-analysis within a reproducible Snakemake workflow.

This pipeline has been used in multiple studies, including:

> Humphrey, J., Brophy, E., Kosoy, R. et al. Long-read RNA sequencing atlas of human microglia isoforms elucidates disease-associated genetic regulation of splicing. Nat Genet 57, 604–615 (2025). https://doi.org/10.1038/s41588-025-02099-0

> Jang, B., BP, K., Tokolyi, A. et al. A meta-analysis of single-nucleus expression quantitative trait loci linking genetic risk to brain disorders. Nat Genet 58, 737–747 (2026). https://doi.org/10.1038/s41588-026-02541-x

---

# Installation

Clone the repository

```bash
git clone https://github.com/RajLabMSSM/mmQTL-pipeline.git
cd mmQTL-pipeline
```

## Software

Load required software modules (or configure your environment)

- snakemake ≥6.9
- plink2/2.3
- gcta/1.94.1
- bcftools/1.9
- tabix/0.2.6
- R/4.2.0
- Python/3.6.8

---

# Pipeline overview

For each dataset the pipeline performs:

1. Sample matching
2. Genotype preprocessing
3. Cross-cohort genotype harmonization
4. Phenotype normalization
5. PEER factor & Covariate regression
7. Phenotype harmonization across datasets
8. GRM construction
9. mmQTL association testing per cohort
10. Meta-analysis across cohorts
11. Result collation

---

# Input data for the meta-analysis

User must assemble a  **dataset key** (see example) created by the user in the format of a dataframe. The dataset key contains one row per dataset: **1. Dataset, 2. sampleKey, 3. VCF, 4. Phenotype file, 5. Covariate file**; where the column provides the path to each of these files. 

Each row corresponds to one cohort.

Example:

| dataset | sample_key | genotypes | phenotypes | PEER | counts |
|----------|------------|-----------|------------|------|--------|
| ROSMAP_DLPFC | sample_key.tsv | genotype.vcf.gz | gene_tpm.tsv.gz | 15 | gene_counts.tsv.gz |
| Mayo_TCX | sample_key.tsv | genotype.vcf.gz | gene_tpm.tsv.gz | 20 | gene_counts.tsv.gz |

Optional columns

- covariates

See `example/example_data_key.tsv`

---

## Input data per cohort:

The inputs required by the user are a *sampleKey* matching RNA samples to DNA samples, a *VCF file* containing the genotype information, a *phenotype matrix* and a *gene annotation file* (GENCODE or similar). 

The *sample key* must have two columns named **participant_id** listing the names of the genotype IDs and **sample_id** listing the name of the phenotype sample IDs.

The *genotypes* must be in VCF format with all chromsomes in a single file.

The pipeline automatically

- converts to PLINK
- removes multiallelic variants
- removes blacklist regions
- removes variants with excessive missingness
- removes singleton variants across cohorts
- filters rare variants using user-defined minor allele frequency thresholds
- optional variant extraction
- creates GRMs

The *phenotype matrix* must be a tab-separated table with one feature per row, the first column should be the unique identifier of the feature, followed by each sample ID in the dataset. The matrix should contain normalised values for each sample. See below for types of phenotypes that mmQTL can work with.

Gene metadata can be generated using

```
scripts/get_gene_meta_from_gtf.R
```

The *PEER factor* should be a number of PEER factors to be regressed out of the phenotype matrix in each dataset. Ideally this has already been determined using the [QTL-mapping-pipeline](https://github.com/RajLabMSSM/QTL-mapping-pipeline)

*counts* For mapping gene expression QTLs, TPM values are used for removing lowly expressed genes, and read counts are used as the phenotype. For all other QTL types a counts column is not necessary.

The user may optionally provide a table of known covariates to be regressed from the phenotype matrix in addition to the inferred PEER factors. The covariate file should be provided in **long format**, where each row corresponds to a covariate and each column corresponds to a sample. The first column should contain the covariate name, followed by one column per sample matching the sample IDs in the phenotype matrix.

Example:

| ID | Sample1 | Sample2 | Sample3 |
|----|--------:|--------:|--------:|
| age | 72 | 65 | 81 |
| biological_sex | 1 | 0 | 1 |
| astrocyte | 0.12 | 0.08 | 0.31 |
| microglia | 0.15 | 0.08 | 0.09 |

Both categorical and continuous covariates are supported, including demographic variables (e.g. age and biological sex), technical covariates (e.g. RIN, batch, sequencing metrics), and estimated cell-type proportions.

To enable known covariates, set

```yaml
known_covars: True
```

and include a `covariates` column in the dataset key pointing to the corresponding covariate file for each dataset. The pipeline will merge the supplied covariates with the inferred PEER factors prior to covariate regression.

---

The pipeline supports both *cis* and *trans* association testing.

## *cis* QTL

```yaml
QTL_type: cis
```

Features are tested within the specified *cis*-window. The pipeline harmonizes phenotype matrices across cohorts before association testing.

---

## *trans* QTL

```yaml
QTL_type: trans
```

Additional preprocessing is performed.

The pipeline

- constructs harmonized *trans* phenotype metadata
- prepares *trans*-only phenotype matrices
- optional extraction of a predefined set of variants (primarily for *trans*-QTL analyses)

```yaml
variantsToExtract = "trans_variants_to_extract.txt"
```

A text file with one variant per line pre-selected for **trans**-QTL mapping. This substantially reduces computational burden for targeted *trans*-QTL analyses.

---

# Conditional QTL mapping

The pipeline also supports mapping multiple independent QTL signals per feature.

```yaml
QTL_number: 1
```

returns only the primary association.

Setting

```yaml
QTL_number: 2
```

maps the primary and secondary independent associations.

Similarly,

```yaml
QTL_number: 5
```

returns up to five independent QTL signals. If set to greater than 5, it will only return 5. 

The pipeline automatically

- runs mmQTL for each requested conditional peak
- collates peak-specific results
- produces separate summary files

Example outputs

```
ROSMAP_peak_1_top_assoc.tsv.gz

ROSMAP_peak_2_top_assoc.tsv.gz

ROSMAP_peak_3_top_assoc.tsv.gz
```

and corresponding genome-wide association files

```
*_peak_1_full_assoc.tsv.gz

*_peak_2_full_assoc.tsv.gz
```

---

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

---

## Meta-analysis approaches

Association statistics from individual cohorts are combined using both fixed- and random-effects meta-analysis. The user can choose between a standard Wald Random association `scrips/MMQTL26_DAK` or the Buhm/Han random effect model implementation `scrips/MMQTL26a`

---

## Wrapper script to create all metadata!

`prepare_qtl_metadata.lsf` creates all QTL metadata files based on your GTF and list of junctions

---

# Configuration

The pipeline is controlled through `config.yaml`. Values not supplied by the user fall back to the defaults shown below.

## Required settings

| Parameter | Default | Description |
|------------|---------|-------------|
| `dataKey` | — | Dataset key file with one row per cohort and paths to cohort-specific inputs |
| `dataCode` | `Default` | Output prefix used to name result files |
| `phenoMeta` | `""` | Phenotype metadata file for the main analysis |
| `QTL_type` | `cis` | QTL mode: *cis* or *trans* |

## Optional settings

| Parameter | Default | Description |
|------------|---------|-------------|
| `phenoMetaTrans` | same as `phenoMeta` | Phenotype metadata used for *trans* analysis |
| `group` | `False` | Transcript usage mode; if `True`, transcripts are divided by gene TPM |
| `leafcutter` | `False` | Enable LeafCutter mode |
| `SUPPA` | `False` | Enable SUPPA mode |
| `SUPPA_events` | `""` | File listing SUPPA events |
| `edqtl` | `False` | Enable RNA-editing QTL mode |
| `variantsToExtract` | `"/dev/null"` | Optional list of variants to retain before association testing |
| `minDatasets` | `2` | Minimum number of datasets a feature must be present in for meta-analysis |
| `phenoThreshold` | `1` | TPM threshold used for phenotype filtering |
| `phenoFraction` | `0.5` | Fraction of samples that must exceed `phenoThreshold` |
| `known_covars` | `False` | Merge user-supplied covariates with PEER factors |
| `GTF` | `""` | GTF annotation file used for metadata generation |
| `outFolder` | `results/` | Output directory |
| `mmQTL_bin` | `/sc/arion/projects/bigbrain/MMQTL26a` | Path to the mmQTL binary |
| `crossmap_file` | `"/dev/null"` | Crossmap file used during collation |
| `snp_to_feature_file` | `"/dev/null"` | SNP-to-feature mapping file used during collation |
| `threads` | `4` | Number of threads for mmQTL association testing |
| `threads_collate` | `4` | Number of threads for result collation |
| `R_version` | `R/4.2.0` | R module version |
| `PLINK_version` | `plink2/2.3` | PLINK2 module version |
| `TABIX_version` | `tabix/0.2.6` | Tabix module version |
| `BCFTOOLS_version` | `bcftools/1.9` | BCFtools module version |
| `GCTA_version` | `gcta/1.94.1` | GCTA module version |

## Notes

- `QTL_number` is capped at `5` by the workflow.
- If `phenoMetaTrans` is not set, the pipeline uses `phenoMeta`.
- If `known_covars: True`, the `dataKey` should include a `covariates` column with the path to each cohort’s known covariate file.
- If `variantsToExtract` is left as `"/dev/null"`, all available variants are retained.
- A minimal example config is provided in `config.yaml`.

---

# Output structure

```
results/
    genotypes/
        filtered PLINK files
        GRMs
    <dataset>/
        normalized phenotypes
        PEER covariates
        regressed phenotypes
    mmQTL/
        phenotype_metadata.tsv
        phenotype_metadata_trans.tsv (if *trans* mode)
        pheno_list.txt
        chr*_geno_list.txt
        *_peak_1_top_assoc.tsv.gz
        *_peak_1_full_assoc.tsv.gz
        *_peak_2_full_assoc.tsv.gz
        *_peak_3_full_assoc.tsv.gz
```

---

# Parallelization

> Association testing is automatically split into chromosome chunks.
> Chunk sizes are proportional to chromosome length, improving load balancing across compute nodes.
> QTL mapping is parallelized per rule based on number of threads provided.
> Sorting, and merging is also parallelized based on thread and results generated after completion!

---

# A note from the developers

Like any research software, this pipeline is a work in progress. While it has been used successfully in multiple large-scale studies, we cannot guarantee that every edge case has been anticipated.

If you run into unexpected behavior, please don't hesitate to open a GitHub Issue. Contributions, suggestions, feature requests, and pull requests are always welcome. We hope this pipeline saves you time! :)
