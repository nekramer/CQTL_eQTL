# CQTL_eQTL

The pipelines used for processing eQTLs, reQTLs, and colocalizations for the CQTL project.

## Workflow

### Condition-separated eQTLs

1. Clone this entire repo into working directory:
    ```bash
    git clone https://github.com/nekramer/CQTL_eQTL.git
    ```
2. Make comma-separated `samplesheet.csv`, `donorSamplesheet.csv`, and `dnaSamplesheet.csv` giving sample
sequencing information, various experimental covariates, and various donor covariates. The scripts `makeSamplesheet.R`,
`makeDonorSamplesheet.R`, and `makeDNAsamplesheet.R` were used to get specific subsets of internal sequencing
sample sheets.

3. Update `config/config_QTLtools_eQTL.yaml` with relevant software versions, genome assemblies, and other processing
parameters.

4. Map eQTLs with:

    ```bash
    sbatch run_QTLtools_eQTL
    ```
This command will launch all required steps necessary to prepare RNA-seq and genotyping data, get covariates,
and run QTLtools to map nominal eQTLs and eQTLs with a permutation pass. It can process QTLs with varying numbers 
of PEER factors in parallel.

Note: Based on qc results, certain donors have been omitted from eQTL analyses and are noted in `omitted_donors.txt`.

### response eQTLs

1. Update `config/config_reQTL.yaml` prior to eQTL processing for automatic addition of necessary file paths.
This will result in `config/config_reQTL_final.yaml` to be used in reQTL processing. Otherwise, update 
`config/config_reQTL_final.yaml` directly with final covariate inclusions, necessary data file paths, and software 
versions.

2. Map response eQTLs from standard eQTLs with:

    ```bash
    sbatch run_responseQTL
    ```
This processing gets rsIDs and LD buddies with rsIDs for final sets of standard eQTL leads, format appropriate data
for linear modeling in R, and tests the significance of the `genotype:Condition` interaction term for lead FN-f 
variant/eGene pairs (using `lme4` and `ANOVA` in R). Further filtering of reQTL signals (i.e. only found in FN-f
treatment condition, genotype group sample sizes, difference in PBS/FN-f effect sizes) is done in 
downstream analyses.

### Colocalization with GWAS

1. Update `config/config_colocalization.yaml` with necessary data paths/file prefixes and software versions.
`eqtlN` is the number of donors in the analysis to be inputted to `coloc`.

2. Run PBS eQTL and FN-f response and non-response colocalization with osteoarthritis (OA) GWAS with:

    ```bash
    sbatch run_colocalization
    ```
This workflow obtains MAFs for variants and runs datasets through coloc. It only tests for signal colocalization
if the lead eQTL variant and lead GWAS variant are at least in moderate LD (r2 > 0.5).

## Boer et al. 2021 GWAS

Colocalizations were performed with [Boer et al. 2021](https://www.cell.com/cell/fulltext/S0092-8674(21)00941-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421009417%3Fshowall%3Dtrue) 
GWAS, which includes data for 11 defined phenotypes encompassing major sites of OA. The scripts used to process
lead variants and summary statistics, including cleaning, lift over, rsID parsing, and LD buddy R2 calculations
are included in the `/GWAS` folder. LD buddy calculations were performed using 1000 Genomes, with results included 
for LD with all 1000G samples and European-only 1000G samples. rsIDs were obtained by position and allele matching 
against parsed, .csv versions of dbSNP references.
