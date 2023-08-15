# CQTL_eQTL

The pipeline used for processing eQTLs and reQTLs for the CQTL project.

## Workflow
1. Clone this entire repo into working directory:
    ```bash
    git clone https://github.com/nekramer/CQTL_eQTL.git
    ```
2. Make comma-separated `samplesheet.csv`, `donorSamplesheet.csv`, and `dnaSamplesheet.csv` giving sample
sequencing information, various experimental covariates, and various donor covariates. The scripts `makeSamplesheet.R`,
`makeDonorSamplesheet.R`, and `makeDNAsamplesheet.R` were used to get specific subsets of internal sequencing
sample sheets.

3. Map eQTLs with:

    ```bash
    sbatch run_QTLtools_eQTL
    ```
This command will launch all required steps necessary to prepare RNA-seq and genotyping data, get covariates,
and run QTLtools to map nominal eQTLs and eQTLs with a permutation pass.                                
