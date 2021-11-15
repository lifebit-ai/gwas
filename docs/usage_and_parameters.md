# Usage

In order to use this pipeline, you can run the following example:

**Binary**

```bash
nextflow run main.nf \
  --grm_plink_input "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
  --pheno_data "<pheno_file>" \
  --trait_type "binary" \
  --vcfs_list "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/vcfs.csv"
```

**Quantitative**

```bash
nextflow run main.nf \
  --grm_plink_input "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
  --pheno_data "<pheno_file>" \
  --trait_type "quantitative" \
  --vcfs_list "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/vcfs.csv"
```

# Parameters

## **ESSENTIAL**

- **--vcfs_list** : path/url to CSV file containing chr chunk information, path to aggregated VCFs, VCFs index.
- **--grm_plink_input** : path/url to folder that contains bed, bim, fam files for aggregated VCFs.
- **--pheno_data** : path/url to file that contains phenotypic information about cohort to be analysed in plink phenotypic format (.phe).

## **Binary**

- **--trait_type** : Should be set as 'binary'.

## **Quantitative**

- **--trait_type** : Should be set to 'quantitative'.


## **Optional**

- **--q_filter** : Minimum allele frequency filter for selecting sites.
- **--thres_m** : Minimum threshold for missingess.
- **--hwe_threshold** : Significance threshold for Hardy-Weinberg Equilibrium. Default = 1e-5
- **--hwe_test** : Type of test done for Hardy-Weinberge Equilibrium. Default = 'midp' which stands for the mid-p adjustmenttends (plink's recommended option). Favours filter of variants with more missing data.

- **--ld_window_size** Default = 50. Window size for LD-pruning.
- **--ld_step_size** Default = 10. Step size for LD-pruning.
- **--ld_r2_threshold** Default = 0.1. R2 correlation threshold for LD-pruning.

- **--saige_step1_extra_flags** : Additional flags for SAIGE, they should be formatted as "--LOCO=FALSE".
- **--outdir** : Output directory for results.
- **--gwas_cat** : Path to GWAS catalog CSV file. Defaults to 's3://lifebit-featured-datasets/projects/gel/gel-gwas/gwascat.csv'.
- **--output_tag** : Prefix to identify output files.
- **--top_n_sites** : Minimum number of top sites to be included in output.
- **--max_top_n_sites** : Maximum number of top sites to be included in output.
- **--saige_filename_pattern** : File pattern specifically for SAIGE files.