# GEL GWAS

![](bin/covid_1_manhattan.png)

## Example usage
Read more about parameters [here](https://github.com/lifebit-ai/gel-gwas/blob/marcos-integrates-cb-output/docs/usage_and_parameters.md)
```bash
nextflow run main.nf \
  --grm_plink_input "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
  --phenofile "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/d602bec4b31d5d75f74f1dbb408bd392db57bdb6/cohort_data_phenos.csv" \
  --metadata "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/d602bec4b31d5d75f74f1dbb408bd392db57bdb6/metadata.csv" \
  --continuous_var_aggregation "mean" \
  --continuous_var_transformation "zscore" \
  --pheno_col "Specimen type" \
  --design_mode 'case_vs_control_contrast' \
  --case_group "NOSE" \
  --trait_type "binary" \
  --vcfs_list "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/vcfs.csv"
```
