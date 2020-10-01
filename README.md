# GEL GWAS

![](bin/covid_1_manhattan.png)

## Example usage
```bash
nextflow run main.nf \
  --plinkFile "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
  --cohort_browser_phenofile "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/d602bec4b31d5d75f74f1dbb408bd392db57bdb6/cohort_data_phenos.csv" \
  --input_meta_data "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/d602bec4b31d5d75f74f1dbb408bd392db57bdb6/metadata.csv" \
  --continuous_var_transformation "mean" \
  --phenoCol "Specimen type" \
  --mode 'case_vs_control_contrast' \
  --case_group "NOSE" \
  --vcfsList "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/vcfs.csv"
```
