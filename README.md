# GEL GWAS

![](bin/covid_1_manhattan.png)

## Example usage
Read more about parameters [here](https://github.com/lifebit-ai/gel-gwas/docs/usage_and_parameters.md)

```bash
nextflow run main.nf \
  --grm_plink_input "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
  --pheno_data "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/d602bec4b31d5d75f74f1dbb408bd392db57bdb6/cohort_data_phenos.csv" \
  --trait_type "binary" \
  --vcfs_list "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/vcfs.csv" \
```
