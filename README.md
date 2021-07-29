# GEL GWAS

![](bin/covid_1_manhattan.png)

## Example usage
Read more about parameters [here](https://github.com/lifebit-ai/gel-gwas/docs/usage_and_parameters.md)

```bash
nextflow run main.nf \
  --grm_plink_input "s3://lifebit-featured-datasets/pipelines/simulate/ukb-simulated-results/simulated_hapgen-100000ind-updated.merged.*{bed,bim,fam}" \
  --pheno_data "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/traits_design_matrix_control_all_case_2.phe" \
  --trait_type "binary" \
  --vcfs_list "s3://lifebit-featured-datasets/pipelines/simulate/ukb-simulated-results/vcfs_ukbio.csv"
```
