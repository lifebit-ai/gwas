# GWAS

![](bin/covid_1_manhattan.png)

## Example usage
Read more about parameters [here](https://github.com/lifebit-ai/gwas/blob/dev/docs/usage_and_parameters.md)

Example run using VCF files as input:

```bash
nextflow run main.nf \
  --grm_plink_input "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
  --pheno_data "s3://lifebit-featured-datasets/projects/gel/gel-gwas/cb_binary_pheno.phe" \
  --trait_type "binary" \
  --genotype_files_list "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/vcfs.csv" \
```
