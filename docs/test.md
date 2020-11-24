# Testing modes

## 1. **Testing mode**

This features allows the user to use any real or simulated CB pheno data with testing genotypic data. It will take the `Platekey_in_aggregate_VCF` column and match it to the genotypic data `Platekeys` 

`test_data_munging.R` takes the following inputs:
- **--input_file** : Pheno data file that will be wrangled.
- **--ids_columns** : Name of column containing platekeys. This column with be switched to a consecution of `1_1`, `2_2`, `3_3` and so on to match testing genotypic data.
- **--outprefix** : Prefix of output files

It will output the following files:

- GWAS ready file: Removes any diagnosis code, HPO or ICD columns that are not meant to be run on a GWAS.
- PheWAS ready file: It doesn't remove diagnosis code
- All in: For is the same as PheWAS ready file, just in case future changes happen on the GWAS file. 

Subsequent steps in gel-gwas will take the GWAS ready file and ignore the others. 

This steps should then send the resulting file and its metadata to the processes to run a GWAS. 

Note: It doesn't assume that sex comes from any column in particular. If present it will be kept as a covariate for GWAS and transformed as any other categorical variable.
