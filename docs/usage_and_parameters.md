# Usage

In order to use this pipeline, you can run the following example:

```bash
nextflow run main.nf \
  --plink_file "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
  --phenofile "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/d602bec4b31d5d75f74f1dbb408bd392db57bdb6/cohort_data_phenos.csv" \
  --metadata "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/d602bec4b31d5d75f74f1dbb408bd392db57bdb6/metadata.csv" \
  --continuous_var_aggregation "mean" \
  --continuous_var_transformation "zscore" \
  --pheno_col "Specimen type" \
  --mode 'case_vs_control_contrast' \
  --case_group "NOSE" \
  --trait_type "binary" \
  --vcfs_list "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/vcfs.csv"
```

# Parameters

**ESSENTIAL**
- **--vcfs_list** : path/url to CSV file containing chr chunk information, path to aggregated VCFs, VCFs index.
- **--plink_file** : path/url to S3 bucket that contains bed, bim, fam files for aggregated VCFs.
- **--pheno_col** : String with the name of the phenotypic column to be used as trait. Note for CB users, it must match the name of the column, i.e. 'Specimen type'.
- **--phenofile** : path/url to file that contains phenotypic information about cohort to be analysed.
- **--metadata** : path/url to file that contains metadata from phenotypic information.
- **--query** : Under development. It will allow the user to query through multiple instances & arrays of the same phenotype.
- **--programme** : Under development. It will allow the user to select between cohorts.
**Binary**
- **--trait_type** : Should be set as 'binary'.
- **--design_mode** : String containing the design matrix configuration to be used. Allows the user to select between the following scenarios:

|| Value | Description | Needs | Added value |
|--|--|--|--|--|
| Scenario 1 | case_vs_control_contrast | User wants to run on a particular case but wants to use all the rest of cases as controls. | Subset the particular case group and select all the remaining individuals as control | Find significant associations exclusive to the case group you are interested |
| Scenario 2 | case_vs_group_contrasts | User wants to run on a particular case but wants to compare to each of the other cases as controls independently | Subset case vs each group as control | Find associations that are different to an specific group |
| Scenario 3 | all_vs_all | User doesn't have a particular group in mind and wants to run an exploration on the phenotype | All vs All approach | Allows for exploration or assumptions free analysis |

- **case_group** : String containing name of the case group selected for contrasts.
**Quantitative**
- **trait_type** : Should be set to 'quantitative'.
**Optional**
- **continuous_var_transformation** : Transforms continuous variables using 'log', 'log2', 'zscores' or 'None'.
- **continuous_var_aggregation** : Defines how to aggregate different measurements. Choose between 'max', 'min', 'median' or 'mean'.
- **q_filter** : Minimum allele frequency filter for selecting sites.
- **thres_m** : Minimum threshold for missingess.
- **thres_HWE** : Minimum threshold for Hardy-Weinberg Equilibrium
- **plink_keep_pheno** : Space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis.
- **saige_step1_extra_flags** : Additional flags for SAIGE, they should be formatted as "--LOCO=FALSE".
- **outdir** : Output directory for results.
- **gwas_cat** : Path to GWAS catalog CSV file. Defaults to 's3://lifebit-featured-datasets/projects/gel/gel-gwas/gwascat.csv'.
- **output_tag** : Prefix to identify output files.
- **top_n_sites** : Minimum number of top sites to be included in output.
- **max_top_n_sites** : Maximum number of top sites to be included in output.
- **saige_filename_pattern** : File pattern specifically for SAIGE files.


