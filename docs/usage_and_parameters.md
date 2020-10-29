# Usage

In order to use this pipeline, you can run the following example:
**Binary**
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

**Quantitative**
```bash
nextflow run main.nf \
  --grm_plink_input "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
  --phenofile "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/d602bec4b31d5d75f74f1dbb408bd392db57bdb6/cohort_data_phenos.csv" \
  --metadata "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/d602bec4b31d5d75f74f1dbb408bd392db57bdb6/metadata.csv" \
  --continuous_var_aggregation "mean" \
  --continuous_var_transformation "log" \
  --pheno_col "Height (HCM)" \
  --trait_type "quantitative" \
  --vcfs_list "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/vcfs.csv"
```

**LDSC - Genetic correlation with an external GWAS summary stats**
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
  --vcfs_list "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/vcfs.csv" \
  --post_analysis "genetic_correlation_h2" \
  --gwas_summary "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/e4f8fc5bd62c70ef38c6cedfdfaa6d087f586054/gwas_summary_qt.csv"
```

**LDSC - Genetic correlation with GWAS catalogue summary stats**
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
  --vcfs_list "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/vcfs.csv" \
  --post_analysis "genetic_correlation_h2" \
  --gwas_cat_study_id "GCST004420-ci" \
  --gwas_cat_study_size = 3631 \
  --gwas_catalogue_ftp = "https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/projects/gel/prs/ftp_locations_harmonized.csv" \
  --hapmap3_snplist = "https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/projects/gel/gel-gwas/assets/w_hm3.snplist" \
  --ld_scores_tar_bz2 = "https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/eur_w_ld_chr.tar.bz2" \
```

**LDSC - Heritability**

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
  --vcfs_list "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/vcfs.csv" \
  --post_analysis "heritability" 
```

# Parameters

## **ESSENTIAL**

- **--vcfs_list** : path/url to CSV file containing chr chunk information, path to aggregated VCFs, VCFs index.
- **--grm_plink_input** : path/url to folder that contains bed, bim, fam files for aggregated VCFs.
- **--pheno_col** : String with the name of the phenotypic column to be used as trait. Note for CB users, it must match the name of the column, i.e. 'Specimen type'.
- **--phenofile** : path/url to file that contains phenotypic information about cohort to be analysed.
- **--metadata** : path/url to file that contains metadata from phenotypic information.
- **--query** : Under development. It will allow the user to query through multiple instances & arrays of the same phenotype.
- **--programme** : Under development. It will allow the user to select between cohorts.

## **Binary**

- **--trait_type** : Should be set as 'binary'.
- **--design_mode** : String containing the design matrix configuration to be used. Allows the user to select between the following scenarios:

|| Value | Description | Needs | Added value |
|--|--|--|--|--|
| Scenario 1 | case_vs_control_contrast | User wants to run on a particular case but wants to use all the rest of cases as controls. | Subset the particular case group and select all the remaining individuals as control | Find significant associations exclusive to the case group you are interested |
| Scenario 2 | case_vs_group_contrasts | User wants to run on a particular case but wants to compare to each of the other cases as controls independently | Subset case vs each group as control | Find associations that are different to an specific group |
| Scenario 3 | all_vs_all | User doesn't have a particular group in mind and wants to run an exploration on the phenotype | All vs All approach | Allows for exploration or assumptions free analysis |

- **--case_group** : String containing name of the case group selected for contrasts.

## **Quantitative**

- **--trait_type** : Should be set to 'quantitative'.


## **Optional**

- **--continuous_var_transformation** : Transforms continuous variables using 'log', 'log10', 'log2', 'zscores' or 'None'.
- **--continuous_var_aggregation** : Defines how to aggregate different measurements. Choose between 'max', 'min', 'median' or 'mean'.
- **--q_filter** : Minimum allele frequency filter for selecting sites.
- **--thres_m** : Minimum threshold for missingess.
- **--thres_HWE** : Minimum threshold for Hardy-Weinberg Equilibrium
- **--plink_keep_pheno** : Space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis.
- **--saige_step1_extra_flags** : Additional flags for SAIGE, they should be formatted as "--LOCO=FALSE".
- **--outdir** : Output directory for results.
- **--gwas_cat** : Path to GWAS catalog CSV file. Defaults to 's3://lifebit-featured-datasets/projects/gel/gel-gwas/gwascat.csv'.
- **--output_tag** : Prefix to identify output files.
- **--top_n_sites** : Minimum number of top sites to be included in output.
- **--max_top_n_sites** : Maximum number of top sites to be included in output.
- **--saige_filename_pattern** : File pattern specifically for SAIGE files.

## LDSC
- **--post_analysis** : String with `genetic_correlation_h2` or `heritability` for running genetic correlation analysis or heritability after GWAS.
- **--gwas_summary** : Path/URL to external gwas summary statistics to run genetic correlation analysis between cohort of interest and external GWAS summary statistics. The following column names and format (can also be comma-separated instead of whitespace-separated) are required to ensure it works:

```
snpid hg18chr bp a1 a2 or se pval info ngt CEUaf
rs3131972	1	742584	A	G	1.092	0.0817	0.2819	0.694	0	0.16055
rs3131969	1	744045	A	G	1.087	0.0781	0.2855	0.939	0	0.133028
rs3131967	1	744197	T	C	1.093	0.0835	0.2859	0.869	0	.
rs1048488	1	750775	T	C	0.9158	0.0817	0.2817	0.694	0	0.836449
rs12562034	1	758311	A	G	0.9391	0.0807	0.4362	0.977	0	0.0925926
rs4040617	1	769185	A	G	0.9205	0.0777	0.2864	0.98	0	0.87156
rs28576697	1	860508	T	C	1.079	0.2305	0.7423	0.123	0	0.74537
rs1110052	1	863421	T	G	1.088	0.2209	0.702	0.137	0	0.752294
rs7523549	1	869180	T	C	1.823	0.8756	0.4929	0.13	0	0.0137615
``` 
- **--gwas_cat_study_id** : String containing study ID for GWAS study in GWAS catalogue.
- **--gwas_cat_study_size** : Integer contaning sample size for GWAS study.
- **--gwas_catalogue_ftp** : String https link or path to where gwas catalogue download is provided. Defaults to `https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/projects/gel/prs/ftp_locations_harmonized.csv`.
- **--external_gwas_tag** : String containing tag to be used to identify external GWAS resource for LDSC.
- **--hapmap3_snplist** : String containing link to HapMap3 SNP list. i.e. "https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/projects/gel/gel-gwas/assets/w_hm3.snplist"
- **--ld_scores_tar_bz2** : String containing link to LD Scores `*.tar.bz2` file. i.e. "https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/eur_w_ld_chr.tar.bz2"

## Testing mode
- **--testing** : String that allows the user to use synthetic genotypic data with real pheno data for developing and testing purposes. i.e. "True"