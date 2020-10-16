
# Colocalisation analysis integration with pheWAS pipeline

PheWAS output must be combined with GWAS in order to compare results and establish if a phenotype of interest and a trait of interest share common variants.

## 1. Information about the method

Vignette: https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html#org192e16b:

### 1.1 - Inputs
- **--post_analysis** : String containing `"coloc"` to specify that colocalization analysis should be executed.
- **--gwas_input** : CSV file containing summary statistics from a GWAS study. This file must contain the following columns:

    1. `SNPID` : Contains IDs of SNPs (i.e "rs000001").
    2. `p.value` : Contains significance levels for each one of the SNPs tested.
    3. `N` : How many individuals with the SNP are present.
    4. `N.Cases` : In `gwas_trait_type = "binary"` how many individuals in the case group.
    5. `N.Controls`: In `gwas_trait_type = "binary"` how many individuals in the control group.

- **gwas_trait_type** : String containing the type of trait analysed in GWAS: 'binary' or 'quantitative'


### 1.2 - Outputs
`colocalization` folder containing:
- `coloc_results.csv` CSV containing Posterior probabilities resulting from colocalization analysis between each phenotype and the respective GWAS study.
- `coloc_heatmap.png` Heatmap of the Posterior probabilities for all phenotypes with description of the group.


## 2. Q&A / Considerations
1. It has to work with individual VCFs
2. Should users only use pheWAS with ICD10 columns?
3. Can I run with multi-level phenotypes? -> Test and if not adapt design_matrix.R
4. Does the report needs updating?
5. Not all samples might have any ICD10 code at all. 

## 3. Tasks
- **(1) Make it compatible with SAIGE GWAS output**
    - [ ] Check which stats are needed.
    - [ ] Check that format of table is compatible
    - [ ] Create script to transform SAIGE output into ldsc.py input
    - [ ] Check that `ldsc.py` runs with the data provided from SAIGE
- **(2) Build script to run quantitative and binary coloc analysis**
    - [x] Create function to run across all phenotypes.
    - [x] Produce final table & plot
- **(3) Modify report to include colocalization analysis when selected**
    - [x] Add new table
    - [x] Add new plot
- **(4) Add necessary nextflow processes**
    - [ ] Create switch for turning the analysis on and off 
    - [ ] Bifurcate report into two process depending on scenarios
- **(5) Update nextflow.config and Dockerfile**
   - [ ] Create new image for ldsc
- **(6) Test on CloudOS**

## 4. Tests

### Iter2

https://cloudos.lifebit.ai/app/jobs/5f878fc28f81710113099d68 :tada:

## 5. Future Work
- **(1) Add GWAS catalogue**
- **(2) Improve plot, and explore other plotting options**
- **(3) Add other types of data for colocalization (eQTLs, pQTLs)**
- **(4) Add scatterplots for comparison between genetic variants of interest between datasets**




