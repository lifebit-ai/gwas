
# CB integration with GWAS pipeline

## 1. **Ingestion**

Two `.csv` files:
- **CB phenotypic data:** Contains columns selected by user with all respective measurements
- **CB metadata for phenotypes:** Contains information about the columns

These two files need to be passed to the pipeline in order to make it work.

## 2. **Aggregation & Transformation**

Adds a script that takes the phenotypic data and metadata associated and performs the following tasks:
- Cleans the data from files with missing genotypic data
- Reads column by column selected by the user and applies a corresponding aggregation multiple measurements and additionally transform the data if needed. Currently, is compatible with:
  - Categorical (multi-level or not) -> Selects the first measurement until querying is allowed. Adds unknown label for NA, but they are not selected for contrast.
  - Integer/Continuous -> Applies aggregation using mean, min, max, avg across measurements (instances and arrays) and transformation (log, log2, Z-score, None)
  - Dates and Time -> transforms it into YYYYMMDD integer which can be used as covariates
  - Text -> Removes free text
- Generates `.phe` file for GWAS

Assumes that sex comes from `participant phenotypic sex` -> This behaviour will be change in the future.

## 3. **Multiple design matrices**

  Given a categorical phenotype, explores all the potential combinations of interest for users in terms of contrasts to be run. Different scenarios to consider. By convention 1=case, 0=control

### 3.1 - Possible scenarios

|| Description | Needs | Added value |
|--|--|--|--|
| Scenario 1 | User wants to run on a particular case but wants to use all the rest of cases as controls. | Subset the particular case group and select all the remaining individuals as control | Find significant associations exclusive to the case group you are interested |
| Scenario 2 | User wants to run on a particular case but wants to compare to each of the other cases as controls independently | Subset case vs each group as control | Find associations that are different to an specific group |
| Scenario 3 | User doesn't have a particular group in mind and wants to run an exploration on the phenotype | All vs All approach | Allows for exploration or assumptions free analysis |


# Future work

- [ ] **(1) Add querying for data in instances**


- [ ] **(2) Implement multiple contrast reporting and handling of results**


- [ ] **(3) GWAS post-analysis**
   - [ ] Add MAGMA & friends to run pathway-analysis
   - [ ] Expand on dataviz: 
      - [ ] Add IntAssoPlot https://www.frontiersin.org/articles/10.3389/fgene.2020.00260/full
      - [ ] Expand report with modules being run.
      - [ ] Add interactivity when exploring results.
   - [ ] PRS within pipeline
   - [ ] LDSC within pipeline
   - [ ] Annovar within pipeline
   - [ ] GWAMA/METAL within pipeline
