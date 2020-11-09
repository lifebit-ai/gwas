# Simulate pheno data & testing modes

## 1. **Simulate mode**

This option enables the user to simulate their own cohort of participants and test it with the testdata provided for development. This feature is meant to allow CI tests to run on independent simulated data and ensure that the pipeline's further development doesn't have any further dependencies. The following scenarios are considered:

|| Description | Needs | Flags | 
|--|--|--|--|
| Scenario 1 | Simulate from config file, no real pheno data present. | Config file similar to the example below |  **--sim_config_file** <path to config.yml> |
| Scenario 2 | User has real pheno data that can be used to generate new samples | Config file similar to example below, only having the number of samples to be generated. Real phenotypic data to upsample from | **--sim_config_file** <path to config. yml> **--sim_pheno_data <path to pheno data.csv>**|

Example of config file:
```
params:    
    n_samples: 100
    col_params:
        cat_col:
            n_instances: 1
            n_arrays: 4
            type: "Categorical"
            values: ['A', 'B', 'NA']
            prob: [0.4, 0.4, 0.2]
            n_missing_col: 2
        int_col:
            n_instances: 2
            n_arrays: 4
            type: "Integer"
            mean: 1
            sd: 4
            distribution: 'normal'
            min: 0
            max: 2
            frac_missing: 0.33
            n_missing_col: 3
            positive_only: True
        cont_col:
            n_instances: 2
            n_arrays: 4
            type: "Continuous"
            mean: 1
            sd: 4
            distribution: 'normal'
            frac_missing: 0.33
            n_missing_col: 3
            positive_only: True
        date_col:
            n_instances: 2
            n_arrays: 4
            type: "Date"
            starting_date: "1954-12-01"
            end_date: "2020-08-01"
            frac_missing: 0.33
            n_missing_col: 3
    

```

It can be possible to provide a JSON instead, but refactor of the script might be necessary to make it work with this option.

`simulate_phenodata.R` takes the following inputs:
- **--pheno_data** : Path to pheno data file. It can be use to generate more samples by upsampling with replacement. (Optional)
- **--config_file** : Path to yml file that will be use to generate as many samples as required or generate the whole cohort data out of the blue. (Required)
- **--outprefix** :

It produces two files:
- **<outprefix>_pheno_data.csv** : Emulates the pheno data that the CB run analysis or export button produces.
- **<outprefix>_pheno_metadata.csv** : Emulates a simple version metadata data that the CB run analysis or export button produces and it's used to process CB outputs. It generates only the 3 columns required (id of the column, name of the column, datatype of the column).


- **query.json** : Is **not generated** yet, can be included in future updates to test query feature as well.

## 2. **Testing mode**

This features allows the user to use any real or simulated CB pheno data with testing genotypic data. It will take the `Platekey_in_aggregate_VCF` column and match it to the genotypic data `Platekeys` 

`test_data_munging.R` takes the following inputs:
- **--input_file** : Pheno data file that will be wrangled.
- **--ids_columns** : Name of column containing platekeys. This column with be switched to a consecution of `1_1`, `2_2`, `3_3` and so on to match testing genotypic data.
- **--outprefix** : Prefix of output files

It will output the following files:

- GWAS ready file: Removes any diagnosis code columns that are not meant to be run on a GWAS.
- PheWAS ready file: It doesn't remove diagnosis code
- All in: For is the same as PheWAS ready file, just in case future changes happen on the GWAS file. 

Subsequent steps in gel-gwas will take the GWAS ready file and ignore the others. 

This steps should then send the resulting file and its metadata to the processes to run a GWAS. 

Note: It doesn't assume that sex comes from any column in particular. If present it will be kept as a covariate for GWAS and transformed as any other categorical variable.


## Considerations

1. Running simulation and testing is absolutely compatible. Platekeys are generated during simulation, but to ensure they match testing data they can be use in conjuntion with the testing step to couple it to synthetic data.
2. Simulations gives freedom of testing and development to this pipe allowing users to implement the same data model the CB uses in order to add new features and improvements to this pipeline, without depending on the CB itself. 