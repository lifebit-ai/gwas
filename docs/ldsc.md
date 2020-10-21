# Genetic Correlation and Heritability within GWAS pipeline

Adds the option of running heritability in your GWAS summary statistics as well as computing the genetic correlation between your trait of interest and a second trait with gwas summary statistics  

## 1. Information about the method

Vignette: https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation

### 1.1 - Inputs
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




### 1.2 - Outputs
`genetic_correlation` folder containing:
- `<filename>.log` Txt file with the results of the analysis between trait being analysed in the GWAS pipeline and the respective external GWAS study.

`heritability` folder containing:
- `<filename>.log` Txt file with the results of the heritability analysis.

## 2. Q&A / Considerations
1. Results in test might look weird because I used the same testdata with different traits. 

## 3. Tasks
- **(1) Make it compatible with SAIGE GWAS output**
    - [x] Check which stats are needed.
    - [x] Check that format of columns is compatible
    - [x] Write script to reformat file
- **(2) Build processes for heritability**
    - [x] Create process to munge statistics to the required format for LDSC.
    - [x] Test the process and see how well it connects
- **(3) Modify report to include LDSC analysis when selected**
    - [x] Dump content of log file into the report as no tables are produced
- **(4) Add necessary nextflow processes**
    - [x] Create switch for turning the analysis on and off 
    - [x] Bifurcate report into two process depending on scenarios
- **(6) Update nextflow.config and Dockerfile**
- **(7) Test on CloudOS**

## 4. Tests

### Genetic correlation & Heritability - iter 1

Genetic correlation
Heritability

## 5. Future Work
- **(1) Add GWAS catalogue**
- **(2) Find ways to better structure the results from the log files**
- **(3) Explore other visualizations for this analysis**

