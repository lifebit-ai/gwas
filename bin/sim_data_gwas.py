# Script to prepare GWAS input files using output simulated data in VCF format from hapgen_simulate pipeline.

# ## Prepare simulated data for GWAS

import pandas as pd
import numpy as np

# ### Create a design file with vcf files

import glob

df = pd.DataFrame(columns=['chr','file','index'])
s3_location = "s3://lifebit-user-data-444144d0-32be-49f1-b866-a81a61b1f19b/deploit/teams/5c6d3e9bd954e800b23f8c62/users/60be2cf0303ee601a69e01a3/projects/619b69ddfb26cf01dcede60a/notebook-sessions/619b6a07fb26cf01dcede610/"
for vcf in glob.glob('*.vcf.gz'):
    print(vcf)
    chrom = vcf.split('_')[1].lstrip('chr')
    
    df = df.append({'chr':chrom,'file':s3_location+vcf,'index':s3_location+vcf+'.csi'},ignore_index=True)
    

df.to_csv('vcfs_60K_design_file.csv',index=False)

# ### Simulate phenotype
# Run e.g. bcftools query -l x.vcf.gz > 60K_samples.txt to get list of sample IDs
samples = pd.read_csv('60K_samples.txt',sep='\t',header=None)

samples.head()

# randomly assign 10000 cases, rest are controls
cases = samples.sample(3000).index
samples.loc[cases,'phenotype'] = 1
samples['phenotype'] = samples['phenotype'].replace(np.nan,0)
samples['phenotype'] = samples['phenotype'].astype(int)

#set.seed(123)
samples['age'] = pd.Series(np.random.randint(35, 96, size=len(samples)))

samples['sex'] = pd.Series(np.random.randint(1,3, size=len(samples)))

samples['IID'] = samples[0]
samples['FID'] = samples[0]


samples[['FID','IID','phenotype','age','sex']].to_csv('pheno_cov.tsv',sep='\t',index=False)