#!/usr/bin/env nextflow
/*
========================================================================================
                         lifebit-ai/gel-gwas
========================================================================================
 lifebit-ai/gel-gwas GWAS pipeline built for Genomics England using SAIGE
 #### Homepage / Documentation
 https://github.com/lifebit-ai/gel-gwas
----------------------------------------------------------------------------------------
*/

// Channel setup
Channel
  .fromPath(params.phenoFile)
  .ifEmpty { exit 1, "Pheno file not found: ${params.phenoFile}" }
  .set { phenoCh }
Channel
  .fromFilePairs("${params.plinkFile}.{bed,bim,fam}",size:3, flat : true)
  .ifEmpty { exit 1, "PLINK files not found: ${params.plinkFile}" }
  .set { plinkCh }

/*--------------------------------------------------
  GWAS Analysis 1 with SAIGE - Fit the null mixed-model
---------------------------------------------------*/

process fit_null_glmm {
  tag "$name"
  publishDir params.outdir, mode: 'copy'

  input:
  set val(plink_GRM_snps), file(bed), file(bim), file(fam) from plinkCh
  each file(phenoFile) from phenoCh

  output:
  file "*" into fit_null_glmm_results

  script:
  // Rscript bin/step1_fitNULLGLMM.R\
  //     --plinkFile=testdata/sampleA \
  //     --phenoFile=testdata/new_sample.phe \
  //     --phenoCol=PHE \
  //     --sampleIDColinphenoFile=IID \
  //     --traitType=binary        \
  //     --outputPrefix=step1_out \
  //     --outputPrefix_varRatio=step1_ \
  //     --nThreads=2 \
  //     --LOCO=FALSE \
  //     --IsOverwriteVarianceRatioFile=TRUE
  """
  step1_fitNULLGLMM.R \
    --plinkFile=${plink_GRM_snps} \
    --phenoFile=${phenoFile} \
    --phenoCol=${params.phenoCol} \
    --sampleIDColinphenoFile=IID \
    --traitType=${params.traitType} \
    --outputPrefix=step1_${params.phenoCol}_out \
    --outputPrefix_varRatio=step1_${params.phenoCol} \
    --nThreads=${task.cpus} ${saigeStep1ExtraFlags}
  """
}