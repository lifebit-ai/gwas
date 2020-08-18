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

/*--------------------------------------------------
  Channel setup
---------------------------------------------------*/
Channel
  .fromPath(params.phenoFile)
  .ifEmpty { exit 1, "Pheno file not found: ${params.phenoFile}" }
  .into { phenoCh; phenoCh_gwas_filtering}
Channel
  .fromFilePairs("${params.plinkFile}",size:3, flat : true)
  .ifEmpty { exit 1, "PLINK files not found: ${params.plinkFile}" }
  .set { plinkCh }
Channel
  .fromPath(params.plink_keep_pheno)
  .set {plink_keep_pheno_ch}
Channel
  .fromPath(params.vcfsList)
  .ifEmpty { exit 1, "Cannot find CSV VCFs file : ${params.vcfsList}" }
  .splitCsv(skip:1)
  .map { chr, vcf, index -> [file(vcf).simpleName, chr, file(vcf), file(index)] }
  .set { vcfsCh }

/*--------------------------------------------------
  Pre-GWAS filtering - download, filter and convert VCFs
---------------------------------------------------*/

process gwas_filtering {
  tag "$name"
  publishDir "${params.outdir}/gwas_filtering", mode: 'copy'

  input:
  set val(name), val(chr), file(vcf), file(index) from vcfsCh
  each file(phenofile) from phenoCh_gwas_filtering
  each file(plink_keep_file) from plink_keep_pheno_ch

  output:
  set val(name), val(chr), file("${name}.filtered_final.vcf.gz"), file("${name}.filtered_final.vcf.gz.csi") into filteredVcfsCh
  file("${name}_filtered.{bed,bim,fam}") into plinkTestCh

  script:
  // TODO: (High priority) Only extract needed individuals from VCF files with `bcftools -S samples.txt` - get from samples file?
  // TODO: (Not required) `bcftools -T sites_to_extract.txt`
  // Optional parameters
  extra_plink_filter_missingness_options = params.plink_keep_pheno != "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/nodata" ? "--keep ${plink_keep_file}" : ""
  """
  # Download, filter and convert (bcf or vcf.gz) -> vcf.gz
  bcftools view -q ${params.qFilter} $vcf -Oz -o ${name}_filtered.vcf.gz
  bcftools index ${name}_filtered.vcf.gz
 
  # Create PLINK binary from vcf.gz
  plink2 \
    --make-bed \
    --set-missing-var-ids @:#,\\\$r,\\\$a \
    --vcf ${name}_filtered.vcf.gz \
    --out ${name}_filtered \
    --vcf-half-call m \
    --double-id \
    --set-hh-missing \
    --new-id-max-allele-len 60 missing

  #Filter missingness
  plink \
    --bfile ${name}_filtered \
    --pheno $phenofile \
    --pheno-name ${params.phenoCol} \
    --allow-no-sex \
    --test-missing midp \
    --out ${name} \
    --1 \
    --keep-allele-order \
    ${extra_plink_filter_missingness_options}

  awk '\$5 < ${params.thres_m} {print}' ${name}.missing > ${name}.missing_FAIL

  #Filter HWE
  plink \
    --bfile ${name}_filtered \
    --pheno $phenofile \
    --pheno-name ${params.phenoCol} \
    --allow-no-sex \
    --hwe ${params.thres_HWE} midp \
    --out ${name}.misHWEfiltered \
    --make-just-bim \
    --exclude ${name}.missing_FAIL \
    --1 \
    --keep-allele-order \
    ${extra_plink_filter_missingness_options}

  bcftools view ${name}_filtered.vcf.gz | awk -F '\\t' 'NR==FNR{c[\$1\$4\$6\$5]++;next}; c[\$1\$2\$4\$5] > 0' ${name}.misHWEfiltered.bim - | bgzip > ${name}.filtered_temp.vcf.gz
  bcftools view -h ${name}_filtered.vcf.gz -Oz -o ${name}_filtered.header.vcf.gz
  cat ${name}_filtered.header.vcf.gz ${name}.filtered_temp.vcf.gz > ${name}.filtered_final.vcf.gz
  bcftools index ${name}.filtered_final.vcf.gz
  """
}

/*--------------------------------------------------
  GWAS Analysis 1 with SAIGE - Fit the null mixed-model
---------------------------------------------------*/

process gwas_1_fit_null_glmm {
  tag "$plink_GRM_snps"
  publishDir "${params.outdir}/gwas_1_fit_null_glmm", mode: 'copy'

  input:
  set val(plink_GRM_snps), file(bed), file(bim), file(fam) from plinkCh
  each file(phenoFile) from phenoCh

  output:
  file "*" into fit_null_glmm_results
  file ("step1_${params.phenoCol}_out.rda") into rdaCh
  file ("step1_${params.phenoCol}.varianceRatio.txt") into varianceRatioCh

  script:
  """
  step1_fitNULLGLMM.R \
    --plinkFile=${plink_GRM_snps} \
    --phenoFile=${phenoFile} \
    --phenoCol=${params.phenoCol} \
    --sampleIDColinphenoFile=IID \
    --traitType=${params.traitType} \
    --outputPrefix=step1_${params.phenoCol}_out \
    --outputPrefix_varRatio=step1_${params.phenoCol} \
    --nThreads=${task.cpus} ${params.saigeStep1ExtraFlags}
  """
}

/*--------------------------------------------------
  GWAS Analysis 2 with SAIGE - Perform mixed-model association testing
---------------------------------------------------*/

process gwas_2_spa_tests {
  tag "$name"
  publishDir "${params.outdir}/gwas_2_spa_tests", mode: 'copy'

  input:
  set val(name), val(chr), file(vcf), file(index) from filteredVcfsCh
  each file(rda) from rdaCh
  each file(varianceRatio) from varianceRatioCh

  output:
  file "*" into results

  script:
  """
  step2_SPAtests.R \
    --vcfFile=${vcf} \
    --vcfFileIndex=${index} \
    --vcfField=GT \
    --chrom=${chr} \
    --minMAC=20 \
    --sampleFile=day0_covid.samples \
    --GMMATmodelFile=${rda} \
    --varianceRatioFile=${varianceRatio} \
    --SAIGEOutputFile=${params.phenoCol}.${name}.SAIGE.gwas.txt \
    --numLinesOutput=2 \
    --IsOutputAFinCaseCtrl=TRUE \
    --IsDropMissingDosages=FALSE \
    --IsOutputNinCaseCtrl=TRUE \
    --IsOutputHetHomCountsinCaseCtrl=TRUE
  """
}

// TODO: Post-GWAS filtering
// TODO: Creation of post-GRM plink files (bcftools) - only needs to be run once
