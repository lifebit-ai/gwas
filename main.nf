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
if (params.trait_type == 'binary') {
  inv_normalisation = 'FALSE'
}
else if (params.trait_type == 'quantitative') {
  inv_normalisation = 'TRUE'
}
else exit 1, "Trait type is not recognised. Please check input for --trait_type parameter."


ch_pheno = params.pheno_data ? Channel.value(file(params.pheno_data)) : Channel.empty()
(phenoCh_gwas_filtering, phenoCh) = ch_pheno.into(2)

Channel
  .fromFilePairs("${params.grm_plink_input}",size:3, flat : true)
  .ifEmpty { exit 1, "PLINK files not found: ${params.grm_plink_input}.\nPlease specify a valid --grm_plink_input value. eg. testdata/*.{bed,bim,fam}" }
  .set { plinkCh }
Channel
  .fromPath(params.plink_keep_pheno)
  .set {plink_keep_pheno_ch}
Channel
  .fromPath(params.vcfs_list)
  .ifEmpty { exit 1, "Cannot find CSV VCFs file : ${params.vcfs_list}" }
  .splitCsv(skip:1)
  .map { chr, vcf, index -> [file(vcf).simpleName, chr, file(vcf), file(index)] }
  .set { vcfsCh }
Channel
  .fromPath(params.gwas_cat)
  .ifEmpty { exit 1, "Cannot find GWAS catalogue CSV  file : ${params.gwas_cat}" }
  .set { ch_gwas_cat }



  /*--------------------------------------------------
  Pre-GWAS filtering - download, filter and convert VCFs
  ---------------------------------------------------*/
if (params.trait_type == 'binary'){
  process gwas_filtering_bin {
    tag "$name"
    publishDir "${params.outdir}/gwas_filtering", mode: 'copy'

    input:
    set val(name), val(chr), file(vcf), file(index) from vcfsCh
    each file(phe_file) from phenoCh_gwas_filtering
    each file(plink_keep_file) from plink_keep_pheno_ch

    output:
    set val(name), val(chr), file("${name}.filtered_final.vcf.gz"), file("${name}.filtered_final.vcf.gz.csi") into filteredVcfsCh
    
    script:
    // TODO: (High priority) Only extract needed individuals from VCF files with `bcftools -S samples.txt` - get from samples file?
    // TODO: (Not required) `bcftools -T sites_to_extract.txt`
    // Optional parameters
    extra_plink_filter_missingness_options = params.plink_keep_pheno != "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/nodata" ? "--keep ${plink_keep_file}" : ""
    """
    # Download, filter and convert (bcf or vcf.gz) -> vcf.gz
    bcftools view -q ${params.q_filter} $vcf -Oz -o ${name}_filtered.vcf.gz
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
      --pheno $phe_file \
      --pheno-name PHE \
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
      --pheno $phe_file \
      --pheno-name PHE \
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
}


if (params.trait_type != 'binary') {
  process gwas_filtering_qt {
  tag "$name"
  publishDir "${params.outdir}/gwas_filtering", mode: 'copy'

  input:
  set val(name), val(chr), file(vcf), file(index) from vcfsCh
  each file(phe_file) from phenoCh_gwas_filtering
  each file(plink_keep_file) from plink_keep_pheno_ch

  output:
  set val(name), val(chr), file("${name}.filtered_final.vcf.gz"), file("${name}.filtered_final.vcf.gz.csi") into filteredVcfsCh

  script:
  // TODO: (High priority) Only extract needed individuals from VCF files with `bcftools -S samples.txt` - get from samples file?
  // TODO: (Not required) `bcftools -T sites_to_extract.txt`
  // Optional parameters
  extra_plink_filter_missingness_options = params.plink_keep_pheno != "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/nodata" ? "--keep ${plink_keep_file}" : ""
  """
  # Download, filter and convert (bcf or vcf.gz) -> vcf.gz
  bcftools view -q ${params.q_filter} $vcf -Oz -o ${name}_filtered.vcf.gz
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

  #Filter HWE
  plink \
    --bfile ${name}_filtered \
    --pheno $phe_file \
    --pheno-name PHE \
    --allow-no-sex \
    --hwe ${params.thres_HWE} midp \
    --out ${name}.misHWEfiltered \
    --make-just-bim \
    --1 \
    --keep-allele-order \
    ${extra_plink_filter_missingness_options}

  bcftools view ${name}_filtered.vcf.gz | awk -F '\\t' 'NR==FNR{c[\$1\$4\$6\$5]++;next}; c[\$1\$2\$4\$5] > 0' ${name}.misHWEfiltered.bim - | bgzip > ${name}.filtered_temp.vcf.gz
  bcftools view -h ${name}_filtered.vcf.gz -Oz -o ${name}_filtered.header.vcf.gz
  cat ${name}_filtered.header.vcf.gz ${name}.filtered_temp.vcf.gz > ${name}.filtered_final.vcf.gz
  bcftools index ${name}.filtered_final.vcf.gz
  """
  }
}


/*--------------------------------------------------
  GWAS Analysis 1 with SAIGE - Fit the null mixed-model
---------------------------------------------------*/


process gwas_1_fit_null_glmm_bin {
  tag "$plink_grm_snps"
  publishDir "${params.outdir}/gwas_1_fit_null_glmm", mode: 'copy'

  input:
  set val(plink_grm_snps), file(bed), file(bim), file(fam) from plinkCh
  each file(phenoFile) from phenoCh

  output:
  file "*" into fit_null_glmm_results
  file ("step1_${phenoFile.baseName}_out.rda") into rdaCh
  file ("step1_${phenoFile.baseName}.varianceRatio.txt") into varianceRatioCh

  script:
  """
  step1_fitNULLGLMM.R \
    --plinkFile=${plink_grm_snps} \
    --phenoFile="${phenoFile}" \
    --phenoCol="PHE" \
    --invNormalize=${inv_normalisation} \
    --traitType=binary       \
    --sampleIDColinphenoFile=IID \
    --outputPrefix="step1_${phenoFile.baseName}_out" \
    --outputPrefix_varRatio="step1_${phenoFile.baseName}" \
    --nThreads=${task.cpus} ${params.saige_step1_extra_flags}
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
  file("*.SAIGE.gwas.txt") into ch_saige_output

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
    --SAIGEOutputFile="step2_SPAtests.${name}.SAIGE.gwas.txt" \
    --numLinesOutput=2 \
    --IsOutputAFinCaseCtrl=TRUE \
    --IsDropMissingDosages=FALSE \
    --IsOutputNinCaseCtrl=TRUE \
    --IsOutputHetHomCountsinCaseCtrl=TRUE
  """
}

/*--------------------------------------------------
  GWAS Analysis 2 with SAIGE - Combine SAIGE outputs
---------------------------------------------------*/

process prepare_files {
  tag "preparation_files"
  publishDir "${params.outdir}/MultiQC/", mode: 'copy'

  input:
  file(saige_output) from ch_saige_output.collect()

  output:
  set file("*top_n.csv"), file("*${params.output_tag}.csv") into ch_report_input

  script:

  """

  # creates 2 .csv files, saige_results_<params.output_tag>.csv, saige_results_top_n.csv
  concat_chroms.R \
    --saige_output_name='saige_results' \
    --filename_pattern='${params.saige_filename_pattern}' \
    --output_tag='${params.output_tag}' \
    --top_n_sites=${params.top_n_sites} \
    --max_top_n_sites=${params.max_top_n_sites}
  """
}

process create_report {
tag "report"
publishDir "${params.outdir}/MultiQC/", mode: 'copy'

input:
file(saige_output) from ch_report_input.collect()
file(gwas_cat) from ch_gwas_cat

output:
file "multiqc_report.html" into ch_report_outputs
set file("*png"), file("*ipynb"), file("*csv") into ch_report_outputs_all

script:

"""
cp /opt/bin/* .

# creates gwascat_subset.csv
subset_gwascat.R \
  --saige_output='saige_results_${params.output_tag}.csv' \
  --gwas_cat='${gwas_cat}'

# creates <params.output_tag>_manhattan.png with analysis.csv as input
manhattan.R \
  --saige_output='saige_results_${params.output_tag}.csv' \
  --output_tag='${params.output_tag}'

# creates <params.output_tag>_qqplot_ci.png with analysis.csv as input
qqplot.R \
  --saige_output='saige_results_${params.output_tag}.csv' \
  --output_tag='${params.output_tag}'

# Generates the report
Rscript -e "rmarkdown::render('gwas_report.Rmd', params = list(manhattan='${params.output_tag}_manhattan.png',qqplot='${params.output_tag}_qqplot_ci.png', gwascat='gwascat_subset.csv', saige_results='saige_results_top_n.csv', trait_type='${params.trait_type}'))"
mv gwas_report.html multiqc_report.html

# Generates the ipynb
jupytext --to ipynb gwas_report.Rmd
"""
}

