#!/usr/bin/env nextflow
/*
========================================================================================
                         lifebit-ai/gwas
========================================================================================
 lifebit-ai/gwas GWAS pipeline using SAIGE linear mixed model approach for association testing
 #### Homepage / Documentation
 https://github.com/lifebit-ai/gwas
----------------------------------------------------------------------------------------
*/


/*---------------------------------------------------
  Define and show header with all params information 
-----------------------------------------------------*/

// Header log info

def summary = [:]

if (workflow.revision) summary['Pipeline Release'] = workflow.revision

summary['Max Resources']                  = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
summary['Output dir']                     = params.outdir
summary['Launch dir']                     = workflow.launchDir
summary['Working dir']                    = workflow.workDir
summary['Script dir']                     = workflow.projectDir
summary['User']                           = workflow.userName

summary['vcfs_list']                      = params.vcfs_list
summary['grm_plink_input']                = params.grm_plink_input
summary['pheno_data']                     = params.pheno_data
summary['covariate_cols']                 = params.covariate_cols

summary['q_filter']                       = params.q_filter
summary['miss_test_p_threshold']          = params.miss_test_p_threshold
summary['variant_missingness']            = params.variant_missingness
summary['hwe_threshold']                  = params.hwe_threshold
summary['trait_type']                     = params.trait_type
summary['saige_step1_extra_flags']        = params.saige_step1_extra_flags
summary['gwas_cat']                       = params.gwas_cat
summary['output_tag']                     = params.output_tag
summary['top_n_sites']                    = params.top_n_sites
summary['max_top_n_sites']                = params.max_top_n_sites
summary['saige_filename_pattern']         = params.saige_filename_pattern

summary['outdir']                         = params.outdir

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"


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


params.output_tag ? Channel.value(params.output_tag).into {ch_output_tag_report; ch_output_tag } : Channel.value(params.phenotype_colname).into {ch_output_tag_report; ch_output_tag }


ch_pheno = params.pheno_data ? Channel.value(file(params.pheno_data)) : Channel.empty()
(phenoCh_gwas_filtering, ch_pheno_for_saige, phenoCh, ch_pheno_vcf2plink) = ch_pheno.into(4)
ch_covariate_cols = params.covariate_cols ? Channel.value(params.covariate_cols) : "null"

if (params.grm_plink_input) {
  Channel
  .fromFilePairs("${params.grm_plink_input}",size:3, flat : true)
  .ifEmpty { exit 1, "PLINK files not found: ${params.grm_plink_input}.\nPlease specify a valid --grm_plink_input value. eg. testdata/*.{bed,bim,fam}" }
  .set { external_ch_plink_pruned }
}
Channel
  .fromPath(params.vcfs_list)
  .ifEmpty { exit 1, "Cannot find CSV VCFs file : ${params.vcfs_list}" }
  .splitCsv(skip:1)
  .map { chr, vcf, index -> [file(vcf).simpleName, chr, file(vcf), file(index)] }
  .into { vcfsCh; inputVcfCh}
Channel
  .fromPath(params.high_LD_long_range_regions)
  .ifEmpty { exit 1, "Cannot find file containing long-range LD regions for exclusion : ${params.high_LD_long_range_regions}" }
  .set { ch_high_ld_regions }
  
Channel
  .fromPath(params.gwas_cat)
  .ifEmpty { exit 1, "Cannot find GWAS catalogue CSV  file : ${params.gwas_cat}" }
  .set { ch_gwas_cat }



  /*--------------------------------------------------
  Pre-GWAS filtering - download, filter and convert VCFs
  ---------------------------------------------------*/
process vcf2plink {
  tag "$name"
  publishDir "${params.outdir}/gwas_filtering", mode: 'copy'

  input:
  set val(name), val(chr), file(vcf), file(index) from inputVcfCh
  each file(phe_file) from ch_pheno_vcf2plink

  output:
  set val(name), val(chr), file('*.bed'), file('*.bim'), file('*.fam') into filteredPlinkCh

  script:
  """
  # Download, filter and convert (bcf or vcf.gz) -> vcf.gz
  tail -n +2 ${phe_file}| cut -f2 > samples.txt
  bcftools view -S samples.txt $vcf -Oz -o ${name}_downsampled.vcf.gz
  bcftools view -q ${params.q_filter} ${name}_downsampled.vcf.gz -Oz -o ${name}_filtered.vcf.gz
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
"""   
}

process filter_missingness {
  tag "$name"
  publishDir "${params.outdir}/filter_miss", mode: 'copy'
  input:
  set val(name), val(chr), file(bed), file(bim), file(fam) from filteredPlinkCh
  each file(phe_file) from phenoCh

  output:
  set val(name), val(chr), file('*_miss_filtered.bed'), file('*_miss_filtered.bim'), file('*_miss_filtered.fam') into ch_plink_for_hwe
  
  script:
  if ( params.trait_type == "binary" )
    """
   plink \
     --bfile ${name}_filtered \
     --pheno $phe_file \
     --pheno-name ${params.phenotype_colname} \
     --allow-no-sex \
     --test-missing midp \
     --out ${name} \
     --1 \
     --keep-allele-order \

   awk '\$5 < ${params.miss_test_p_threshold} {print \$2 }' ${name}.missing > ${name}.missing_FAIL 

   plink --bfile ${name}_filtered \
     --keep-allele-order \
     --allow-no-sex \
     --exclude ${name}.missing_FAIL \
     --make-bed \
     --out ${name}_miss_filtered
   """
else if ( params.trait_type == "quantitative" )
  """
     plink \
     --bfile ${name}_filtered \
     --allow-no-sex \
     --missing \
     --out ${name} \
     --1 \
     --keep-allele-order

    awk '\$5 > ${params.variant_missingness} {print \$2 }' ${name}.lmiss > ${name}.missing_FAIL

    plink --bfile ${name}_filtered \
     --keep-allele-order \
     --allow-no-sex \
     --exclude ${name}.missing_FAIL \
     --make-bed \
     --out ${name}_miss_filtered
  """

}

process calculate_hwe {
  tag "$name"
  publishDir "${params.outdir}/filter_miss", mode: 'copy'

  input:
  set val(name), val(chr), file(bed), file(bim), file(fam) from ch_plink_for_hwe
  each file(phe_file) from phenoCh_gwas_filtering

  output:
  set val(name), val(chr), file("${name}.filtered_final.vcf.gz"), file("${name}.filtered_final.vcf.gz.csi") into filteredVcfsCh, ch_filtered_vcfs_for_pruning
  file("${name}.misHWEfiltered*") into ch_filtered_plink
  
  script:
  """
  plink \
    --bfile ${name}_miss_filtered \
    --pheno $phe_file \
    --pheno-name ${params.phenotype_colname} \
    --allow-no-sex \
    --hwe ${params.hwe_threshold} ${params.hwe_test} \
    --out ${name}.misHWEfiltered \
    --make-bed \
    --1 \
    --keep-allele-order

  plink \
    --bfile ${name}.misHWEfiltered  \
    --keep-allele-order \
    --recode vcf-iid bgz \
    --out ${name}_filtered_vcf

  bcftools view ${name}_filtered_vcf.vcf.gz | awk -F '\\t' 'NR==FNR{c[\$1\$4\$6\$5]++;next}; c[\$1\$2\$4\$5] > 0' ${name}.misHWEfiltered.bim - | bgzip > ${name}.filtered_temp.vcf.gz
  bcftools view -h ${name}_filtered_vcf.vcf.gz -Oz -o ${name}_filtered.header.vcf.gz
  cat ${name}_filtered.header.vcf.gz ${name}.filtered_temp.vcf.gz > ${name}.filtered_final.vcf.gz
  bcftools index ${name}.filtered_final.vcf.gz
  """
}


process merge_plink {
    tag "merge_plink"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    file("*") from ch_filtered_plink.collect()

    output:
    set file('merged.bed'), file('merged.bim'), file('merged.fam') into ch_plink_merged

    when: !params.grm_plink_input

    script:
    """
    ls *.bed > bed.txt
    ls *.bim > bim.txt
    ls *.fam > fam.txt
    paste bed.txt bim.txt fam.txt > merge.temp.list
    tail -n +2 merge.temp.list > merge.list
    bed_file=\$(head -n1 merge.temp.list | cut -f1)
    bed_prefix=`echo "\${bed_file%.*}"`
    plink --keep-allele-order \
    --bfile \${bed_prefix} \
    --merge-list merge.list \
    --allow-no-sex \
    --make-bed \
    --out merged

    """
}


process ld_prune {
    tag "LD-prune plink set"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set file(bed), file(bim), file(fam) from ch_plink_merged
    file(long_range_ld_regions) from ch_high_ld_regions

    output:
    set val('merged_pruned'), file('merged_pruned.bed'), file('merged_pruned.bim'), file('merged_pruned.fam') into ch_plink_pruned

    when: !params.grm_plink_input

    script:
    """
    plink \
    --bfile merged \
    --keep-allele-order \
    --indep-pairwise ${params.ld_window_size} ${params.ld_step_size} ${params.ld_r2_threshold} \
    --exclude range ${long_range_ld_regions} \
    --allow-no-sex \
    --out merged
    plink \
    --bfile merged \
    --keep-allele-order \
    --extract merged.prune.in \
    --make-bed \
    --allow-no-sex \
    --out merged_pruned 
    """
}

/*--------------------------------------------------
  GWAS Analysis 1 with SAIGE - Fit the null mixed-model
---------------------------------------------------*/
ch_plink_input_for_grm = params.grm_plink_input ? external_ch_plink_pruned : ch_plink_pruned

process gwas_1_fit_null_glmm {
  tag "$plink_grm_snps"
  label 'saige'
  publishDir "${params.outdir}/gwas_1_fit_null_glmm", mode: 'copy'

  input:
  set val(plink_grm_prefix), file(pruned_bed), file(pruned_bim), file(pruned_fam) from ch_plink_input_for_grm
  each file(phenoFile) from ch_pheno_for_saige
  val(cov_columns) from ch_covariate_cols
    
  output:
  file "*" into fit_null_glmm_results
  file ("step1_${phenoFile.baseName}_out.rda") into rdaCh
  file ("step1_${phenoFile.baseName}.varianceRatio.txt") into varianceRatioCh

  script:
  cov_columns_arg = params.covariate_cols ? "--covarColList=${cov_columns}" : ""

  """
  step1_fitNULLGLMM.R \
    --plinkFile=${plink_grm_prefix} \
    --phenoFile="${phenoFile}" \
    ${cov_columns_arg} \
    --phenoCol=${params.phenotype_colname} \
    --invNormalize=${inv_normalisation} \
    --traitType=${params.trait_type}       \
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
  label 'saige'
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
  val(output_tag) from ch_output_tag

  output:
  set file("*top_n.csv"), file("*${final_output_tag}.csv") into ch_report_input

  script:
  final_output_tag = output_tag.toString().toLowerCase().replaceAll(" ", "_")

  """

  # creates 2 .csv files, saige_results_<output_tag>.csv, saige_results_top_n.csv
  concat_chroms.R \
    --saige_output_name='saige_results' \
    --filename_pattern='${params.saige_filename_pattern}' \
    --output_tag='${final_output_tag}' \
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
val(output_tag) from ch_output_tag_report

output:
file "multiqc_report.html" into ch_report_outputs
set file("*png"), file("*ipynb"), file("*csv") into ch_report_outputs_all

script:
final_output_tag = output_tag.toString().toLowerCase().replaceAll(" ", "_")
"""
cp /opt/bin/* .

# creates gwascat_subset.csv
subset_gwascat.R \
  --saige_output='saige_results_${final_output_tag}.csv' \
  --gwas_cat='${gwas_cat}'

# creates <params.output_tag>_manhattan.png with analysis.csv as input
manhattan.R \
  --saige_output='saige_results_${final_output_tag}.csv' \
  --output_tag='${final_output_tag}'

# creates <params.output_tag>_qqplot_ci.png with analysis.csv as input
qqplot.R \
  --saige_output='saige_results_${final_output_tag}.csv' \
  --output_tag='${final_output_tag}'

# Generates the report
Rscript -e "rmarkdown::render('gwas_report.Rmd', params = list(manhattan='${final_output_tag}_manhattan.png',qqplot='${final_output_tag}_qqplot_ci.png', gwascat='gwascat_subset.csv', saige_results='saige_results_top_n.csv', trait_type='${params.trait_type}'))"
mv gwas_report.html multiqc_report.html

# Generates the ipynb
jupytext --to ipynb gwas_report.Rmd
"""
}

