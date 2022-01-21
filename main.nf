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

summary['saige']                          = params.saige
summary['bolt_lmm']                       = params.bolt_lmm
summary['regenie']                        = params.regenie

summary['genotype_files_list']            = params.genotype_files_list
summary['genotype_format']                = params.genotype_format
summary['grm_plink_input']                = params.grm_plink_input
summary['pheno_data']                     = params.pheno_data
summary['covariate_cols']                 = params.covariate_cols

summary['input_folder_location']          = params.input_folder_location
summary['file_pattern']                   = params.file_pattern
summary['file_suffix']                    = params.file_suffix
summary['index_suffix']                   = params.index_suffix
summary['number_of_files_to_process']     = params.number_of_files_to_process

summary['q_filter']                       = params.q_filter
summary['miss_test_p_threshold']          = params.miss_test_p_threshold
summary['variant_missingness']            = params.variant_missingness
summary['hwe_threshold']                  = params.hwe_threshold
summary['trait_type']                     = params.trait_type
summary['saige_step1_extra_flags']        = params.saige_step1_extra_flags
summary['saige_analysis_type']            = params.saige_analysis_type
summary['gwas_cat']                       = params.gwas_cat
summary['output_tag']                     = params.output_tag
summary['top_n_sites']                    = params.top_n_sites
summary['max_top_n_sites']                = params.max_top_n_sites
summary['saige_filename_pattern']         = params.saige_filename_pattern

summary['outdir']                         = params.outdir

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"


def extractInt( String input ) {
  return input.replaceAll("[^0-9]", "").toInteger()
}

def get_chromosome( file ) {
    // using RegEx to extract chromosome number from file name
    regexpPE = /(?:chr)[a-zA-Z0-9]+/
    (file =~ regexpPE)[0].replaceAll('chr','')
    
}

def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
        log.warn "Unknown parameter: ${it}"
        return false
    }
    return true
}

def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
}

def defineSaigeAnalysisTypeList() {
    return [
        'additive',
        'recessive',
        'dominant'

    ]
}

def defineFormatList() {
    return [
      'bgen',
      'vcf'
    ]
}
/*--------------------------------------------------
  Channel setup
---------------------------------------------------*/
if (params.input_folder_location) {
Channel.fromPath("${params.input_folder_location}/**${params.file_pattern}*.{${params.file_suffix},${params.index_suffix}}")
       .map { it -> [ get_chromosome(file(it).simpleName.minus(".${params.index_suffix}").minus(".${params.file_suffix}")), "s3:/"+it] }
       .groupTuple(by:0)
       .map { chr, files_pair -> [ chr, files_pair[0], files_pair[1] ] }
       .map { chr, vcf, index -> [ file(vcf).simpleName, chr, file(vcf), file(index) ] }
       .take( params.number_of_files_to_process )
       .set { inputVcfCh }
}

projectDir = workflow.projectDir

Channel
  .fromPath("${projectDir}/bin/concat_covariates.R",  type: 'file', followLinks: false)
  .set { ch_concat_covariates_r }

if (params.trait_type == 'binary') {
  inv_normalisation = 'FALSE'
}
else if (params.trait_type == 'quantitative') {
  inv_normalisation = 'TRUE'
}
else exit 1, "Trait type is not recognised. Please check input for --trait_type parameter."


params.output_tag ? Channel.value(params.output_tag).into {ch_output_tag_report; ch_output_tag } : Channel.value(params.phenotype_colname).into {ch_output_tag_report; ch_output_tag }

saige_analysis_list = defineSaigeAnalysisTypeList()
saige_analysis = params.saige_analysis_type
if (saige_analysis.contains(',')) exit 1, 'You can choose only one analysis model, see --help for more information'
if (!checkParameterExistence(saige_analysis, saige_analysis_list)) exit 1, "Unknown analysis mode: ${saige_analysis}, see --help for more information"


ch_pheno = params.pheno_data ? Channel.value(file(params.pheno_data)) : Channel.empty()
(phenoCh_gwas_filtering, ch_pheno_pca, ch_pheno_for_saige, ch_pheno_for_regenie, ch_pheno_for_bolt_lmm, phenoCh, ch_pheno_vcf2plink, ch_pheno_bgen) = ch_pheno.into(8)
ch_covariate_cols = params.covariate_cols ? Channel.value(params.covariate_cols) : "null"

if (params.grm_plink_input) {
  Channel
  .fromFilePairs("${params.grm_plink_input}",size:3, flat : true)
  .ifEmpty { exit 1, "PLINK files not found: ${params.grm_plink_input}.\nPlease specify a valid --grm_plink_input value. eg. testdata/*.{bed,bim,fam}" }
  .into { external_ch_plink_pruned_saige; external_ch_plink_pruned_bolt_lmm;   external_ch_plink_pruned_pca }
}
if (params.genotype_files_list && params.genotype_format == 'vcf') {
Channel
  .fromPath(params.genotype_files_list)
  .ifEmpty { exit 1, "Cannot find CSV VCFs file : ${params.genotype_files_list}" }
  .splitCsv(skip:1)
  .map { chr, vcf, index -> [file(vcf).simpleName, chr, file(vcf), file(index)] }
  .set { inputVcfCh }
}
else if (params.genotype_files_list && params.genotype_format == 'bgen') {
  Channel
  .fromPath(params.genotype_files_list)
  .ifEmpty { exit 1, "Cannot find CSV file containing paths to .bgen/.sample files: ${params.genotype_files_list}" }
  .splitCsv(skip:1)
  .map { chr, bgen, bgi_index -> [file(bgen).simpleName, chr, file(bgen), file(bgi_index)] }
  .set { ch_input_bgen }

  Channel
  .fromPath(params.bgen_sample_file)
  .set { ch_bgen_sample_file }

}
else if (!params.genotype_files_list && !params.input_folder_location) {
  exit 1, "File containing paths to genotype files not specified. Please specify a .csv file with paths using --genotype_files_list parameter, or a input s3 path using --input_folder_location parameter."
}
else if (params.genotype_format != 'vcf' && params.genotype_format != 'bgen') {
  exit 1, "Genotype format not supported. Please choose out of valid formats."

}

if (!params.phenotype_colname) {
  exit 1, "Phenotype column name has to be specified."
}

Channel
  .fromPath(params.high_LD_long_range_regions)
  .ifEmpty { exit 1, "Cannot find file containing long-range LD regions for exclusion : ${params.high_LD_long_range_regions}" }
  .set { ch_high_ld_regions }
  
Channel
  .fromPath(params.gwas_cat)
  .ifEmpty { exit 1, "Cannot find GWAS catalogue CSV  file : ${params.gwas_cat}" }
  .set { ch_gwas_cat }

if (params.run_pca) {
  Channel
      .from( 1..params.number_pcs )
      .flatMap { it -> "PC$it" }
      .toList()
      .set { ch_pca_cols }
} else {
  ch_pca_cols = Channel.fromList(['null'])
}

if (params.ld_scores) {
  Channel
      .fromPath(params.ld_scores)
      .ifEmpty { exit 1, "Cannot find file containing LD scores : ${params.ld_scores}" }
      .set { ch_ld_scores }

}
    
  /*--------------------------------------------------
  Pre-GWAS filtering - download, filter and convert VCFs
  ---------------------------------------------------*/
if (params.genotype_format == 'vcf') {
  process vcf2plink {
    tag "$name"
    label 'high_memory'
    publishDir "${params.outdir}/gwas_filtering", mode: 'copy'

    input:
    set val(name), val(chr), file(vcf), file(index) from inputVcfCh
    each file(phe_file) from ch_pheno_vcf2plink

    output:
    set val(name), val(chr), file('*.bed'), file('*.bim'), file('*.fam') into filteredPlinkCh

    script:
    plink_memory = extractInt(task.memory.toString()) * 1000
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
      --memory ${plink_memory} \
      --double-id \
      --set-hh-missing \
      --new-id-max-allele-len 60 missing
  """   
  }
}
else if (params.genotype_format == 'bgen') {
  process bgen2plink {
    tag "$name"
    label "mid_memory"
    publishDir "${params.outdir}/gwas_filtering", mode: 'copy'

    input:
    set val(name), val(chr), file(bgen), file(index) from ch_input_bgen
    each file(phe_file) from ch_pheno_bgen
    each file(sample_file) from ch_bgen_sample_file

    output:
    set val(name), val(chr), file('*.bed'), file('*.bim'), file('*.fam') into filteredPlinkCh

    script:
    plink_memory = extractInt(task.memory.toString()) * 1000
    """
    # Create a sample ID file for --keep
    tail -n +2 ${phe_file}| cut -f1,2 > samples.txt

    # Create PLINK binary from vcf.gz
    plink2 \
      --make-bed \
      --bgen ${bgen}\
      --out ${name}_filtered \
      --maf ${params.maf_filter} \
      --double-id \
      --memory ${plink_memory} \
      --keep samples.txt \
      --sample ${sample_file}
    """ 
  }
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
  plink_memory = extractInt(task.memory.toString()) * 1000
  if ( params.trait_type == "binary" )
    """
   plink \
     --bfile ${name}_filtered \
     --pheno $phe_file \
     --pheno-name ${params.phenotype_colname} \
     --allow-no-sex \
     --test-missing midp \
     --memory ${plink_memory} \
     --out ${name} \
     --1 \
     --keep-allele-order \

   awk '\$5 < ${params.miss_test_p_threshold} {print \$2 }' ${name}.missing > ${name}.missing_FAIL 

   plink --bfile ${name}_filtered \
     --keep-allele-order \
     --allow-no-sex \
     --exclude ${name}.missing_FAIL \
     --memory ${plink_memory} \
     --make-bed \
     --out ${name}_miss_filtered
   """
else if ( params.trait_type == "quantitative" )
  """
     plink \
     --bfile ${name}_filtered \
     --allow-no-sex \
     --missing \
     --memory ${plink_memory} \
     --out ${name} \
     --1 \
     --keep-allele-order

    awk '\$5 > ${params.variant_missingness} {print \$2 }' ${name}.lmiss > ${name}.missing_FAIL

    plink --bfile ${name}_filtered \
     --keep-allele-order \
     --allow-no-sex \
     --exclude ${name}.missing_FAIL \
     --memory ${plink_memory} \
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
  set val(name), val(chr), file("${name}.filtered_final.vcf.gz"), file("${name}.filtered_final.vcf.gz.csi") into filteredVcfsCh_saige, filteredVcfsCh_bolt_lmm, ch_filtered_vcfs_for_pruning
  file("${name}.misHWEfiltered*") into ch_filtered_plink
  
  script:
  plink_memory = extractInt(task.memory.toString()) * 1000
  """
  plink \
    --bfile ${name}_miss_filtered \
    --pheno $phe_file \
    --pheno-name ${params.phenotype_colname} \
    --memory ${plink_memory} \
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
    --memory ${plink_memory} \
    --out ${name}_filtered_vcf

  bcftools view ${name}_filtered_vcf.vcf.gz | awk -F '\\t' 'NR==FNR{c[\$1\$4\$6\$5]++;next}; c[\$1\$2\$4\$5] > 0' ${name}.misHWEfiltered.bim - | bgzip > ${name}.filtered_temp.vcf.gz
  bcftools view -h ${name}_filtered_vcf.vcf.gz -Oz -o ${name}_filtered.header.vcf.gz
  cat ${name}_filtered.header.vcf.gz ${name}.filtered_temp.vcf.gz > ${name}.filtered_final.vcf.gz
  bcftools index ${name}.filtered_final.vcf.gz
  """
}

if (params.bolt_lmm || !params.grm_plink_input) {

  process merge_plink {
      tag "merge_plink"

      publishDir "${params.outdir}", mode: 'copy'

      input:
      file("*") from ch_filtered_plink.collect()

      output:
      set file('merged.bed'), file('merged.bim'), file('merged.fam') into ch_plink_merged, ch_plink_merged_bgen



      script:
      plink_memory = extractInt(task.memory.toString()) * 1000
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
      --memory ${plink_memory} \
      --make-bed \
      --out merged

      """
  }
}

if (!params.grm_plink_input) {

  process ld_prune {
      tag "LD-prune plink set"
      publishDir "${params.outdir}", mode: 'copy'

      input:
      set file(bed), file(bim), file(fam) from ch_plink_merged
      file(long_range_ld_regions) from ch_high_ld_regions

      output:
      set val('merged_pruned'), file('merged_pruned.bed'), file('merged_pruned.bim'), file('merged_pruned.fam') into ( ch_plink_pruned_saige, ch_plink_pruned_bolt_lmm, ch_plink_pruned_pca)

      script:
      plink_memory = extractInt(task.memory.toString()) * 1000
      """
      plink \
      --bfile merged \
      --keep-allele-order \
      --indep-pairwise ${params.ld_window_size} ${params.ld_step_size} ${params.ld_r2_threshold} \
      --exclude range ${long_range_ld_regions} \
      --allow-no-sex \
      --memory ${plink_memory} \
      --out merged
      plink \
      --bfile merged \
      --keep-allele-order \
      --extract merged.prune.in \
      --make-bed \
      --memory ${plink_memory} \
      --allow-no-sex \
      --out merged_pruned 
      """
  }
}
ch_plink_pruned_for_pca = params.grm_plink_input ? external_ch_plink_pruned_pca: ch_plink_pruned_pca

process run_pca {
    tag "Run PCA"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set val(plink_prefix), file(bed), file(bim), file(fam) from ch_plink_pruned_for_pca
    each file(phenotype_file) from ch_pheno_pca
    file(concat_covariates_script) from ch_concat_covariates_r

    output:
    set file('pca_results.eigenvec'), file('pca_results.eigenval') into ch_pca_files
    file('covariates_with_PCs.tsv') into (ch_full_covariate_file_saige, ch_full_covariate_file_bolt_lmm, ch_full_covariate_file_regenie)

    when: params.run_pca

    script:
    """
    if [ \$(wc -l ${bim} | cut -d " " -f1) -lt 220 ]; then
        echo "Error: PCA requires data to contain at least 220 variants."
        exit 1
    fi
    if [ \$(wc -l ${fam} | cut -d " " -f1) -gt 5000 ]; then
    echo "### It seems that you are using a large sample size. \
    The randomised algorithm originally implemented for <Galinsky KJ et al. 2016> \
    will be used to perform the PC calculations." 

    plink2 --bfile ${plink_prefix} --pca ${params.number_pcs} approx --out pca_results
    fi
    if [ \$(wc -l ${fam} | cut -d " " -f1) -lt 5001 ]; then
    plink2 --bfile ${plink_prefix} --pca ${params.number_pcs} --out pca_results
    fi

    concat_covariates.R \
    --phenotype_colname=${params.phenotype_colname} \
    --pcs_file=pca_results.eigenvec \
    --phenotype_file=${phenotype_file} \
    --output_file=covariates_with_PCs.tsv

    """
}

/*--------------------------------------------------
  GWAS Analysis 1 with SAIGE - Fit the null mixed-model
---------------------------------------------------*/
ch_plink_input_for_grm_saige = params.grm_plink_input ? external_ch_plink_pruned_saige : ch_plink_pruned_saige
ch_plink_input_for_grm_bolt_lmm = params.grm_plink_input ? external_ch_plink_pruned_bolt_lmm : ch_plink_pruned_bolt_lmm
ch_covariate_file_for_saige = params.run_pca ? ch_full_covariate_file_saige : ch_pheno_for_saige
ch_covariate_file_for_bolt_lmm = params.run_pca ? ch_full_covariate_file_bolt_lmm : ch_pheno_for_bolt_lmm
ch_covariate_file_for_regenie = params.run_pca ? ch_full_covariate_file_regenie : ch_pheno_for_regenie


if (params.bolt_lmm || params.regenie) {

  process convert2bgen   {
    tag "convert2bgen"
    label "convert2bgen"
    publishDir "${params.outdir}/bgen", mode: 'copy'

    input:
    set file(bed), file(bim), file(fam) from ch_plink_merged_bgen

    output:
    set file('merged_bgen.bgen'), file('merged_bgen.sample') into ch_merged_bgen_bolt_lmm, ch_merged_bgen_regenie_step1, ch_merged_bgen_regenie_step2

    script: 
    """
    plink2 --bfile ${bed.baseName} --export bgen-1.2 bits=8 --out merged_bgen
    """
  }
}

if (params.regenie) {

/*--------------------------------------------------
  GWAS using REGENIE
---------------------------------------------------*/

  process regenie_step1_fit_model {
    tag "regenie_step1_fit_model"
    label 'regenie'
    publishDir "${params.outdir}/regenie", mode: 'copy'

    input:
    set file(bgen), file(sample_file) from ch_merged_bgen_regenie_step1
    file(pheno_covariates) from ch_full_covariate_file_regenie

    output:
    set file("*.loco"), file("*_pred.list") into ch_regenie_step1_pred
    file "covariates.txt" into ch_regenie_cov
    file "pheno.txt" into ch_regenie_pheno

    script:
    covariates = params.covariate_cols ? "--covarColList ${params.covariate_cols}" : ''
    """
    sed -e '1s/^.//' ${pheno_covariates}| sed 's/\t/ /g' > full_pheno_covariates.txt
    pheno_col=`awk -v RS=' ' '/${params.phenotype_colname}/{print NR; exit}' full_pheno_covariates.txt`
    cut -d' ' -f1,2,\$pheno_col full_pheno_covariates.txt > pheno.txt
    cut -d' ' --complement -f\$pheno_col full_pheno_covariates.txt > covariates.txt

    regenie \
     --step 1 \
      --bgen ${bgen} \
      --covarFile covariates.txt \
      $covariates \
      --phenoFile pheno.txt \
      --bsize 100 \
      --threads ${task.cpus} \
      ${params.trait_type == "binary" ? '--bt' : ''} \
      --lowmem \
      --lowmem-prefix tmp_rg \
      --out regenie_fit_out
    """

  }

  process regenie_step2_association_testing {
    tag "regenie_step2_association_tests"
    label 'regenie'
    publishDir "${params.outdir}/regenie", mode: 'copy'
    stageInMode 'copy'

    input:
    set file(loco), file(pred) from ch_regenie_step1_pred
    set file(bgen), file(sample_file) from ch_merged_bgen_regenie_step1
    file(covariates_file) from ch_regenie_cov
    file(pheno) from ch_regenie_pheno

    output:
    file "*" into ch_regenie_step2_assoc

    script:
    covariates = params.covariate_cols ? "--covarColList ${params.covariate_cols}" : ''
    """
    mv ${loco} tmp.txt
    cp tmp.txt regenie_fit_out_1.loco
    rm tmp.txt
    mv ${pred} tmp.txt
    cp tmp.txt regenie_fit_out_pred.list
    regenie \
      --step 2 \
      --bgen ${bgen} \
      --covarFile ${covariates_file} \
      --phenoFile ${pheno} \
      --bsize 200 \
      ${params.trait_type == "binary" ? '--bt' : ''} \
      --minMAC ${params.regenie_min_mac} \
      --minINFO ${params.regenie_min_imputation_score} \
      --firth --approx \
      --pThresh 0.01 \
      --gz \
      --ignore-pred \
      --out regenie_firth

    """

  }


}

/*--------------------------------------------------
  GWAS using BOLT-LMM
---------------------------------------------------*/
if (params.bolt_lmm) {
  process run_bolt_lmm {
  tag "$name"
  label 'bolt_lmm'
  publishDir "${params.outdir}/bolt_lmm", mode: 'copy'

  input:
  set val(plink_grm_prefix), file(pruned_bed), file(pruned_bim), file(pruned_fam) from ch_plink_input_for_grm_bolt_lmm
  each file(full_covariate_file) from ch_covariate_file_for_bolt_lmm
  set val(name), val(chr), file(vcf), file(index) from filteredVcfsCh_bolt_lmm
  file(ld_scores) from ch_ld_scores
  set file(bgen), file(sample_file) from ch_merged_bgen_bolt_lmm


  output:
  file "*" into ch_bolt_lmm_results

  script:
  """
  sed -e '1s/^.//' ${full_covariate_file}| sed 's/\t/ /g' > pheno_covariates.txt


  bolt --bfile=${plink_grm_prefix} \
       --phenoFile=pheno_covariates.txt \
       --phenoCol=${params.phenotype_colname} \
       --LDscoresFile=${ld_scores} \
       --qCovarCol=${params.bolt_lmm_quant_covariates} \
       --verboseStats \
       --lmm \
       --bgenFile=${bgen} \
       --sampleFile=${sample_file} \
       --LDscoresMatchBp \
       --covarFile=pheno_covariates.txt \
      --covarCol=${params.bolt_lmm_categ_covariates} \
      --statsFileBgenSnps=bgen_snps_stats.gz \
      --statsFile=stats.tab 

  """

}
}

if (params.saige) {
process gwas_1_fit_null_glmm_with_pcs {
  tag "$plink_grm_snps"
  label 'saige'
  publishDir "${params.outdir}/gwas_1_fit_null_glmm", mode: 'copy'

  input:
  set val(plink_grm_prefix), file(pruned_bed), file(pruned_bim), file(pruned_fam) from ch_plink_input_for_grm_saige
  each file(full_covariate_file) from ch_covariate_file_for_saige
  val(cov_columns) from ch_covariate_cols
  val(pc_columns) from ch_pca_cols.collect()
    
  output:
  file "*" into fit_null_glmm_results
  file ("step1_${params.phenotype_colname}_out.rda") into rdaCh
  file ("step1_${params.phenotype_colname}.varianceRatio.txt") into varianceRatioCh

  script:
  pc_cols = pc_columns.join(',')
  covariate_columns = params.run_pca ? cov_columns + ',' + pc_cols: cov_columns
  cov_columns_arg = params.covariate_cols ? "--covarColList=${covariate_columns}" : ""

  """
  step1_fitNULLGLMM.R \
    --plinkFile=${plink_grm_prefix} \
    --phenoFile="${full_covariate_file}" \
    ${cov_columns_arg} \
    --phenoCol=${params.phenotype_colname} \
    --invNormalize=${inv_normalisation} \
    --traitType=${params.trait_type}       \
    --sampleIDColinphenoFile=IID \
    --outputPrefix="step1_${params.phenotype_colname}_out" \
    --outputPrefix_varRatio="step1_${params.phenotype_colname}" \
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
  set val(name), val(chr), file(vcf), file(index) from filteredVcfsCh_saige
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
    --analysisType=${params.saige_analysis_type} \
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
  --output_tag='${final_output_tag}' \
  --p_value_cutoff=${params.p_significance_threshold}

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
}

