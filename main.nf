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
  Initial checks
--------------------------------------------------*/

if (params.gwas_summary && params.gwas_cat_study_id) {
  exit 1, "You have provided both external gwas summary statistics and GWAS catalogue summary statistics: only one cohort can be used. \
  \nPlease use only one of these inputs (i.e. only --gwas_summary or only --gwas_cat_study_id)."
}

/*--------------------------------------------------
  Channel setup
---------------------------------------------------*/
ch_hapmap3_snplist =  params.hapmap3_snplist ? Channel.value(file(params.hapmap3_snplist)) :  "null"
ch_ld_scores_tar_bz2 =  params.ld_scores_tar_bz2 ? Channel.value(file(params.ld_scores_tar_bz2)) :  "null"
ch_query =  params.query ? Channel.value(file(params.query)) : "None"
ch_pheno_data = params.pheno_data ? Channel.value(file(params.pheno_data)) : Channel.empty()
ch_pheno_metadata = params.pheno_metadata ? Channel.value(file(params.pheno_metadata)) : Channel.empty()
ch_gwas_summary = params.gwas_summary ? Channel.value(file(params.gwas_summary)) : Channel.empty()

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
  Testing mode: 
  Change platekeys by testing data platekeys
---------------------------------------------------*/
if (params.pheno_data && params.testing){
  ch_pheno_data.into{ch_pheno_data_test}
  process switch_platekeys {
    tag "$name"
    publishDir "${params.outdir}/switch_platekeys", mode: 'copy'
    input:
    file(pheno_data) from ch_pheno_data_test

    output:
    file("${params.output_tag}_gwas.csv") into ch_pheno_data_test2
    file("${params.output_tag}_phewas.csv") into ch_pheno_data_phewas
    file("${params.output_tag}_IDs.csv") into ch_conversion_platekeys

    script:
    """
    test_data_munging.R --input_file "${pheno_data}" \
                        --ids_column "${params.test_ids_column}" \
                        --outprefix "${params.output_tag}"
    """
  }
  if (params.query){
    process transforms_cb_output_testing {
    tag "$name"
    publishDir "${params.outdir}/design_matrix", mode: 'copy'

    input:
    file(pheno_data) from ch_pheno_data_test2
    file(pheno_metadata) from ch_pheno_metadata
    file(query_file) from ch_query

    output:
    file("${params.output_tag}_.phe") into ch_transform_cb
    file("*.json") into ch_encoding_json
    file("*.csv") into ch_encoding_csv

    script:
    """
    cp /opt/bin/* .

    mkdir -p ${params.outdir}/design_matrix
    
    transform_cb_output.R --input_cb_data "$pheno_data" \
                          --input_meta_data "$pheno_metadata" \
                          --phenoCol "${params.pheno_col}" \
                          --query_file "${query_file}" \
                          --continuous_var_transformation "${params.continuous_var_transformation}" \
                          --continuous_var_aggregation "${params.continuous_var_aggregation}" \
                          --outdir "." \
                          --output_tag "${params.output_tag}"
    """
   }
  }
  if (!params.query){
    process transforms_cb_output_testing {
    tag "$name"
    publishDir "${params.outdir}/design_matrix", mode: 'copy'

    input:
    file(pheno_data) from ch_pheno_data_test2
    file(pheno_metadata) from ch_pheno_metadata

    output:
    file("${params.output_tag}_.phe") into ch_transform_cb
    file("*.json") into ch_encoding_json
    file("*.csv") into ch_encoding_csv

    script:
    """
    cp /opt/bin/* .

    mkdir -p ${params.outdir}/design_matrix
    
    transform_cb_output.R --input_cb_data "$pheno_data" \
                          --input_meta_data "$pheno_metadata" \
                          --phenoCol "${params.pheno_col}" \
                          --query_file "${ch_query}" \
                          --continuous_var_transformation "${params.continuous_var_transformation}" \
                          --continuous_var_aggregation "${params.continuous_var_aggregation}" \
                          --outdir "." \
                          --output_tag "${params.output_tag}"
    """
   }
  }
}


/*--------------------------------------------------
  Ingest output from CB
---------------------------------------------------*/
if (params.pheno_data && !params.testing){
  process transforms_cb_output {
    tag "$name"
    publishDir "${params.outdir}/design_matrix", mode: 'copy'

    input:
    file(pheno_data) from ch_pheno_data
    file(pheno_metadata) from ch_pheno_metadata
    val(query_file) from ch_query

    output:
    file("${params.output_tag}_.phe") into ch_transform_cb
    file("*.json") into ch_encoding_json
    file("*.csv") into ch_encoding_csv

    script:
    """
    cp /opt/bin/* .

    mkdir -p ${params.outdir}/design_matrix
    
    transform_cb_output.R --input_cb_data "$pheno_data" \
                          --input_meta_data "$pheno_metadata" \
                          --phenoCol "${params.pheno_col}" \
                          --query_file "${query_file}" \
                          --continuous_var_transformation "${params.continuous_var_transformation}" \
                          --continuous_var_aggregation "${params.continuous_var_aggregation}" \
                          --outdir "." \
                          --output_tag "${params.output_tag}"
    """
  }
}
//TODO: Check this later and finish it with the processes 
if (params.trait_type == 'binary' && params.pheno_data && params.case_group && params.design_mode != 'all_contrasts') {
  process add_design_matrix_case_group {
    tag "$name"
    publishDir "${params.outdir}/contrasts", mode: 'copy'

    input:
    file(pheFile) from ch_transform_cb
    file(json) from ch_encoding_json

    output:
    file("${params.output_tag}_design_matrix_control_*.phe") into phenoCh_gwas_filtering

    script:
    """
    cp /opt/bin/* .

    mkdir -p ${params.outdir}/contrasts

    create_design.R --input_file ${pheFile} \
                    --mode "${params.design_mode}" \
                    --case_group "${params.case_group}" \
                    --outdir . \
                    --output_tag ${params.output_tag} \
                    --phenoCol "${params.pheno_col}"
                      
    """
  }
}

if (params.trait_type == 'binary' && params.pheno_data && params.design_mode == 'all_contrasts') {
  process add_design_matrix_all{
    tag "$name"
    publishDir "${params.outdir}/contrasts", mode: 'copy'

    input:
    file(pheFile) from ch_transform_cb
    file(json) from ch_encoding_json

    output:
    file("${params.output_tag}_design_matrix_control_*.phe") into phenoCh_gwas_filtering

    script:
    """
    cp /opt/bin/* .

    mkdir -p ${params.outdir}/contrasts

    create_design.R --input_file ${pheFile} \
                    --mode "${params.design_mode}" \
                    --outdir . \
                    --output_tag ${params.output_tag} \
                    --phenoCol "${params.pheno_col}"
                      
    """
  }
}

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
  ch_transform_cb.into{phenoCh_gwas_filtering}
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
// Create channel for this process
phenoCh_gwas_filtering.into{phenoCh}

if (params.trait_type == 'binary'){
  process gwas_1_fit_null_glmm_bin {
    tag "$plink_grm_snps"
    publishDir "${params.outdir}/gwas_1_fit_null_glmm", mode: 'copy'

    input:
    set val(plink_grm_snps), file(bed), file(bim), file(fam) from plinkCh
    each file(phenoFile) from phenoCh

    output:
    file "*" into fit_null_glmm_results
    file ("step1_${params.pheno_col.replaceAll(/\s/,'_').replaceAll(/\(|\)/, '')}_out.rda") into rdaCh
    file ("step1_${params.pheno_col.replaceAll(/\s/,'_').replaceAll(/\(|\)/, '')}.varianceRatio.txt") into varianceRatioCh

    script:
    """
    step1_fitNULLGLMM.R \
      --plinkFile=${plink_grm_snps} \
      --phenoFile="${phenoFile}" \
      --phenoCol="PHE" \
      --traitType=binary       \
      --sampleIDColinphenoFile=IID \
      --outputPrefix="step1_${params.pheno_col.replaceAll(/\s/,'_').replaceAll(/\(|\)/, '')}_out" \
      --outputPrefix_varRatio="step1_${params.pheno_col.replaceAll(/\s/,'_').replaceAll(/\(|\)/, '')}" \
      --nThreads=${task.cpus} ${params.saige_step1_extra_flags}
    """
  }

}

if (params.trait_type == 'quantitative'){
  process gwas_1_fit_null_glmm_qt {
    tag "$plink_grm_snps"
    publishDir "${params.outdir}/gwas_1_fit_null_glmm", mode: 'copy'

    input:
    set val(grm_plink_input), file(bed), file(bim), file(fam) from plinkCh
    each file(phenoFile) from phenoCh

    output:
    file "*" into fit_null_glmm_results
    file ("step1_${params.pheno_col.replaceAll(/\s/,'_').replaceAll(/\(|\)/, '')}_out.rda") into rdaCh
    file ("step1_${params.pheno_col.replaceAll(/\s/,'_').replaceAll(/\(|\)/, '')}.varianceRatio.txt") into varianceRatioCh

    script:
    """
    step1_fitNULLGLMM.R \
      --plinkFile=${grm_plink_input} \
      --phenoFile="${phenoFile}" \
      --phenoCol="PHE" \
      --traitType=quantitative       \
	    --invNormalize=TRUE	\
      --sampleIDColinphenoFile=IID \
      --outputPrefix="step1_${params.pheno_col.replaceAll(/\s/,'_').replaceAll(/\(|\)/, '')}_out" \
      --outputPrefix_varRatio="step1_${params.pheno_col.replaceAll(/\s/,'_').replaceAll(/\(|\)/, '')}" \
      --nThreads=${task.cpus} ${params.saige_step1_extra_flags}
    """
  }

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
    --SAIGEOutputFile="${params.pheno_col.replaceAll(/\s/,'_').replaceAll(/\(|\)/, '')}.${name}.SAIGE.gwas.txt" \
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
  set file("*top_n.csv"), file("*${params.output_tag}.csv") into (ch_ldsc_input, ch_report_input)

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
/*--------------------------------------------------
  LDSC - Genetic correlation and heritability
---------------------------------------------------*/
if (params.post_analysis == 'heritability' || params.post_analysis == 'genetic_correlation_h2'){
  process prepare_files_ldsc {
    tag "preparation_files"
    publishDir "${params.outdir}/ldsc_inputs/", mode: 'copy'

    input:
    set file("*top_n.csv"), file(summary_stats) from ch_ldsc_input

    output:
    file("${params.output_tag}_transformed_gwas_stats.txt") into ch_ldsc_input2

    script:

    """
    convert_output.R \
      --gwas_stats "$summary_stats" \
      --output_tag ${params.output_tag}
    """
  }
  process munge_saige_output {
    tag "munge_saige_output"
    publishDir "${params.outdir}/ldsc_inputs/", mode: 'copy'

    input:
    file(saige_summary_stats) from ch_ldsc_input2
    file(hapmap3_snplist) from ch_hapmap3_snplist

    output:
    file("${params.output_tag}_ldsc.sumstats.gz") into ch_saige_ldsc

    script:

    """
    munge_sumstats.py --sumstats $saige_summary_stats \
                      --out "${params.output_tag}_ldsc" \
                      --merge-alleles $hapmap3_snplist \
                      --a1 Allele1 \
                      --a2 Allele2 \
                      --signed-sumstats Tstat,0 \
                      --p p.value \
                      --snp SNPID \
                      --info inputationInfo
    """
  }
}

if (params.post_analysis == 'heritability'){

  process heritability {
    tag "heritability"
    publishDir "${params.outdir}/heritability/", mode: 'copy'

    input:
    file(saige_output) from ch_saige_ldsc
    file(ld_scores_tar_bz2) from ch_ld_scores_tar_bz2

    output:
    file("${params.output_tag}_h2.log") into ch_ldsc_report_input

    script:
    """
    tar -xvjf ${ld_scores_tar_bz2}

    ldsc.py \
      --h2 $saige_output \
      --ref-ld-chr ${ld_scores_tar_bz2.simpleName}/ \
      --w-ld-chr ${ld_scores_tar_bz2.simpleName}/ \
      --out ${params.output_tag}_h2
    """
  }
}

if (params.post_analysis == 'genetic_correlation_h2' && params.gwas_summary){
  process prepare_gwas_summary_ldsc {
    tag "preparation_gwas_summary_ldsc"
    publishDir "${params.outdir}/ldsc_inputs/", mode: 'copy'

    input:
    val(gwas_summary_file) from ch_gwas_summary

    output:
    file("${params.external_gwas_tag}_transformed_gwas_stats.txt") into ch_gwas_summary_ldsc

    script:

    """
    convert_output.R \
      --gwas_stats "$gwas_summary_file" \
      --output_tag "${params.external_gwas_tag}"
    """
  }
  //* Munge gwas stats

  process munge_gwas_summary {
    tag "munge_gwas_summary"
    publishDir "${params.outdir}/ldsc_inputs/", mode: 'copy'

    input:
    file(summary_stats) from ch_gwas_summary_ldsc
    file(hapmap3_snplist) from ch_hapmap3_snplist

    output:
    file("${params.external_gwas_tag}_gwas_summary.sumstats.gz") into ch_gwas_summary_ldsc2

    script:

    """
    munge_sumstats.py \
          --sumstats "$summary_stats" \
          --out "${params.external_gwas_tag}_gwas_summary" \
          --merge-alleles $hapmap3_snplist
    """
  }

  //* Run genetic correlation
  process genetic_correlation_h2 {
    tag "genetic_correlation_h2"
    publishDir "${params.outdir}/genetic_correlation/", mode: 'copy'

    input:
    file(gwas_summary_ldsc) from ch_gwas_summary_ldsc2
    file(saige_ldsc) from ch_saige_ldsc
    file(ld_scores_tar_bz2) from ch_ld_scores_tar_bz2

    output:
    file("${params.output_tag}_genetic_correlation.log") into ch_ldsc_report_input

    script:

    """
    tar -xvjf ${ld_scores_tar_bz2}

    ldsc.py \
          --rg $saige_ldsc,$gwas_summary_ldsc \
          --ref-ld-chr ${ld_scores_tar_bz2.simpleName}/ \
          --w-ld-chr ${ld_scores_tar_bz2.simpleName}/ \
          --out ${params.output_tag}_genetic_correlation \
          --no-intercept
    """
  }

}

 //* gwas catalogue

if (params.post_analysis == 'genetic_correlation_h2' && params.gwas_cat_study_id){

  gwas_catalogue_ftp_ch = Channel.fromPath(params.gwas_catalogue_ftp, checkIfExists: true)
    .ifEmpty { exit 1, "GWAS catalogue ftp locations not found: ${params.gwas_catalogue_ftp}" }
    .splitCsv(header: true)
    .map { row -> tuple(row.study_accession, row.ftp_link_harmonised_summary_stats) }
    .filter{ it[0] == params.gwas_cat_study_id}
    .ifEmpty { exit 1, "The GWAS study accession number you provided does not come as a harmonized dataset that can be used as a base cohort "}
    .flatten()
    .last()

  process download_gwas_catalogue {
    label "high_memory"
    publishDir "${params.outdir}/GWAS_cat", mode: "copy"
    
    input:
    val(ftp_link) from gwas_catalogue_ftp_ch
    
    output:
    file("*.h.tsv*") into downloaded_gwas_catalogue_ch
    
    script:
    """
    wget ${ftp_link}
    """
  }

  process transform_gwas_catalogue {
    label "high_memory"
    publishDir "${params.outdir}/GWAS_cat", mode: "copy"
    
    input:
    file gwas_catalogue_base from downloaded_gwas_catalogue_ch
    
    output:
    file("${params.gwas_cat_study_id}.data") into transformed_base_ch
    
    script:
    """
    transform_gwas_catalogue.R --input_gwas_cat "${gwas_catalogue_base}" \
                               --outprefix "${params.gwas_cat_study_id}"
    """
    }

  
  //* Munge gwas cat stats

  process munge_gwas_cat_summary {
    tag "munge_gwas_summary"
    publishDir "${params.outdir}/ldsc_inputs/", mode: 'copy'

    input:
    file(summary_stats) from transformed_base_ch
    file(hapmap3_snplist) from ch_hapmap3_snplist

    output:
    file("${params.gwas_cat_study_id}_gwas_summary.sumstats.gz") into ch_gwas_summary_ldsc2

    script:

    """
    munge_sumstats.py \
          --sumstats "$summary_stats" \
          --out "${params.gwas_cat_study_id}_gwas_summary" \
          --merge-alleles $hapmap3_snplist \
          --signed-sumstats BETA,0 \
          --N ${params.gwas_cat_study_size}
    """
  }

  //* Run genetic correlation
  process genetic_correlation_h2_gwas_cat {
    tag "genetic_correlation_h2"
    publishDir "${params.outdir}/genetic_correlation/", mode: 'copy'

    input:
    file(gwas_summary_ldsc) from ch_gwas_summary_ldsc2
    file(saige_ldsc) from ch_saige_ldsc
    file(ld_scores_tar_bz2) from ch_ld_scores_tar_bz2

    output:
    file("${params.output_tag}_genetic_correlation.log") into ch_ldsc_report_input

    script:
    """
    tar -xvjf ${ld_scores_tar_bz2}

    ldsc.py \
          --rg $saige_ldsc,$gwas_summary_ldsc \
          --ref-ld-chr ${ld_scores_tar_bz2.simpleName}/ \
          --w-ld-chr ${ld_scores_tar_bz2.simpleName}/ \
          --out ${params.output_tag}_genetic_correlation \
          --no-intercept
    """
  }

}

if(params.post_analysis){
  process create_report_ldsc {
    tag "report"
    publishDir "${params.outdir}/MultiQC/", mode: 'copy'

    input:
    file(saige_output) from ch_report_input.collect()
    file(gwas_cat) from ch_gwas_cat
    file(ldsc_log) from ch_ldsc_report_input

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
    Rscript -e "rmarkdown::render('gwas_report.Rmd', params = list(manhattan='${params.output_tag}_manhattan.png',qqplot='${params.output_tag}_qqplot_ci.png', gwascat='gwascat_subset.csv', saige_results='saige_results_top_n.csv', trait_type='${params.trait_type}', ldsc_log='$ldsc_log'))"
    mv gwas_report.html multiqc_report.html

    # Generates the ipynb
    jupytext --to ipynb gwas_report.Rmd
    """
  }
}
if(!params.post_analysis){
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
}
