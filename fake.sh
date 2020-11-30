Command executed:

  # Download, filter and convert (bcf or vcf.gz) -> vcf.gz
  bcftools view -q 0.005:minor sampleA_chr10.vcf.gz -Oz -o sampleA_chr10_filtered.vcf.gz
  bcftools index sampleA_chr10_filtered.vcf.gz
  
  # Create PLINK binary from vcf.gz
  plink2       --make-bed       --set-missing-var-ids @:#,\$r,\$a       --vcf sampleA_chr10_filtered.vcf.gz       --out sampleA_chr10_filtered       --vcf-half-call m       --double-id       --set-hh-missing       --new-id-max-allele-len 1000 missing       --output-chr  chrM
  
  #Filter missingness
  plink       --bfile sampleA_chr10_filtered       --pheno covid_1_design_matrix_control_all_case_1.phe       --pheno-name PHE       --allow-no-sex       --test-missing midp       --out sampleA_chr10       --1       --keep-allele-order              --output-chr chrM
  
  awk -v thresm=${params.thres_m} '$5 < thresm {print}'  sampleA_chr10.missing > sampleA_chr10.missing_FAIL 
  
  #Filter HWE
  plink       --bfile sampleA_chr10_filtered       --pheno covid_1_design_matrix_control_all_case_1.phe       --pheno-name PHE       --allow-no-sex       --hwe 1e-5 midp       --out sampleA_chr10.misHWEfiltered       --make-just-bim       --exclude sampleA_chr10.missing_FAIL       --1       --keep-allele-order              --output-chr chrM
  
  bcftools view sampleA_chr10_filtered.vcf.gz | awk -F '\t' 'NR==FNR{c[$1$4$6$5]++;next}; c[$1$2$4$5] > 0' sampleA_chr10.misHWEfiltered.bim - | bgzip > sampleA_chr10.filtered_temp.vcf.gz
  bcftools view -h sampleA_chr10_filtered.vcf.gz -Oz -o sampleA_chr10_filtered.header.vcf.gz
  cat sampleA_chr10_filtered.header.vcf.gz sampleA_chr10.filtered_temp.vcf.gz > sampleA_chr10.filtered_final.vcf.gz
  bcftools index sampleA_chr10.filtered_final.vcf.gz
