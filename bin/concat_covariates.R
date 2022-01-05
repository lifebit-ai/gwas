#!/usr/bin/env Rscript

############################## ARGUMENTS SECTION #############################
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no all arguments passed or help needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) | (length(args) == 1) ) {
  cat("
      The helper R Script concat_covariates.R

      Mandatory arguments:
        --pcs_file         - Path to file containing eigenvectors from PCA step (output of plink --pca command).

        --phenotype_file        - Path to file containing phenotype of interest and covariates. Supplied via --pheno_data to the pipeline.

        --output_file           - Path to resultant file with phenotype, covariates and PCs combined.

         --help                    - helpful documentation.


     Usage:

          The typical command for running the script is as follows:

          ./concat_covariates.R --pcs_file='pca_results.eigenvec' --phenotype_file='pheno_and_covariates.tsv' --output_file='full_covariates_pheno_files.tsv'

     Output:

      Returns a single .tsv file {saige_output_name}_{output_tag}.csv with covariates, phenotype of interest and PCs combined into one file.

      See ./concat_covariates.R --help for more details.

      \n")

  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")

argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL
rm(argsL)


## Give some value to optional arguments if not provided
# if(is.null(args$filename_pattern)) {args$filename_pattern = ".SAIGE.gwas.txt"} else {args$filename_pattern=as.character(args$filename_pattern)}
# if(is.null(args$saige_output_name)) {args$saige_output_name = "saige_results"} else {args$saige_output_name=as.character(args$saige_output_name)}
# if(is.null(args$top_n_sites)) {args$top_n_sites = 200} else {args$top_n_sites=as.numeric(args$top_n_sites)}
# if(is.null(args$max_top_n_sites)) {args$max_top_n_sites = 1000} else {args$max_top_n_sites=as.numeric(args$max_top_n_sites)}

############################## LIBRARIES SECTION #############################

suppressWarnings(suppressMessages(library(plyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(snakecase)))
suppressWarnings(suppressMessages(library(dplyr)))

# ######################### VARIABLES REASSIGNMENT SECTION ###############################

# Facilitates testing and protects from wh-spaces, irregular chars

pcs_file <- as.character(args$pcs_file)
phenotype_file <- as.character(args$phenotype_file)
output_file <- as.character(args$output_file)


cat("\n")
cat("ARGUMENTS SUMMARY")
cat("\n")
cat("pcs_file         : ", pcs_file         ,"\n",sep="")
cat("phenotype_file   : ", phenotype_file   ,"\n",sep="")
cat("output_file  : ", output_file  ,"\n",sep="")


# ############################### SCRIPT SECTION ###############################

principal_components = fread(pcs_file)
covariates_file = fread(phenotype_file)

output_df <- dplyr::left_join(principal_components, covariates_file, by='IID') %>% mutate(across(where(is.numeric), coalesce, -9))
output_df <- output_df[output_df$phenotype != -9]
data.table::fwrite(output_df, "covariates_with_PCs.tsv", sep = "\t")



