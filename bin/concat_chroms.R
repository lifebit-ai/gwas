#!/usr/bin/env Rscript

############################## ARGUMENTS SECTION #############################
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no all arguments passed or help needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) | (length(args) == 1) ) {
  cat("
      The helper R Script concat_chroms.R

      Mandatory arguments:
        --output_tag=value         - A string in single quotes used as the identifier of the analysis
                                     in the output files name. Don't use whitespaces or irregular characters.
                                     The name will be converted to snakecase (eg. snake_case)

         --help                    - you are reading it

      Optional arguments:
        --filename_pattern=value  - The pattern of the SAIGE generated individual results files,
                                    used in list.files() for finding the files to concatenate.
                                    Default: '.SAIGE.gwas.txt'

        --saige_output_name=path   - The desired output identifier of the concatenated SAIGE output file with all the individual regions.
                                     Default: 'saige_results'
                                     NOTE:
                                     The name will be converted to snakecase (eg. snake_case)
                                     The final file name will be created by combining the output_tag value and this optional argument.
                                     eg.
                                     output_tag = 'breast_cancer_EAS_cohort'
                                     saige_output_name = 'saige_results',
                                     then the final file name will be:
                                     'saige_results_breast_cancer_EAS_cohort.csv'

        --top_n_sites=int          - The top N sites from the SAIGE results to be displayed in the report.
                                     The ranking is by ascending p-value.
                                     Default: 200
                                     NOTE:
                                     The maximum limit for this option is set by the parameter max_top_n_sites

        --max_top_n_sites=int      - The maximum allowed top N sites from the SAIGE results to be displayed in the report.
                                     Default: 1000
                                     NOTE:
                                     This is the upper bound for the values that top_n_sites can take.

     Usage:

          The typical command for running the script is as follows:

          ./concat_chroms.R --output_tag='summary_statistics'

     Output:

      Returns one .csv file {saige_output_name}_{output_tag}.csv, with the concatenated results
      across all tested genomic regions in the GWAS.

      See ./concat_chroms.R --help for more details.

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
if(is.null(args$filename_pattern)) {args$filename_pattern = ".SAIGE.gwas.txt"} else {args$filename_pattern=as.character(args$filename_pattern)}
if(is.null(args$saige_output_name)) {args$saige_output_name = "saige_results"} else {args$saige_output_name=as.character(args$saige_output_name)}
if(is.null(args$top_n_sites)) {args$top_n_sites = 200} else {args$top_n_sites=as.numeric(args$top_n_sites)}
if(is.null(args$max_top_n_sites)) {args$max_top_n_sites = 1000} else {args$max_top_n_sites=as.numeric(args$max_top_n_sites)}

############################## LIBRARIES SECTION #############################

suppressWarnings(suppressMessages(library(plyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(snakecase)))

# ######################### VARIABLES REASSIGNMENT SECTION ###############################

# Facilitates testing and protects from wh-spaces, irregular chars

# required
output_tag         <- as.character(args$output_tag)
# optional
filename_pattern   <- as.character(args$filename_pattern)
saige_output_name  <- snakecase::to_snake_case(as.character(args$saige_output_name))
max_top_n_sites    <- as.numeric(args$max_top_n_sites)
top_n_sites        <- ifelse(as.numeric(args$top_n_sites) > max_top_n_sites, max_top_n_sites, as.numeric(args$top_n_sites))

cat("\n")
cat("ARGUMENTS SUMMARY")
cat("\n")
cat("output_tag         : ", output_tag         ,"\n",sep="")
cat("filename_pattern   : ", filename_pattern   ,"\n",sep="")
cat("saige_output_name  : ", saige_output_name  ,"\n",sep="")
cat("top_n_sites        : ", top_n_sites        ,"\n",sep="")

# ############################### SCRIPT SECTION ###############################

paths <- list.files(".", pattern = filename_pattern, full.names = TRUE)
list_of_dfs <- lapply(paths,data.table::fread)

if (length(paths) == 1){
  saige_results <- data.frame(list_of_dfs[[1]])
} else {
  saige_results <- plyr::rbind.fill(list_of_dfs)
}

saige_results['p.value'][is.na(saige_results['p.value'])] = 0.999999
saige_results_sorted_topN <- saige_results[order(saige_results[['p.value']]),][1:top_n_sites,]
data.table::fwrite(saige_results, paste0(saige_output_name, "_", output_tag, ".csv"), sep = ",")
data.table::fwrite(saige_results_sorted_topN, "saige_results_top_n.csv", sep = ",")
