#!/usr/bin/env Rscript

############################## ARGUMENTS SECTION #############################
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no all arguments passed or help needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) | (length(args) == 1) ) {
  cat("
      The helper R Script subset_gwascat.R

      Mandatory arguments:
          --saige_output=path       - The path to the SAIGE output file.
                                      NOTE:
                                      If you have performed the analysis per chromosome,
                                      concatenate outputs across all chromosomes first.

         --gwas_cat=path            - The path to the NHGRI-EBI Catalog of published genome-wide association studies
                                      provided as a .csv file.
                                      NOTE:
                                      This can be retrieved using the bioconductor R package 'gwascat':
                                      gwascat <- as.data.frame(gwascat::makeCurrentGwascat())
                                      data.table::fwrite(gwascat, 'gwascat.csv', sep = ',')

         --help                     - you are reading it

     Usage:

          The typical command for running the script is as follows:

          ./subset_gwascat.R --saige_output='saige_results.csv' --gwas_cat='gwascat.csv'

     Output:

      Returns the GWAS Catalogue information as a .csv file of the top N ranked sites based on the SAIGE p.values column.

      \n")

  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")

argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL
rm(argsL)

############################## LIBRARIES SECTION #############################

suppressWarnings(suppressMessages(library(data.table)))

# ######################### VARIABLES REASSIGNMENT SECTION ###############################

# Facilitates testing and protects from wh-spaces, irregular chars

# required
saige_output     <- as.character(args$saige_output)
gwas_cat         <- as.character(args$gwas_cat)

cat("\n")
cat("ARGUMENTS SUMMARY")
cat("\n")
cat("saige_output     : ", saige_output     ,"\n",sep="")
cat("gwas_cat         : ", gwas_cat         ,"\n",sep="")

# ############################### SCRIPT SECTION ###############################

a <- as.data.frame(data.table::fread(saige_output))
to_keep <- a[["SNPID"]]
b <- as.data.frame(data.table::fread(gwas_cat))
c <- b[ b$SNPS %in% to_keep, ]
data.table::fwrite(c, file = "gwascat_subset.csv", sep = ",")
