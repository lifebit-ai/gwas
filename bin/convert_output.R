#!/usr/bin/env Rscript



####################
# Import libraries #
####################

suppressPackageStartupMessages({
library(optparse)
library(data.table)
library(tidyverse)
    })

options(warn=-1)

##########################################################
# Parse arguments                                        
##########################################################

option_list = list(
  make_option(c("--gwas_stats"), action="store", default='saige_output.csv', type='character',
              help="String containing input GWAS summary statistics."),
  make_option(c("--output_tag"), action="store", default='CB', type='character',
              help="String containing the prefix to be used in the output files")
)

args = parse_args(OptionParser(option_list=option_list))

gwas_stats                = args$gwas_stats
outprefix                 = args$output_tag

data = fread(gwas_stats) %>%
       write.table(paste0(outprefix, '_transformed_gwas_stats.txt'), sep=' ', row.names=FALSE, quote=FALSE)