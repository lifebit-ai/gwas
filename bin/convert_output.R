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
 
)

args = parse_args(OptionParser(option_list=option_list))

gwas_stats                = args$gwas_stats

data = fread(gwas_stats) %>%
       write.table(paste0('transformed_', gwas_stats), sep=' ', row.names=FALSE, quotes=FALSE)