#!/usr/bin/env Rscript

###########################
# Import libraries
###########################
suppressPackageStartupMessages({
library(data.table)
library(tidyverse)
library(optparse)
    })

###########################
# Functions
###########################

###########################
# CLI & Arguments parsing
###########################

option_list = list(
  make_option(c("--input_file"), action="store", default='data/sample_multilevel.phe', type='character',
              help="String containing input Real Pheno file."),
  make_option(c("--ids_column"), action="store", default='None', type='character',
              help="String containing column key for IDs."),
  make_option(c("--outprefix"), action="store", default='', type='character',
              help="String containing output simulated prefixes.")
)

args = parse_args(OptionParser(option_list=option_list))

# Args to variables
input_file              = args$input_file
ids_column              = args$ids_column
outprefix               = paste0(args$outprefix, "_")


###########################
# Data import
###########################

data = fread(input_file) %>% as.data.frame()

###########################
# Transform ID columns
###########################

old_ID = data[[ids_column]]
new_ID = paste0(1:length(old_ID),'_', 1:length(old_ID))

conversion_table = data.frame(old=old_ID, new=new_ID)
colnames(conversion_table) = c(ids_column, "new")
data[ids_column] = conversion_table$new

###########################
# save conversion table for cohort data
###########################
conversion_table %>% write.csv(paste0(outprefix,'_IDs.csv'), row.names=FALSE, quote=FALSE)
###########################
# save GWAS data
###########################
mask = colnames(data)[str_detect(colnames(data) %>% str_to_lower(), 'diagnosis_codes|icd|hpo')]
data %>% select(-all_of(mask)) %>% write.csv(paste0(outprefix,'_gwas.csv'), quote=FALSE, row.names=FALSE)

###########################
# save PheWAS data
###########################
data %>% write.csv(paste0(outprefix, '_phewas.csv'), quote=FALSE, row.names=FALSE)
###########################
# save original data
###########################
data %>% write.csv(paste0(outprefix, '_original.csv'), quote=FALSE, row.names=FALSE)