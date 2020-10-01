#!/usr/bin/env Rscript

###########################
# Import libraries
###########################
suppressPackageStartupMessages({
library(data.table)
library(tidyr)
library(dplyr)
library(optparse)
    })

###########################
# Functions
###########################

build_phenoFiles = function(row, data, design, out_path){
  control = design$control[row]
  case = design$case[row]
  # subsets the control case samples
  design_matrix = data[data$PHE %in% c(control, case),]
  #Apply control (0) vs case (1) convention
  design_matrix$PHE = case_when(design_matrix$PHE == control ~ 0,
                                design_matrix$PHE == case ~ 1,
                                TRUE ~ 0
  )
  # Writes it into file and returns the dataframe in case future steps are added
  write.table(design_matrix, paste0(out_path,'design_matrix_control_', control,'_case_',case,'.phe'), sep='\t',  quote=FALSE, row.names=FALSE)
  return(design_matrix)
}
###########################
# CLI & Arguments parsing
###########################

option_list = list(
  make_option(c("--input_file"), action="store", default='data/sample_multilevel.phe', type='character',
              help="String containing input phenoFile."),
  make_option(c("--case_group"), action="store", default='None', type='character',
              help="String containing case group from the desired phenotypic column."),
  make_option(c("--outprefix"), action="store", default='', type='character',
              help="String containing output phenoFile prefixes."),
  make_option(c("--outdir"), action="store", default='.', type='character',
              help="String containing output phenoFile dir (ie. results).")
)

args = parse_args(OptionParser(option_list=option_list))

# Args to variables
input_file              = args$input_file
case_group              = args$case_group
outprefix               = paste0(args$outprefix, "_")
outdir                  = sub("/$","",args$outdir)

system(paste0("mkdir -p ", outdir), intern=T)

out_path = paste0(outdir, "/", outprefix)

###########################
# Data import & creation of files
###########################

data = fread(input_file)

# In case group not specified it will run all combinations, 
# otherwise it will subset for the case group and generate only those files

if (case_group != 'None'){
case_group = as.integer(case_group)
design = crossing(control = data$PHE, case = data$PHE) %>%
            filter((case == case_group) & !(case == control))
design_list = sapply(1:nrow(design), function(x) build_phenoFiles(x, data, design, out_path), simplify=FALSE)

data$PHE = case_when(data$PHE == case_group ~ 1,
                     TRUE ~ 0)
write.table(data, paste0(out_path,'design_matrix_control_all','_case_',case_group,'.phe'), sep='\t',  quote=FALSE, row.names=FALSE)
  

}

if (case_group == 'None'){
design = crossing(control = data$PHE, case = data$PHE) %>%
            filter(!(case == control))
# Removes inverse duplicates (ie. 1-0, 0-1)
design = design[!duplicated(apply(design,1,function(x) paste(sort(x),collapse=''))), ]
design_list = sapply(1:nrow(design), function(x) build_phenoFiles(x, data, design, out_path), simplify=FALSE)
}