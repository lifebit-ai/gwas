#!/usr/bin/env Rscript

###########################
# Import libraries
###########################
suppressPackageStartupMessages({
library(data.table)
library(tidyverse)
library(optparse)
library(yaml)
    })

###########################
# CLI & Arguments parsing
###########################

option_list = list(
  make_option(c("--pheno_data"), action="store", default='None', type='character',
              help="String containing path/URL to input pheno data use as template for sampling."),
  make_option(c("--config_file"), action="store", default='None', type='character',
              help="String containing path/URL to input pheno data use as template for sampling."),
  make_option(c("--outprefix"), action="store", default='sim_1', type='character',
              help="String containing output simulated prefixes.")
)

args = parse_args(OptionParser(option_list=option_list))

# Args to variables
pheno_data              = args$pheno_data
config_file             = args$config_file
outprefix               = args$outprefix


###########################
# Load config
###########################
config = read_yaml(config_file)$params

###########################
# Resample from a template file
###########################

if (pheno_data != 'None') {
    pheno_df = fread(pheno_data)
    sym_data = pheno_df %>% sample_n(config$n_samples, replacement=T)
    write.csv(sym_data, paste0(outprefix,'_pheno_data.csv'),  quote=FALSE, row.names=FALSE)
}

##########################################
# Generate new samples according to config
##########################################
simulate_pheno = function(config, col_names){
    #create column names for instances and arrays
    combinations = expand.grid(col_names, 
                                0:(config[['col_params']][[col_names]][['n_instances']]-1),
                                0:(config[['col_params']][[col_names]][['n_arrays']]-1))
    combinations = paste0(combinations$Var1,
                            '-', 
                            combinations$Var2,
                            '.', 
                            combinations$Var3)
    # Categorical
    if (config[['col_params']][[col_names]][['type']] == 'Categorical'){
        # Sample categorical data from the values in the config
        sym_cols = sapply(combinations,
                          function(x) sample(config[['col_params']][[col_names]][['values']], 
                                 config[['n_samples']],
                                 replace=T,
                                 prob=config[['col_params']][[col_names]][['prob']]))
        # Get all NA for n random columns
        col_to_na = sample(colnames(sym_cols), config[['col_params']][[col_names]][['n_missing_col']])
        sym_cols[sym_cols == 'NA'] = NA
        sym_cols[, col_to_na] = NA
        return(sym_cols %>% as.tibble())
        
    }
    if (config[['col_params']][[col_names]][['type']] == 'Integer'){
        #Sample from a normal distribution
        if (config[['col_params']][[col_names]][['distribution']] == 'normal'){
            sym_cols = sapply(combinations,
                              function(x) rnorm(config[['n_samples']],
                                    config[['col_params']][[col_names]][['mean']],
                                    config[['col_params']][[col_names]][['sd']]))
        }
        #Sample from a uniform distribution
        if (config[['col_params']][[col_names]][['distribution']] == 'uniform'){
            sym_cols = sapply(combinations,
                              function(x) runif(config[['n_samples']],
                                    config[['col_params']][[col_names]][['min']],
                                    config[['col_params']][[col_names]][['max']]))
        }
        #If only want positive values, use absolute values for negatives
        if (config[['col_params']][[col_names]][['positive_only']]) {
            sym_cols = abs(sym_cols)
        }
        #Transform columns in NA
        col_to_na = sample(colnames(sym_cols), config[['col_params']][[col_names]][['n_missing_col']])
        sym_cols[, col_to_na] = NA
        mask = sapply(combinations,
        function(x) sample(c(TRUE, FALSE), config[['n_samples']], replace=T, prob=c(config[['col_params']][[col_names]][['frac_missing']], 1-config[['col_params']][[col_names]][['frac_missing']]))
        )
        sym_cols[mask] = NA
        return(floor(sym_cols) %>% as.tibble())
    }
    if (config[['col_params']][[col_names]][['type']] == 'Continuous'){
        #Repeat process with continuous
        if (config[['col_params']][[col_names]][['distribution']] == 'normal'){
            sym_cols = sapply(combinations,
                              function(x) rnorm(config[['n_samples']],
                                    config[['col_params']][[col_names]][['mean']],
                                    config[['col_params']][[col_names]][['sd']]))
        }
        if (config[['col_params']][[col_names]][['distribution']] == 'uniform'){
            sym_cols = sapply(combinations,
                              function(x) runif(config[['n_samples']],
                                    config[['col_params']][[col_names]][['min']],
                                    config[['col_params']][[col_names]][['max']]))
        }
        
        if (config[['col_params']][[col_names]][['positive_only']]) {
            sym_cols = abs(sym_cols)
        }

        col_to_na = sample(colnames(sym_cols), config[['col_params']][[col_names]][['n_missing_col']])
        sym_cols[, col_to_na] = NA
        mask = sapply(combinations,
        function(x) sample(c(TRUE, FALSE), config[['n_samples']], replace=T, prob=c(config[['col_params']][[col_names]][['frac_missing']], 1-config[['col_params']][[col_names]][['frac_missing']]))
        )
        sym_cols[mask] = NA
        return(sym_cols %>% as.tibble())        
    }
    if (config[['col_params']][[col_names]][['type']] == 'Date'){
        #For dates, generate sequence of dates between starting and end date and sample the vector
        sym_cols = sapply(combinations,
                          function(x) sample(seq(as.Date(config[['col_params']][[col_names]][['starting_date']]),
                                                 as.Date(config[['col_params']][[col_names]][['end_date']]), 
                                                 by="day"), 
                                size=config[['n_samples']], 
                                replace=T), simplify=FALSE) %>% as.tibble()
        #Introduce NAs
        col_to_na = sample(colnames(sym_cols), config[['col_params']][[col_names]][['n_missing_col']])
        sym_cols[, col_to_na] = NA
        mask = sapply(combinations,
        function(x) sample(c(TRUE, FALSE), config[['n_samples']], replace=T, prob=c(config[['col_params']][[col_names]][['frac_missing']], 1-config[['col_params']][[col_names]][['frac_missing']]))
        )
        sym_cols[mask] = NA
        return(sym_cols)  
    }
}
if (pheno_data == 'None'){
    #Run everything and add fake platekey IDs
    sym_data = sapply(names(config[['col_params']]), function(x) simulate_pheno(config, x), simplify=FALSE) %>% bind_cols()
    sym_data['Platekey_in_aggregate_VCF-0.0'] = 1:(config[['n_samples']])
    sym_data['i'] = 1:(config[['n_samples']])
    write.csv(sym_data, paste0(outprefix, '_pheno_data.csv'), quote=FALSE, row.names=FALSE)

}
# Create synthetic metadata file
value_type = sapply(names(config[['col_params']]), function(x) config[['col_params']][[x]][['type']]) %>% as.vector()

metadata = data.frame('id'=1:length(names(config[['col_params']])), valueType=value_type, name=names(config[['col_params']])) 
platekeys_ids = data.frame(id=c(length(names(config[['col_params']]))+1, length(names(config[['col_params']]))+2),
                           valueType=c("Categorical", "Categorical"),
                           name=c("Platekey_in_aggregate_VCF", "i"))

rbind(metadata, platekeys_ids) %>% write.csv(paste0(outprefix, '_pheno_metadata.csv'), quote=FALSE, row.names=FALSE)