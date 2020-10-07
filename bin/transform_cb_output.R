#!/usr/bin/env Rscript



####################
# Import libraries #
####################

suppressPackageStartupMessages({
library(optparse)
library(data.table)
library(tidyverse)
library(jsonlite)
    })

options(warn=-1)

##########################################################
# Parse arguments                                        
##########################################################

option_list = list(
  make_option(c("--input_cb_data"), action="store", default='data/cohort_data_phenos_v4.csv', type='character',
              help="String containing input Cohort Browser data."),
  make_option(c("--input_meta_data"), action="store", default='assets/Metadata phenotypes - Mapping file.csv', type='character',
              help="String containing input metadata for columns in Cohort Browser output."),
  make_option(c("--phenoCol"), action="store", default='None', type='character',
              help="String representing phenotype that will be used for GWAS comparison(s)."),
  make_option(c("--continuous_var_transformation"), action="store", default='log', type='character',
              help="String representing the type of transformation desired for integer and continuous input data"),
  make_option(c("--continuous_var_aggregation"), action="store", default='mean', type='character',
              help="String representing the type of aggregation desired for input data"),
  make_option(c("--outdir"), action="store", default='.', type='character',
              help="String containing the output directory"),
  make_option(c("--outprefix"), action="store", default='CB', type='character',
              help="String containing the prefix to be used in the output files")
 
)

args = parse_args(OptionParser(option_list=option_list))

input_cb_data                 = args$input_cb_data
input_meta_data               = args$input_meta_data
phenoCol                      = args$phenoCol
aggregation                   = args$continuous_var_aggregation
transformation                = args$continuous_var_transformation
outprefix                     = paste0(args$outprefix, "_")
outdir                        = sub("/$","",args$outdir)

system(paste0("mkdir -p ", outdir), intern=T)

out_path = paste0(outdir, "/", outprefix)



if (!(aggregation %in% c('mean', 'max', 'min', 'median'))){
    stop('Selected transformation for continuous variables not supported.')
}

##########################################################
# Import cohort browser (cb) data and contrast phenotype 
##########################################################

cb_data = fread(input_cb_data) %>% as.tibble

# Remove columns full of NAs (empty string in CSV)
cb_data = cb_data %>% select_if(~!all(is.na(.)))

##################################################
# Keep only participants for which we have a VCF #
##################################################

cb_data = cb_data %>% filter(!`Platekey in aggregate VCF-0.0`== "")

################################
# Re-encode cb_data phenotypes #
################################

# Trim suffix that denotes multiple entries of columns and replace spaces by "-"
#colnames(cb_data) = colnames(cb_data) %>% str_replace("-[^-]+$", "")
colnames(cb_data) = colnames(cb_data) %>% 
        str_replace_all(" ", "_") %>% 
        str_replace_all("\\(|\\)","") %>%
        str_to_lower()

# Use phenotype metadata (data dictionary) to determine the type of each phenotype -> This will be given by CB
pheno_dictionary = fread(input_meta_data) %>%
        as.tibble # Change by metadata input var
pheno_dictionary$'Field Name' = str_replace_all(pheno_dictionary$'Field Name'," ", "_") %>%
        str_replace_all('\\(|\\)',"") %>% 
        str_to_lower()

#Compress multiple measures into a single measurement


encode_pheno_values = function(column, data, pheno_dictionary, transformation, aggregation){
    
    #Clean column name
    pheno_cols = data[, str_detect(colnames(data), column)]

    pheno_dtype = filter(pheno_dictionary, str_detect(pheno_dictionary$`Field Name`, column)) %>% 
            pull(`FieldID Type`)
    ################################
    # Individual ID                #
    ################################
    if (column == "individual_id"){

        pheno_cols = data[[column]]
        return(as.vector(pheno_cols))
    }
    ################################
    # Categorical               #
    ################################
    if (str_detect(pheno_dtype, "Categorical") == TRUE){
        if (str_detect(column, 'platekey')){
            pheno_cols = pheno_cols[[1]] %>% as.vector
            return(pheno_cols)
        }
        # Fill the gaps and get list of unique values
        pheno_cols[pheno_cols == ''] = "UNKNOWN"
        pheno_values = pheno_cols %>% unlist() %>% sort() %>% unique()
        # Decide aggregation behaviour for samples with paired measures
        if (dim(pheno_cols)[2] > 1) {
            # Arbitrary : get the first column
            # Adds variable called query match that is specific for the column 
            pheno_cols = apply(pheno_cols, 1, function(x) x[1])
        }
        # Encode unique values and create mapping list
        encoding = as.list(1:length(pheno_values))
        names(encoding) = pheno_values
        # Store .json with encoding mappings, will be used later on.
        encoding_json = toJSON(encoding,keep_vec_names=TRUE)
        write(encoding_json, file = file.path(column, ".json", fsep = ""))
        #Use mapping list on aggregated columns to get
        encoded_col = lapply(pheno_cols, function(x) encoding[x]) %>% unlist() %>% as.vector
        return(encoded_col)
    }
    ################################
    # Year of Birth                #
    ################################   
    if ((str_detect(column,"birth") == TRUE)){
        # Transform year of birth into age
        current_year = format(Sys.time(), "%Y") %>% as.integer
        age = current_year - data[[column]] %>% as.vector
        return(age)
    }
    ################################
    # Integers and Continuous      #
    ################################ 
    if (str_detect(pheno_dtype, 'Integer|Continuous')){
        
        # pick transformation function - tried a case_when but it seems... 
        # ...I cannot make it give back functions
        if (aggregation == 'mean'){
            aggregation = function(x) mean(x)
        } else if (aggregation == 'median') {
            aggregation = function(x) median(x)
        } else if (aggregation == 'max') {
            aggregation = function(x) max(x)
        } else if (aggregation){
            aggregation = function(x) min(x)
        }

        #Apply aggregation & transformation
        ## Get unique sets of measurements
        if (dim(pheno_cols)[2] > 1){
            #Finds group of instances
            sets_measures = str_extract(colnames(pheno_cols), "-[:digit:]") %>% unique()
            ## Group by the same group of arrays
            ##Merge arrays per instances
            pheno_cols = sapply(sets_measures, function(value) apply(pheno_cols[, str_detect(colnames(pheno_cols), value)], 1, function(x) aggregation(x)))
            #Group by instances
            pheno_cols = apply(pheno_cols, 1, function(x) aggregation(x))
        }else{
            pheno_cols = lapply(pheno_cols, function(x) aggregation(x))
        }
        pheno_cols = pheno_cols %>% as.vector

        if (transformation == 'log'){
            pheno_cols = log(pheno_cols)
        }
        if (transformation == 'log2') {
            pheno_cols = log2(pheno_cols)
        } 
        if (transformation == 'zscore') {
            pheno_cols = (pheno_cols - mean(pheno_cols)) / sd(pheno_cols)
        }
        if (transformation == 'None'){
            pheno_cols = pheno_cols
        }

        return(pheno_cols)

    }
    ################################
    # Dates                        #
    ################################ 
    if (str_detect(pheno_dtype, 'Time|Date')){
        # Transform - turns it into a big integer
        # Fill empty gaps with current date
        pheno_cols[pheno_cols == ''] = format(Sys.time(), "%d/%m/%Y")
        ## Multiple array support
        if (dim(pheno_cols)[2] > 1) {
            # Turns the dates into a big integer
            pheno_cols = apply(pheno_cols, 1, function(x) format(as.Date(x, "%d/%m/%Y"), "%Y%m%d") %>% as.integer)
            # Aggregate - gets the first column - arbitrary
            pheno_cols = apply(pheno_cols, 1, function(x) x[1])
        }else{
            # If only one array, applies directly the transformation
            pheno_cols = lapply(pheno_cols, function(x) format(as.Date(x, "%d/%m/%Y"), "%Y%m%d") %>% as.integer) %>% as.vector
        }
        return(pheno_cols[[1]])
    }
    ################################
    # Free text              #
    ################################ 
    if (str_detect(pheno_dtype, 'Text')){
        ## Sets text to NA
        return(rep(NA, dim(pheno_cols)[1]))
    }

}

# Run across all columns
# encode_pheno_values('specimen_type', cb_data, pheno_dictionary, transformation)
columns_to_transform = colnames(cb_data) %>%
        str_replace("-[^-]+$", "") %>%
        unique
cb_data_transformed = sapply(columns_to_transform, function(x) encode_pheno_values(x, cb_data, pheno_dictionary, transformation, aggregation), simplify=FALSE) %>% as.data.frame

#####################
# Make final output #
#####################

#TODO: Add more covariates
column_to_PHE = phenoCol %>% str_replace('\\(|\\)',"") %>% 
        str_replace("-[^-]+$", "") %>% 
        str_replace(' ','_') %>% 
        str_to_lower
cb_data_transformed = as_tibble(cb_data_transformed)

##Build the .phe file format
cb_data_transformed$FID = cb_data_transformed[['platekey_in_aggregate_vcf']]
cb_data_transformed$IID = cb_data_transformed[['platekey_in_aggregate_vcf']]
cb_data_transformed$PAT = 0
cb_data_transformed$MAT = 0
# This should be provided either by default from the CB output or as an argument or calculated from the VCF data
cb_data_transformed$PHE = cb_data_transformed[[column_to_PHE]]

cb_data_transformed[['individual_id']] = NULL


##################################################
# Write .phe file                                #
##################################################
cb_data_transformed = cb_data_transformed %>% select(FID, IID, PAT, MAT, PHE, everything(), -!!as.symbol(column_to_PHE), -`platekey`, -`platekey_in_aggregate_vcf`)
### FID, IID this has to be the platekey metadata -> Agg VCF columns. 
write.table(cb_data_transformed, paste0(out_path,'.phe'), sep='\t',  quote=FALSE, row.names=FALSE)



