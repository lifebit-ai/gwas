#!/usr/bin/env Rscript



####################
# Import libraries #
####################

suppressPackageStartupMessages({
library(optparse)
library(data.table)
library(tidyverse)
library(jsonlite)
library(snakecase)
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
  make_option(c("--query_file"), action="store", default='None', type='character',
              help="String containing path/URL to query file."),  
  make_option(c("--phenoCol"), action="store", default='None', type='character',
              help="String representing phenotype that will be used for GWAS comparison(s)."),
  make_option(c("--continuous_var_transformation"), action="store", default='log', type='character',
              help="String representing the type of transformation desired for integer and continuous input data"),
  make_option(c("--continuous_var_aggregation"), action="store", default='mean', type='character',
              help="String representing the type of aggregation desired for input data"),
  make_option(c("--outdir"), action="store", default='.', type='character',
              help="String containing the output directory"),
  make_option(c("--output_tag"), action="store", default='CB', type='character',
              help="String containing the prefix to be used in the output files")
 
)

args = parse_args(OptionParser(option_list=option_list))

input_cb_data                 = args$input_cb_data
input_meta_data               = args$input_meta_data
query_file                    = args$query_file
phenoCol                      = args$phenoCol
aggregation                   = args$continuous_var_aggregation
transformation                = args$continuous_var_transformation
outprefix                     = paste0(args$output_tag, "_")
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


################################
# Re-encode cb_data phenotypes #
################################

# Trim suffix that denotes multiple entries of columns and replace spaces by "-"
#colnames(cb_data) = colnames(cb_data) %>% str_replace("-[^-]+$", "")
colnames(cb_data) = colnames(cb_data) %>% 
        str_replace_all(" ", "_") %>% 
        str_replace_all("\\(","") %>%
        str_replace_all("\\)","") %>%
        str_to_lower()

# Use phenotype metadata (data dictionary) to determine the type of each phenotype -> This will be given by CB
pheno_dictionary = fread(input_meta_data) %>%
        as.tibble # Change by metadata input var
pheno_dictionary$`name` = str_replace_all(pheno_dictionary$`name`," ", "_") %>%
        str_replace_all("\\(","") %>%
        str_replace_all("\\)","") %>% 
        str_to_lower()

##########################################################
# Read query file and prepare list of lists
##########################################################

if (query_file != 'None'){
    query_df = fromJSON(query_file, flatten=T)$search[, c('values','column.id')]
    query_df = left_join(query_df, pheno_dictionary, by = c('column.id' = 'id')) %>% select(values, name)
}
if (query_file == 'None'){
    query_df = 'None'
}
##################################################
# Keep only participants for which we have a VCF #
##################################################

cb_data = cb_data %>% filter(!`platekey_in_aggregate_vcf-0.0`== "") %>% select(-i)
#Compress multiple measures into a single measurement

##################################################
# Functions
##################################################

encode_pheno_values = function(column, data, pheno_dictionary, transformation, aggregation, query_df = 'None'){
    
    #Clean column name
    pheno_cols = data[, str_detect(colnames(data), column)]

    pheno_dtype = filter(pheno_dictionary, str_detect(pheno_dictionary$`name`, column)) %>% 
            pull(`valueType`)
    
    
    ################################
    # Individual ID                #
    ################################
    if (column == "individual_id|i"){

        pheno_cols = data[[column]]
        return(as.vector(pheno_cols))
    }
    ################################
    # Categorical                  #
    ################################
    if (str_detect(pheno_dtype, "Categorical") == TRUE){
        if (str_detect(column, 'platekey')){
            pheno_cols = pheno_cols[[1]] %>% as.vector()
            return(pheno_cols)
        }
        # Fill the gaps and get list of unique values
        pheno_cols[pheno_cols == ''] = "UNKNOWN"
        pheno_cols[pheno_cols == NA] = "UNKNOWN"
        pheno_values = pheno_cols %>% unlist() %>% sort() %>% unique()
        # Decide aggregation behaviour for samples with paired measures
        condition = dim(pheno_cols)[2]
        if (condition > 1 & query_df == 'None') {
            # Arbitrary : get the first column
            # Adds variable called query match that is specific for the column 
            pheno_cols = apply(pheno_cols, 1, function(x) x[1])

        }
        if (condition > 1 & query_df != 'None') {
            if (sum(str_detect(query_df$`name`, column)) > 0){
                # identify rows with queried values
                query_values = query_df[str_detect(query_df$name, column),]
                query_mask = apply(pheno_cols, 2, function(x) x %in% query_values$values)
                # get the values that are in the query
                values = sapply(1:dim(pheno_cols)[1], function(x) pheno_cols[x, query_mask[x,]])
                # get the first entry for the list of values queried for each row with queried values 
                # and get a list of the first value encounter in each row
                # There shouldn't be empty rows in this variables because the query was used to generate the cohort
                values = as.vector(unlist(apply(rbind(values), 2, function(x) unlist(x)[1])))
                # Substitute the original first column by the new column of first encountered values
                pheno_cols[[1]] = values
                #Select the first column of values
                pheno_cols = apply(pheno_cols, 1, function(x) x[1])
            }
            #Cannot make it flat because if query == 'None' it would break the pipeline
            if(sum(str_detect(query_df$`name`, column)) == 0){
                #when the column is not on the query file, just apply the standard filter -> take the first column
                pheno_cols = apply(pheno_cols, 1, function(x) x[1])
            }
            
        }    
        # Encode unique values and create mapping list
        encoding = as.list(1:length(pheno_values))
        names(encoding) = pheno_values
        # Store .json & csv with encoding mappings, will be used later on.
        #csv
        encoding_csv = data.frame(code = 1:length(pheno_values),
                                  original = pheno_values)
        write.csv(encoding_csv, file.path(column, ".csv", fsep = ""), quote=TRUE, row.names=FALSE)
        #json
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
        age = current_year - pheno_cols %>% as.vector()
        return(age)
    }
    ################################
    # Integers and Continuous      #
    ################################ 
    if (str_detect(pheno_dtype, 'Integer|Continuous')){
        
        # pick transformation function - tried a case_when but it seems... 
        # ...I cannot make it give back functions
        if (aggregation == 'mean'){
            aggregation_fun = function(x) mean(x, na.rm=TRUE)
        }
        if (aggregation == 'median') {
            aggregation_fun = function(x) median(x, na.rm=TRUE)
        }
        if (aggregation == 'max') {
            aggregation_fun = function(x) max(x, na.rm=TRUE)
        }
        if (aggregation == 'min'){
            aggregation_fun = function(x) min(x, na.rm=TRUE)
        }

        #Apply aggregation & transformation
        ## Get unique sets of measurements
        if (dim(pheno_cols)[2] > 1){
            #Finds group of instances
            sets_measures = str_extract(colnames(pheno_cols), "-[:digit:]") %>% unique()
            ## Group by the same group of arrays
            ##Merge arrays per instances
            pheno_cols = sapply(sets_measures, function(value) apply(pheno_cols[, str_detect(colnames(pheno_cols), value)], 1, function(x) aggregation_fun(x)))
            #Group by instances
            pheno_cols = apply(pheno_cols, 1, function(x) aggregation_fun(x))
        }
        if (is.vector(pheno_cols) && length(dim(pheno_cols)) == 1) {
            pheno_cols = lapply(pheno_cols, function(x) aggregation_fun(x))
        }
        pheno_cols = pheno_cols[[1]] %>% as.vector()
        
        if (str_detect(column, 'pc[0-9]')){
            transformation = 'None'
        }

        if (transformation == 'log'){
            pheno_cols = log(pheno_cols)
        }
        if (transformation == 'log10'){
            pheno_cols = log(pheno_cols, 10)
        }
        if (transformation == 'log2') {
            pheno_cols = log2(pheno_cols)
        } 
        #Deals with sd of vectors with only 1 non-NA value
        if (transformation == 'zscore' & sum(!is.na(pheno_cols)) < 2 ) {
            pheno_cols = pheno_cols
        }
        if (transformation == 'zscore' & sum(!is.na(pheno_cols)) >= 2 ) {
            pheno_cols = (pheno_cols - mean(pheno_cols, na.rm=TRUE)) / sd(pheno_cols, na.rm=TRUE)
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
        }
        if (is.vector(pheno_cols) && length(dim(pheno_cols)) == 1) {
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
        pheno_cols = rep(NA, dim(pheno_cols)[1])
        return(pheno_cols)
    }
}
#####################
# Run encoding
#####################

# Run across all columns
# encode_pheno_values('specimen_type', cb_data, pheno_dictionary, transformation)
columns_to_transform = colnames(cb_data) %>%
        str_replace("-[^-]+$", "") %>%
        unique
cb_data_transformed = sapply(columns_to_transform, function(x) encode_pheno_values(x, cb_data, pheno_dictionary, transformation, aggregation, query_df), simplify=FALSE) %>% as.data.frame

#####################
# Make final output #
#####################

#TODO: Add more covariates
column_to_PHE = phenoCol %>% str_replace_all("\\(","") %>%
        str_replace_all("\\)","") %>%
        str_replace_all("-[^-]+$", "") %>% 
        str_replace_all(' ','_') %>% 
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
cb_data_transformed = cb_data_transformed %>% select(FID, IID, PAT, MAT, PHE, everything(), -!!as.symbol(column_to_PHE), -`platekey_in_aggregate_vcf`)
### FID, IID this has to be the platekey metadata -> Agg VCF columns. 
write.table(cb_data_transformed, paste0(out_path,'.phe'), sep='\t',  quote=FALSE, row.names=FALSE)



