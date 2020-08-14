#!/usr/bin/env Rscript

############################## ARGUMENTS SECTION #############################
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no all arguments passed or help needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) | (length(args) == 1) ) {
  cat("
      The helper R Script manhattan.R
      Mandatory arguments:
          --saige_output=path        - The path to the SAIGE output file.
                                       NOTE:
                                       If you have performed the analysis per chromosome,
                                       concatenate outputs across all chromosomes first.

          --output_tag=value         - A string in single quotes used as the identifier of the analysis
                                       in the output files name. Don't use whitesspaces or irregular characters.
                                       The name will be converted to snakecase (eg. snake_case)
         --help                      - you are reading it
         
      Optionnal arguments:

          --title=chr                 - A string in single quotes used as the title of the plot
                                        Default: 'Manhattan plot with rsid annotation (p-value cuttof: $p_value_cutoff, gwas id: $output_tag)'
          --width=val                 - A value to denote width of plot (base R)
                                        Default: 2200
          --height=val                - A value to denote height of plot (base R)
                                        Default: 1400
          --units=chr                 - A value to denote units of plot (base R)
                                        Default: 'px'
          --res=val                   - A value to denote res of plot (base R)
                                        Default: 300
          --type=chr                  - A value to denote type of plot (base R)
                                        Default: 'cairo'
                                        NOTE: If unsure what you should use you can check your system's 
                                        default with the command getOption('bitmapType')
          --cex=val                   - A value to denote the cex for the manhattan plot
                                        Default: 0.6
          --cex_axis=val              - A value to denote the cex.axis for the manhattan plot
                                        Default: 0.6
          --suggestive_line=boolean   - Boolean. If set to FALSE the suggestive line
                                        in the manhattan plot is omitted.
                                        Default: TRUE
          --manhattan_colour_1=chr    - A string in single quotes, can be either base R colour name or hexcode.
                                        Default: '#4e4b4c' (Genomics England dark gray)
          --manhattan_colour_2=chr    - A string in single quotes, can be either base R colour name or hexcode.
                                        Default: '#4dc5ce' (Genomics England teal)
          --p_value_cutoff=val        - Significance cuttof to be used for annotation of snps.
                                        Default: 0.01
          --bp_col                    - The name of the column with the position information.
                                        Default: 'POS' (SAIGE)
          --snp_col                   - The name of the column with the rsid.
                                        Default: 'SNPID' (SAIGE)
          --p_val_col                 - The name of the column with the p-values.
                                        Default: 'p.value' (SAIGE)              

     Usage:
      
          The typical command for running the script is as follows:
    
          ./manhattan.R --saige_output='saige_out.txt' --output_tag='covid_2'
      
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
if(is.null(args$p_value_cutoff)) {args$p_value_cutoff = 0.01 } else {args$p_value_cutoff=as.numeric(args$p_value_cutoff)}
if(is.null(args$title)) {args$title  = paste0(title = "Manhattan plot with rsid annotation (p-value < ", as.character(args$p_value_cutoff), ", gwas id: ", args$output_tag, ")")} else {args$title=as.character(args$title)}
if(is.null(args$width)) {args$width = 2200} else {args$width=as.numeric(args$width)}
if(is.null(args$height)) {args$height = 1400} else {args$height=as.numeric(args$height)}
if(is.null(args$units)) {args$units = "px"} else {args$units=as.character(args$units)}
if(is.null(args$type)) {args$type = "cairo" }  else {args$type=as.character(args$type)}
if(is.null(args$res)) {args$res = 300} else {args$res=as.character(args$res)}
if(is.null(args$cex)) {args$cex = 0.6} else {args$cex=as.character(args$cex)}
if(is.null(args$cex_axis)) {args$cex_axis = 0.6} else {args$cex_axis=as.character(args$cex_axis)}
if(is.null(args$suggestive_line)) {args$suggestive_line  = TRUE} else {args$suggestive_line=as.logical(args$suggestive_line)}
if(is.null(args$manhattan_colour_1)) {args$manhattan_colour_1 = "#4e4b4c"} else {args$manhattan_colour_1=as.character(args$manhattan_colour_1)}
if(is.null(args$manhattan_colour_2)) {args$manhattan_colour_2 = "#4dc5ce"} else {args$manhattan_colour_2=as.character(args$manhattan_colour_2)}
if(is.null(args$bp_col)) {args$bp_col = "POS"} else {args$bp_col=as.character(args$bp_col)}
if(is.null(args$snp_col)) {args$snp_col = "SNPID"} else {args$snp_col=as.character(args$snp_col)}
if(is.null(args$p_val_col)) {args$p_val_col = "p.value"} else {args$p_val_col=as.character(args$p_val_col)}

############################## LIBRARIES SECTION #############################

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(qqman)))
suppressWarnings(suppressMessages(library(snakecase)))


# ######################### VARIABLES REASSIGNMENT SECTION ###############################

# Facilitates testing and protects from wh-spaces, irregular chars

# required
saige_output       <- args$saige_output
output_tag         <- snakecase::to_snake_case(args$output_tag)

# optional
title              <- args$title
width              <- args$width
height             <- args$height
units              <- args$units
type               <- args$type
res                <- args$res
cex                <- args$cex
cex_axis           <- args$cex_axis
suggestive_line    <- args$suggestive_line
manhattan_colour_1 <- args$manhattan_colour_1
manhattan_colour_2 <- args$manhattan_colour_2
p_value_cutoff     <- args$p_value_cutoff
bp_col             <- args$bp_col
snp_col            <- args$snp_col
p_val_col          <- args$p_val_col

cat("\n")
cat("ARGUMENTS SUMMARY")
cat("\n")
cat("saige_output       : ", saige_output       ,"\n",sep="")
cat("output_tag         : ", output_tag         ,"\n",sep="")
cat("p_value_cutoff     : ", p_value_cutoff     ,"\n",sep="")
cat("bp_col             : ",  bp_col            ,"\n",sep="")
cat("snp_col            : ",  snp_col           ,"\n",sep="")
cat("p_val_col          : ",  p_val_col         ,"\n",sep="")
cat("suggestive_line    : ", suggestive_line    ,"\n",sep="")
cat("manhattan_colour_1 : ", manhattan_colour_1 ,"\n",sep="")
cat("manhattan_colour_2 : ", manhattan_colour_2 ,"\n",sep="")
cat("width              : ", width              ,"\n",sep="")
cat("height             : ", height             ,"\n",sep="")
cat("units              : ", units              ,"\n",sep="")
cat("type               : ", type               ,"\n",sep="")
cat("res                : ", res                ,"\n",sep="")
cat("cex                : ", cex                ,"\n",sep="")
cat("cex_axis           : ", cex_axis           ,"\n",sep="")
cat("suggestive_line    : ", suggestive_line    ,"\n",sep="")
cat("manhattan_colour_1 : ", manhattan_colour_1 ,"\n",sep="")
cat("manhattan_colour_2 : ", manhattan_colour_2 ,"\n",sep="")

# ############################### SCRIPT SECTION ###############################

#Plotting manhattan plots

analysis <- data.table::fread(saige_output)

png(filename = paste0(output_tag, "_", "manhattan" , ".png") ,
    width    = width   ,
    height   = height  ,
    units    = units  ,
    res      = res,
    type     = type)

manhattan_plot <- qqman::manhattan(analysis,
                                   bp   = bp_col,
                                   p    = p_val_col,
                                   snp  = snp_col,
                                   main = title,
                                   col  = c(manhattan_colour_1, manhattan_colour_2),
                                   highlight      =  analysis[[snp_col]] [analysis[[p_val_col]] < p_value_cutoff],
                                   annotatePval   = p_value_cutoff,
                                   suggestiveline = suggestive_line)

dev.off()

