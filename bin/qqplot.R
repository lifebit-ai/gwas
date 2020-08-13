#!/usr/bin/env Rscript

############################## ARGUMENTS SECTION #############################
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no all arguments passed or help needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) | (length(args) == 1) ) {
  cat("
      The helper R Script qqplot.R
      Mandatory arguments:
          --saige_output=path        - The path to the SAIGE output file
          --output_tag=value         - A string in single quotes used as the identifier of the analysis
                                       in the output files name. Don't use whitesspaces or irregular characters.
                                       The name will be converted to snakecase (eg. snake_case)
         --help                      - you are reading it
         
      Optionnal arguments:

          --width=val                 - A value to denote width of plot (base R)
                                        Default: 12
          --height=val                - A value to denote height of plot (base R)
                                        Default: 6
          --units=chr                 - A value to denote units of plot (base R)
                                        Default: 'in'
          --res=val                   - A value to denote res of plot (base R)
                                        Default: 300
          --type=chr                  - A value to denote type of plot (base R)
                                        Default: getOption('bitmapType')

     Usage:
      
          The typical command for running the script is as follows:
    
          ./qqplot.R --saige_output='saige_out.txt' --output_tag='covid_2'
      
      WARNING : here put all the things the user has to know
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
if(is.null(args$width)) {args$width = 12} else {args$width=as.numeric(args$width)}
if(is.null(args$height)) {args$height = 6} else {args$height=as.numeric(args$height)}
if(is.null(args$units)) {args$units = "in"} else {args$units=as.character(args$units)}
if(is.null(args$type)) {args$type = getOption("bitmapType") }  else {args$type=as.character(args$type)}
if(is.null(args$res)) {args$res = 300} else {args$res=as.character(args$res)}

############################## LIBRARIES SECTION #############################

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(GWASTools)))
suppressWarnings(suppressMessages(library(snakecase)))


# ######################### VARIABLES REASSIGNMENT SECTION ###############################

# Facilitates testing and protects from wh-spaces, irregular chars

# required
saige_output      <- args$saige_output
output_tag        <- snakecase::to_snake_case(args$output_tag)

# optional
width             <- args$width
height            <- args$height
units             <- args$units
type              <- args$type
res               <- args$res

cat("\n")
cat("ARGUMENTS SUMMARY")
cat("\n")
cat("saige_output      : ", saige_output,     "\n",sep="")
cat("output_tag        : ", output_tag,       "\n",sep="")
cat("saige_output      : ", saige_output,     "\n",sep="")
cat("output_tag        : ", output_tag,       "\n",sep="")
cat("saige_output      : ", saige_output,     "\n",sep="")
cat("output_tag        : ", output_tag,       "\n",sep="")
cat("saige_output      : ", saige_output,     "\n",sep="")
cat("output_tag        : ", output_tag,       "\n",sep="")

# ############################### FUNCTIONS SECTION ###############################

lambda<-function(pvalues){
    chisq <- qchisq(1-pvalues,1)
    lambda <- median(chisq,na.rm=TRUE)/qchisq(0.5,1)
    return(lambda)
}

# ############################### SCRIPT SECTION ###############################

#Plotting qqplots

analysis <- data.table::fread(saige_output)

colnames(analysis)[which(colnames(analysis)=="POS")]="BP"
colnames(analysis)[which(colnames(analysis)=="p.value")]="P"

png(filename =  paste0(output_tag,"_", "qqplot_ci" , ".png") ,
    width  = width,
    height = height,
    units  = units,
    bg     = bg,
    res    = res,
    type   = type)


GWASTools::qqPlot(analysis$P, main=paste0(analysis_tag,",","lambda=",round(lambda(analysis$P),2)))
dev.off()

#check whether we can just save the R objects and also give the option for exporting through a future airlock
