#!/usr/bin/env Rscript

library(plyr)
library(data.table)

paths <- list.files(".", pattern = "txt", full.names = TRUE)
list_of_dfs <- lapply(paths,data.table::fread)
analysis <- plyr::rbind.fill(list_of_dfs)
data.table::fwrite(analysis, "analysis.csv", sep = ",")