
#!/usr/bin/env Rscript

library(data.table)

a <- as.data.frame(data.table::fread("analysis.csv"))
to_keep <- a[["SNPID"]]
b <- as.data.frame(data.table::fread("gwascat.csv"))
c <- b[ b$SNPS %in% to_keep, ]
data.table::fread(c, file = "gwascat_subset.csv", sep = ",")