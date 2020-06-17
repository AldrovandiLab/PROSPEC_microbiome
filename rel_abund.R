#! /usr/bin/Rscript

# convert an OTU/ASV table to relative abundance
source("utils.R")

args <- commandArgs(T)
if (length(args) < 2) {
	cat("USAGE: ./rel_abund.R count_table_file out_file\n")
	q()
}

tab_file <- args[1]
out_file <- args[2]

data <- read.table(tab_file, header=T, as.is=T, sep="\t", row.names=1, comment.char="", skip=0, quote="")
out <- normalizeByCols(data)
tmp <- colnames(out)
out$Taxon <- rownames(out)
out <- out[,c("Taxon", tmp)]

write.table(out, file=out_file, row.names=F, col.names=T, quote=F, sep="\t")

