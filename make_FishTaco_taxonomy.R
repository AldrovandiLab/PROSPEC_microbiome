#!/usr/bin/Rscript


library(useful)
library(stringi)

args <- commandArgs(T)
infile <- args[1]
outfile <- args[2]

data <- read.table(infile, header=T, as.is=T, sep="\t", row.names=1)

taxa <- rownames(data)
tax_headers <- c("d", "p", "c", "o", "f", "g", "s")


df <- as.data.frame(matrix(NA, nrow=length(taxa), ncol=length(tax_headers)))
colnames(df) <- tax_headers
rownames(df) <- taxa

for (str in taxa) {
	arr <- unlist(stri_split_fixed(str, "|"))
	tmp <- data.frame(prefix=unlist(lapply(arr, function(x) unlist(stri_split_fixed(x, "__"))[1])), value=unlist(lapply(arr, function(x) unlist(stri_split_fixed(x, "__"))[2])))
	rownames(tmp) <- tmp$prefix
	df[str, tax_headers] <- sprintf("%s__%s", tax_headers, as.character(tmp[tax_headers, "value"]))
}

write.table(df, file=outfile, row.names=T, col.names=F, sep="\t", quote=F)
