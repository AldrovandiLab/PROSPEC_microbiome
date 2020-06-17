#!/usr/bin/Rscript


library(useful)
source("utils.R")


for (gr in c("fibrosis", "steatosis")) {
	data <- read.table(sprintf("/Lab_Share/Carolyn_Yanavich/PICRUSt/metagenome_contributions.L0.relabund.%s.MUSiCC_corrected.Pathway.txt", gr), header=T,as.is=T, sep="\t", row.names=1)

	# convert to relative abundance
	data <- normalizeByCols(data)

	# filter
	inds <- rowMeans(data) > 0.001
	data <- data[inds,]

	# reorder and add 'Pathway' column back
	sel <- colnames(data)
	data$Pathway <- rownames(data)
	out <- data[,c("Pathway", sel)]

	write.table(out, file=sprintf("/Lab_Share/Carolyn_Yanavich/PICRUSt/metagenome_contributions.L0.relabund.%s.MUSiCC_corrected.Pathway.filtered.txt", gr), quote=F, sep="\t", row.names=F, col.names=T)
}

####
# after running compute_differential_abundance.py
for (gr in c("fibrosis", "steatosis")) {
	res=read.table(sprintf("/Lab_Share/Carolyn_Yanavich/PICRUSt/DA_functions.%s.tab", gr),header=T,as.is=T,sep="\t",row.names=1)
	res$padj=p.adjust(res$pval, method="fdr")
	res=res[order(res$pval),]

	# store list of DA pathways at padj<0.2
	out <- rownames(subset(res, padj<0.2))
	write.table(out, file=sprintf("/Lab_Share/Carolyn_Yanavich/PICRUSt/function_abundance_KEGG.MUSiCC_corrected.Pathway.filtered.%s.selected_DA_pathways.txt", gr), quote=F, sep="\t", row.names=F, col.names=F)
}


## hack the singLogP column to get reverse enrichment
#res$singLogP <- -res$singLogP
#res$StatValue <- -res$StatValue
#a <- res$meanCases
#res$meanCases <- res$meanControls
#res$meanControls <- a
#sel <- setdiff(colnames(res), "padj")
#res$Function <- rownames(res)
#res <- res[,c("Function", sel)]
#colnames(res) <- c("Function", "meanCases", "meanControls", "StatValue", "pval", "singLogP", "Bonf", "FDR-0.01", "FDR-0.05", "FDR-0.1", "Metadata")
#write.table(res, file="/Lab_Share/Kathy_NAFLD/combined/FishTaco/DA_functions.reversed.tab", quote=F, row.names=F, col.names=T, sep="\t")





