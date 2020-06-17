
source("/Lab_Share/fanli/code/Carolyn_Yanavich/utils.R")


# generate taxon abundance and functional abundance files

ta <- read.table("/Lab_Share/Carolyn_Yanavich/071018/humann2/fishtaco/merged_abundance_table.kraken_raw.Complete.L7.txt", header=T, as.is=T, sep="\t", quote="", comment.char="", row.names=1)
fa <- read.table("/Lab_Share/Carolyn_Yanavich/071018/humann2/fishtaco/merged_genefamilies_KEGG.txt", header=T, as.is=T, sep="\t", quote="", comment.char="", row.names=1)
colnames(fa) <- gsub(".cat.filtered2_Abundance.RPKs", "", colnames(fa))

ta <- normalizeByCols(ta)
fa <- normalizeByCols(fa)
ta$Taxon <- rownames(ta)
fa$Function <- rownames(fa)
ta <- ta[,c("Taxon", setdiff(colnames(ta), "Taxon"))]
fa <- fa[,c("Function", setdiff(colnames(fa), "Function"))]

#colnames(ta)[1] <- "Taxon"
#colnames(fa)[1] <- "Function"

for (gr in c("fibrosis", "steatosis")) {
	labels <- read.table(sprintf("/Lab_Share/Carolyn_Yanavich/071018/humann2/fishtaco/labels_%s.txt", gr), header=T, as.is=T, sep="\t")
	sel <- labels$Sample
	
	# convert to relative abundance
	ta.sel <- ta[, sel]
	rownames(ta.sel) <- ta$Taxon
	# filter at 1% in at least 10% of samples
	ta.sel <- ta.sel[which(rowSums(ta.sel > 0.01) >= ceiling(ncol(ta.sel)/10)),]
	ta.sel$Taxon <- rownames(ta.sel)
	ta.sel <- ta.sel[,c("Taxon", setdiff(colnames(ta.sel), "Taxon"))]
	write.table(ta.sel, file=sprintf("/Lab_Share/Carolyn_Yanavich/071018/humann2/fishtaco/merged_abundance_table.kraken_raw.Complete.L7.relabund.%s.txt", gr), quote=F, row.names=F, col.names=T, sep="\t")
	
	# convert to relative abundance
	fa.sel <- fa[, sel]
	rownames(fa.sel) <- fa$Function
	# filter at 0.5% in at least 10% of samples
	fa.sel <- fa.sel[which(rowSums(fa.sel > 0.001) >= ceiling(ncol(fa.sel)/10)),]
	fa.sel$Function <- rownames(fa.sel)
	fa.sel <- fa.sel[,c("Function", setdiff(colnames(fa.sel), "Function"))]
	write.table(fa.sel, file=sprintf("/Lab_Share/Carolyn_Yanavich/071018/humann2/fishtaco/merged_genefamilies_KEGG.%s.txt", gr), quote=F, row.names=F, col.names=T, sep="\t")
}



## write files for FishTaco
#genefamilies <- read.table("/Lab_Share/Kathy_NAFLD/Shotgun/merged_genefamilies_KEGG.txt", header=T, as.is=T, sep="\t", row.names=1, comment.char="")
#colnames(genefamilies) <- gsub("PCMP_", "", gsub(".cat.filtered2_Abundance.RPKs", "", colnames(genefamilies)))
#genefamilies <- genefamilies[, colnames(data.norm)]
#inds_to_remove <- which(rowSums(genefamilies) < 1000)
#genefamilies <- genefamilies[setdiff(1:nrow(genefamilies), inds_to_remove),]
#inds_to_remove <- which(rowSums(genefamilies>100) < 10)
#genefamilies <- genefamilies[setdiff(1:nrow(genefamilies), inds_to_remove),]
#genefamilies.norm <- normalizeByCols(genefamilies)
#write.table(genefamilies.norm, file="/Lab_Share/Kathy_NAFLD/combined/FishTaco/function_abundance_KEGG.txt", quote=F, sep="\t", row.names=T, col.names=T)
