

mapping <- read.table("/Lab_Share/Carolyn_Yanavich/CY_Mapping_with_metadata.080817.txt", header=T, as.is=T, sep="\t", row.names=1, comment.char="")
mapping.sel <- subset(mapping, SampleType == "RectalSwab")

# generate sample labels for FishTaco
tmp <- mapping.sel[, c("Group"), drop=F]; colnames(tmp) <- "Label"; tmp$Sample <- rownames(tmp); tmp <- tmp[,c("Sample", "Label")]
for (gr in c("fibrosis", "steatosis")) {
	out <- subset(tmp, Label %in% c("normal", gr))
	out$Label <- ifelse(as.character(out$Label) == "normal", 0, 1)
	write.table(out, file=sprintf("/Lab_Share/Carolyn_Yanavich/PICRUSt/sample_labels.%s.txt", gr), quote=F, sep="\t", row.names=F, col.names=T)
}

# generate taxon abundance and functional abundance files
ta <- read.table("/Lab_Share/Carolyn_Yanavich/fastq/sample_counts.relabund.tab", header=T, as.is=T, sep="\t")
fa <- read.table("/Lab_Share/Carolyn_Yanavich/fastq/metagenome_contributions.L0.relabund.txt", header=T, as.is=T, sep="\t")
colnames(fa)[1] <- "Function"
for (gr in c("fibrosis", "steatosis")) {
	sel <- rownames(subset(mapping.sel, Group %in% c("normal", gr)))
	out <- ta[, c("Taxon", sel)]
	write.table(out, file=sprintf("/Lab_Share/Carolyn_Yanavich/PICRUSt/sample_counts.relabund.%s.tab", gr), quote=F, row.names=F, col.names=T, sep="\t")
	out <- fa[, c("Function", sel)]
	write.table(out, file=sprintf("/Lab_Share/Carolyn_Yanavich/PICRUSt/metagenome_contributions.L0.relabund.%s.txt", gr), quote=F, row.names=F, col.names=T, sep="\t")
}

# generate FishTaco taxonomy file (for visualization)

