# tabulate demographics table

library(useful)


mapping_fn <- "/Lab_Share/Carolyn_Yanavich/CY_Mapping_with_metadata.013019.txt"
mapping <- read.table(mapping_fn, header=T, as.is=T, sep="\t", comment.char="")
mapping <- subset(mapping, SampleType=="RectalSwab")
rownames(mapping) <- mapping$pr_pront

additional <- read.table("/Lab_Share/Carolyn_Yanavich/Copy of EnglishPROSPEC_swab_FINAL_for_Carolyn_Nov2017.txt", header=T, as.is=T, sep="\t", row.names=1)

mapping <- merge(mapping, additional, by="row.names"); rownames(mapping) <- mapping$Row.names
mapping$Group <- factor(mapping$Group, levels=c("Normal", "Steatosis", "Fibrosis"))



# lab values
lab_vars <- c("ALT", "AST", "tap", "GGT", "Albumin", "Glucose", "Choletsterol", "LDL", "HDL", "Triglycerides")
for (lv in lab_vars) {
	print(aggregate(as.formula(sprintf("%s ~ Group", lv)), mapping, summary))
	print(kruskal.test(as.formula(sprintf("%s ~ Group", lv)), mapping))
}

lab_vars <- c("Statin.use", "Triglyceride.med")
for (lv in lab_vars) {
	tab <- table(mapping[,c("Group", lv)])
	print(tab)
	print(fisher.test(tab))
	
}
