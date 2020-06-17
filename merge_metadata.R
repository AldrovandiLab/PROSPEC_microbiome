#!/usr/bin/Rscript

library(useful)

#mapping <- read.table("/Lab_Share/Carolyn_Yanavich/CY_Mapping.080817.txt", header=T, as.is=T, sep="\t", comment.char="")
#metadata <- read.table("/Lab_Share/Carolyn_Yanavich/metadata.012418.txt", header=T, as.is=T, sep="\t", comment.char="")
#addtl <- read.table("/Lab_Share/Carolyn_Yanavich/metadata.012518.txt", header=T, as.is=T, sep="\t", comment.char="")

### specimen collection section
## make an ID to link to mapping sheet
#mapping$pr_pront <- gsub("^([A-Z]*)", "", mapping$X.SampleID)
#inds <- match(mapping$pr_pront, metadata$pr_pront) # from mapping -> metadata

### grab extra t_arv_atual values
##metadata[match(addtl$pr_pront, metadata$pr_pront), "t_arv_atual2"] <- addtl$t_arv_atual

## get variables of interest
#mapping$Group <- metadata[inds, "group"]
#mapping$Sex <- ifelse(metadata[inds, "Sex"]==0, "F", "M")
#mapping$Age <- metadata[inds, "Age..Years."]
#mapping$BMI <- metadata[inds, "BMI"]
#mapping$Nadir_CD4 <- metadata[inds, "nadir_cd4_final"]
#mapping$Current_CD4 <- metadata[inds, "cd4_final"]
#mapping$Diabetes <- ifelse(metadata[inds, "Diabetes"]==0, "No", "Yes")
#mapping$Hypertension <- ifelse(metadata[inds, "Hypertension"]==0, "No", "Yes")
#mapping$CoronaryArteryDisease <- ifelse(metadata[inds, "Coronary.Artery.Disease"]==0, "No", "Yes") # all negative in the set of 82
#mapping$Dyslipidemia <- ifelse(metadata[inds, "Dislypidemia"]==0, "No", "Yes")
#mapping$Etiol <- metadata[inds, "etiol"] # only a single value != 1
#mapping$ARV <- metadata[inds, "t_total_art"]
#mapping$bb_cum <- ifelse(metadata[inds, "bb_cum"]==0, "TDF", "AZT")

#colnames(mapping)[1] <- "#SampleID"
#write.table(mapping, file="/Lab_Share/Carolyn_Yanavich/CY_Mapping_with_metadata.080817.txt", quote=F, row.names=F, col.names=T, sep="\t")


# add diet and race data to mapping
mapping <- read.table("/Lab_Share/Carolyn_Yanavich/CY_Mapping_with_metadata.121018.txt", header=T, as.is=T, sep="\t", comment.char="", quote="")
diet <- read.table("/Lab_Share/Carolyn_Yanavich/dietary.txt", header=T, as.is=T, sep="\t", row.names=1)
sel <- intersect(rownames(diet), mapping$pr_pront)
diet <- diet[sel,]
for (mvar in c("total_calories", "pct_carbohydrates", "pct_protein", "pct_fat", "cholesterol", "total_fiber")) {
	mapping[, mvar] <- diet[mapping$pr_pront, mvar]
}

race <- read.table("/Lab_Share/Carolyn_Yanavich/race.txt", header=T, as.is=T, sep="\t", row.names=1)
mapping$Race <- race[mapping$pr_pront, "Race"]


write.table(mapping, file="/Lab_Share/Carolyn_Yanavich/CY_Mapping_with_metadata.013019.txt", quote=F, sep="\t", row.names=F, col.names=T)



