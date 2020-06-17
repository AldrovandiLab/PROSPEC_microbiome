#!/usr/bin/Rscript

library(ggplot2)
library(ape)
library(plyr)
library(reshape2)
library(cluster)
library(RColorBrewer)
library(phyloseq)
library(grid)
library(gridExtra)
library(gplots)
library(vegan)
library(irr)
library(useful)
library(pscl)
library(parallel)
library(poLCA)
library(igraph)
library(randomForest)
library(ROCR)
library(stringi)
library(missForest)
library(DESeq2)
library(tableone)

source("/Lab_Share/fanli/code/Carolyn_Yanavich/utils.R")
source("/Lab_Share/fanli/code/Carolyn_Yanavich/mcc.R")

distance_metrics <- c("bray", "jaccard", "jsd")
alpha_metrics <- c("Chao1", "Shannon", "Simpson", "Observed")

## Full analysis of all PROSPEC shotgun samples run at QIAGEN (073018)

#########################################################################################################
## read in mapping files
out_txt <- sprintf("/Lab_Share/Carolyn_Yanavich/071018/phyloseq/phyloseq_output.%s.%s.txt", "PROSPEC_shotgun", format(Sys.Date(), "%m%d%y"))
out_pdf <- sprintf("/Lab_Share/Carolyn_Yanavich/071018/phyloseq/phyloseq_output.%s.%s.pdf", "PROSPEC_shotgun", format(Sys.Date(), "%m%d%y"))

mapping <- read.table("/Lab_Share/Carolyn_Yanavich/CY_Mapping_with_metadata.013019.txt", header=T, as.is=T, sep="\t", comment.char="", row.names=1)
mapping <- subset(mapping, Project=="PROSPEC")

## color tables
color_table <- read.table("/Lab_Share/fanli/code/Core.16S/taxa_coloring.Genus.050818.txt", header=T, as.is=T, sep="\t", comment.char="")
coloring <- color_table$Color
names(coloring) <- color_table$Genus
ordering <- rev(names(coloring))
cols <- colorRampPalette(c("white", "red"), space = "rgb")

ordering.genus <- color_table
ordering.genus$Phylum <- factor(ordering.genus$Phylum)
ordering.genus$Class <- factor(ordering.genus$Class)
inds=order(ordering.genus$Phylum, ordering.genus$Class, ordering.genus$Genus)
ordering.genus <- ordering.genus[inds,]
ordering.genus$Genus <- factor(ordering.genus$Genus, levels=unique(ordering.genus$Genus))

coloring.group <- c("#888888", "orange", "red"); names(coloring.group) <- c("normal", "steatosis", "fibrosis")
dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")

##########################################################################################################################
## Shotgun analysis
data <- read.table("/Lab_Share/Carolyn_Yanavich/071018/merged_abundance_table.kraken_raw.Complete.L7.txt", header=T, as.is=T, sep="\t", comment.char="", row.names=1, quote="")
colnames(data) <- gsub(".Complete.kraken.L7", "", colnames(data))
colnames(data) <- gsub("^X", "", colnames(data))
colnames(data) <- unlist(lapply(colnames(data), function(x) unlist(strsplit(x, "_"))[1]))
rownames(data) <- gsub("\\|", ";", rownames(data))
rownames(data) <- gsub("d__", "k__", rownames(data))

ranklist <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxlist <- lapply(rownames(data), function(x) parse_taxonomy_qiime(x))
taxa <- matrix(NA, nrow=nrow(data), ncol=length(ranklist)); colnames(taxa) <- ranklist
for (i in 1:length(taxlist)) {
	taxa[i, names(taxlist[[i]])] <- taxlist[[i]]
}
tt <- tax_table(taxa); rownames(tt) <- rownames(data)
rownames(mapping) <- sprintf("PROSPEC.%s", mapping$ProjectCode)
sel <- intersect(rownames(mapping), colnames(data))
data <- data[, sel]
mapping <- mapping[sel,]
ps <- phyloseq(otu_table(data, taxa_are_rows=T), tt, sample_data(mapping))
data <- as.data.frame(otu_table(ps))

pdf(out_pdf, width=12)

## read counts
df <- melt(sample_sums(ps)); df$SampleID <- rownames(df); df <- df[order(df$value),]; df$SampleID <- factor(df$SampleID, levels=df$SampleID)
df$log10value <- log10(df$value)
read_counts <- df # store for later use
p <- ggplot(df, aes(x=SampleID, y=value)) + geom_bar(stat="identity") + theme_classic() + ggtitle("Read counts") + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8))
print(p)
p <- ggplot(df, aes(x=SampleID, y=log10value)) + geom_bar(stat="identity") + theme_classic() + ggtitle("Read counts (log10)") + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8))
print(p)

## taxa barplots of all samples (QC check)
ps.relative <- transform_sample_counts(ps, function(x) x / sum(x) )
ps.sel <- ps.relative
ordi <- ordinate(ps.sel, method = "PCoA", distance = "bray")
ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))
otu.filt <- normalizeByCols(as.data.frame(otu_table(ps.sel)))
otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level="Genus")
agg <- aggregate(. ~ Genus, otu.filt, sum)
genera <- agg$Genus
agg <- agg[,-1]
agg <- normalizeByCols(agg)
inds_to_grey <- which(rowMeans(agg)<0.005)
genera[inds_to_grey] <- "Other"
agg$Genus <- genera
df <- melt(agg, variable.name="SampleID")
df2 <- aggregate(as.formula("value ~ Genus + SampleID"), df, sum)
df2$SampleID <- as.character(df2$SampleID)
df2$SampleIDfactor <- factor(df2$SampleID, levels=ordering.pc1)
p <- ggplot(df2, aes_string(x="SampleID", y="value", fill="Genus", order="Genus")) + geom_bar(stat="identity", position="stack") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring) + ggtitle(sprintf("Taxa barplots all samples (QC)")) + guides(col = guide_legend(ncol = 3)) + ylim(c(-.1, 1.01))
print(p)
p <- ggplot(df2, aes_string(x="SampleIDfactor", y="value", fill="Genus", order="Genus")) + geom_bar(stat="identity", position="stack") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring) + ggtitle(sprintf("Taxa barplots all samples (QC)")) + guides(col = guide_legend(ncol = 3)) + ylim(c(-.1, 1.01))
print(p)

##########################################################################################################################
## Trim down to final dataset of PROSPEC samples for analysis
## format metadata
metadata_variables <- read.table("/Lab_Share/Carolyn_Yanavich/metadata_types.PROSPEC.txt", header=T, as.is=T, sep="\t", row.names=1)
sel <- intersect(rownames(metadata_variables), colnames(mapping))
metadata_variables <- metadata_variables[sel,, drop=F]
mapping.sel <- mapping[rownames(sample_data(ps)), sel]
## fix column types
for (mvar in rownames(metadata_variables)) {
	if (metadata_variables[mvar, "type"] == "factor") {
		mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
		if (metadata_variables[mvar, "baseline"] != "") {
			mapping.sel[,mvar] <- relevel(mapping.sel[,mvar], metadata_variables[mvar, "baseline"])
		}
	} else if (metadata_variables[mvar, "type"] == "numeric") {
		mapping.sel[,mvar] <- as.numeric(as.character(mapping.sel[,mvar]))
	}
}
## impute missing values
set.seed(sum(dim(mapping.sel)))
imputed <- missForest(mapping.sel)$ximp
mapping.sel <- imputed 
sample_data(ps) <- mapping.sel

## demographics tables
demo_vars <- c("Age", "BMI", "Sex", "Hypertension", "Diabetes", "Dyslipidemia", "ARV", "bb_cum", "bb_curr", "Nadir_CD4", "Current_CD4", "total_calories", "pct_carbohydrates", "pct_protein", "pct_fat", "cholesterol", "total_fiber", "Race")
nonnormal_vars <- c("Age", "BMI", "Sex", "Hypertension", "Diabetes", "Dyslipidemia", "ARV", "bb_cum", "bb_curr", "Nadir_CD4", "Current_CD4")
write.table(print(CreateTableOne(vars=demo_vars, strata=c("Group"), data=mapping.sel, smd=T), nonnormal=nonnormal_vars), file="/Lab_Share/Carolyn_Yanavich/071018/phyloseq/Table_1.Group.txt", quote=F, sep="\t", row.names=T, col.names=T)
mapping.tmp <- subset(mapping.sel, Group %in% c("Normal", "Steatosis")); mapping.tmp$Group <- droplevels(mapping.tmp$Group)
write.table(print(CreateTableOne(vars=demo_vars, data=mapping.tmp, smd=T), nonnormal=nonnormal_vars), file="/Lab_Share/Carolyn_Yanavich/071018/phyloseq/Table_1.Normal_Steatosis.txt", quote=F, sep="\t", row.names=T, col.names=T)
mapping.tmp <- subset(mapping.sel, Group %in% c("Normal", "Fibrosis")); mapping.tmp$Group <- droplevels(mapping.tmp$Group)
write.table(print(CreateTableOne(vars=demo_vars, strata=c("Group"), data=mapping.tmp, smd=T), nonnormal=nonnormal_vars), file="/Lab_Share/Carolyn_Yanavich/071018/phyloseq/Table_1.Normal_Fibrosis.txt", quote=F, sep="\t", row.names=T, col.names=T)
mapping.tmp <- subset(mapping.sel, Group %in% c("Steatosis", "Fibrosis")); mapping.tmp$Group <- droplevels(mapping.tmp$Group)
write.table(print(CreateTableOne(vars=demo_vars, strata=c("Group"), data=mapping.tmp, smd=T), nonnormal=nonnormal_vars), file="/Lab_Share/Carolyn_Yanavich/071018/phyloseq/Table_1.Steatosis_Fibrosis.txt", quote=F, sep="\t", row.names=T, col.names=T)

## make relative and rarefied objects
ps.relative <- transform_sample_counts(ps, function(x) x / sum(x) )
ps.rarefied <- rarefy_even_depth(ps, sample.size=5018, rngseed=nsamples(ps))

## distance matrices
dm <- list()
dm[["PROSPEC"]] <- list()
for (distance_metric in distance_metrics) {
	dm[["PROSPEC"]][[distance_metric]] <- as.matrix(phyloseq::distance(ps.relative, method=distance_metric))
}

## thresholds for association/permutation tests
nsamps_threshold <- 0.01 # relative abundance needed to call a sample positive
filt_threshold <- 0.1 # fraction of samples that need to be positive to keep an OTU for association testing
nperm <- 100000
siglevel <- 0.05

##########################################################################################################################
## analysis

ps.sel <- ps.relative
mapping.sel <- as(sample_data(ps.sel), "data.frame")
mapping.sel$SampleID <- rownames(mapping.sel)

## PCoA
for (distance_metric in distance_metrics) {
	ordi <- ordinate(ps.sel, method = "PCoA", distance = distance_metric)
	p <- plot_ordination(ps.sel, ordi, "samples", color = "Group") + theme_classic() + ggtitle(distance_metric) + theme_classic() + stat_ellipse(type="t")
	print(p)
}
## PERMANOVA
sink(out_txt, append=F)
print("PERMANOVA on PROSPEC samples")
for (distance_metric in distance_metrics) {
	print(distance_metric)
	form <- as.formula(sprintf("as.dist(dm[[\"PROSPEC\"]][[distance_metric]]) ~ %s", paste(rownames(subset(metadata_variables, useForPERMANOVA=="yes")), collapse="+")))
	res <- adonis(form , data=as(sample_data(ps.relative), "data.frame"), permutations=999)
	print(res)
}
sink()

## taxa barplots
# order by PC1 (Bray-Curtis)
ordi <- ordinate(ps.sel, method = "PCoA", distance = "bray")
ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))
otu.filt <- as.data.frame(otu_table(ps.sel))
otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level="Genus")
agg <- aggregate(. ~ Genus, otu.filt, sum)
genera <- agg$Genus
agg <- agg[,-1]
agg <- sweep(agg, 2, colSums(agg), "/")
genera[which(rowMeans(agg)<0.01)] <- "Other"
agg$Genus <- genera
df <- melt(agg, variable.name="SampleID")
agg <- aggregate(value~Genus+SampleID, df, sum)
agg$SampleID <- as.character(agg$SampleID)
agg$SampleIDfactor <- factor(agg$SampleID, levels=ordering.pc1)
agg$Genus <- factor(agg$Genus, levels=levels(ordering.genus$Genus))
# taxa barplot representation
df.SampleIDstr <- unique(agg[,c("SampleID", "SampleIDfactor")])
df.SampleIDstr$Group <- as.character(mapping.sel[df.SampleIDstr$SampleID, "Group"])
p <- ggplot(agg, aes(x=SampleIDfactor, y=value, fill=Genus, order=Genus)) + geom_bar(stat="identity", position="stack") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=4)) + ggtitle(sprintf("%s taxa summary (L5, ordered by PC1)", "PROSPEC")) + scale_fill_manual(values=coloring) + ylim(c(-.1, 1.01)) + annotate("rect", xmin = as.numeric(df.SampleIDstr$SampleIDfactor)-0.5, xmax = as.numeric(df.SampleIDstr$SampleIDfactor)+0.5, ymin = -0.04, ymax = -0.02, fill=coloring.group[df.SampleIDstr$Group])
print(p)
agg$Group <- mapping.sel[agg$SampleID, "Group"]
p <- ggplot(agg, aes(x=SampleIDfactor, y=value, fill=Genus, order=Genus)) + geom_bar(stat="identity", position="stack") + facet_wrap(~Group, scales="free_x") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=4)) + ggtitle(sprintf("%s taxa summary (L5, ordered by PC1)", "PROSPEC")) + scale_fill_manual(values=coloring) + ylim(c(-.1, 1.01))
print(p)

## alpha diversity by Group
adiv <- estimate_richness(ps.rarefied, measures=alpha_metrics)
adiv$SampleID <- rownames(adiv)
adiv <- merge(adiv, mapping.sel, by="row.names"); rownames(adiv) <- adiv$SampleID
sink(out_txt, append=T)
for (mvar in c("Group")) {
	plotlist <- list()
	for (alpha_metric in alpha_metrics) {
		if (nlevels(adiv[,mvar]) > 2) {
			test <- kruskal.test(as.formula(sprintf("%s ~ %s", alpha_metric, mvar)), adiv)
		}
		else {
			test <- wilcox.test(as.formula(sprintf("%s ~ %s", alpha_metric, mvar)), adiv)
		}
		p <- ggplot(adiv, aes_string(x=mvar, y=alpha_metric)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("%s by %s (%s p=%.4g)", alpha_metric, mvar, test$method, test$p.value)) + theme(title=element_text(size=8))
		plotlist[[length(plotlist)+1]] <- p
		# GLM
		mod <- glm(as.formula(sprintf("%s ~ Group + %s", alpha_metric, paste(rownames(subset(metadata_variables, useAsConfounder=="yes")), collapse="+"))), family="gaussian", adiv)
		print(sprintf("GLM with %s", alpha_metric))
		print(summary(mod))
	}
	multiplot(plotlist=plotlist, cols=2, rows=2)
}
sink()

## beta diversity by Group
pairs <- as.data.frame(t(combn(mapping.sel$SampleID, 2))); colnames(pairs) <- c("s1", "s2"); pairs$s1 <- as.character(pairs$s1); pairs$s2 <- as.character(pairs$s2)
pairs$g1 <- mapping.sel[pairs$s1, "Group"]; pairs$g2 <- mapping.sel[pairs$s2, "Group"]
pairs$Group <- factor(ifelse(pairs$g1==pairs$g2, "Within", "Between"))
for (distance_metric in distance_metrics) {
	pairs[, distance_metric] <- dm[["PROSPEC"]][[distance_metric]][cbind(pairs$s1, pairs$s2)]
	test <- kruskal.test(as.formula(sprintf("%s ~ %s", distance_metric, "Group")), pairs)
	p <- ggplot(pairs, aes_string(x="Group", y=distance_metric)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("%s by %s (%s p=%.4g)", distance_metric, "Group", test$method, test$p.value)) + theme(title=element_text(size=8))
	print(p)
}


## DESeq2
ps.deseq2 <- ps
res <- {}
for (mvar in c("Group")) {
	dds <- phyloseq_to_deseq2(ps.deseq2, as.formula(sprintf("~ %s", mvar)))
	data <- as.data.frame(otu_table(ps))
	data.rel <- normalizeByCols(data)
	ftk <- names(which(unlist(apply(data.rel, 1, function(x) length(which(x>=nsamps_threshold)))) >= ceiling(filt_threshold*ncol(data))))
	dds <- dds[ftk,]
	dds <- DESeq(dds)
	for (rn in setdiff(resultsNames(dds), "Intercept")) {
		tmp <- results(dds, name=rn)
		tmp$taxa <- unlist(lapply(rownames(tmp), function(x) {
			tmp <- unlist(strsplit(x, ";"))
			tmp[length(tmp)]
		}))
		tmp$dir <- ifelse(tmp$padj < siglevel, ifelse(sign(tmp$log2FoldChange)==1, "up", "down"), "NS")
		tmp$comparison <- rn
		tmp$metadata_variable <- mvar
		res <- rbind(res, tmp)
	}
}
res <- as.data.frame(res)
res$padj <- p.adjust(res$pvalue, method="fdr")
res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$log2FoldChange)==1, "up", "down"), "NS")
write.table(res, file="/Lab_Share/Carolyn_Yanavich/071018/phyloseq/Shotgun_DESeq2.Species.txt", quote=F, sep="\t", row.names=T, col.names=T)
# forest plots for each variable
res <- read.table("/Lab_Share/Carolyn_Yanavich/071018/phyloseq/Shotgun_DESeq2.Species.txt", header=T, row.names=1, sep="\t", as.is=T, quote="", comment.char="")
for (comp in unique(res$comparison)) {
	resSig <- subset(res, padj<siglevel & comparison==comp); resSig <- resSig[order(resSig$log2FoldChange),]
	resSig$taxa <- factor(resSig$taxa, levels=resSig$taxa)
	p <- ggplot(resSig, aes(x=taxa, y=log2FoldChange)) + geom_bar(stat="identity", fill="#aaaaaa") + geom_text(aes(label=taxa), y=0, size=2, hjust=0.5) + coord_flip() + theme_classic() + ggtitle(sprintf("DESeq2 hits (%s)", comp)) + theme(axis.text.y=element_blank())
	print(p)
	lims <- max(abs(resSig$log2FoldChange) + abs(resSig$lfcSE))*1.0
	p <- ggplot(resSig, aes(x=taxa, y=log2FoldChange, color=dir)) + geom_point() + geom_errorbar(aes(x=taxa, ymin=log2FoldChange-lfcSE, max=log2FoldChange+lfcSE), width=0.2) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("%s (%s)", "DESeq2 hits", comp)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
	print(p)
}
# heatmap of all results
df <- res
df2 <- dcast(df, taxa ~ comparison, value.var="log2FoldChange"); rownames(df2) <- df2$taxa; df2 <- df2[,-1]
df2.sig <- dcast(df, taxa ~ comparison, value.var="padj"); rownames(df2.sig) <- df2.sig$taxa; df2.sig <- df2.sig[,-1]
df2.sigstr <- matrix("", nrow=nrow(df2.sig), ncol=ncol(df2.sig))
df2.sigstr[which(df2.sig < siglevel, arr.ind=T)] <- "*"
#df2[which(df2.sig >= siglevel, arr.ind=T)] <- 0
heatmap.2(as.matrix(df2), Rowv=T, Colv=T, dendrogram="both", trace="none", col=colorRampPalette(c("blue", "white", "red"))(150), margins=c(6,24), cexCol=0.8, cexRow=0.4, cellnote=df2.sigstr, notecex=0.8, notecol="black", main=sprintf("log2FC values of all taxa"))
# heatmap of VSD values
mvar <- "Group"
dds <- phyloseq_to_deseq2(ps.deseq2, as.formula(sprintf("~ %s", mvar)))
data <- as.data.frame(otu_table(ps))
data.rel <- normalizeByCols(data)
ftk <- names(which(unlist(apply(data.rel, 1, function(x) length(which(x>=nsamps_threshold)))) >= ceiling(filt_threshold*ncol(data))))
dds <- dds[ftk,]
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind=T)
df <- assay(vsd)
colcolors <- coloring.group[mapping.sel[colnames(df), "Group"]]
heatmap.2(as.matrix(df), Rowv=T, Colv=T, dendrogram="both", trace="none", col=colorRampPalette(c("white", "blue"))(150), ColSideColors=colcolors, margins=c(6,24), cexCol=0.8, cexRow=0.4, main=sprintf("VSD of all taxa"))



dev.off()


