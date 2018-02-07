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
library(parallel)
library(useful)
library(pscl)
library(MASS)
library(boot)
library(igraph)
library(ggfortify)
library(metagenomeSeq)
library(stringi)
library(car)
library(randomForest)
library(ROCR)
library(Hmisc)
library(missForest)

source("utils.R")
source("mcc.R")
indir <- "/Lab_Share/Carolyn_Yanavich"

distance_metrics <- c("bray", "jaccard", "jsd")
alpha_metrics <- c("Chao1", "Shannon", "Simpson", "Observed")
icc_levels <- c(5,6,7)
names(icc_levels) <- c("Family", "Genus", "Species")
ncores <- 16
cols.sig <- c("black", "red", "grey"); names(cols.sig) <- c("ns", "sig", "NA")

#################################################################
## handoff to phyloseq
seqtab.nochim <- readRDS(sprintf("%s/data/merged_seqtab.nochim.rds", indir))
dada2_fn <- sprintf("%s/data/DADA2.RData", indir)
blast_fn <- sprintf("%s/data/BLAST_results.parsed.txt", indir)
mapping_fn <- sprintf("%s/mapping/CY_Mapping_with_metadata.080817.txt", indir)
out_txt <- sprintf("%s/phyloseq/phyloseq_output.%s.%s.txt", indir, "PROSPEC", format(Sys.Date(), "%m%d%y"))
out_pdf <- sprintf("%s/phyloseq/phyloseq_output.%s.%s.pdf", indir, "PROSPEC", format(Sys.Date(), "%m%d%y"))

## load all samples run to first examine controls
load(dada2_fn)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
mapping <- read.table(mapping_fn, header=T, sep="\t", comment.char="", row.names=1, as.is=T)
sel <- intersect(rownames(mapping), rownames(seqtab.nochim))
mapping <- mapping[sel,]
seqtab.nochim <- seqtab.nochim[sel,]
mapping$NumReadsOTUTable <- rowSums(seqtab.nochim)[rownames(mapping)]
mapping$SampleIDstr <- sprintf("%s (%d)", rownames(mapping), rowSums(seqtab.nochim)[rownames(mapping)])

### load BLAST taxa
#taxa <- read.table(blast_fn, header=F, as.is=T, sep="\t", row.names=1, quote="")
#inds_with_taxa <- as.numeric(gsub("seq", "", rownames(taxa))) # some sequences had no BLAST hits, so exclude these from the phyloseq object
#seqtab.nochim <- seqtab.nochim[, inds_with_taxa]
#blast_taxa <- do.call(rbind, lapply(taxa$V2, function(x) unlist(stri_split_fixed(x, ";"))))
#rownames(blast_taxa) <- colnames(seqtab.nochim); colnames(blast_taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#ps <- phyloseq(otu_table(t(seqtab.nochim), taxa_are_rows=TRUE), sample_data(mapping), tax_table(blast_taxa))
ps <- phyloseq(otu_table(t(seqtab.nochim), taxa_are_rows=TRUE), sample_data(mapping), tax_table(taxa))
set.seed(prod(dim(seqtab.nochim)))

color_table <- read.table(sprintf("%s/taxa_coloring.Jeff_022616.txt", indir), header=T, as.is=T, sep="\t", comment.char="")
color_table$Family <- gsub("f__", "", color_table$Family)
coloring <- color_table$Color
names(coloring) <- color_table$Family
ordering <- rev(names(coloring))
cols <- colorRampPalette(c("white", "red"), space = "rgb")
coloring.group <- c("black", "blue", "green"); names(coloring.group) <- c("normal", "fibrosis", "steatosis")
coloring.sampletype <- brewer.pal(4, "Set1"); names(coloring.sampletype) <- c("BacterialMock", "NegativeControl", "PCRBlank", "Stool")
coloring.bmkamount <- colorRampPalette(c("gray70", "black"), space = "rgb")(7)
coloring.bmkamount2 <- colorRampPalette(c("white", "black"), space = "rgb")(7)

ordering.family <- color_table
ordering.family$Class <- factor(ordering.family$Class, levels=unique(ordering.family$Class))
inds=order(ordering.family$Class, ordering.family$Family)
ordering.family <- ordering.family[inds,]
ordering.family$Family <- factor(ordering.family$Family, levels=unique(ordering.family$Family))

pdf(out_pdf, width=12)

##################################################################################
## trim to PROSPEC
ps <- subset_samples(ps, Project=="PROSPEC")
taxa_to_remove <- names(which(taxa_sums(ps)==0))
ps <- prune_taxa(setdiff(taxa_names(ps), taxa_to_remove), ps)
## remove empty samples
empty_samples <- names(which(sample_sums(ps)==0))
ps <- prune_samples(setdiff(sample_names(ps), empty_samples), ps)
mapping <- as(sample_data(ps), "data.frame")

## control samples
ps.sel <- subset_samples(ps, SampleType %in% c("NegativeControl", "PCRBlank"))
mapping.sel <- as(sample_data(ps.sel), "data.frame")
otu.filt <- normalizeByCols(as.data.frame(otu_table(ps.sel)))
otu.filt$Family <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level="Family")
agg <- aggregate(. ~ Family, otu.filt, sum)
families <- agg$Family
agg <- agg[,-1]
agg <- normalizeByCols(agg)
families[which(rowMeans(agg)<0.01)] <- "Other"
agg$Family <- families
df <- melt(agg, variable.name="SampleID")
df2 <- aggregate(as.formula("value ~ Family + SampleID"), df, sum)
df2$SampleType <- mapping.sel[as.character(df2$SampleID), "SampleType"]
df2$SampleIDstr <- mapping.sel[as.character(df2$SampleID), "SampleIDstr"]
# ordered by NumReadsOTUTable
ordering <- mapping.sel[order(mapping.sel$NumReadsOTUTable, decreasing=T), "SampleIDstr"]
df2$SampleIDstr <- factor(df2$SampleIDstr, levels=ordering)
p <- ggplot(df2, aes_string(x="SampleIDstr", y="value", fill="Family", order="Family")) + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.01)) + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring) + ggtitle(sprintf("Negative.Control (ordered by NumReadsOTUTable)")) + guides(col = guide_legend(ncol = 3))
print(p)

### top 10 OTUs for each control
otu_table <- as.matrix(as.data.frame(otu_table(ps.sel)))
#for (sid in rownames(mapping.sel)[order(mapping.sel$NumReadsOTUTable, decreasing=T)]) {
#	counts <- sort(otu_table[,sid], decreasing=T)
#	names(counts) <- paste(names(counts), tax_table(ps.sel)[names(counts), "Species"], sep="\n")
#	df <- melt(counts[1:10]); df$taxstr <- names(counts)[1:10]; df$taxstr <- factor(df$taxstr, levels=df$taxstr)
#	p <- ggplot(df, aes(x=taxstr, y=value)) + geom_bar(stat="identity") + theme_classic() + ggtitle(sprintf("Major contaminants in %s (%d reads)", sid, sum(counts))) + theme(axis.text.x=element_text(size=4))
#	print(p)
#}

## distribution of read counts by blank vs. sample for each SV
negids <- rownames(subset(mapping.sel, SampleType %in% c("NegativeControl", "PCRBlank")))
rs <- rowSums(otu_table[,negids])
rs.true <- rowSums(otu_table(ps)[, setdiff(sample_names(ps), negids)])
pct_blank <- 100* (rs / (rs + rs.true))
hist(pct_blank, breaks=40)
otus_to_exclude <- names(which(pct_blank > 10))

##################################################################################
## store metadata variables
metadata_variables <- read.table(sprintf("%s/mapping/metadata_types.PROSPEC.txt", indir), header=T, as.is=T, sep="\t", row.names=1)
sel <- intersect(rownames(metadata_variables), colnames(mapping))
metadata_variables <- metadata_variables[sel,, drop=F]
mapping.sel <- mapping[rownames(sample_data(ps)), sel]
# fix column types
for (mvar in rownames(metadata_variables)) {
	if (metadata_variables[mvar, "type"] == "factor") {
		mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
		if (!(is.na(metadata_variables[mvar, "baseline"])) && metadata_variables[mvar, "baseline"] != "") {
			mapping.sel[,mvar] <- relevel(mapping.sel[,mvar], metadata_variables[mvar, "baseline"])
		}
	} else if (metadata_variables[mvar, "type"] == "numeric") {
		mapping.sel[,mvar] <- as.numeric(as.character(mapping.sel[,mvar]))
	} else if (metadata_variables[mvar, "type"] == "date") {
		mapping.sel[,mvar] <- as.Date(sprintf("%06d", mapping.sel[,mvar]), format="%m%d%y")
		mapping.sel[,mvar] <- factor(as.character(mapping.sel[,mvar]), levels=as.character(unique(sort(mapping.sel[,mvar]))))
	}
}
sample_data(ps) <- mapping.sel


##################################################################################
## un-rarefied/filtered data for overall analysis including blanks
ps.relative <- transform_sample_counts(ps, function(x) x / sum(x) )

## order by PC1 (Bray-Curtis)
ordi <- ordinate(ps.relative, method = "PCoA", distance = "bray")
ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))

## PCoA
# ordination - all samples
for (distance_metric in distance_metrics) {
	ordi <- ordinate(ps.relative, method = "PCoA", distance = distance_metric)
	p <- plot_ordination(ps.relative, ordi, "samples", color = "SampleType") + theme_classic() + ggtitle(distance_metric)
	print(p)
	p <- plot_ordination(ps.relative, ordi, "samples", color = "Group") + theme_classic() + ggtitle(distance_metric)
	print(p)
}

## overall taxa barplots
#ps.sel <- subset_samples(ps.relative, SampleType=="BacterialMock")
ps.sel <- ps.relative
otu.filt <- as.data.frame(otu_table(ps.sel))
otu.filt$Family <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level="Family")
otu.filt$Family[which(is.na(otu.filt$Family))] <- "Other"
otu.filt$Family[which(otu.filt$Family=="")] <- "Other"
agg <- aggregate(. ~ Family, otu.filt, sum)
families <- agg$Family
agg <- agg[,-1]
agg <- sweep(agg, 2, colSums(agg), "/")
families[which(rowMeans(agg)<0.01)] <- "Other"
agg$Family <- families
df <- melt(agg, variable.name="SampleID")
agg <- aggregate(value~Family+SampleID, df, sum)
agg$SampleID <- as.character(agg$SampleID)
agg$SampleIDfactor <- factor(agg$SampleID, levels=ordering.pc1)
agg$Family <- factor(agg$Family, levels=levels(ordering.family$Family))
# taxa barplot representation
df.SampleIDstr <- unique(agg[,c("SampleID", "SampleIDfactor")])
df.SampleIDstr$SampleType <- as.character(mapping[df.SampleIDstr$SampleID, "SampleType"])
p <- ggplot(agg, aes(x=SampleIDfactor, y=value, fill=Family, order=Family)) + geom_bar(stat="identity", position="stack") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=2)) + ggtitle(sprintf("%s taxa summary (L5, ordered by PC1)", "Overall")) + scale_fill_manual(values=coloring) + ylim(c(-.1, 1.01))
print(p)

## use contam-filtered, pruned phyloseq objects for all further analysis
ps.relative <- prune_taxa(setdiff(taxa_names(ps.relative), otus_to_exclude), ps.relative)
ps <- prune_taxa(setdiff(taxa_names(ps), otus_to_exclude), ps)

# plot read counts
df <- melt(sample_sums(ps)); df$SampleID <- rownames(df); df <- df[order(df$value),]; df$SampleID <- factor(df$SampleID, levels=df$SampleID)
df$SampleType <- ifelse(mapping.sel[as.character(df$SampleID), "SampleType"] %in% c("PCRBlank", "NegativeControl"), "Blank", "Sample")
p <- ggplot(df, aes(x=SampleID, y=value, fill=SampleType)) + geom_bar(stat="identity") + theme_classic() + ggtitle("Read counts after contaminant SV removal") + scale_fill_manual(values=c("red", "black"))
print(p)

# make rarefied and relative ps objects
ps.rarefied <- rarefy_even_depth(ps, sample.size=27402, rngseed=nsamples(ps))
ps <- prune_samples(sample_names(ps.rarefied), ps)
ps.relative <- transform_sample_counts(ps, function(x) x / sum(x) )

##################################################################################
# trim to PROSPEC samples w/o blanks and controls, make rarefied and relative ps objects
psPROSPEC <- subset_samples(ps, SampleType=="RectalSwab")
psPROSPEC.relative <- subset_samples(ps.relative, SampleType=="RectalSwab")
psPROSPEC.rarefied <- subset_samples(ps.rarefied, SampleType=="RectalSwab")

# impute missing metadata values
mapping.sel <- as(sample_data(psPROSPEC), "data.frame")
set.seed(sum(dim(mapping.sel)))
imputed <- missForest(mapping.sel)$ximp
#mapping.sel[which(is.na(mapping.sel$Age)), "Age"] <- mean(mapping.sel$Age, na.rm=T)
#mapping.sel[which(is.na(mapping.sel$BMI)), "BMI"] <- mean(mapping.sel$BMI, na.rm=T)
sample_data(psPROSPEC) <- imputed
sample_data(psPROSPEC.relative) <- imputed
sample_data(psPROSPEC.rarefied) <- imputed

## distance matrices
dm <- list()
dm[["all"]] <- list()
dm[["PROSPEC"]] <- list()
for (distance_metric in distance_metrics) {
	dm[["all"]][[distance_metric]] <- as.matrix(distance(ps.relative, method=distance_metric))
	dm[["PROSPEC"]][[distance_metric]] <- as.matrix(distance(psPROSPEC.relative, method=distance_metric))
}

## thresholds for association/permutation tests
nsamps_threshold <- 10 # number of reads to call a sample positive
filt_threshold <- 0.1 # fraction of samples that need to be positive to keep an OTU for association testing
nperm <- 100000

## PERMANOVA
sink(out_txt, append=F)
print("PERMANOVA on PROSPEC samples")
for (distance_metric in distance_metrics) {
	print(distance_metric)
	form <- as.formula(sprintf("as.dist(dm[[\"PROSPEC\"]][[distance_metric]]) ~ %s", paste(rownames(subset(metadata_variables, useForPERMANOVA=="yes")), collapse="+")))
	res <- adonis(form , data=as(sample_data(psPROSPEC.relative), "data.frame"), permutations=999)
	print(res)
}
sink()

##########################################################################################
### Analysis of PROSPEC samples
ps.sel <- psPROSPEC.relative
mapping.sel <- as(sample_data(ps.sel), "data.frame")
mapping.sel$SampleID <- rownames(mapping.sel)
## PCoA
for (distance_metric in distance_metrics) {
	ordi <- ordinate(ps.sel, method = "PCoA", distance = distance_metric)
	p <- plot_ordination(ps.sel, ordi, "samples", color = "Group") + theme_classic() + ggtitle(distance_metric) + theme_classic()
	print(p)
}

## taxa barplots
# order by PC1 (Bray-Curtis)
ordi <- ordinate(ps.sel, method = "PCoA", distance = "bray")
ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))
otu.filt <- as.data.frame(otu_table(ps.sel))
otu.filt$Family <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level="Family")
#	otu.filt$Family[which(is.na(otu.filt$Family))] <- "NA"
agg <- aggregate(. ~ Family, otu.filt, sum)
families <- agg$Family
agg <- agg[,-1]
agg <- sweep(agg, 2, colSums(agg), "/")
families[which(rowMeans(agg)<0.01)] <- "Other"
agg$Family <- families
df <- melt(agg, variable.name="SampleID")
agg <- aggregate(value~Family+SampleID, df, sum)
agg$SampleID <- as.character(agg$SampleID)
agg$SampleIDfactor <- factor(agg$SampleID, levels=ordering.pc1)
agg$Family <- factor(agg$Family, levels=levels(ordering.family$Family))
# taxa barplot representation
df.SampleIDstr <- unique(agg[,c("SampleID", "SampleIDfactor")])
df.SampleIDstr$Group <- as.character(mapping.sel[df.SampleIDstr$SampleID, "Group"])
p <- ggplot(agg, aes(x=SampleIDfactor, y=value, fill=Family, order=Family)) + geom_bar(stat="identity", position="stack") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=4)) + ggtitle(sprintf("%s taxa summary (L5, ordered by PC1)", "PROSPEC")) + scale_fill_manual(values=coloring) + ylim(c(-.1, 1.01)) + annotate("rect", xmin = as.numeric(df.SampleIDstr$SampleIDfactor)-0.5, xmax = as.numeric(df.SampleIDstr$SampleIDfactor)+0.5, ymin = -0.04, ymax = -0.02, fill=coloring.group[df.SampleIDstr$Group])
print(p)
agg$Group <- mapping.sel[agg$SampleID, "Group"]
p <- ggplot(agg, aes(x=SampleIDfactor, y=value, fill=Family, order=Family)) + geom_bar(stat="identity", position="stack") + facet_wrap(~Group, scales="free_x") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=4)) + ggtitle(sprintf("%s taxa summary (L5, ordered by PC1)", "PROSPEC")) + scale_fill_manual(values=coloring) + ylim(c(-.1, 1.01))
print(p)

## alpha diversity by Group
adiv <- estimate_richness(psPROSPEC.rarefied, measures=alpha_metrics)
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

### taxa significance - OTU level (ZINB)
#otu.filt <- as.data.frame(otu_table(psPROSPEC.rarefied))
#agg <- otu.filt
#ftk <- names(which(rowSums(agg)>=100))
#ftk2 <- names(which(unlist(apply(agg, 1, function(x) length(which(x>=nsamps_threshold)))) > ceiling(filt_threshold*ncol(agg))))
#agg <- agg[intersect(ftk,ftk2),] # filter out less than 100 reads to reduce testing burden and improve model fit
#agg$OTU <- rownames(agg)
#out <- mclapply(agg$OTU, function(f) {
#	df <- melt(agg[f,]); colnames(df) <- c("OTU", "SampleID", "value"); df$SampleID <- as.character(df$SampleID)
#	tmp <- {}
#	print(f)
#	for (mvar in rownames(subset(metadata_variables, useForNB=="yes"))) {
#		#print(sprintf("%s %s", f, mvar))
#		df2 <- df
#		df2[, mvar] <- mapping.sel[df2$SampleID, mvar]
#		tt <- try(m <- zeroinfl(as.formula(sprintf("value ~ %s | 1", mvar)), data = df2, dist = "negbin", EM = FALSE, maxit=100), silent=T) # using EM=TRUE causes certain models to hang...
#		if (class(tt) == "zeroinfl") {
#			coef <- summary(m)$coefficients$count # rows are [(Intercept), comparisons, Log(theta)], columns are [Estimate, SE, Z, pval]
#			tmp <- rbind(tmp, cbind(f, paste(tax_table(psPROSPEC.rarefied)[f,1:7], collapse=","), mvar, class(mapping.sel[,mvar]), "ZINB", rownames(coef), coef))
#		} else if (class(tt) == "try-error") {
#			tt <- try(m <- glm.nb(as.formula(sprintf("value ~ %s", mvar)), data = df2), silent=T)
#			if (class(tt)[1] == "negbin") { 
#				coef <- summary(m)$coefficients # rows are [(Intercept), comparisons], columns are [Estimate, SE, Z, pval]
#				tmp <- rbind(tmp, cbind(f, paste(tax_table(psPROSPEC.rarefied)[f,1:7], collapse=","), mvar, class(mapping.sel[,mvar]), "NB", rownames(coef), coef))
#			}
#		}
#		#print(sprintf("finished %s %s", f, mvar))
#	}
#	print(sprintf("finished %s", f))
#	tmp
#}, mc.cores=16)
#res <- as.data.frame(do.call(rbind, out))
#colnames(res) <- c("OTU", "taxonomy", "metadata_variable", "class", "method", "coefficient", "Estimate", "SE", "Z", "pvalue")
#res <- subset(res, !coefficient %in% c("(Intercept)", "Log(theta)")) # remove intercept and log-theta terms before FDR correction
#res$pvalue <- as.numeric(as.character(res$pvalue))
#res$padj <- p.adjust(res$pvalue, method="fdr")
#res <- res[order(res$padj, decreasing=F),]
#resSig <- subset(res, padj<0.05)
#write.table(res, file=sprintf("%s/phyloseq/taxa_significance.L8.txt", indir), quote=F, sep="\t", row.names=F, col.names=T)
#write.table(resSig, file=sprintf("%s/phyloseq/taxa_significance.L8_sighits.txt", indir), quote=F, sep="\t", row.names=F, col.names=T)

## taxa significance - Species level (ZINB)
otu.filt <- as.data.frame(otu_table(psPROSPEC.rarefied))
otu.filt$Species <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psPROSPEC.rarefied), level="Species")
agg <- aggregate(. ~ Species, otu.filt, sum)
species <- agg$Species
agg <- agg[,-1]
rownames(agg) <- species
ftk <- names(which(rowSums(agg)>=100))
ftk2 <- names(which(unlist(apply(agg, 1, function(x) length(which(x>=nsamps_threshold)))) > ceiling(filt_threshold*ncol(agg))))
agg <- agg[intersect(ftk,ftk2),] # filter out less than 100 reads to reduce testing burden and improve model fit
# zero-inflated negative binomial or negative binomial regression
agg$Species <- rownames(agg)
res2 <- {}
out <- mclapply(agg$Species, function(f) {
	df <- melt(agg[f,]); colnames(df) <- c("Species", "SampleID", "value"); df$SampleID <- as.character(df$SampleID)
	tmp <- {}
	print(f)
	nsamps_detected <- length(which(df$value>=nsamps_threshold))
	for (mvar in rownames(subset(metadata_variables, useForNB=="yes"))) {
		#print(sprintf("%s %s", f, mvar))
		df2 <- df
		df2[, mvar] <- mapping.sel[df2$SampleID, mvar]
		for (confounder in rownames(subset(metadata_variables, useAsConfounder=="yes"))) {
			df2[, confounder] <- mapping.sel[df2$SampleID, confounder]
		}
		tt <- try(m <- zeroinfl(as.formula(sprintf("value ~ %s + %s | 1", mvar, paste(rownames(subset(metadata_variables, useAsConfounder=="yes")), collapse="+"))), data = df2, dist = "negbin", EM = FALSE, maxit=100), silent=T) # using EM=TRUE causes certain models to hang...
		if (class(tt) == "zeroinfl") {
			coef <- summary(m)$coefficients$count # rows are [(Intercept), comparisons, Log(theta)], columns are [Estimate, SE, Z, pval]
			tmp <- rbind(tmp, cbind(f, f, nsamps_detected, mvar, class(mapping.sel[,mvar]), "ZINB", rownames(coef), coef))
		} else if (class(tt) == "try-error") {
			tt <- try(m <- glm.nb(as.formula(sprintf("value ~ %s + %s", mvar, paste(rownames(subset(metadata_variables, useAsConfounder=="yes")), collapse="+"))), data = df2), silent=T)
			if (class(tt)[1] == "negbin") { 
				coef <- summary(m)$coefficients # rows are [(Intercept), comparisons], columns are [Estimate, SE, Z, pval]
				tmp <- rbind(tmp, cbind(f, f, nsamps_detected, mvar, class(mapping.sel[,mvar]), "NB", rownames(coef), coef))
			}
		}
		#print(sprintf("finished %s %s", f, mvar))
	}
	print(sprintf("finished %s", f))
	tmp
}, mc.cores=16)
res <- as.data.frame(do.call(rbind, out))
colnames(res) <- c("OTU", "taxonomy", "nsamps_detected", "metadata_variable", "class", "method", "coefficient", "Estimate", "SE", "Z", "pvalue")
# remove intercept, log-theta, and confounder terms before FDR correction
to_remove <- colSums(do.call(rbind, lapply(c("(Intercept)", "Log(theta)", rownames(subset(metadata_variables, useAsConfounder=="yes"))), function(conf) {grepl(conf, res$coefficient, fixed=T)}))) > 0
res <- res[!to_remove,]
res$pvalue <- as.numeric(as.character(res$pvalue))
# p-value adjust
res$padj <- p.adjust(res$pvalue, method="fdr")
res <- res[order(res$padj, decreasing=F),]
write.table(res, file=sprintf("%s/phyloseq/taxa_significance.L7.txt", indir), quote=F, sep="\t", row.names=F, col.names=T)

## taxa significance - Genus level (ZINB)
otu.filt <- as.data.frame(otu_table(psPROSPEC.rarefied))
otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psPROSPEC.rarefied), level="Genus")
agg <- aggregate(. ~ Genus, otu.filt, sum)
genera <- agg$Genus
agg <- agg[,-1]
rownames(agg) <- genera
ftk <- names(which(rowSums(agg)>=100))
ftk2 <- names(which(unlist(apply(agg, 1, function(x) length(which(x>=nsamps_threshold)))) > ceiling(filt_threshold*ncol(agg))))
agg <- agg[intersect(ftk,ftk2),] # filter out less than 100 reads to reduce testing burden and improve model fit
# zero-inflated negative binomial or negative binomial regression
agg$Genus <- rownames(agg)
res2 <- {}
out <- mclapply(agg$Genus, function(f) {
	df <- melt(agg[f,]); colnames(df) <- c("Genus", "SampleID", "value"); df$SampleID <- as.character(df$SampleID)
	tmp <- {}
	print(f)
	nsamps_detected <- length(which(df$value>=nsamps_threshold))
	for (mvar in rownames(subset(metadata_variables, useForNB=="yes"))) {
		#print(sprintf("%s %s", f, mvar))
		df2 <- df
		df2[, mvar] <- mapping.sel[df2$SampleID, mvar]
		for (confounder in rownames(subset(metadata_variables, useAsConfounder=="yes"))) {
			df2[, confounder] <- mapping.sel[df2$SampleID, confounder]
		}
		tt <- try(m <- zeroinfl(as.formula(sprintf("value ~ %s + %s | 1", mvar, paste(rownames(subset(metadata_variables, useAsConfounder=="yes")), collapse="+"))), data = df2, dist = "negbin", EM = FALSE, maxit=100), silent=T) # using EM=TRUE causes certain models to hang...
		if (class(tt) == "zeroinfl") {
			coef <- summary(m)$coefficients$count # rows are [(Intercept), comparisons, Log(theta)], columns are [Estimate, SE, Z, pval]
			tmp <- rbind(tmp, cbind(f, f, nsamps_detected, mvar, class(mapping.sel[,mvar]), "ZINB", rownames(coef), coef))
		} else if (class(tt) == "try-error") {
			tt <- try(m <- glm.nb(as.formula(sprintf("value ~ %s + %s", mvar, paste(rownames(subset(metadata_variables, useAsConfounder=="yes")), collapse="+"))), data = df2), silent=T)
			if (class(tt)[1] == "negbin") { 
				coef <- summary(m)$coefficients # rows are [(Intercept), comparisons], columns are [Estimate, SE, Z, pval]
				tmp <- rbind(tmp, cbind(f, f, nsamps_detected, mvar, class(mapping.sel[,mvar]), "NB", rownames(coef), coef))
			}
		}
		#print(sprintf("finished %s %s", f, mvar))
	}
	print(sprintf("finished %s", f))
	tmp
}, mc.cores=16)
res <- as.data.frame(do.call(rbind, out))
colnames(res) <- c("OTU", "taxonomy", "nsamps_detected", "metadata_variable", "class", "method", "coefficient", "Estimate", "SE", "Z", "pvalue")
# remove intercept, log-theta, and confounder terms before FDR correction
to_remove <- colSums(do.call(rbind, lapply(c("(Intercept)", "Log(theta)", rownames(subset(metadata_variables, useAsConfounder=="yes"))), function(conf) {grepl(conf, res$coefficient, fixed=T)}))) > 0
res <- res[!to_remove,]
res$pvalue <- as.numeric(as.character(res$pvalue))
# p-value adjust
res$padj <- p.adjust(res$pvalue, method="fdr")
res <- res[order(res$padj, decreasing=F),]
write.table(res, file=sprintf("%s/phyloseq/taxa_significance.L6.txt", indir), quote=F, sep="\t", row.names=F, col.names=T)

## taxa significance - Family level (ZINB)
otu.filt <- as.data.frame(otu_table(psPROSPEC.rarefied))
otu.filt$Family <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psPROSPEC.rarefied), level="Family")
agg <- aggregate(. ~ Family, otu.filt, sum)
families <- agg$Family
agg <- agg[,-1]
rownames(agg) <- families
ftk <- names(which(rowSums(agg)>=100))
ftk2 <- names(which(unlist(apply(agg, 1, function(x) length(which(x>=nsamps_threshold)))) > ceiling(filt_threshold*ncol(agg))))
agg <- agg[intersect(ftk,ftk2),] # filter out less than 100 reads to reduce testing burden and improve model fit
# zero-inflated negative binomial or negative binomial regression
agg$Family <- rownames(agg)
res2 <- {}
out <- mclapply(agg$Family, function(f) {
	df <- melt(agg[f,]); colnames(df) <- c("Family", "SampleID", "value"); df$SampleID <- as.character(df$SampleID)
	tmp <- {}
	print(f)
	nsamps_detected <- length(which(df$value>=nsamps_threshold))
	for (mvar in rownames(subset(metadata_variables, useForNB=="yes"))) {
		#print(sprintf("%s %s", f, mvar))
		df2 <- df
		df2[, mvar] <- mapping.sel[df2$SampleID, mvar]
		for (confounder in rownames(subset(metadata_variables, useAsConfounder=="yes"))) {
			df2[, confounder] <- mapping.sel[df2$SampleID, confounder]
		}
		tt <- try(m <- zeroinfl(as.formula(sprintf("value ~ %s + %s | 1", mvar, paste(rownames(subset(metadata_variables, useAsConfounder=="yes")), collapse="+"))), data = df2, dist = "negbin", EM = FALSE, maxit=100), silent=T) # using EM=TRUE causes certain models to hang...
		if (class(tt) == "zeroinfl") {
			coef <- summary(m)$coefficients$count # rows are [(Intercept), comparisons, Log(theta)], columns are [Estimate, SE, Z, pval]
			tmp <- rbind(tmp, cbind(f, f, nsamps_detected, mvar, class(mapping.sel[,mvar]), "ZINB", rownames(coef), coef))
		} else if (class(tt) == "try-error") {
			tt <- try(m <- glm.nb(as.formula(sprintf("value ~ %s + %s", mvar, paste(rownames(subset(metadata_variables, useAsConfounder=="yes")), collapse="+"))), data = df2), silent=T)
			if (class(tt)[1] == "negbin") { 
				coef <- summary(m)$coefficients # rows are [(Intercept), comparisons], columns are [Estimate, SE, Z, pval]
				tmp <- rbind(tmp, cbind(f, f, nsamps_detected, mvar, class(mapping.sel[,mvar]), "NB", rownames(coef), coef))
			}
		}
		#print(sprintf("finished %s %s", f, mvar))
	}
	print(sprintf("finished %s", f))
	tmp
}, mc.cores=16)
res <- as.data.frame(do.call(rbind, out))
colnames(res) <- c("OTU", "taxonomy", "nsamps_detected", "metadata_variable", "class", "method", "coefficient", "Estimate", "SE", "Z", "pvalue")
# remove intercept, log-theta, and confounder terms before FDR correction
to_remove <- colSums(do.call(rbind, lapply(c("(Intercept)", "Log(theta)", rownames(subset(metadata_variables, useAsConfounder=="yes"))), function(conf) {grepl(conf, res$coefficient, fixed=T)}))) > 0
res <- res[!to_remove,]
res$pvalue <- as.numeric(as.character(res$pvalue))
# p-value adjust
res$padj <- p.adjust(res$pvalue, method="fdr")
res <- res[order(res$padj, decreasing=F),]
write.table(res, file=sprintf("%s/phyloseq/taxa_significance.L5.txt", indir), quote=F, sep="\t", row.names=F, col.names=T)

taxlevel_lookup <- c("L5", "L6", "L7")
names(taxlevel_lookup) <- c("Family", "Genus", "Species")

## ZINB results
for (level in c("Family", "Genus")) {
	taxlevel <- taxlevel_lookup[level]
	res <- read.table(sprintf("%s/phyloseq/taxa_significance.%s.txt", indir, taxlevel), header=T, sep="\t", as.is=T, quote="")
	res$Estimate <- as.numeric(res$Estimate); res$EstimateWeight <- ifelse(is.na(res$Estimate), 0.5, res$Estimate)
	res$EstimateColor <- ifelse(res$Estimate>0, "red", "blue"); res$EstimateColor[which(is.na(res$EstimateColor))] <- "black"
	resSig <- subset(res, padj<0.1)
#	df <- dcast(resSig, OTU ~ coefficient, value.var="Z"); rownames(df) <- df$OTU; df <- df[,-1,drop=F]; df <- t(as.matrix(df))
#	df2 <- df; df2[which(is.na(df2))] <- 0
#	heatmap.2(df, Rowv=F, Colv=F, dendrogram="none", trace="none", col=colorRampPalette(c("blue", "white","red"))(150), margins=c(10,10), cexCol=0.8, cexRow=0.8, sepwidth=c(0.02,0.02), colsep=1:ncol(df), rowsep=1:nrow(df), sepcolor="grey")
#	heatmap.2(df2, Rowv=T, Colv=F, dendrogram="row", trace="none", col=colorRampPalette(c("blue", "white","red"))(150), margins=c(10,10), cexCol=0.8, cexRow=0.8, sepwidth=c(0.02,0.02), colsep=1:ncol(df), rowsep=1:nrow(df), sepcolor="grey")
#	heatmap.2(df2, Rowv=F, Colv=T, dendrogram="column", trace="none", col=colorRampPalette(c("blue", "white","red"))(150), margins=c(10,10), cexCol=0.8, cexRow=0.8, sepwidth=c(0.02,0.02), colsep=1:ncol(df), rowsep=1:nrow(df), sepcolor="grey")
#	heatmap.2(df2, Rowv=T, Colv=T, dendrogram="both", trace="none", col=colorRampPalette(c("blue", "white","red"))(150), margins=c(10,10), cexCol=0.8, cexRow=0.8, sepwidth=c(0.02,0.02), colsep=1:ncol(df), rowsep=1:nrow(df), sepcolor="grey")
	
	# also draw violin plots of distributions
	ps.sel <- psPROSPEC.relative
	otu.filt <- as.data.frame(otu_table(ps.sel))
  otu.filt.abundance <- normalizeByCols(otu.filt)
  otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level=level)
  otu.filt.abundance[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level=level)
  # rename Prevotella_6, etc -> Prevotella
  otu.filt[[level]] <- gsub("_\\d$", "", otu.filt[[level]])
  otu.filt.abundance[[level]] <- gsub("_\\d$", "", otu.filt.abundance[[level]])
  agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
  agg.abundance <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt.abundance, sum)
  lvl <- agg[[level]]
  lvl.abundance <- agg.abundance[[level]]
  agg <- agg[,-1]
  agg.abundance <- agg.abundance[,-1]
  rownames(agg) <- lvl
  rownames(agg.abundance) <- lvl.abundance
  agg <- agg[rownames(agg.abundance),]
  agg.melt <- melt(as.matrix(agg), as.is=T); colnames(agg.melt) <- c("taxa", "SampleID", "value")
	agg.melt$Group <- droplevels(mapping.sel[agg.melt$SampleID, "Group"])

	df <- subset(agg.melt, taxa %in% resSig$OTU)
	df$taxa <- factor(df$taxa, levels=unique(resSig$OTU))
	p <- ggplot(df, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap(~taxa, scales="free", ncol=3, drop=F) + theme_classic() + ggtitle(sprintf("Rel. abund. of ZINB taxa (%s)", level)) + coord_flip()
	print(p)
}

## randomForest classification of Group (one-vs-one setup, X vs Normal baseline)
set.seed(011818)
ps.sel <- psPROSPEC.relative
mapping.sel <- as(sample_data(ps.sel), "data.frame")
for (level in c("Species", "Genus")){
  otu.filt <- as.data.frame(otu_table(ps.sel))
  otu.filt.abundance <- normalizeByCols(otu.filt)
  otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level=level)
  otu.filt.abundance[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level=level)
  # rename Prevotella_6, etc -> Prevotella
  otu.filt[[level]] <- gsub("_\\d$", "", otu.filt[[level]])
  otu.filt.abundance[[level]] <- gsub("_\\d$", "", otu.filt.abundance[[level]])
  agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
  agg.abundance <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt.abundance, sum)
  lvl <- agg[[level]]
  lvl.abundance <- agg.abundance[[level]]
  agg <- agg[,-1]
  agg.abundance <- agg.abundance[,-1]
  rownames(agg) <- lvl
  rownames(agg.abundance) <- lvl.abundance
  agg.abundance <- agg.abundance[which(rowSums(agg.abundance > 0.01) >= ceiling(ncol(agg.abundance)/10)),]
  agg <- agg[rownames(agg.abundance),]
  agg.melt.stored <- melt(as.matrix(agg), as.is=T); colnames(agg.melt.stored) <- c("taxa", "SampleID", "value") # used later for violin plots
  res.mean <- {}; res.sd <- {}
	for (mvar_level in setdiff(levels(mapping.sel[, "Group"]), "normal")) {
		sel <- sample_names(subset_samples(ps.relative, Group %in% c(mvar_level, "normal")))
		data.sel <- as.data.frame(t(agg[,sel]))
		for (confounder in rownames(subset(metadata_variables, useAsConfounder=="yes"))) {
			data.sel[, confounder] <- unlist(sample_data(ps.sel)[rownames(data.sel), confounder])
		}
		agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), rownames(subset(metadata_variables, useAsConfounder=="yes" & type=="factor")))]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "taxa", "value")
		response <- droplevels(mapping.sel[sel, "Group"]); names(response) <- sel

#		# after running for the first time, COMMENT OUT THIS BLOCK ##
#		num_iter <- 100
#		ncores <- 20
#		out <- mclapply(1:num_iter, function (dummy) {
#				importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
#			}, mc.cores=ncores )
#		collated.importance <- do.call(cbind, out)
#		out <- mclapply(1:num_iter, function (dummy) {
#				rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.75)$error.cv
#			}, mc.cores=ncores )
#		collated.cv <- do.call(cbind, out)

#		write.table(collated.importance, file=sprintf("%s/phyloseq/randomForest.%s.%s.importance.txt", indir, mvar_level, level), quote=F, sep="\t", row.names=T, col.names=F)
#		write.table(collated.cv, file=sprintf("%s/phyloseq/randomForest.%s.%s.cv.txt", indir, mvar_level, level), quote=F, sep="\t", row.names=T, col.names=F)
#		## END BLOCK TO COMMENT ##
		
		collated.importance <- read.table(sprintf("%s/phyloseq/randomForest.%s.%s.importance.txt", indir, mvar_level, level), header=F, as.is=T, sep="\t", row.names=1)
		collated.cv <- read.table(sprintf("%s/phyloseq/randomForest.%s.%s.cv.txt", indir, mvar_level, level), header=F, as.is=T, sep="\t", row.names=1)
		importance.mean <- rowMeans(collated.importance)
		importance.sd <- unlist(apply(collated.importance, 1, sd))
		cv.mean <- rowMeans(collated.cv)
		cv.sd <- unlist(apply(collated.cv, 1, sd))
		res.mean <- rbind(res.mean, importance.mean)
		res.sd <- rbind(res.sd, importance.sd)
		inds <- order(importance.mean, decreasing=T)
#		inds <- inds[1:as.numeric(names(which.min(cv.mean)))] # edit as appropriate
		inds <- inds[which(importance.mean[inds]>0.001)]
		
		## after running for the first time, COMMENT OUT THIS BLOCK ##
		# using a sparse model with N predictors
#		sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
#		save(sparseRF, file=sprintf("%s/phyloseq/randomForest.%s.%s.model", indir, mvar_level, level))
		load(sprintf("%s/phyloseq/randomForest.%s.%s.model", indir, mvar_level, level))
		# accuracy of final sparseRF model
		pred <- predict(sparseRF, type="prob")
		pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
		confusion_matrix <- table(pred_df[, c("true", "predicted")])
		accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
		vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
		mccvalue <- mcc(vec.pred, vec.true)
		df <- confusion_matrix
		p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s) (accuracy = %.2f%%, MCC = %.4f)", mvar_level, level, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
		print(p)
		write.table(confusion_matrix, file=sprintf("%s/phyloseq/randomForest.%s.%s.confusion_matrix.txt", indir, mvar_level, level), quote=F, sep="\t", row.names=T, col.names=T)
		# ROC analysis
		if (nlevels(response)==2) {
			pred <- predict(sparseRF, type="prob")
			pred2 <- prediction(pred[,2], ordered(response))
			perf <- performance(pred2, "tpr", "fpr")
			perf.auc <- performance(pred2, "auc")
			print(plot(perf, main=sprintf("ROC %s %s (sparseRF final model)", mvar_level, level)) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g", unlist(perf.auc@y.values))))
		}
		## END BLOCK TO COMMENT ##
		
		load(sprintf("%s/phyloseq/randomForest.%s.%s.model", indir, mvar_level, level))
		# plotting - per-group sparse model
		df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
		colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
		print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s", mvar_level, level)))
		
		# plotting - per-group variables
		df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
		res <- read.table(sprintf("%s/phyloseq/taxa_significance.%s.txt", indir, taxlevel_lookup[level]), header=T, sep="\t", as.is=T, quote="")
		res <- subset(res, coefficient == paste("Group", mvar_level, sep="")); rownames(res) <- res$OTU
		df$Z <- res[as.character(df$OTU), "Z"]; df$padj <- res[as.character(df$OTU), "padj"]
		df$OTU_string <- sprintf("%s (Z=%.4g, padj=%.4g)", df$OTU, df$Z, df$padj)
		df$sig <- ifelse(df$padj < 0.1, "sig", "ns"); df$sig[which(is.na(df$sig))] <- "ns"
		print(ggplot(df, aes(x=OTU, y=importance, label=OTU)) + geom_bar(position=position_dodge(), stat="identity", fill="#999999") + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string, color=sig), size=3, hjust=0) + ggtitle(sprintf("%s explanatory %s", mvar_level, level)) + scale_color_manual(values=cols.sig))
		# violin plots of relative abundance
		agg.melt <- agg.melt.stored
		agg.melt$Group <- droplevels(mapping.sel[agg.melt$SampleID, "Group"])
		agg.melt <- subset(agg.melt, taxa %in% rownames(df) & Group %in% c(mvar_level, "normal"))
		agg.melt$taxa <- factor(agg.melt$taxa, levels=rownames(df))
		p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap(~taxa, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF taxa (%s, %s)", mvar_level, level)) + coord_flip()
		print(p)
	}
	# plotting - all groups importance values
	rownames(res.mean) <- setdiff(levels(mapping.sel[, "Group"]), "normal"); rownames(res.sd) <- rownames(res.mean)
	df <- res.mean; df[which(df<0.001, arr.ind=T)] <- 0
	heatmap.2(df, Rowv=T, Colv=T, dendrogram="both", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(10,10), cexCol=1, cexRow=1, main=sprintf("RF importance values (%s)", level))
	
}

## randomForest classification of Group (multiclass)
set.seed(011918)
ps.sel <- psPROSPEC.relative
mapping.sel <- as(sample_data(ps.sel), "data.frame")
for (level in c("Species", "Genus")){
  otu.filt <- as.data.frame(otu_table(ps.sel))
  otu.filt.abundance <- normalizeByCols(otu.filt)
  otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level=level)
  otu.filt.abundance[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level=level)
  # rename Prevotella_6, etc -> Prevotella
  otu.filt[[level]] <- gsub("_\\d$", "", otu.filt[[level]])
  otu.filt.abundance[[level]] <- gsub("_\\d$", "", otu.filt.abundance[[level]])
  agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
  agg.abundance <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt.abundance, sum)
  lvl <- agg[[level]]
  lvl.abundance <- agg.abundance[[level]]
  agg <- agg[,-1]
  agg.abundance <- agg.abundance[,-1]
  rownames(agg) <- lvl
  rownames(agg.abundance) <- lvl.abundance
  agg.abundance <- agg.abundance[which(rowSums(agg.abundance > 0.01) >= ceiling(ncol(agg.abundance)/10)),]
  agg <- agg[rownames(agg.abundance),]
  agg.melt.stored <- melt(as.matrix(agg), as.is=T); colnames(agg.melt.stored) <- c("taxa", "SampleID", "value") # used later for violin plots
  res.mean <- {}; res.sd <- {}
	data.sel <- as.data.frame(t(agg))
	for (confounder in rownames(subset(metadata_variables, useAsConfounder=="yes"))) {
		data.sel[, confounder] <- unlist(sample_data(ps.sel)[rownames(data.sel), confounder])
	}
	agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), rownames(subset(metadata_variables, useAsConfounder=="yes" & type=="factor")))]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "taxa", "value")
	response <- droplevels(mapping.sel[, "Group"]); names(response) <- rownames(mapping.sel)

#	# after running for the first time, COMMENT OUT THIS BLOCK ##
#	num_iter <- 100
#	ncores <- 20
#	out <- mclapply(1:num_iter, function (dummy) {
#			importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
#		}, mc.cores=ncores )
#	collated.importance <- do.call(cbind, out)
#	out <- mclapply(1:num_iter, function (dummy) {
#			rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.75)$error.cv
#		}, mc.cores=ncores )
#	collated.cv <- do.call(cbind, out)

#	write.table(collated.importance, file=sprintf("%s/phyloseq/randomForest.%s.%s.importance.txt", indir, "multiclass", level), quote=F, sep="\t", row.names=T, col.names=F)
#	write.table(collated.cv, file=sprintf("%s/phyloseq/randomForest.%s.%s.cv.txt", indir, "multiclass", level), quote=F, sep="\t", row.names=T, col.names=F)
#	## END BLOCK TO COMMENT ##
	
	collated.importance <- read.table(sprintf("%s/phyloseq/randomForest.%s.%s.importance.txt", indir, "multiclass", level), header=F, as.is=T, sep="\t", row.names=1)
	collated.cv <- read.table(sprintf("%s/phyloseq/randomForest.%s.%s.cv.txt", indir, "multiclass", level), header=F, as.is=T, sep="\t", row.names=1)
	importance.mean <- rowMeans(collated.importance)
	importance.sd <- unlist(apply(collated.importance, 1, sd))
	cv.mean <- rowMeans(collated.cv)
	cv.sd <- unlist(apply(collated.cv, 1, sd))
	res.mean <- rbind(res.mean, importance.mean)
	res.sd <- rbind(res.sd, importance.sd)
	inds <- order(importance.mean, decreasing=T)
#		inds <- inds[1:as.numeric(names(which.min(cv.mean)))] # edit as appropriate
	inds <- inds[which(importance.mean[inds]>0.001)]
	
	## after running for the first time, COMMENT OUT THIS BLOCK ##
	# using a sparse model with N predictors
#	sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
#	save(sparseRF, file=sprintf("%s/phyloseq/randomForest.%s.%s.model", indir, "multiclass", level))
	load(sprintf("%s/phyloseq/randomForest.%s.%s.model", indir, "multiclass", level))
	# accuracy of final sparseRF model
	pred <- predict(sparseRF, type="prob")
	pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
	confusion_matrix <- table(pred_df[, c("true", "predicted")])
	class_errors <- unlist(lapply(levels(mapping.sel$Group), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$Group)
	accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
	vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
	mccvalue <- mcc(vec.pred, vec.true)
	df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
	p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s) (accuracy = %.2f%%, MCC = %.4f)", "multiclass", level, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
	print(p)
	write.table(confusion_matrix, file=sprintf("%s/phyloseq/randomForest.%s.%s.confusion_matrix.txt", indir, "multiclass", level), quote=F, sep="\t", row.names=T, col.names=T)
	# ROC analysis
	if (nlevels(response)==2) {
		pred <- predict(sparseRF, type="prob")
		pred2 <- prediction(pred[,2], ordered(response))
		perf <- performance(pred2, "tpr", "fpr")
		perf.auc <- performance(pred2, "auc")
		print(plot(perf, main=sprintf("ROC %s %s (sparseRF final model)", "multiclass", level)) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g", unlist(perf.auc@y.values))))
	}
	## END BLOCK TO COMMENT ##
	
	load(sprintf("%s/phyloseq/randomForest.%s.%s.model", indir, "multiclass", level))
	# plotting - per-group sparse model
	df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
	colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
	print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s", "multiclass", level)))
	
	# plotting - per-group variables
	df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
	res <- read.table(sprintf("%s/phyloseq/taxa_significance.%s.txt", indir, taxlevel_lookup[level]), header=T, sep="\t", as.is=T, quote="")
	res <- res[grep("Group", res$coefficient),]
	res <- aggregate(padj ~ OTU, res, min); rownames(res) <- res$OTU
	df$padj <- res[as.character(df$OTU), "padj"]
	df$OTU_string <- sprintf("%s (min padj=%.4g)", df$OTU, df$padj)
	df$sig <- ifelse(df$padj < 0.1, "sig", "ns"); df$sig[which(is.na(df$sig))] <- "ns"
	print(ggplot(df, aes(x=OTU, y=importance, label=OTU)) + geom_bar(position=position_dodge(), stat="identity", fill="#999999") + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string, color=sig), size=3, hjust=0) + ggtitle(sprintf("%s explanatory %s", "multiclass", level)) + scale_color_manual(values=cols.sig))
	# violin plots of relative abundance
	agg.melt <- agg.melt.stored
	agg.melt$Group <- droplevels(mapping.sel[agg.melt$SampleID, "Group"])
	agg.melt <- subset(agg.melt, taxa %in% rownames(df))
	agg.melt$taxa <- factor(agg.melt$taxa, levels=rownames(df))
	p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap(~taxa, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF taxa (%s, %s)", "multiclass", level)) + coord_flip()
	print(p)
	
}





dev.off()




