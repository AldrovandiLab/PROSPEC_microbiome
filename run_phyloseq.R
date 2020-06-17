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
library(stringi)
library(car)
library(randomForest)
library(ROCR)
library(Hmisc)
library(missForest)

source("utils.R")
source("mcc.R")
indir <- "data"
output_dir <- "output"

distance_metrics <- c("bray", "jaccard", "jsd")
alpha_metrics <- c("Chao1", "Shannon", "Simpson", "Observed")
icc_levels <- c(5,6,7)
names(icc_levels) <- c("Family", "Genus", "Species")
ncores <- 16
cols.sig <- c("black", "red", "grey"); names(cols.sig) <- c("ns", "sig", "NA")

#################################################################
## handoff to phyloseq
seqtab.nochim <- readRDS(sprintf("%s/merged_seqtab.nochim.rds", indir))
dada2_fn <- sprintf("%s/DADA2.RData", indir)
blast_fn <- sprintf("%s/BLAST_results.parsed.txt", indir)
mapping_fn <- sprintf("mapping/CY_Mapping_with_metadata.013019.txt")
out_txt <- sprintf("%s/phyloseq_output.%s.%s.txt", output_dir, "PROSPEC", format(Sys.Date(), "%m%d%y"))
out_pdf <- sprintf("%s/phyloseq_output.%s.%s.pdf", output_dir, "PROSPEC", format(Sys.Date(), "%m%d%y"))

## load all samples run to first examine controls
load(dada2_fn)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
mapping <- read.table(mapping_fn, header=T, sep="\t", comment.char="", row.names=1, as.is=T, quote="")
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

color_table <- read.table("taxa_coloring.Jeff_022616.txt", header=T, as.is=T, sep="\t", comment.char="")
color_table$Family <- gsub("f__", "", color_table$Family)
coloring <- color_table$Color
names(coloring) <- color_table$Family
ordering <- rev(names(coloring))
cols <- colorRampPalette(c("white", "red"), space = "rgb")
coloring.group <- c("black", "blue", "green"); names(coloring.group) <- c("normal", "fibrosis", "steatosis")
coloring.sampletype <- brewer.pal(4, "Set1"); names(coloring.sampletype) <- c("BacterialMock", "NegativeControl", "PCRBlank", "Stool")
coloring.bmkamount <- colorRampPalette(c("gray70", "black"), space = "rgb")(7)
coloring.bmkamount2 <- colorRampPalette(c("white", "black"), space = "rgb")(7)
dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")


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
metadata_variables <- read.table(sprintf("%s/metadata_types.PROSPEC.txt", "mapping"), header=T, as.is=T, sep="\t", row.names=1)
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
	p <- plot_ordination(ps.relative, ordi, "samples", color = "SampleType") + theme_classic() + ggtitle(distance_metric) + stat_ellipse(type="t")
	print(p)
	p <- plot_ordination(ps.relative, ordi, "samples", color = "Group") + theme_classic() + ggtitle(distance_metric) + stat_ellipse(type="t")
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

# impute missing metadata values (but not diet data)
diet_variables <- c("total_calories", "pct_carbohydrates", "pct_protein", "pct_fat", "cholesterol", "total_fiber")
mapping.sel <- as(sample_data(psPROSPEC), "data.frame")
mapping.sel <- mapping.sel[, setdiff(colnames(mapping.sel), diet_variables)]
set.seed(sum(dim(mapping.sel)))
imputed <- missForest(mapping.sel)$ximp
imputed <- cbind(imputed, as(sample_data(psPROSPEC), "data.frame")[, diet_variables])
#mapping.sel[which(is.na(mapping.sel$Age)), "Age"] <- mean(mapping.sel$Age, na.rm=T)
#mapping.sel[which(is.na(mapping.sel$BMI)), "BMI"] <- mean(mapping.sel$BMI, na.rm=T)
sample_data(psPROSPEC) <- imputed
sample_data(psPROSPEC.relative) <- imputed
sample_data(psPROSPEC.rarefied) <- imputed
mapping.sel <- imputed

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
siglevel <- 0.05

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
	p <- plot_ordination(ps.sel, ordi, "samples", color = "Group") + theme_classic() + ggtitle(distance_metric) + theme_classic() + stat_ellipse(type="t")
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

## taxa significance - OTU level (ZINB)
otu.filt <- as.data.frame(otu_table(psPROSPEC.rarefied))
agg <- otu.filt
ftk <- names(which(rowSums(agg)>=100))
ftk2 <- names(which(unlist(apply(agg, 1, function(x) length(which(x>=nsamps_threshold)))) > ceiling(filt_threshold*ncol(agg))))
agg <- agg[intersect(ftk,ftk2),] # filter out less than 100 reads to reduce testing burden and improve model fit
agg$OTU <- rownames(agg)
out <- mclapply(agg$OTU, function(f) {
	df <- melt(agg[f,]); colnames(df) <- c("OTU", "SampleID", "value"); df$SampleID <- as.character(df$SampleID)
	tmp <- {}
	print(f)
	for (mvar in rownames(subset(metadata_variables, useForNB=="yes"))) {
		#print(sprintf("%s %s", f, mvar))
		df2 <- df
		df2[, mvar] <- mapping.sel[df2$SampleID, mvar]
		tt <- try(m <- zeroinfl(as.formula(sprintf("value ~ %s | 1", mvar)), data = df2, dist = "negbin", EM = FALSE, maxit=100), silent=T) # using EM=TRUE causes certain models to hang...
		if (class(tt) == "zeroinfl") {
			coef <- summary(m)$coefficients$count # rows are [(Intercept), comparisons, Log(theta)], columns are [Estimate, SE, Z, pval]
			tmp <- rbind(tmp, cbind(f, paste(tax_table(psPROSPEC.rarefied)[f,1:7], collapse=","), mvar, class(mapping.sel[,mvar]), "ZINB", rownames(coef), coef))
		} else if (class(tt) == "try-error") {
			tt <- try(m <- glm.nb(as.formula(sprintf("value ~ %s", mvar)), data = df2), silent=T)
			if (class(tt)[1] == "negbin") { 
				coef <- summary(m)$coefficients # rows are [(Intercept), comparisons], columns are [Estimate, SE, Z, pval]
				tmp <- rbind(tmp, cbind(f, paste(tax_table(psPROSPEC.rarefied)[f,1:7], collapse=","), mvar, class(mapping.sel[,mvar]), "NB", rownames(coef), coef))
			}
		}
		#print(sprintf("finished %s %s", f, mvar))
	}
	print(sprintf("finished %s", f))
	tmp
}, mc.cores=16)
res <- as.data.frame(do.call(rbind, out))
colnames(res) <- c("OTU", "taxonomy", "metadata_variable", "class", "method", "coefficient", "Estimate", "SE", "Z", "pvalue")
res <- subset(res, !coefficient %in% c("(Intercept)", "Log(theta)")) # remove intercept and log-theta terms before FDR correction
res$pvalue <- as.numeric(as.character(res$pvalue))
res$padj <- p.adjust(res$pvalue, method="fdr")
res <- res[order(res$padj, decreasing=F),]
resSig <- subset(res, padj<0.05)
write.table(res, file=sprintf("%s/taxa_significance.L8.txt", output_dir), quote=F, sep="\t", row.names=F, col.names=T)
write.table(resSig, file=sprintf("%s/taxa_significance.L8_sighits.txt", output_dir), quote=F, sep="\t", row.names=F, col.names=T)

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
write.table(res, file=sprintf("%s/taxa_significance.L7.txt", output_dir), quote=F, sep="\t", row.names=F, col.names=T)

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
write.table(res, file=sprintf("%s/taxa_significance.L6.txt", output_dir), quote=F, sep="\t", row.names=F, col.names=T)

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
write.table(res, file=sprintf("%s/taxa_significance.L5.txt", output_dir), quote=F, sep="\t", row.names=F, col.names=T)

taxlevel_lookup <- c("L5", "L6", "L7")
names(taxlevel_lookup) <- c("Family", "Genus", "Species")

## ZINB results
for (level in c("Family", "Genus")) {
	taxlevel <- taxlevel_lookup[level]
	res <- read.table(sprintf("%s/taxa_significance.%s.txt", output_dir, taxlevel), header=T, sep="\t", as.is=T, quote="")
	res$Estimate <- as.numeric(res$Estimate); res$EstimateWeight <- ifelse(is.na(res$Estimate), 0.5, res$Estimate)
	res$EstimateColor <- ifelse(res$Estimate>0, "red", "blue"); res$EstimateColor[which(is.na(res$EstimateColor))] <- "black"
	resSig <- subset(res, padj<0.1)
	df <- dcast(resSig, OTU ~ coefficient, value.var="Z"); rownames(df) <- df$OTU; df <- df[,-1,drop=F]; df <- t(as.matrix(df))
	df2 <- df; df2[which(is.na(df2))] <- 0
#	heatmap.2(df, Rowv=F, Colv=F, dendrogram="none", trace="none", col=colorRampPalette(c("blue", "white","red"))(150), margins=c(10,10), cexCol=0.8, cexRow=0.8, sepwidth=c(0.02,0.02), colsep=1:ncol(df), rowsep=1:nrow(df), sepcolor="grey")
#	heatmap.2(df2, Rowv=T, Colv=F, dendrogram="row", trace="none", col=colorRampPalette(c("blue", "white","red"))(150), margins=c(10,10), cexCol=0.8, cexRow=0.8, sepwidth=c(0.02,0.02), colsep=1:ncol(df), rowsep=1:nrow(df), sepcolor="grey")
#	heatmap.2(df2, Rowv=F, Colv=T, dendrogram="column", trace="none", col=colorRampPalette(c("blue", "white","red"))(150), margins=c(10,10), cexCol=0.8, cexRow=0.8, sepwidth=c(0.02,0.02), colsep=1:ncol(df), rowsep=1:nrow(df), sepcolor="grey")
	heatmap.2(df2, Rowv=T, Colv=T, dendrogram="both", trace="none", col=colorRampPalette(c("blue", "white","red"))(150), margins=c(10,10), cexCol=0.8, cexRow=0.8, sepwidth=c(0.02,0.02), colsep=1:ncol(df), rowsep=1:nrow(df), sepcolor="grey")
	
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

#		write.table(collated.importance, file=sprintf("%s/randomForest.%s.%s.importance.txt", output_dir, mvar_level, level), quote=F, sep="\t", row.names=T, col.names=F)
#		write.table(collated.cv, file=sprintf("%s/randomForest.%s.%s.cv.txt", output_dir, mvar_level, level), quote=F, sep="\t", row.names=T, col.names=F)
#		## END BLOCK TO COMMENT ##
		
		collated.importance <- read.table(sprintf("%s/randomForest.%s.%s.importance.txt", output_dir, mvar_level, level), header=F, as.is=T, sep="\t", row.names=1)
		collated.cv <- read.table(sprintf("%s/randomForest.%s.%s.cv.txt", output_dir, mvar_level, level), header=F, as.is=T, sep="\t", row.names=1)
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
#		save(sparseRF, file=sprintf("%s/randomForest.%s.%s.model", output_dir, mvar_level, level))
		load(sprintf("%s/randomForest.%s.%s.model", output_dir, mvar_level, level))
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
		write.table(confusion_matrix, file=sprintf("%s/randomForest.%s.%s.confusion_matrix.txt", output_dir, mvar_level, level), quote=F, sep="\t", row.names=T, col.names=T)
		# ROC analysis
		if (nlevels(response)==2) {
			pred <- predict(sparseRF, type="prob")
			pred2 <- prediction(pred[,2], ordered(response))
			perf <- performance(pred2, "tpr", "fpr")
			perf.auc <- performance(pred2, "auc")
			print(plot(perf, main=sprintf("ROC %s %s (sparseRF final model)", mvar_level, level)) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g", unlist(perf.auc@y.values))))
		}
		## END BLOCK TO COMMENT ##
		
		load(sprintf("%s/randomForest.%s.%s.model", output_dir, mvar_level, level))
		# plotting - per-group sparse model
		df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
		colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
		print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s", mvar_level, level)))
		
		# plotting - per-group variables
		df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
		res <- read.table(sprintf("%s/taxa_significance.%s.txt", output_dir, taxlevel_lookup[level]), header=T, sep="\t", as.is=T, quote="")
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
  taxlist <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level=level)
  otu.filt[[level]] <- taxlist
  otu.filt.abundance[[level]] <- taxlist
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

#	write.table(collated.importance, file=sprintf("%s/randomForest.%s.%s.importance.txt", output_dir, "multiclass", level), quote=F, sep="\t", row.names=T, col.names=F)
#	write.table(collated.cv, file=sprintf("%s/randomForest.%s.%s.cv.txt", output_dir, "multiclass", level), quote=F, sep="\t", row.names=T, col.names=F)
#	## END BLOCK TO COMMENT ##
	
	collated.importance <- read.table(sprintf("%s/randomForest.%s.%s.importance.txt", output_dir, "multiclass", level), header=F, as.is=T, sep="\t", row.names=1)
	collated.cv <- read.table(sprintf("%s/randomForest.%s.%s.cv.txt", output_dir, "multiclass", level), header=F, as.is=T, sep="\t", row.names=1)
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
#	save(sparseRF, file=sprintf("%s/randomForest.%s.%s.model", output_dir, "multiclass", level))
	load(sprintf("%s/randomForest.%s.%s.model", output_dir, "multiclass", level))
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
	write.table(confusion_matrix, file=sprintf("%s/randomForest.%s.%s.confusion_matrix.txt", output_dir, "multiclass", level), quote=F, sep="\t", row.names=T, col.names=T)
	# ROC analysis
	if (nlevels(response)==2) {
		pred <- predict(sparseRF, type="prob")
		pred2 <- prediction(pred[,2], ordered(response))
		perf <- performance(pred2, "tpr", "fpr")
		perf.auc <- performance(pred2, "auc")
		print(plot(perf, main=sprintf("ROC %s %s (sparseRF final model)", "multiclass", level)) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g", unlist(perf.auc@y.values))))
	}
	## END BLOCK TO COMMENT ##
	
	load(sprintf("%s/randomForest.%s.%s.model", output_dir, "multiclass", level))
	# plotting - per-group sparse model
	df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
	colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
	print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s", "multiclass", level)))
	
	# plotting - per-group variables
	df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
	res <- read.table(sprintf("%s/taxa_significance.%s.txt", output_dir, taxlevel_lookup[level]), header=T, sep="\t", as.is=T, quote="")
	res <- res[grep("Group", res$coefficient),]
	res <- aggregate(padj ~ OTU, res, min); rownames(res) <- res$OTU
	df$padj <- res[as.character(df$OTU), "padj"]
	df$OTU_string <- sprintf("%s (min padj=%.4g)", df$OTU, df$padj)
	df$sig <- ifelse(df$padj < 0.1, "sig", "ns"); df$sig[which(is.na(df$sig))] <- "ns"
	print(ggplot(df, aes(x=OTU, y=importance, label=OTU)) + geom_bar(position=position_dodge(), stat="identity", fill="#999999") + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string, color=sig), size=3, hjust=0) + ggtitle(sprintf("%s explanatory %s", "multiclass", level)) + scale_color_manual(values=cols.sig))
	# shading rectangles of importance values
	df.rect <- df
	df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
	df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
	p <- ggplot(df.rect, aes(x=x, y=OTU, color=sig, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s explanatory", "multiclass", level)) + scale_fill_gradient(low="white", high="black") + scale_color_manual(values=cols.sig)
	print(p)
	# violin plots of relative abundance
	agg.melt <- agg.melt.stored
	agg.melt$Group <- droplevels(mapping.sel[agg.melt$SampleID, "Group"])
	agg.melt <- subset(agg.melt, taxa %in% rownames(df))
	agg.melt$taxa <- factor(agg.melt$taxa, levels=rownames(df))
	p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap(~taxa, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF taxa (%s, %s)", "multiclass", level)) + coord_flip()
	print(p)
	agg.melt2 <- subset(agg.melt, taxa %in% taxlist)
	p <- ggplot(agg.melt2, aes(x=taxa, y=value, color=Group)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", "multiclass", level)) + coord_flip()
	print(p)
	
}


## DESeq2 analysis of PICRUSt-predicted metagenomes
library(DESeq2)
for (lvl in c("L0", "L1", "L2", "L3")) {
	metagenome <- read.table(sprintf("data/metagenome_contributions.%s.txt", lvl), header=T, as.is=T, sep="\t", comment.char="", skip=1, quote="", row.names=1)
	metagenome <- metagenome[, rownames(mapping.sel)]
	dds <- DESeqDataSetFromMatrix(metagenome, mapping.sel, design= ~ BMI + Age + Sex + Group)
	#inds <- which(rowSums(counts(dds)>=10) >= 5) # filter to at least 5 samples with 10+ reads
	#dds <- dds[inds,]
	dds <- DESeq(dds)
	for (mvar in c("fibrosis", "steatosis")) {
		res <- results(dds, name=sprintf("Group_%s_vs_normal", mvar))
		resSig <- as.data.frame(subset(res, padj<0.1)); resSig <- resSig[order(resSig$log2FoldChange),]; resSig$taxa <- rownames(resSig); resSig$taxa <- factor(resSig$taxa, levels=resSig$taxa); resSig$hdir <- 1-ceiling(sign(resSig$log2FoldChange)/2)
		if (nrow(resSig) > 0) {
			p <- ggplot(resSig, aes(x=taxa, y=log2FoldChange)) + geom_bar(stat="identity", fill="#aaaaaa") + geom_text(aes(label=taxa), y=0, size=2, hjust=0.5) + coord_flip() + theme_classic() + ggtitle(sprintf("DESeq2 hits (%s, %s)", lvl, mvar)) + theme(axis.text.y=element_blank())
			print(p)
		}
		write.table(res, file=sprintf("%s/PICRUSt_DESeq2.%s.%s.txt", output_dir, lvl, mvar), quote=F, sep="\t", row.names=T, col.names=T)
	}
}



##########################################################################################
## Dietary analysis
mapping.sel <- as(sample_data(psPROSPEC.relative), "data.frame")
subjects_to_exclude <- unique(rownames(which(is.na(mapping.sel), arr.ind=T)))
mapping.sel <- mapping.sel[setdiff(rownames(mapping.sel), subjects_to_exclude),]
psdiet.relative <- prune_samples(rownames(mapping.sel), psPROSPEC.relative)
psdiet.rarefied <- prune_samples(rownames(mapping.sel), psPROSPEC.rarefied)


## distance matrices
dm <- list()
dm[["diet"]] <- list()
for (distance_metric in distance_metrics) {
	dm[["diet"]][[distance_metric]] <- as.matrix(distance(psdiet.relative, method=distance_metric))
}

## PERMANOVA
sink(out_txt, append=T)
print("PERMANOVA on dietary samples")
for (distance_metric in distance_metrics) {
	print(distance_metric)
	form <- as.formula(sprintf("as.dist(dm[[\"diet\"]][[distance_metric]]) ~ %s", paste(diet_variables, collapse="+")))
	res <- adonis(form , data=as(sample_data(psPROSPEC.relative), "data.frame"), permutations=999)
	print(res)
}
sink()

## alpha diversity
adiv <- estimate_richness(psdiet.rarefied, measures=alpha_metrics)
adiv$SampleID <- rownames(adiv)
adiv <- merge(adiv, mapping.sel, by="row.names"); rownames(adiv) <- adiv$SampleID
res <- {}
for (mvar in diet_variables) {
	plotlist <- list()
	for (alpha_metric in alpha_metrics) {
		test <- cor.test(as.formula(sprintf("~ %s + %s", alpha_metric, mvar)), adiv, method="spearman")
		p <- ggplot(adiv, aes_string(x=mvar, y=alpha_metric)) + geom_point() + stat_smooth() + theme_classic() + ggtitle(sprintf("%s by %s (%s p=%.4g)", alpha_metric, mvar, test$method, test$p.value)) + theme(title=element_text(size=8))
		plotlist[[length(plotlist)+1]] <- p
		res <- rbind(res, c(mvar, alpha_metric, test$method, test$estimate, test$p.value))
	}
	multiplot(plotlist=plotlist, cols=2, rows=2)
}
res <- as.data.frame(res)
colnames(res) <- c("variable", "alpha_metric", "method", "estimate", "pval")
res$pval <- as.numeric(as.character(res$pval))
res$padj <- p.adjust(res$pval)
res <- res[order(res$pval),]
write.table(res, file=sprintf("%s/adiv.diet.txt", output_dir), quote=F, sep="\t", row.names=F, col.names=T)

## taxa significance - Genus level (ZINB)
otu.filt <- as.data.frame(otu_table(psdiet.rarefied))
otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psdiet.rarefied), level="Genus")
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
	for (mvar in diet_variables) {
		#print(sprintf("%s %s", f, mvar))
		df2 <- df
		df2[, mvar] <- mapping.sel[df2$SampleID, mvar]
		tt <- try(m <- zeroinfl(as.formula(sprintf("value ~ %s | 1", mvar)), data = df2, dist = "negbin", EM = FALSE, maxit=100), silent=T) # using EM=TRUE causes certain models to hang...
		if (class(tt) == "zeroinfl") {
			coef <- summary(m)$coefficients$count # rows are [(Intercept), comparisons, Log(theta)], columns are [Estimate, SE, Z, pval]
			tmp <- rbind(tmp, cbind(f, f, nsamps_detected, mvar, class(mapping.sel[,mvar]), "ZINB", rownames(coef), coef))
		} else if (class(tt) == "try-error") {
			tt <- try(m <- glm.nb(as.formula(sprintf("value ~ %s", mvar)), data = df2), silent=T)
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
res$Estimate <- as.numeric(as.character(res$Estimate)); res$SE <- as.numeric(as.character(res$SE))
# p-value adjust
res$padj <- p.adjust(res$pvalue, method="fdr")
res <- res[order(res$padj, decreasing=F),]
write.table(res, file=sprintf("%s/taxa_significance_diet.L6.txt", output_dir), quote=F, sep="\t", row.names=F, col.names=T)
# forest plots
for (mv in levels(res$metadata_variable)) {
	res.sel <- subset(res, metadata_variable==mv); res.sel <- res.sel[order(res.sel$Estimate),]
	res.sel$taxonomy <- factor(res.sel$taxonomy, levels=res.sel$taxonomy)
	res.sel$dir <- ifelse(res.sel$padj < siglevel, ifelse(sign(res.sel$Estimate)==1, "up", "down"), "NS")
	lims <- max(abs(res.sel$Estimate) + abs(res.sel$SE))*1.2
	inds_to_remove <- which(is.na(res.sel$Z))
	res.sel <- res.sel[setdiff(1:nrow(res.sel), inds_to_remove),]
	p <- ggplot(res.sel, aes(x=taxonomy, y=Estimate, color=dir)) + geom_point() + geom_errorbar(aes(x=taxonomy, ymin=Estimate-SE, max=Estimate+SE), width=0.2) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("%s (%s)", "ZINB hits", mv)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
	print(p)
}

####################################################################
####################################################################
####################################################################
### Cohesion analysis

# Online script to generate cohesion metrics for a set of samples 
# CMH 06Dec17; cherren@wisc.edu

# User instructions: read in a sample table (in absolute or relative abundance) as object "b".
# If using a custom correlation matrix, read in that matrix at the designated line.
# Run the entire script, and the 4 vectors (2 of connectedness and 2 of cohesion) are generated for each sample at the end.
# Parameters that can be adjusted include pers.cutoff (persistence cutoff for retaining taxa in analysis), 
# iter (number of iterations for the null model), tax.shuffle (whether to use taxon shuffle or row shuffle randomization), 
# and use.custom.cors (whether to use a pre-determined correlation matrix)

####################create necessary functions######################

#find the number of zeroes in a vector
zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)
}

#create function that averages only negative values in a vector
neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals) == 0) n.mean <- 0
  return(n.mean)
}

#create function that averages only positive values in a vector
pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals) == 0) p.mean <- 0
  return(p.mean)
}

###################################################################
###################################################################
### Workflow options ####
###################################################################
###################################################################

## Choose a persistence cutoff (min. fraction of taxon presence) for retaining taxa in the analysis
pers.cutoff <- 0.10
## Decide the number of iterations to run for each taxon. (>= 200 is recommended)
# Larger values of iter mean the script takes longer to run
iter <- 200
## Decide whether to use taxon/column shuffle (tax.shuffle = T) or row shuffle algorithm (tax.shuffle = F)
tax.shuffle <- T
## Option to input your own correlation table
# Note that your correlation table MUST have the same number of taxa as the abundance table. 
# There should be no empty (all zero) taxon vectors in the abundance table. 
# Even if you input your own correlation table, the persistence cutoff will be applied
use.custom.cors <- F
## Number of cores to use for 'mclapply' [set ncores=1 if deterministic results desired]
ncores <- 16

###################################################################
###################################################################

# Read in dataset
## Data should be in a matrix where each row is a sample. 
b <- t(as.data.frame(otu_table(psPROSPEC)))

# Read in custom correlation matrix, if desired. Must set "use.custom.cors" to TRUE
if(use.custom.cors == T) {
  custom.cor.mat <- read.csv("your_path_here.csv", header = T, row.names = 1)
  custom.cor.mat <- as.matrix(custom.cor.mat)
  #Check that correlation matrix and abundance matrix have the same dimension
  print(dim(b)[2] == dim(custom.cor.mat)[2])
}


# Suggested steps to re-format data. At the end of these steps, the data should be in a matrix "c" where there are no
# empty samples or blank taxon columns. 
c <- as.matrix(b)
c <- c[rowSums(c) > 0, colSums(c) > 0]

# Optionally re-order dataset to be in chronological order. Change date format for your data. 
#c <- c[order(as.Date(rownames(c), format = "%m/%d/%Y")), ]

# Save total number of individuals in each sample in the original matrix. This will be 1 if data are in relative abundance, 
# but not if matrix c is count data
rowsums.orig <- rowSums(c)

# Based on persistence cutoff, define a cutoff for the number of zeroes allowed in a taxon's distribution
zero.cutoff <- ceiling(pers.cutoff * dim(c)[1])
  
# Remove taxa that are below the persistence cutoff
d <- c[ , apply(c, 2, zero) < (dim(c)[1]-zero.cutoff) ]
# Remove any samples that no longer have any individuals, due to removing taxa
d <- d[rowSums(d) > 0, ]

#If using custom correlation matrix, need to remove rows/columns corresponding to the taxa below persistence cutoff
if(use.custom.cors == T){
  custom.cor.mat.sub <- custom.cor.mat[apply(c, 2, zero) < (dim(c)[1]-zero.cutoff), apply(c, 2, zero) < (dim(c)[1]-zero.cutoff)]
}

# Create relative abundance matrix.  
rel.d <- d / rowsums.orig
# Optionally, check to see what proportion of the community is retained after cutting out taxa
hist(rowSums(rel.d), breaks=40)

# Create observed correlation matrix
cor.mat.true <- cor(rel.d)

# Create vector to hold median otu-otu correlations for initial otu
med.tax.cors <- vector()

# set a seed
set.seed(sum(dim(rel.d)))

# Run this loop for the null model to get expected pairwise correlations
# Bypass null model if the option to input custom correlation matrix is TRUE
if(use.custom.cors == F) {
	if(tax.shuffle) {
		med.tax.cors <- do.call(cbind, mclapply(1:ncol(rel.d), function(which.taxon) {
		  #create vector to hold correlations from every permutation for each single otu
		  ## perm.cor.vec.mat stands for permuted correlations vector matrix
		  perm.cor.vec.mat <- vector()
		  for(i in 1:iter){
		    # create permuted matrix (with which.taxon column fixed)
		    perm.rel.d <- apply(rel.d, 2, sample)
		    perm.rel.d[, which.taxon] <- rel.d[, which.taxon]
		    rownames(perm.rel.d) <- rownames(rel.d)
		    # Calculate correlation matrix of permuted matrix
		    cor.mat.null <- cor(perm.rel.d)
		    # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
		    perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
		  }
		  # For large datasets, this can be helpful to know how long this loop will run
		  if(which.taxon %% 20 == 0){print(which.taxon)}
		  # Save the median correlations between the focal taxon and all other taxa  
		  apply(perm.cor.vec.mat, 1, median)
		}, mc.cores=ncores))
	} else {
		for(which.taxon in 1:dim(rel.d)[2]){
		med.tax.cors <- do.call(cbind, mclapply(1:ncol(rel.d), function(which.taxon) {
		  
		  #create vector to hold correlations from every permutation for each single otu
		  ## perm.cor.vec.mat stands for permuted correlations vector matrix
		  perm.cor.vec.mat <- vector()
		  
		  for(i in 1:iter){
		    #Create duplicate matrix to shuffle abundances
		    perm.rel.d <- rel.d 
		    
		    #For each taxon
		    for(j in 1:dim(rel.d)[1]){ 
		      which.replace <- which(rel.d[j, ] > 0 ) 
		      # if the focal taxon is greater than zero, take it out of the replacement vector, so the focal abundance stays the same
		      which.replace.nonfocal <- setdiff(which.replace, which.taxon)
		      
		      #Replace the original taxon vector with a vector where the values greater than 0 have been randomly permuted 
		      perm.rel.d[j, which.replace.nonfocal] <- sample(rel.d[ j, which.replace.nonfocal]) 
		    }

		    # Calculate correlation matrix of permuted matrix
		    cor.mat.null <- cor(perm.rel.d)
		    
		    # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
		    perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
		    
		  }
		  # For large datasets, this can be helpful to know how long this loop will run
		  if(which.taxon %% 20 == 0){print(which.taxon)}
		  # Save the median correlations between the focal taxon and all other taxa  
		  apply(perm.cor.vec.mat, 1, median)
		}, mc.cores=ncores))
	 }
	}
}

  
# Save observed minus expected correlations. Use custom correlations if use.custom.cors = TRUE
if(use.custom.cors == T) {
  obs.exp.cors.mat <- custom.cor.mat.sub
  } else {
    obs.exp.cors.mat <- cor.mat.true - med.tax.cors
  }
  
diag(obs.exp.cors.mat) <- 0

#### 
#### Produce desired vectors of connectedness and cohesion 

# Calculate connectedness by averaging positive and negative observed - expected correlations
connectedness.pos <- apply(obs.exp.cors.mat, 2, pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat, 2, neg.mean)

# Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg

cohesion <- cbind(cohesion.pos, cohesion.neg); colnames(cohesion) <- c("Cohesion.Positive", "Cohesion.Negative")
cohesion <- merge(cohesion, mapping.sel, by="row.names")

for (cohesion_var in c("Cohesion.Positive", "Cohesion.Negative")) {
	test <- kruskal.test(as.formula(sprintf("%s ~ Group", cohesion_var)), cohesion)
	p <- ggplot(cohesion, aes_string(x="Group", y=cohesion_var)) + geom_violin() + geom_jitter(position=position_jitter(0.2)) + stat_summary(fun.y=mean, geom="point", shape=23, size=4, fill="red") + ggtitle(sprintf("%s ~ Group (KW p=%.4g)", cohesion_var, test$p.value)) + theme_classic()
	print(p)
}

dev.off()




