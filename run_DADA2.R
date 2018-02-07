#!/usr/bin/Rscript

library(ggplot2)
library(plyr)
library(reshape2)
library(dada2)
library(stringi)

multiplot <- function(..., plotlist=NULL, file, cols=1, rows=1) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
	
	i = 1
	while (i < numPlots) {
		numToPlot <- min(numPlots-i+1, cols*rows)
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(i, i+cols*rows-1), ncol = cols, nrow = rows, byrow=T)
		if (numToPlot==1) {
		  print(plots[[i]])
		} else {
		  # Set up the page
		  grid.newpage()
		  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		  # Make each plot, in the correct location
		  for (j in i:(i+numToPlot-1)) {
		    # Get the i,j matrix positions of the regions that contain this subplot
		    matchidx <- as.data.frame(which(layout == j, arr.ind = TRUE))
		    print(plots[[j]], vp = viewport(layout.pos.row = matchidx$row,
		                                    layout.pos.col = matchidx$col))
		  }
		}
		i <- i+numToPlot
  }
}

## Quality profiles
ncores <- 16
pdf(sprintf("/Lab_Share/Carolyn_Yanavich/phyloseq/qualityProfiles.pdf"))
fnFs <- list.files(sprintf("/Lab_Share/Carolyn_Yanavich/fastq/split.R1/"), pattern=".fastq", full.names=T)
fnRs <- list.files(sprintf("/Lab_Share/Carolyn_Yanavich/fastq/split.R2/"), pattern=".fastq", full.names=T)
samples <- sapply(fnFs, function(x) gsub(".fastq", "", basename(x)))
# plot quality profiles
p <- plotQualityProfile(fnFs, aggregate=T) + ggtitle(sprintf("Quality profile: (R1)"))
print(p)
p <- plotQualityProfile(fnRs, aggregate=T) + ggtitle(sprintf("Quality profile: (R2)"))
print(p)
dev.off()

## set truncation length parameters - these can vary as long as the same amplicon can be reconstructed
## DO NOT use trimLeft
truncation_lengths <- list(c(145,130))
names(truncation_lengths) <- c("080817")

run <- "080817"
ncores=16

## DADA2 inference
# make list of FASTQ files
fnFs <- list.files(sprintf("/Lab_Share/Carolyn_Yanavich/fastq/split.R1/"), pattern=".fastq", full.names=T)
fnRs <- list.files(sprintf("/Lab_Share/Carolyn_Yanavich/fastq/split.R2/"), pattern=".fastq", full.names=T)
samples <- sapply(fnFs, function(x) gsub(".fastq", "", basename(x)))
# filter FASTQs
filt_path <- file.path(sprintf("/Lab_Share/Carolyn_Yanavich/fastq/DADA2_filtered"))
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=truncation_lengths[[run]], maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=8, matchIDs=TRUE, verbose=TRUE)
# dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- samples
names(derepRs) <- samples
# learn error rates from subset of data
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=ncores)
errF <- dadaFs.lrn[[1]]$err_out
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=ncores)
errR <- dadaRs.lrn[[1]]$err_out
# sample inference with learned error rates
dadaFs <- dada(derepFs, err=errF, selfConsist = FALSE, multithread=ncores)
dadaRs <- dada(derepRs, err=errR, selfConsist = FALSE, multithread=ncores)
# plot error rates
pdf(sprintf("/Lab_Share/Carolyn_Yanavich/phyloseq/errorRates.pdf"))
for (i in seq_along(derepFs)) {
	print(plotErrors(dadaFs[[i]], nominalQ=TRUE))
	print(plotErrors(dadaRs[[i]], nominalQ=TRUE))
}
dev.off()
# merge read pairs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# make OTU table
seqtab <- makeSequenceTable(mergers)
save.image(file=sprintf("/Lab_Share/Carolyn_Yanavich/fastq/DADA2_workspace.RData"))
# remove singletons (sequences that only appear in one sample)
#	seqtab <- seqtab[, colSums(seqtab>0)>1]
# in-silico size restriction
table(nchar(getSequences(seqtab)))
seqtab <- seqtab[, nchar(colnames(seqtab)) %in% seq(from=252,to=254)] ### EDIT MANUALLY
#	seqtab <- seqtab[, nchar(colnames(seqtab)) %in% seq(from=232,to=234)] ### EDIT MANUALLY
# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE, multithread=ncores)
saveRDS(seqtab.nochim, sprintf("/Lab_Share/Carolyn_Yanavich/fastq/merged_seqtab.nochim.rds"))

# assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/Lab_Share/DADA2/rdp_train_set_14.fa.gz", multithread=ncores)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
taxa <- addSpecies(taxa, "/Lab_Share/DADA2/rdp_species_assignment_14.fa.gz", verbose=TRUE)
unname(head(taxa))
save(seqtab.nochim, taxa, file="/Lab_Share/Carolyn_Yanavich/fastq/Carolyn_Yanavich_DADA2.RData")

# save for BLAST taxonomy assignment
out <- data.frame(id=sprintf(">seq%06d", 1:ncol(seqtab.nochim)), seq=colnames(seqtab.nochim))
write.table(out, file="/Lab_Share/Carolyn_Yanavich/fastq/DADA2_sequences.fasta", quote=F, sep="\n", row.names=F, col.names=F)





