library(ShortRead)
library(biom) # note: use Joey's biom latest dev version; library(devtools); install_github("joey711/biom");

project_dir <- "/Lab_Share/Carolyn_Yanavich/fastq"

# 1) Make study db
# grab study seqs
seqtab.nochim <- readRDS(sprintf("%s/merged_seqtab.nochim.rds", project_dir))
seqs_study <- colnames(seqtab.nochim)
ids_study <- paste("study", 1:ncol(seqtab.nochim), sep = "_")
# merge db and study seqs
db_out <- data.frame(ids=ids_study,seqs=seqs_study,count=colSums(seqtab.nochim))
fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))
# write study fasta for filtering
writeFasta(fasta, file = sprintf("%s/gg_13_5_study_db.fasta.pre", project_dir))
# filter sequences that diverge from gg_13_5 by 97%
# depending on how well greengenes covers your study sequences, consider reducing 97% to 70 or 50%
system(sprintf("vsearch --usearch_global %s/gg_13_5_study_db.fasta.pre --db /home/fanli/miniconda3/pkgs/qiime-default-reference-0.1.3-py27_0/lib/python2.7/site-packages/qiime_default_reference/gg_13_5_otus/rep_set/97_otus.fasta --matched %s/gg_13_5_study_db.fasta --id 0.97", project_dir, project_dir))
id_filtered <- as.character(id(readFasta(sprintf("%s/gg_13_5_study_db.fasta", project_dir))))
db_out_filt <- db_out[db_out$ids%in%id_filtered,]
seqtab_biom <- t(seqtab.nochim)

# 2) output seq variant count data as biom;
# subset seqtab and output sample count biom
seqtab_biom <- seqtab_biom[rownames(seqtab_biom)%in%db_out_filt$seqs,]
rownames(seqtab_biom) <- db_out_filt[db_out_filt$seqs%in%rownames(seqtab_biom),"ids"]
biom_object <- biom::make_biom(data = seqtab_biom)
biom::write_biom(biom_object, biom_file = sprintf("%s/sample_counts.biom", project_dir))
# create final study db
system(sprintf("cat /home/fanli/miniconda3/pkgs/qiime-default-reference-0.1.3-py27_0/lib/python2.7/site-packages/qiime_default_reference/gg_13_5_otus/rep_set//97_otus.fasta >> %s/gg_13_5_study_db.fasta", project_dir))

