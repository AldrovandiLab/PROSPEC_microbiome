#!/usr/bin/Rscript


library(useful)
library(ggplot2)
library(FishTacoPlot)
library(igraph)

input_dir <- "/Lab_Share/Carolyn_Yanavich/071018/humann2/fishtaco"

## large plots first
out_pdf <- sprintf("%s/fishtaco_plot_LARGE.pdf", input_dir)
pdf(out_pdf, width=36, height=48)
df <- {}

# steatosis UP
prefix <- "fishtaco_out.steatosis_UP"
tax_file <- sprintf("%s/fishtaco_out.steatosis.taxonomy.txt", input_dir)
print(ggplot(df) + geom_blank() + theme_classic() + annotate("text", x=1, y=1, label=sprintf("steatosis UP"), size=16))
p <- MultiFunctionTaxaContributionPlots(input_dir=input_dir, input_prefix=prefix, input_taxa_taxonomy=tax_file, sort_by="list", plot_type="bars", add_predicted_da_markers=TRUE, show_only_diff_abun_taxa=TRUE, color_small_cont_by_phyla=FALSE, add_names_in_bars=TRUE)
print(p)

# steatosis DOWN
prefix <- "fishtaco_out.steatosis_DOWN"
tax_file <- sprintf("%s/fishtaco_out.steatosis.taxonomy.txt", input_dir)
print(ggplot(df) + geom_blank() + theme_classic() + annotate("text", x=1, y=1, label=sprintf("steatosis DOWN"), size=16))
p <- MultiFunctionTaxaContributionPlots(input_dir=input_dir, input_prefix=prefix, input_taxa_taxonomy=tax_file, sort_by="list", plot_type="bars", add_predicted_da_markers=TRUE, show_only_diff_abun_taxa=TRUE, color_small_cont_by_phyla=FALSE, add_names_in_bars=TRUE)
print(p)

dev.off()

## normal size plots
out_pdf <- sprintf("%s/fishtaco_plot.pdf", input_dir)
pdf(out_pdf, width=12)


# steatosis UP
prefix <- "fishtaco_out.steatosis_UP"
tax_file <- sprintf("%s/fishtaco_out.steatosis.taxonomy.txt", input_dir)
print(ggplot(df) + geom_blank() + theme_classic() + annotate("text", x=1, y=1, label=sprintf("steatosis UP"), size=16))
p <- MultiFunctionTaxaContributionPlots(input_dir=input_dir, input_prefix=prefix, input_taxa_taxonomy=tax_file, sort_by="list", plot_type="bars", add_predicted_da_markers=TRUE, show_only_diff_abun_taxa=TRUE, color_small_cont_by_phyla=FALSE)
print(p)

# steatosis DOWN
prefix <- "fishtaco_out.steatosis_DOWN"
tax_file <- sprintf("%s/fishtaco_out.steatosis.taxonomy.txt", input_dir)
print(ggplot(df) + geom_blank() + theme_classic() + annotate("text", x=1, y=1, label=sprintf("steatosis DOWN"), size=16))
p <- MultiFunctionTaxaContributionPlots(input_dir=input_dir, input_prefix=prefix, input_taxa_taxonomy=tax_file, sort_by="list", plot_type="bars", add_predicted_da_markers=TRUE, show_only_diff_abun_taxa=TRUE, color_small_cont_by_phyla=FALSE)
print(p)


dev.off()


