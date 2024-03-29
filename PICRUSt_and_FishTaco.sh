### Carolyn Yanavich project ###

############################################################################
### 1. Generate per-sample FASTQ files for DADA2 and phyloseq analysis
### a. Split FASTQ by barcode
#	INPUT:
#		fastq/${run}/Undetermined_S0_I1_001.fastq.gz
#		fastq/${run}/Undetermined_S0_R1_001.fastq.gz
#		fastq/${run}/Undetermined_S0_R2_001.fastq.gz
#		mapping/${run}.txt
#	EXECUTION:

split_libraries_fastq.py -i /Lab_Share/Carolyn_Yanavich/fastq/Undetermined_S0_R1_001.fastq.gz -b /Lab_Share/Carolyn_Yanavich/fastq/Undetermined_S0_I1_001.fastq.gz -m /Lab_Share/Carolyn_Yanavich/CY_Mapping.080817.txt --store_demultiplexed_fastq -o /Lab_Share/Carolyn_Yanavich/fastq/sl_out.R1
split_libraries_fastq.py -i /Lab_Share/Carolyn_Yanavich/fastq/Undetermined_S0_R2_001.fastq.gz -b /Lab_Share/Carolyn_Yanavich/fastq/Undetermined_S0_I1_001.fastq.gz -m /Lab_Share/Carolyn_Yanavich/CY_Mapping.080817.txt --store_demultiplexed_fastq -o /Lab_Share/Carolyn_Yanavich/fastq/sl_out.R2
split_sequence_file_on_sample_ids.py -i /Lab_Share/Carolyn_Yanavich/fastq/sl_out.R1/seqs.fastq --file_type fastq -o /Lab_Share/Carolyn_Yanavich/fastq/split.R1 &
split_sequence_file_on_sample_ids.py -i /Lab_Share/Carolyn_Yanavich/fastq/sl_out.R2/seqs.fastq --file_type fastq -o /Lab_Share/Carolyn_Yanavich/fastq/split.R2 &

#	OUTPUT:
#		fastq/sl_out.R1, fastq/sl_out.R2
#		fastq/split.R1, fastq/split.R2


### b. Run DADA2 and phyloseq
./run_DADA2.R

blastn -task blastn -query /Lab_Share/Carolyn_Yanavich/fastq/DADA2_sequences.fasta -db /Lab_Share/SILVA/SILVA_128_SSURef_Nr99_tax_silva_trunc -out /Lab_Share/Carolyn_Yanavich/fastq/BLAST_results.txt -num_threads 16 -evalue 0.01 -outfmt 6 -num_alignments 500
ruby parse_blast_output.rb /Lab_Share/Carolyn_Yanavich/fastq/BLAST_results.txt /Lab_Share/SILVA/SILVA_128_SSURef_Nr99_tax_silva_trunc.headers > /Lab_Share/Carolyn_Yanavich/fastq/BLAST_results.parsed.txt


############################################################################
### 2. Analysis using phyloseq
./run_phyloseq.R

############################################################################
### 3. DADA2 -> PICRUSt (https://github.com/vmaffei/dada2_to_picrust)
###	Part 0: Generate gg_16S_counts.tab and gg_ko_counts.tab files (NOTE: only has to be done ONCE ever)
# build gg_16S_counts.tab and gg_ko_counts.tab
PICRUST_DATA_DIR=/share/picrust/picrust-1.1.1/picrust/data
gunzip -c ${PICRUST_DATA_DIR}/16S_13_5_precalculated.tab.gz | sed 's/\#OTU_IDs/taxon_oid/g' > ${PICRUST_DATA_DIR}/gg_16S_counts.tab
gunzip -c ${PICRUST_DATA_DIR}/ko_13_5_precalculated.tab.gz | sed 's/\#OTU_IDs/GenomeID/g' | grep -v metadata_KEGG | cut -f 1-6910 > ${PICRUST_DATA_DIR}/gg_ko_counts.tab
# remove empty lines in gg_ko_counts.tab
sed -i '/^\s*$/d' ${PICRUST_DATA_DIR}/gg_ko_counts.tab
# save KEGG Pathways and Description metadata
echo "" > ${PICRUST_DATA_DIR}/kegg_meta
gunzip -c ${PICRUST_DATA_DIR}/ko_13_5_precalculated.tab.gz | grep metadata_KEGG >> ${PICRUST_DATA_DIR}/kegg_meta

### Part 1. Add study sequences to 97_otus.fasta and tabulate counts (in R)
./DADA2_to_PICRUSt.R

### Part 2. Align seqs and build tree
# align w/ pynast using QIIME scripts; the included options lessen alignment restrictions to prevent alignment failure
# minimum sequence length set by -e
# alignment runtime greatly reduced by parallelization: parallel_align_seqs_pynast.py -O and # of cores
PROJECT_DIR=/Lab_Share/Carolyn_Yanavich/fastq
parallel_align_seqs_pynast.py -e 90 -p 0.1 -i ${PROJECT_DIR}/gg_13_5_study_db.fasta -o ${PROJECT_DIR}/gg_13_5_study_db.fasta.aligned -O 45
# filter alignment with default settings; consider lane filtering by entropy using -e and a low entropy value of ~0.01-0.02
# note: FastTree and/or PICRUSt produce weird errors (segfaults) if -e filters too many lanes
filter_alignment.py -i ${PROJECT_DIR}/gg_13_5_study_db.fasta.aligned/gg_13_5_study_db_aligned.fasta -o ${PROJECT_DIR}/gg_13_5_study_db.fasta.aligned.filtered/
# build tree with fasttree; options are taken from greengenes 13_5 readme notes
# tree building runtime greatly reduced by parallelization: use FastTreeMP w/ same options instead of FastTree
/home/fanli/FastTreeMP -nt -gamma -fastest -no2nd -spr 4 ${PROJECT_DIR}/gg_13_5_study_db.fasta.aligned.filtered/gg_13_5_study_db_aligned_pfiltered.fasta > ${PROJECT_DIR}/study_tree.tree

### Part 3. Create new precalculated files
# format 16S copy number data
mkdir -p ${PROJECT_DIR}/format
mkdir -p ${PROJECT_DIR}/asr
mkdir -p ${PROJECT_DIR}/predict_traits
format_tree_and_trait_table.py -t ${PROJECT_DIR}/study_tree.tree -i ${PICRUST_DATA_DIR}/gg_16S_counts.tab -o ${PROJECT_DIR}/format/16S/
# format kegg IMG data
format_tree_and_trait_table.py -t ${PROJECT_DIR}/study_tree.tree -i ${PICRUST_DATA_DIR}/gg_ko_counts.tab -o ${PROJECT_DIR}/format/KEGG/
# perform ancestral state reconstruction
ancestral_state_reconstruction.py -i ${PROJECT_DIR}/format/16S/trait_table.tab -t ${PROJECT_DIR}/format/16S/pruned_tree.newick -o ${PROJECT_DIR}/asr/16S_asr_counts.tab -c ${PROJECT_DIR}/asr/asr_ci_16S.tab
ancestral_state_reconstruction.py -i ${PROJECT_DIR}/format/KEGG/trait_table.tab -t ${PROJECT_DIR}/format/KEGG/pruned_tree.newick -o ${PROJECT_DIR}/asr/KEGG_asr_counts.tab -c ${PROJECT_DIR}/asr/asr_ci_KEGG.tab
# collect study sequence ids for predict_traits.py -l (greatly reduces runtime)
# convert biom to tsv using biom-format
biom convert -i ${PROJECT_DIR}/sample_counts.biom -o ${PROJECT_DIR}/sample_counts.tab --to-tsv
# predict traits
#predict_traits.py -i ${PROJECT_DIR}/format/16S/trait_table.tab -t ${PROJECT_DIR}/format/16S/reference_tree.newick -r ${PROJECT_DIR}/asr/16S_asr_counts.tab -o ${PROJECT_DIR}/predict_traits/16S_precalculated.tab -a -c ${PROJECT_DIR}/asr/asr_ci_16S.tab -l ${PROJECT_DIR}/sample_counts.tab
predict_traits.py -i ${PROJECT_DIR}/format/16S/trait_table.tab -t ${PROJECT_DIR}/format/16S/reference_tree.newick -r ${PROJECT_DIR}/asr/16S_asr_counts.tab -o ${PROJECT_DIR}/predict_traits/16S_precalculated.tab -l ${PROJECT_DIR}/sample_counts.tab
predict_traits.py -i ${PROJECT_DIR}/format/KEGG/trait_table.tab -t ${PROJECT_DIR}/format/KEGG/reference_tree.newick -r ${PROJECT_DIR}/asr/KEGG_asr_counts.tab -o ${PROJECT_DIR}/predict_traits/ko_precalculated.tab -l ${PROJECT_DIR}/sample_counts.tab
# add KEGG metadata
cp ${PROJECT_DIR}/predict_traits/ko_precalculated.tab ${PROJECT_DIR}/predict_traits/ko_precalculated_no_metadata.tab
cat kegg_meta >> ${PROJECT_DIR}/predict_traits/ko_precalculated.tab

## yay, finally actually run PICRUSt
normalize_by_copy_number.py -i  ${PROJECT_DIR}/sample_counts.biom -o ${PROJECT_DIR}/sample_counts.norm.biom -c ${PROJECT_DIR}/predict_traits/16S_precalculated.tab
predict_metagenomes.py -i ${PROJECT_DIR}/sample_counts.norm.biom -o ${PROJECT_DIR}/metagenome_contributions.L0.biom -c ${PROJECT_DIR}/predict_traits/ko_precalculated.tab
# optional: agglomerate counts by KEGG pathway level
for lvl in 1 2 3 
do
	categorize_by_function.py -i ${PROJECT_DIR}/metagenome_contributions.L0.biom -o ${PROJECT_DIR}/metagenome_contributions.L${lvl}.biom -c KEGG_Pathways -l ${lvl}
	biom convert -i ${PROJECT_DIR}/metagenome_contributions.L${lvl}.biom --to-tsv -o ${PROJECT_DIR}/metagenome_contributions.L${lvl}.txt
done
# convert to text format
biom convert -i ${PROJECT_DIR}/metagenome_contributions.L0.biom --to-tsv -o ${PROJECT_DIR}/metagenome_contributions.L0.txt


############################################################################
### 4. FishTaco with PICRUSt predictions + 16S abundances
PROJECT_DIR=/Lab_Share/Carolyn_Yanavich/fastq

# convert to relative abundance files
./rel_abund.R ${PROJECT_DIR}/sample_counts.tab ${PROJECT_DIR}/sample_counts.relabund.tab
./rel_abund.R ${PROJECT_DIR}/metagenome_contributions.L0.txt ${PROJECT_DIR}/metagenome_contributions.L0.relabund.txt

# generate sample labels files, filter input files for each Group (fibrosis/steatosis)
mkdir -p /Lab_Share/Carolyn_Yanavich/PICRUSt/
./make_FishTaco_input.R

# normalize using MUSiCC
for gr in fibrosis steatosis
do
	run_musicc.py -n -c use_generic -o /Lab_Share/Carolyn_Yanavich/PICRUSt/metagenome_contributions.L0.relabund.${gr}.MUSiCC_corrected.txt /Lab_Share/Carolyn_Yanavich/PICRUSt/metagenome_contributions.L0.relabund.${gr}.txt
done

# map KEGG KOs to pathways prior to differential abundance in FishTaco (NOTE: the map_function_level option does this AFTER DA testing so not helpful)
for gr in fibrosis steatosis
do
	python ~/miniconda3/lib/python3.5/site-packages/fishtaco/compute_pathway_abundance.py -ko /Lab_Share/Carolyn_Yanavich/PICRUSt/metagenome_contributions.L0.relabund.${gr}.MUSiCC_corrected.txt -ko2path ~/miniconda3/lib/python3.5/site-packages/fishtaco/data/KOvsPATHWAY_BACTERIAL_KEGG_2013_07_15.tab --output /Lab_Share/Carolyn_Yanavich/PICRUSt/metagenome_contributions.L0.relabund.${gr}.MUSiCC_corrected.Pathway.txt --output_counts /Lab_Share/Carolyn_Yanavich/PICRUSt/metagenome_contributions.L0.relabund.${gr}.MUSiCC_corrected.Pathway_counts.txt
done

# in R: convert to relative abundance and filter pathways to only those present with mean relative abundance > 0.1%; also manually uses padj<0.2 cutoff to identify DA functions
./make_FishTaco_input2.R

# run DA (for testing and viewing full results)
for gr in fibrosis steatosis
do
	python ~/miniconda3/lib/python3.5/site-packages/fishtaco/compute_differential_abundance.py --class /Lab_Share/Carolyn_Yanavich/PICRUSt/sample_labels.${gr}.txt -o /Lab_Share/Carolyn_Yanavich/PICRUSt/DA_functions.${gr}.tab /Lab_Share/Carolyn_Yanavich/PICRUSt/metagenome_contributions.L0.relabund.${gr}.MUSiCC_corrected.Pathway.filtered.txt
done

# identify DA functions
#./make_FishTaco_input2.R
## NOTE: no DA functions identified in either 'fibrosis' or 'steatosis' group

# run FishTaco
for gr in fibrosis steatosis
do
	run_fishtaco.py -ta /Lab_Share/Carolyn_Yanavich/PICRUSt/sample_counts.relabund.${gr}.tab -fu /Lab_Share/Carolyn_Yanavich/PICRUSt/metagenome_contributions.L0.relabund.${gr}.txt -l /Lab_Share/Carolyn_Yanavich/PICRUSt/sample_labels.${gr}.txt -gc ${PROJECT_DIR}/predict_traits/ko_precalculated_no_metadata.tab -op fishtaco_out.${gr} -log
done


