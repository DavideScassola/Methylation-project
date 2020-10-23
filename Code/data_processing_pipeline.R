setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
source("WGBS_analysis_functions.R", chdir = T)

# fragmets_msr_table
wgbs_file <- "../../MethylationCode/MethylationData/wgbs/ENCFF279HCL_GM12878.Rda"
bed  <- "../../MethylationCode/MethylationData/wgbs/ENCFF279HCL_GM12878.bed.gz"
short_name <- "GM12878"
size <- 1e3
msr_ecdf_file <- "../../MethylationCode/MethylationData/msr_ecdf_1e3.Rda"
produce_and_save_fragments_msr_table(wgbs_file, short_name, size, msr_ecdf_file, na_tolerance = 0.4, minimum_reads=1, methylation_assigner = standard_binaryzer, bed = NA, dir = "../../Rexperiments/")
########

# genes table
short_name <- "HelaS3"
wgbs_file  <- "../../MethylationCode/MethylationData/wgbs/ENCFF157VZT_HeLa_S3.Rda"
genebody_annotation_file <- "../../Rexperiments/detailed_genebody_improved.Rda"
produce_and_save_genes_msr_table(wgbs_file, short_name, genebody_annotation_file, gene_type_filter = NA, na_tolerance = 0.4, no_msr = F, dir = "../../Rexperiments/")
########

# fragments_expression_table
short_name <- "Hela"
wgbs_file <- sum_strands(readRDS(get_file_names(dir = "../../MethylationCode/MethylationData/wgbs/", patterns = short_name, T)))
expression_file <- get_file_names(dir = "../../MethylationCode/MethylationData/rna-seq/", patterns = c(short_name, "poly"), T)
genebody_annotation_file <- "../../Rexperiments/detailed_genebody_improved.Rda"
size <- 1e3
ignore_gene_version = T
produce_and_save_fragments_expression_table(wgbs_file, expression_file, genebody_annotation_file, size, ignore_gene_version, short_name, dir = "../../Rexperiments/", tag = "_protein_coding_", filter_gene_type = "protein_coding")
########


# all
short_names <- c("H1", "K562", "GM12878", "GM23248", "Hela", "endodermal", "lung", "stomach")
short_names <- c("H1", "K562", "GM12878", "GM23248", "Hela", "endodermal", "lung")
size <- 1e3
msr_ecdf_file <- "../../MethylationCode/MethylationData/msr_ecdf_1e3.Rda"
ignore_gene_version = T
dir = "../../Rexperiments/final/"

for(short_name in short_names)
{
wgbs_file <- get_file_names(dir = "../../MethylationCode/MethylationData/wgbs", patterns = short_name, T); cat("\n", wgbs_file)
wgbs_file <- sum_strands(readRDS(wgbs_file))
#expression_file <- get_file_names(dir = "../../MethylationCode/MethylationData/rna-seq", patterns = c(short_name, "poly"), T); cat("\n", expression_file)
#genebody_annotation_file <- "../../Rexperiments/detailed_genebody_improved.Rda"
#produce_and_save_fragments_expression_table(wgbs_file, expression_file, genebody_annotation_file, size, ignore_gene_version, short_name, dir = dir); gc(verbose=F)
#produce_and_save_genes_msr_table(wgbs_file, short_name, genebody_annotation_file, gene_type_filter = NA, na_tolerance = 0.4, no_msr = F, dir = dir); gc(verbose=F)
produce_and_save_fragments_msr_table(wgbs_file, short_name, size, msr_ecdf_file, na_tolerance = 0.4, minimum_reads=1, dir = dir); gc(verbose=F)
}
########
