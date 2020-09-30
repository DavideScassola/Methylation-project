setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
source("WGBS_analysis_functions.R", chdir = T)

# msr_ecdf_1e3 = readRDS("../../MethylationCode/MethylationData/msr_ecdf_1e3.Rda")
# msr_ecdf_1e4 = readRDS("../../MethylationCode/MethylationData/msr_ecdf_1e4.Rda")
# msr_ecdf_1e5 = readRDS("../../MethylationCode/MethylationData/msr_ecdf_1e5.Rda")
# msr_ecdf_1e6 = readRDS("../../MethylationCode/MethylationData/msr_ecdf_1e6.Rda")

# fragmets_msr_table
wgbs_file <- "../../MethylationCode/MethylationData/wgbs/"
new_name  <- "../../MethylationCode/MethylationData/wgbs/ENCFF752NXS_GM23248.rda" # could automate
short_name <- "GM23248"
msr_ecdf_file <- "../../MethylationCode/MethylationData/msr_ecdf_1e3.Rda"
size <- 1e3
produce_and_save_fragments_msr_table(wgbs_file, short_name, size, msr_ecdf_file, na_tolerance = 0.4, minimum_reads=1, methylation_assigner = standard_binaryzer, bed = NA, dir = "../../Rexperiments/")
  

# produce rna table
rna_file = "../../MethylationCode/MethylationData/rna-seq/"
extended <- T
correct_gene_id <- T

# read rna file
rna <- read_rna_file(rna_file, reduced = T, correct_gene_id = correct_gene_id); gc()
genebody_annotation <- readRDS("../../Rexperiments/genebody_improved.Rda"); gc()
extension_tag <- ""
if(extended)
{
  genebody_annotation <- readRDS("../../Rexperiments/human_genes_extended_improved.Rda"); gc()
  if(correct_gene_id)
    genebody_annotation$id <- remove_version_from_gene(genebody_annotation$id)
  extension_tag <- "_extended"
}



rna_table <- make_rna_window_data_frame(wgbs, rna, genebody_annotation, size)
saveRDS(rna_table, file = paste("../../Rexperiments/",short_name, "_rna_table_", size, extension_tag, ".Rda", sep = ""))


# final table test
data_table <- join_rna_and_msr_table(rna_table, msr_table)

rna_table = readRDS("../../Rexperiments/stomach_rna_table_1000_extended.Rda")


########
# genes table
short_name <- "lung_30_female"
wgbs_file  <- "../../MethylationCode/MethylationData/wgbs/ENCFF039JFT_lung_30_female.rda"
genebody_annotation_file <- "../../Rexperiments/detailed_genebody_improved.Rda"
produce_and_save_genes_msr_table(wgbs_file, short_name, genebody_annotation_file, gene_type_filter = NA, na_tolerance = 0.3, no_msr = F, dir = "../../Rexperiments/")
  