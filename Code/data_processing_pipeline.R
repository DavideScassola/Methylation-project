setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
source("WGBS_analysis_functions.R", chdir = T)
msr_ecdf_1e3 = readRDS("../../MethylationCode/MethylationData/msr_ecdf_1e3.Rda")
msr_ecdf_1e4 = readRDS("../../MethylationCode/MethylationData/msr_ecdf_1e4.Rda")
msr_ecdf_1e5 = readRDS("../../MethylationCode/MethylationData/msr_ecdf_1e5.Rda")
msr_ecdf_1e6 = readRDS("../../MethylationCode/MethylationData/msr_ecdf_1e6.Rda")


wgbs_file <- "../../MethylationCode/MethylationData/wgbs/stomach.Rda"
new_name  <- "../../MethylationCode/MethylationData/wgbs/stomach.Rda" # could automate
short_name <- "stomach"


# convert bed file into RDA
rda_convert(wgbs_file, new_name)

# read rda file
wgbs <- sum_strands(readRDS(new_name))

# produce msr table
size <- 1e3
msr_ecdf_ref <- msr_ecdf_1e3
methylation_assigner <- standard_binaryzer
na_tolerance <- 0.1
minimum_reads <- 1
rr_table <- total_spatial_experiment(c(new_name), c(size), c(F,T), c(short_name), methylation_assigner, na_tolerance, F, minimum_reads)
start = rr_table[[1]][[1]]$data$fragments_infos_array[,1]
msr_density = rr_table[[1]][[1]]$data$fragments_infos_array[,2]
true_density = sapply(start, function(x){mean(wgbs$prop[x:(x+size)], na.rm = T)})/100
msr = rr_table[[1]][[1]]$data$fragments_infos_array[,3]
inverted_msr = rr_table[[2]][[1]]$data$fragments_infos_array[,3]
sig <- significance_measure(msr, msr_density, msr_ecdf_ref, inverted = F)
inverted_sig <- significance_measure(inverted_msr, msr_density, msr_ecdf_ref, inverted = T)
median_function <- extract_ecdf_function(msr_ecdf_ref, 0.5)
residual <- msr-median_function(msr_density)
inverted_residual <- inverted_msr-median_function(1-msr_density)

# save msr tables
msr_table = data.frame(start,msr_density,true_density,msr,inverted_msr, sig, inverted_sig, residual, inverted_residual)
saveRDS(msr_table, file = paste("../../Rexperiments/",short_name, "_msr_table_", size, ".Rda", sep = ""))



# produce rna table
rna_file = "../../MethylationCode/MethylationData/rna-seq/ENCFF860DPP_left_ventricle_male_34.tsv"

# read rna file
rna <- read_rna_file(rna_file, reduced = F)
genebody_annotation <- readRDS("../../Rexperiments/genebody_improved.Rda")
rna_table <- make_rna_window_data_frame(wgbs, rna, genebody_annotation, size)
saveRDS(rna_table, file = paste("../../Rexperiments/",short_name, "_rna_table_", size, ".Rda", sep = ""))
  

# final table test
data_table <- join_rna_and_msr_table(rna_table, msr_table)



