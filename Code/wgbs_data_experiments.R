
#################################################################################################
setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
source("WGBS_analysis_functions.R", chdir = T)

#source_directory = "./Scrivania/Tesi/MethylationCode/"
#setwd(source_directory)
file_h1 = "../../MethylationCode/MethylationData/wgbs/ENCFF601NBW_H1_cell_line.bed.gz"
data_h1 <- read_ENCODE_bed(file_h1, verbose = T)

file_stomach = "../../MethylationCode/MethylationData/wgbs/ENCFF844EFX_stomach_man_51.bed.gz"
data_stomach <- read_ENCODE_bed(file_stomach, verbose = T)

file_lung = "../../MethylationData/wgbs/ENCFF039JFT_lung.bed.gz"
data_lung <- read_ENCODE_bed(file_lung, verbose = T)


#######################################################################################################


load(file = "../../Rexperiments/total_exp.Rdata")
load(file = "../../Rexperiments/total_exp_fake.Rdata")
plot(total_exp$`H1_inverted:_FALSE`[[1]]$data$fragments_infos_array[,2], pch = "+")


total_exp$`H1_inverted:_FALSE`[[1]]$data$fragments_infos_array[,2]











#######################################################################################################################
file_h1 = "../../MethylationCode/MethylationData/wgbs/ENCFF601NBW_H1_cell_line.bed.gz"
file_stomach = "../../MethylationCode/MethylationData/wgbs/ENCFF844EFX_stomach_man_51.bed.gz"

files = c(file_h1, file_stomach)
names = c("H1","stomach")
sizes = c(1e3, 1e4, 1e5, 1e6)
inversion = c(F,T)

system.time(
  total_exp <- total_spatial_experiment(files, sizes, inversion, names, methylation_assigner = stochastic_binaryzer, na_tolerance = 0.1, fake_data = F)
)

system.time(
  total_exp_fake <- total_spatial_experiment(files, sizes, inversion, names, methylation_assigner = stochastic_binaryzer, na_tolerance = 0.1, fake_data = T)
)

save(total_exp_fake, file = "total_exp_fake.Rdata")
save(total_exp, file = "total_exp.Rdata")
#######################################################################################################################


#################################################################
files = c(file_h1, file_stomach)
files = c(file_h1)
names = c("H1","stomach")
names = c("H1")
sizes = c(1e4, 1e5, 1e6, 1e7)
sizes = c(1e4, 1e5)
chromosome = "chrY"


data_stomach = read_ENCODE_bed(file_stomach)
pos <- get_methylation_positions(data_stomach, "chr1", sum_strands, stochastic_binaryzer, replace_no_reads_entries)
subpos = subset_positions(pos, 0, 1e5)

a = spatial_MSR_experiment_by_chromosome(subpos, 1e4, F, 10)
  
system.time(
  total_exp_chr1 <- total_spatial_experiment_by_chromosome(files, sizes, chromosome, names, stochastic_binaryzer, fake_data=F, minimum_bin_size = 200, mc = T)
)

system.time(
  total_exp_chr1_fake <- total_spatial_experiment_by_chromosome(files, sizes, chromosome, names, stochastic_binaryzer, fake_data=T, minimum_bin_size = 100)
)

save(total_exp_chr1, file = "total_exp_chr1.Rdata")
save(total_exp_fake_chr1, file = "total_exp_fake_chr1.Rdata")
#################################################################



binary <- get_methylation_CpG_binary_vector(data_h1,strands_handler = sum_strands, methylation_assigner = stochastic_binaryzer, missing_read_handler = keep_nas)
remove(data_h1)
gc()



spatial_exp_CpG_List_H1_1e4 = spatial_MSR_experiment_CpG_list(binary, invert=F, window_size=1e4, fake_data=F, na_tolerance=0.1, meth_assigner = stochastic_binaryzer)
spatial_exp_CpG_List_H1_1e5 = spatial_MSR_experiment_CpG_list(binary, invert=F, window_size=1e5, fake_data=F, na_tolerance=0.1, meth_assigner = stochastic_binaryzer)
spatial_exp_CpG_List_H1_1e5_fake = spatial_MSR_experiment_CpG_list(binary, invert=F, window_size=1e5, fake_data=T, na_tolerance=0.1, meth_assigner = stochastic_binaryzer)
spatial_exp_CpG_List_H1_1e4_fake = spatial_MSR_experiment_CpG_list(binary, invert=F, window_size=1e4, fake_data=T, na_tolerance=0.1, meth_assigner = stochastic_binaryzer)

spatial_exp_CpG_List_H1_1e4_inverted = spatial_MSR_experiment_CpG_list(binary, invert=F, window_size=1e4, fake_data=F, na_tolerance=0.1, meth_assigner = stochastic_binaryzer)
spatial_exp_CpG_List_H1_1e5 = spatial_MSR_experiment_CpG_list(binary, invert=F, window_size=1e5, fake_data=F, na_tolerance=0.1, meth_assigner = stochastic_binaryzer)
spatial_exp_CpG_List_H1_1e5_fake = spatial_MSR_experiment_CpG_list(binary, invert=F, window_size=1e5, fake_data=T, na_tolerance=0.1, meth_assigner = stochastic_binaryzer)
spatial_exp_CpG_List_H1_1e4_fake = spatial_MSR_experiment_CpG_list(binary, invert=F, window_size=1e4, fake_data=T, na_tolerance=0.1, meth_assigner = stochastic_binaryzer)

spatial_exp_CpG_List_H1_1e5$fragments_infos_array
lines(spatial_exp_CpG_List_H1_1e5_fake$fragments_infos_array[,3], lty = 3, col=2)
compare_resolution_relevance_plot(spatial_exp_CpG_List_H1_1e5$rr_list[6])

#######################################################################################################

























# On Genome
pos_h1_chr1 <- get_methylation_positions(data_h1, chromosome = "chr1", strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = replace_no_reads_entries)
rr_pos_h1_chr1 <- genome_MSR(pos_h1_chr1, minimum_bin_size = 20)
rr_plots(List(rr_pos_h1_chr1))

# On Genome but setting all CpGs to 1
CpG_positions_chr1 <- binary_nucleotides_pattern_positions("chr1", "CG", Genome = BSgenome.Hsapiens.UCSC.hg38)
rr_CpG_positions_on_Genome_chr1 <- calculate_relevance_resolution_vector_ignoring_nas(CpG_positions_chr1, minimum_bin_size = 20)
rr_plots(List(rr_CpG_positions_on_Genome, rr), legend_names = c("rr_CpG_positions_on_Genome", "rr"))

# On CpG List
binary_CpG_h1 <- get_methylation_CpG_binary_vector(data_h1, chromosome = "chr1", strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = keep_nas)
rr_binary_CpG_h1 <- calculate_relevance_resolution_vector(binary_CpG_h1, na_tolerance = 0.1, na_values_handler = replace_nas_hybrid_stochastic)
binary_CpG_stomach <- get_methylation_CpG_binary_vector(data_stomach, chromosome = "chr1", strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = keep_nas)
rr_binary_CpG_stomach <- calculate_relevance_resolution_vector(binary_CpG_stomach, na_tolerance = 0.1, na_values_handler = replace_nas_hybrid_stochastic)
rr_plots(list(rr_binary_CpG_h1, rr_binary_CpG_stomach), legend_names = c("h1","stomach"), title = "MSR curve on chr1, only CpG list")

##############################################################################
chr2_exp = methylation_experiment_by_chromosome(List(data_h1, data_stomach),names = c("h1", "stomach"),chromosome = "chr2")
chr2_exp$plotter()

whole_CpGlist_exp = methylation_experiment_CpGlist(List(data_h1, data_stomach, data_lung),names = c("h1", "stomach", "lung"), na_tolerance = 0.1)
whole_CpGlist_exp$plotter()

rr_50_threshold <- methylation_experiment_CpGlist(List(data_h1, data_stomach),names = c("h1", "stomach"), na_tolerance = 0.1, methylation_assigner = standard_binaryzer)
rr_adaptive_threshold <- methylation_experiment_CpGlist(List(data_h1, data_stomach),names = c("h1", "stomach"), na_tolerance = 0.1, methylation_assigner = rate_preserving_threshold_binaryzer)
rr_stochastic <- methylation_experiment_CpGlist(List(data_h1, data_stomach),names = c("h1", "stomach"), na_tolerance = 0.1, methylation_assigner = stochastic_binaryzer)

rr_50_threshold_inv <- methylation_experiment_CpGlist(List(data_h1, data_stomach),names = c("h1", "stomach"), na_tolerance = 0.1, methylation_assigner = standard_binaryzer, invert = T)
rr_adaptive_threshold_inv <- methylation_experiment_CpGlist(List(data_h1, data_stomach),names = c("h1", "stomach"), na_tolerance = 0.1, methylation_assigner = rate_preserving_threshold_binaryzer, invert = T)
rr_stochastic_inv <- methylation_experiment_CpGlist(List(data_h1, data_stomach),names = c("h1", "stomach"), na_tolerance = 0.1, methylation_assigner = stochastic_binaryzer, invert = T)
#################################################################################








size_list <- c(29e6, 5e6, 1e6, 3e5, 1e5, 1e4, 1e3)
rr_scales_stochastic <- different_scales_experiment(List(data_h1,data_stomach), size_list, names = c("h1", "stomach"), na_tolerance = 0.1, methylation_assigner = stochastic_binaryzer, invert =F)
#rr_scales_stochastic_no_tolerance <- different_scales_experiment(List(data_h1,data_stomach), size_list, names = c("h1", "stomach"), na_tolerance = 0, methylation_assigner = stochastic_binaryzer, invert = F)
rr_plots(rr_scales_stochastic[[1]], legend_names = size_list, title = "different scales experiment, stochastic assignment, h1 cells, CpG list")


size_list_genomewide <- c(29e6, 5e6, 1e6, 3e5, 1e5, 1e4, 1e3)*100



size_list_rid <- c(1e8, 1e7)
offset <- 1e7
rr_scales_stochastic_genomewide <- different_scales_experiment_genomewide(List(data_h1,data_stomach), size_list_rid, offset = offset, names = c("h1", "stomach"),chromosome = "chr1", methylation_assigner = stochastic_binaryzer)
rr_plots(rr_scales_stochastic_genomewide[[1]], legend_names = size_list_rid, title = "different scales experiment, stochastic assignment, chr1, h1 cells, on genome, starting from base 1e7")


### FRAGMENTS EXPERIMENTS ##########################################################################

data_list = List(data_h1, data_stomach)
remove(data_h1, data_stomach)
gc()
undersample = 20
names = c("H1, stomach")

rr_base = methylation_experiment_CpGlist(data_list,methylation_assigner = stochastic_binaryzer, invert=T,na_tolerance = 0.1)
base_area_h1      = MSR_area(rr_base$rr_list[[1]])
base_area_stomach = MSR_area(rr_base$rr_list[[2]])

sizes <- c(1e3, 1e4, 1e5, 1e6)

rr_fragments_1e3 <- different_positions_scale_CpG_list_experiment(data_list, invert = T, undersample = undersample, size=1e3, na_tolerance=0.1, strands_handler = sum_strands, methylation_assigner = stochastic_binaryzer, missing_read_handler = keep_nas)
gc()
save(rr_fragments_1e3, file = "rr_fragments_1e3.Rdata")
rr_fragments_1e4 <- different_positions_scale_CpG_list_experiment(data_list, invert = T, undersample = undersample, size=1e4, na_tolerance=0.1, strands_handler = sum_strands, methylation_assigner = stochastic_binaryzer, missing_read_handler = keep_nas)
gc()
save(rr_fragments_1e4, file = "rr_fragments_1e4.Rdata")
rr_fragments_1e5 <- different_positions_scale_CpG_list_experiment(data_list, invert = T, undersample = undersample, size=1e5, na_tolerance=0.1, strands_handler = sum_strands, methylation_assigner = stochastic_binaryzer, missing_read_handler = keep_nas)
gc()
save(rr_fragments_1e5, file = "rr_fragments_1e5.Rdata")

rr_fragments_1e6 <- different_positions_scale_CpG_list_experiment(data_list, invert = T, undersample = undersample, size=1e6, na_tolerance=0.1, strands_handler = sum_strands, methylation_assigner = stochastic_binaryzer, missing_read_handler = keep_nas)
gc()
save(rr_fragments_1e6, file = "rr_fragments_1e6.Rdata")

rr_fragments_1e7 <- different_positions_scale_CpG_list_experiment(data_list, invert = T, undersample = 0, size=1e7, na_tolerance=0.1, strands_handler = sum_strands, methylation_assigner = stochastic_binaryzer, missing_read_handler = keep_nas)
gc()
save(rr_fragments_1e7, file = "rr_fragments_1e7.Rdata")

  
gc()

areas_h1_1e3 = sapply(rr_fragments_1e3[[1]], function(x) MSR_area(x))
areas_h1_1e4 = sapply(rr_fragments_1e4[[1]], function(x) MSR_area(x))
areas_h1_1e5 = sapply(rr_fragments_1e5[[1]], function(x) MSR_area(x))
areas_h1_1e6 = sapply(rr_fragments_1e6[[1]], function(x) MSR_area(x))
areas_h1_1e7 = sapply(rr_fragments_1e7[[1]], function(x) MSR_area(x))

areas_stomach_1e3 = sapply(rr_fragments_1e3[[2]], function(x) MSR_area(x))
areas_stomach_1e4 = sapply(rr_fragments_1e4[[2]], function(x) MSR_area(x))
areas_stomach_1e5 = sapply(rr_fragments_1e5[[2]], function(x) MSR_area(x))
areas_stomach_1e6 = sapply(rr_fragments_1e6[[2]], function(x) MSR_area(x))
areas_stomach_1e7 = sapply(rr_fragments_1e7[[2]], function(x) MSR_area(x))

boxplot.default(areas_h1_1e3, areas_h1_1e4, areas_h1_1e5, areas_h1_1e6, areas_h1_1e7, base_area_h1, names = c(1e3,1e4,1e5,1e6,1e7, "whole"))
title("MSR area calculated for 20 random different fragments (h1 cells, CpG list, stochastic assignment, not methylated)")

boxplot.default(areas_stomach_1e3, areas_stomach_1e4, areas_stomach_1e5, areas_stomach_1e6, areas_stomach_1e7, base_area_stomach, names = c(1e3,1e4,1e5,1e6,1e7, "whole"))
title("MSR area calculated for 20 random different fragments (stomach cells, CpG list, stochastic assignment, not methylated)")

boxplot.default(areas_h1_1e3, areas_stomach_1e3, areas_h1_1e4, areas_stomach_1e4, areas_h1_1e5, areas_stomach_1e5, areas_h1_1e6, areas_stomach_1e6, areas_h1_1e7, areas_stomach_1e7, base_area_h1, base_area_stomach, names = c("H1:1e3","stomach:1e3","H1:1e4","stomach:1e4","H1:1e5","stomach:1e5","H1:1e6","stomach:1e6","H1:1e7","stomach:1e7", "H1:whole","stomach:whole"))
title("MSR area calculated for 20 random different fragments (CpG list, stochastic assignment, not methylated)")




areas_stomach = sapply(rr_fragments[[2]], function(x) MSR_area(x))

boxplot.default(areas_stomach, areas_h1, names = c("stomach", "H1"))
title("MSR area calculated for 20 different fragments of size 1e3 (CpG list, stochastic assignment)")



rr_fragments <- different_positions_scale_CpG_list_experiment(data_list, size, names, na_tolerance=0.1, strands_handler = sum_strands, methylation_assigner = stochastic_binaryzer, missing_read_handler = keep_nas, invert = F, undersample = undersample)
compare_resolution_relevance_plot(rr_fragments[[1]], title = "MSR plot for 29 different fragments of 1e6 sites of stomach cells (CpG list, stochastic assignment, not meth)")

areas_h1 = sapply(rr_fragments[[1]], function(x) MSR_area(x))
areas_stomach = sapply(rr_fragments[[2]], function(x) MSR_area(x))

#plot(areas_stomach, areas_h1)
mean(areas_h1, na.rm=T)
mean(areas_stomach, na.rm=T)

plot(areas_h1)
boxplot.default(areas_stomach, areas_h1, names = c("stomach", "H1"))
title("MSR area calculated for 20 different fragments of size 1e3 (CpG list, stochastic assignment)")



#################################################################################



### FRAGMENTS EXPERIMENTS ##########################################################################



file_h1 = "../../MethylationCode/MethylationData/wgbs/ENCFF601NBW_H1_cell_line.bed.gz"
file_stomach = "../../MethylationCode/MethylationData/wgbs/ENCFF844EFX_stomach_man_51.bed.gz"

data_list = List(read_ENCODE_bed(file_h1, verbose = T), read_ENCODE_bed(file_stomach, verbose = T))
size = 1e8
rr_fragments_stomach_1e8 <- different_positions_scale_by_chromosome_experiment(d, size, invert = F, undersample = 0, chromosome = "chr1", minimum_bin_size = 30, strands_handler = sum_strands, methylation_assigner = stochastic_binaryzer)

par(mfrow=c(3,2))
compare_resolution_relevance_plot(rr_fragments_stomach_1e4, title = "fragment size: 1e4")
compare_resolution_relevance_plot(rr_fragments_stomach_1e5, title = "fragment size: 1e5")
compare_resolution_relevance_plot(rr_fragments_stomach_1e6, title = "fragment size: 1e6")
compare_resolution_relevance_plot(rr_fragments_stomach_1e7, title = "fragment size: 1e7")
compare_resolution_relevance_plot(rr_fragments_stomach_1e8, title = "fragment size: 1e8")
compare_resolution_relevance_plot(List(rr_whole_stomach), title = "whole")
title("MSR curves for different fragments for stomach cells (on chr1, stochastic assignment)")

rr_whole_stomach = genome_MSR(get_methylation_positions(d, "chr1", strands_handler = sum_strands, methylation_assigner = stochastic_binaryzer, missing_read_handler = replace_no_reads_entries),minimum_bin_size = 30)
  
areas_stomach_1e4 = sapply(rr_fragments_stomach_1e4, function(x) MSR_area(x))
areas_stomach_1e5 = sapply(rr_fragments_stomach_1e5, function(x) MSR_area(x))
areas_stomach_1e6 = sapply(rr_fragments_stomach_1e6, function(x) MSR_area(x))
areas_stomach_1e7 = sapply(rr_fragments_stomach_1e7, function(x) MSR_area(x))
areas_stomach_1e8 = sapply(rr_fragments_stomach_1e8, function(x) MSR_area(x))
area_stomach_whole = MSR_area(genome_MSR(get_methylation_positions(d, "chr1", strands_handler = sum_strands, methylation_assigner = stochastic_binaryzer, missing_read_handler = replace_no_reads_entries),minimum_bin_size = 30))

boxplot.default(areas_stomach_1e4, areas_stomach_1e5, areas_stomach_1e6, areas_stomach_1e7, areas_stomach_1e8, area_stomach_whole, names = c(1e4,1e5,1e6,1e7,1e8,"whole"))
title("MSR area calculated for different fragments for stomach cells (on chr1, stochastic assignment)")

plot(areas_stomach, areas_h1)
plot(areas_stomach)
plot(areas_h1)
boxplot.default(areas_stomach, areas_h1, names = c("stomach", "H1"))
title("MSR area calculated for x different fragments of size 1e5 (CpG list, stochastic assignment)")
##########################################################################################################

load(file = "../exps/rr_fragments_1e5.Rdata")
load(file = "../exps/rr_fragments_1e6.Rdata")
load(file = "../exps/rr_fragments_1e7.Rdata")

areas_h1 = sapply(rr_fragments_1e5[[1]], function(x) MSR_area(x))
areas_stomach = sapply(rr_fragments[[2]], function(x) MSR_area(x))

#################################################################################

setwd("./Scrivania/Tesi/MethylationCode/")
directory <- "MethylationData/binary_rate/"
cell_files <- list.files(directory,pattern="(_converted.Rda)$")
cell_names <- sub("_converted.Rda","", cell_files)
methylation_vector <- as.logical(readRDS("MethylationData/binary_rate/GSM3436261_O1_TA_Hi_10_converted.Rda"))

#################################################################################
