#################################################################################

# read the data from file

# count missing data
# count Hemymethylation

# preparing the data for MSR analysis:
    # whole genome or just chromosome (need to know human chromosome size)
    # how to treat strands
        # pick just + or -
        # sum + and - counts
        # do an or
    # how to treat methylation proportions
    # how to treat missing values
    # contiguous CpG or whole genome


# calculate MSR (probably with preliminar binning)
# plot MSR

#########################################################################################

#################################################################################################
source("WGBS_analysis_functions.R", chdir = T)

#source_directory = "./Scrivania/Tesi/MethylationCode/"
#setwd(source_directory)
file_h1 = "../../MethylationCode/MethylationData/wgbs/ENCFF601NBW_H1_cell_line.bed.gz"

data_h1 <- read_ENCODE_bed(file_h1, verbose = T)
h1_prop <- meth_proportion(data_h1, minimum_reads = 2)

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
rr_plots(list(rr_binary_CpG_h1))

##############################################################################

file_stomach = "../../MethylationCode/MethylationData/wgbs/ENCFF844EFX_stomach_man_51.bed.gz"
data_stomach <- read_ENCODE_bed(file_stomach, verbose = T)
stomach_prop <- meth_proportion(data_stomach, minimum_reads = 2)

# On Genome
pos_stomach_chr1 <- get_methylation_positions(data_stomach, chromosome = "chr1", strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = replace_no_reads_entries)
rr_pos_stomach_chr1 <- genome_MSR(pos_stomach_chr1, minimum_bin_size = 20)
rr_plots(List(rr_pos_stomach_chr1))

# On CpG List
methylation_binary_vector_stomach <- get_methylation_CpG_binary_vector(data_stomach, chromosome = "chr1", strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = keep_nas)
rr_CpG_stomach <- calculate_relevance_resolution_vector(methylation_binary_vector_stomach, na_tolerance = 0.1, na_values_handler = replace_nas_hybrid_stochastic)
rr_plots(list(rr_CpG_stomach))


rr_plots(List(rr_pos_h1_chr1, rr_pos_stomach_chr1, rr_CpG_positions_on_Genome_chr1), legend = c("H1", "Stomach", "all CpG methylated"), title = "MSR curve on chromosome 1")
rr_plots(List(rr_binary_CpG_h1, rr_CpG_stomach), legend = c("H1", "Stomach"), title = "MSR curve on chromosome 1, only CpG list")


file_lung = "MethylationData/wgbs/ENCFF039JFT_lung.bed.gz"

data_lung <- read_ENCODE_bed(file_lung, verbose = T)
lung_prop <- meth_proportion(data_lung, minimum_reads = 2)



chr2_exp = methylation_experiment_by_chromosome(List(data_h1, data_stomach),names = c("h1", "stomach"),chromosome = "chr2")
chr2_exp$plotter()

whole_CpGlist_exp = methylation_experiment_CpGlist(List(data_h1, data_stomach, data_lung),names = c("h1", "stomach", "lung"), na_tolerance = 0.1)
whole_CpGlist_exp$plotter()


# stomach_props <- sapply(chromosomes(), function(c) meth_proportion_chromosome(data_stomach, 1,c))

rr_50_threshold <- methylation_experiment_CpGlist(List(data_h1, data_stomach),names = c("h1", "stomach"), na_tolerance = 0.1, methylation_assigner = standard_binaryzer)
rr_adaptive_threshold <- methylation_experiment_CpGlist(List(data_h1, data_stomach),names = c("h1", "stomach"), na_tolerance = 0.1, methylation_assigner = rate_preserving_threshold_binaryzer)
rr_stochastic <- methylation_experiment_CpGlist(List(data_h1, data_stomach),names = c("h1", "stomach"), na_tolerance = 0.1, methylation_assigner = stochastic_binaryzer)


rr_50_threshold_inv <- methylation_experiment_CpGlist(List(data_h1, data_stomach),names = c("h1", "stomach"), na_tolerance = 0.1, methylation_assigner = standard_binaryzer, invert = T)
rr_adaptive_threshold_inv <- methylation_experiment_CpGlist(List(data_h1, data_stomach),names = c("h1", "stomach"), na_tolerance = 0.1, methylation_assigner = rate_preserving_threshold_binaryzer, invert = T)
rr_stochastic_inv <- methylation_experiment_CpGlist(List(data_h1, data_stomach),names = c("h1", "stomach"), na_tolerance = 0.1, methylation_assigner = stochastic_binaryzer, invert = T)
#################################################################################

setwd("./Scrivania/Tesi/MethylationCode/")
directory <- "MethylationData/binary_rate/"
cell_files <- list.files(directory,pattern="(_converted.Rda)$")
cell_names <- sub("_converted.Rda","", cell_files)
methylation_vector <- as.logical(readRDS("MethylationData/binary_rate/GSM3436261_O1_TA_Hi_10_converted.Rda"))

#################################################################################
