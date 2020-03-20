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
setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
source("WGBS_analysis_functions.R", chdir = T)

#source_directory = "./Scrivania/Tesi/MethylationCode/"
#setwd(source_directory)
file_h1 = "../../MethylationCode/MethylationData/wgbs/ENCFF601NBW_H1_cell_line.bed.gz"
data_h1 <- read_ENCODE_bed(file_h1, verbose = T)
h1_prop <- meth_proportion(data_h1, minimum_reads = 2)

file_stomach = "../../MethylationCode/MethylationData/wgbs/ENCFF844EFX_stomach_man_51.bed.gz"
data_stomach <- read_ENCODE_bed(file_stomach, verbose = T)
stomach_prop <- meth_proportion(data_stomach, minimum_reads = 2)

file_lung = "../../MethylationData/wgbs/ENCFF039JFT_lung.bed.gz"
data_lung <- read_ENCODE_bed(file_lung, verbose = T)
lung_prop <- meth_proportion(data_lung, minimum_reads = 2)

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

find_few_na_window <- function(binary_vector, size, rate = 0.1)
{
  start <- 1
  l <- length(binary_vector)
  max_nas <- rate*size
  while(start<(l-size))
  {
    nas <- sum(is.na(binary_vector[start:(start+size)]))
    #cat("start:", start, " nas: ", nas)
    if(nas>max_nas) start = start + (nas-max_nas)
    else return(start)
  }
  
  cat("\nnot found")
  return(0)
}



different_scales_experiment <- function(data_list, size_list, names, na_tolerance, strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = keep_nas, invert = F)
{

  rr_list = lapply(data_list, function(d) 
  {
    binary <- get_methylation_CpG_binary_vector(d, chromosome = "all", strands_handler = strands_handler, methylation_assigner = methylation_assigner, missing_read_handler = missing_read_handler)
    rr_dim_list = lapply(size_list, function(s)
    {
      start <- find_few_na_window(binary, s, rate = 0.05)
      M = sum(binary[(1+start):(s+start)], na.rm = T)
      cat("M:", M, " M/size", M/s, "\n")
      calculate_relevance_resolution_vector(binary[(1+start):(s+start)], na_tolerance = na_tolerance, na_values_handler = replace_nas_hybrid_stochastic, invert = invert)
    })
    return(rr_dim_list)
  })
  
  return(rr_list)
}


different_scales_experiment_genomewide <- function(data_list, size_list, offset=0, chromosome, names, na_tolerance, strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = replace_no_reads_entries, invert = F)
{
  
    rr_list = lapply(data_list, function(d) 
    {
      pos <- get_methylation_positions(d, chromosome, strands_handler, methylation_assigner, missing_read_handler)
      print(max(pos))
      rr_dim_list = mclapply(size_list, mc.cores = 1, function(s)
      {
        new_pos = (pos[pos>offset & pos<(offset+s)])
        cat("M: ", length(new_pos), "M/size:", length(new_pos)/s, "\n")
        genome_MSR(new_pos, minimum_bin_size = 20, invert = F, verbose = F)
      })
      
    })
  
  return(rr_list)
}

different_positions_scale_CpG_list_experiment <- function(data_list, size, names, na_tolerance, strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = keep_nas, invert = F, undersample = 0)
{

  rr_list = lapply(data_list, function(d) 
  {
    binary <- get_methylation_CpG_binary_vector(d, chromosome = "all", strands_handler = strands_handler, methylation_assigner = methylation_assigner, missing_read_handler = missing_read_handler)
    
    l <- length(binary)
    fragments <- floor(l/size)
    start_list <- ((0:(fragments-1))*size)+1
    cat("fragments: ", fragments, "\n")
    
    if(undersample!=0)
      start_list <- sample(start_list, undersample, replace=F)
    
    rr_fragments_list = lapply(start_list, function(s)
    {
      calculate_relevance_resolution_vector(binary[(s):(size+s)], na_tolerance = na_tolerance, na_values_handler = replace_nas_hybrid_stochastic, invert = invert)
    })
    return(rr_fragments_list)
  })
  
  return(rr_list)
}




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
size = 1e6
names = c("H1, stomach")
rr_fragments <- different_positions_scale_CpG_list_experiment(data_list, size, names, na_tolerance=0.1, strands_handler = sum_strands, methylation_assigner = stochastic_binaryzer, missing_read_handler = keep_nas, invert = T, undersample = 0)
compare_resolution_relevance_plot(rr_fragments[[1]], title = "MSR plot for 29 different fragments of 1e6 sites of stomach cells (CpG list, stochastic assignment, not meth)")

areas_h1 = sapply(rr_fragments[[1]], function(x) MSR_area(x))
areas_stomach = sapply(rr_fragments[[2]], function(x) MSR_area(x))

plot(areas_stomach, areas_h1)
plot(areas_stomach)
mean(areas_stomach)
mean(areas_h1, na.rm=T)
plot(areas_h1)
boxplot.default(areas_stomach, areas_h1, names = c("stomach", "H1"))
title("MSR area calculated for 290 different fragments of size 1e6 (CpG list, stochastic assignment, not methilated)")



#################################################################################


### FRAGMENTS EXPERIMENTS ##########################################################################

different_positions_scale_by_chromosome_experiment <- function(data_list, size, names, chromosome, strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = replace_no_reads_entries, invert = F, undersample = 0)
{
  
  rr_list = lapply(data_list, function(d) 
  {
    pos <- get_methylation_positions(d, chromosome, strands_handler, methylation_assigner, missing_read_handler)
    
    l <- max(pos)
    fragments <- floor(l/size)
    start_list <- ((0:(fragments-1))*size)+1
    cat("fragments: ", fragments, "\n")
    
    if(undersample!=0)
      start_list <- sample(start_list, undersample, replace=F)
    
    rr_fragments_list = mclapply(start_list, mc.cores = 18, function(s)
    {
      cat("start: ", s, "end: ", size+s, "max: ", l, "\n")
      genome_MSR(pos[pos>s & pos<=(s+size)], invert = invert,minimum_bin_size = 20 )
    })
    return(rr_fragments_list)
  })
  
  return(rr_list)
}

file_h1 = "../../MethylationCode/MethylationData/wgbs/ENCFF601NBW_H1_cell_line.bed.gz"
file_stomach = "../../MethylationCode/MethylationData/wgbs/ENCFF844EFX_stomach_man_51.bed.gz"

data_list = List(read_ENCODE_bed(file_h1, verbose = T), read_ENCODE_bed(file_stomach, verbose = T))
size = 1e7
names = c("H1, stomach")
rr_fragments <- different_positions_scale_by_chromosome_experiment(data_list, size, names, invert = F, undersample = 0, chromosome = "chr1",strands_handler = sum_strands, methylation_assigner = stochastic_binaryzer)
compare_resolution_relevance_plot(rr_fragments[[1]], title = "MSR plot for x different fragments of 1e6 sites of stomach cells (CpG list, stochastic assignment)")

areas_h1 = sapply(rr_fragments[[1]], function(x) MSR_area(x))
areas_stomach = sapply(rr_fragments[[2]], function(x) MSR_area(x))

plot(areas_stomach, areas_h1)
plot(areas_stomach)
plot(areas_h1)
boxplot.default(areas_stomach, areas_h1, names = c("stomach", "H1"))
title("MSR area calculated for x different fragments of size 1e5 (CpG list, stochastic assignment)")



#################################################################################

setwd("./Scrivania/Tesi/MethylationCode/")
directory <- "MethylationData/binary_rate/"
cell_files <- list.files(directory,pattern="(_converted.Rda)$")
cell_names <- sub("_converted.Rda","", cell_files)
methylation_vector <- as.logical(readRDS("MethylationData/binary_rate/GSM3436261_O1_TA_Hi_10_converted.Rda"))

#################################################################################
