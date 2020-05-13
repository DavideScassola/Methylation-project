
#################################################################################################
setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
source("WGBS_analysis_functions.R", chdir = T)

file_h1 <- "../../MethylationCode/MethylationData/wgbs/H1.Rda"
data_H1 <- readRDS(file_h1)

file_stomach <- "../../MethylationCode/MethylationData/wgbs/stomach.Rda"
data_stomach <- readRDS(file_stomach)

file_K562 <- "../../MethylationCode/MethylationData/wgbs/K562.Rda"
data_K562 <- readRDS(file_K562)

file_HeLa <- "../../MethylationCoùde/MethylationData/wgbs/HeLa_S3.Rda"
data_HeLa <- readRDS(file_HeLa)

#######################################################################################################

####################################################
# MSR EXPERIMENTS ON DIFFERENT WINDOWS
load(file = "../../Rexperiments/total_exp.Rdata")
load(file = "../../Rexperiments/total_exp_fake.Rdata")
load(file = "../../Rexperiments/total_exp_chr1.Rdata")
load(file = "../../Rexperiments/total_exp_chr1_fake.Rdata")
load(file = "../../Rexperiments/CG_exp.Rdata")
#load(file = "../../Rexperiments/CG_exp_small.Rdata")
####################################################

####################################################
# MSR ECDFS
#msr_ecdf_1e3_micro = readRDS("../../MethylationCode/MethylationData/msr_ecdf_1e3_micro.Rda")
msr_ecdf_1e3 = readRDS("../../MethylationCode/MethylationData/msr_ecdf_1e3.Rda")
msr_ecdf_1e4 = readRDS("../../MethylationCode/MethylationData/msr_ecdf_1e4.Rda")
msr_ecdf_1e5 = readRDS("../../MethylationCode/MethylationData/msr_ecdf_1e5.Rda")
msr_ecdf_1e6 = readRDS("../../MethylationCode/MethylationData/msr_ecdf_1e6.Rda")
msr_ecdf = List(msr_ecdf_1e3, msr_ecdf_1e4, msr_ecdf_1e5, msr_ecdf_1e6)
####################################################

data_K562 = sum_strands(data_K562)
data_H1 = sum_strands(data_H1)
data_HeLa = sum_strands(data_HeLa)
data_stomach = sum_strands(data_stomach)
gc(full=T)

compare_meth_sites_histograms <- function(names, min_reads = 1, ...)
{
  i = 1
  for( d in list(...))
  {
    p = d$prop[d$reads>=min_reads]
    hist(p, main = paste(names[i], ": methylation proportion on single sites,", "mean: ", round(mean(p, na.rm = T), 2)), probability = T, ylim = c(0, 0.12))
    i = i + 1
    par(ask=TRUE)
  }
  
  par(ask=FALSE)

}

########################################################

box_comparison <- function(exp, explanation, density = F, names = c("1e3", "1e4", "1e5", "1e6"))
{
  i = 3
  ylab = "msr"
  if(density)
    {
    i = 2
    ylab = "density"
  }
  
  msr_1e3 = exp[[1]]$data$fragments_infos_array[,i]
  msr_1e4 = exp[[2]]$data$fragments_infos_array[,i]
  msr_1e5 = exp[[3]]$data$fragments_infos_array[,i]
  msr_1e6 = exp[[4]]$data$fragments_infos_array[,i]
  
  boxplot(msr_1e3, msr_1e4, msr_1e5, msr_1e6, names = names , xlab = "window size", ylab = ylab)
  title(explanation)
}


density_MSR_correlation <- function(exp, explanation, fix = T, windows = c("1e3", "1e4", "1e5", "1e6"))
{
  
  for(i in 1:length(windows))
  {
    density= exp[[i]]$data$fragments_infos_array[,2]
    msr = exp[[i]]$data$fragments_infos_array[,3]
    if(fix)
    {
      plot(density,msr, xlim = c(0,1), ylim = c(0,0.32), col = alpha(1,0.5))
    }
    else
    {
      plot(density,msr, col = alpha(1,0.5))
    }
    
    grid()
    
    
    alpha = 1e-5
    cls = c(2,3)
    add_msr_confidence_line(msr_ecdf[[i]], confidence = 1-alpha/2, col = cls[1])
    add_msr_confidence_line(msr_ecdf[[i]], confidence = alpha/2, col = cls[1])
    add_msr_confidence_line(msr_ecdf[[i]], confidence = 0.5, col = cls[2])
    axis(at = linspace(0,1,21), side = 1, labels = F)
    axis(at = (linspace(0,0.33,34)), side = 2, labels = F)
    legend("bottomleft", legend=c(paste(1-alpha, "confidence interval"), "median"), col=cls, lty = 2, title = "comparison with random data", y.intersp = 0.5,  cex=1 )

    title(paste(explanation, windows[i], "  corr:", round(cor.test(density, msr)$estimate,2)))
    par(ask=TRUE)
    print(cor.test(density,msr))
  }
  par(ask=FALSE)
}


show_zone <- function(exp, n)
{
  #autocor(exp$data$fragments_infos_array[,3],1)
  plot(exp$data$fragments_infos_array[1:n,3], col = alpha(1,0.2), ylab = "msr", xlab = "slice number", main = sprintf("%s cell, window size: %d, inverted: %s", exp$name, exp$window_size, exp$inverted))
  lines(exp$data$fragments_infos_array[1:n,3], col = alpha(1,0.7))
  par(ask=TRUE)
  
  plot(1-exp$data$fragments_infos_array[1:n,2],  col = alpha(1,0.2), ylab = "density", xlab = "slice number", main = sprintf("%s cell, window size: %d", exp$name, exp$window_size))
  lines(1-exp$data$fragments_infos_array[1:n,2], col= alpha(1,0.7))
  par(ask=FALSE)
  
  cat("\nmean msr: ", mean(exp$data$fragments_infos_array[,3], na.rm=T))
  cat("\nmean density:", mean(1-exp$data$fragments_infos_array[,2], na.rm=T))
}

significance_measure <- function(exp_datas, msr_ecdfs, i)
{
  msr_ecdf = msr_ecdfs[[i]]
  exp_data = exp_datas[[i]]
  density= exp_data$data$fragments_infos_array[,2]
  msr = exp_data$data$fragments_infos_array[,3]
  
  l = length(density)
  zero = function(x) {0}
  
  if(msr_ecdf[[1]]$prop==0) msr_ecdf[[1]]$cdf = zero
  if(msr_ecdf[[length(msr_ecdf)]]$prop==1) msr_ecdf[[length(msr_ecdf)]]$cdf = zero
  
  sapply( 1:l ,function(n)
  {
    cat(n, " ")
    if(is.na(density[n]) || is.na(msr[n])) return(NA)
    general_msr_cdf(msr_ecdf, density[n], msr[n])
  })
  
  #lm(msr~density)$residuals
}

experiment_residuals <- function(exp_datas, i)
{
  exp_data = exp_datas[[i]]
  density= exp_data$data$fragments_infos_array[,2]
  msr = exp_data$data$fragments_infos_array[,3]
  lm(msr~density)$residuals
}

################################################################################################

get_sig_measures <- function(exp)
{
  lapply(1:4, function(i)
  {
    significance_measure(exp, msr_ecdf,i)
  })
}

H1_inverted_FALSE_sig <- get_sig_measures(total_exp$`H1_inverted:_FALSE`)
H1_inverted_TRUE_sig <- get_sig_measures(total_exp$`H1_inverted:_TRUE`)
stomach_inverted_FALSE_sig <- get_sig_measures(total_exp$`stomach_inverted:_FALSE`)
stomach_inverted_TRUE_sig <- get_sig_measures(total_exp$`stomach_inverted:_TRUE`)

get_experiment_table <- function(i)
{
  total_exp$`H1_inverted:_FALSE`[[i]]
  total_exp$`H1_inverted:_TRUE`[[i]]
  total_exp$`stomach_inverted:_FALSE`[[i]]
  total_exp$`stomach_inverted:_TRUE`[[i]]
  
  expriment_table = List(start=, H1_density=, H1_msr=, H1_inverted_msr=, H1_sig=, H1_inverted_sig=, H1_residual=, H1_inverted_residual=, stomach_density=, stomach_msr=, stomach_inverted_msr=, stomach_sig=, stomach_inverted_sig=, stomach_residual=, stomach_inverted_residual=)
}

################################################################################################

par(mfrow=c(1,1))
density_MSR_correlation(total_exp$`H1_inverted:_FALSE`, "H1", fix = T)
box_comparison(total_exp$`H1_inverted:_FALSE`, explanation = "H1", density = T)
box_comparison(total_exp$`H1_inverted:_FALSE`, explanation = "H1", density = F)

density_MSR_correlation(total_exp$`H1_inverted:_TRUE`, "H1, inverted")
box_comparison(total_exp$`H1_inverted:_TRUE`, explanation = "H1, inverted", density = T)
box_comparison(total_exp$`H1_inverted:_TRUE`, explanation = "H1, inverted", density = F)

density_MSR_correlation(total_exp$`stomach_inverted:_FALSE`, "stomach")
box_comparison(total_exp$`stomach_inverted:_FALSE`, explanation = "stomach", density = T)
box_comparison(total_exp$`stomach_inverted:_FALSE`, explanation = "stomach", density = F)

density_MSR_correlation(total_exp$`stomach_inverted:_TRUE`, "stomach, inverted")
box_comparison(total_exp$`stomach_inverted:_TRUE`, explanation = "stomach, inverted", density = T)
box_comparison(total_exp$`stomach_inverted:_TRUE`, explanation = "stomach, inverted", density = F)


##### fake data
density_MSR_correlation(total_exp_fake$`stomach_inverted:_TRUE`, "stomach, inverted")
box_comparison(total_exp_fake$`stomach_inverted:_TRUE`, explanation = "stomach, inverted fake", density = T)
box_comparison(total_exp_fake$`stomach_inverted:_TRUE`, explanation = "stomach, inverted fake", density = F)
######

density_MSR_correlation(CG_exp$CG, "CG", windows = c("1e4", "1e5", "1e6", "1e7"), fix = F)
box_comparison(CG_exp$CG, explanation = "CG", names = c("1e4", "1e5", "1e6", "1e7"), density = T)
box_comparison(CG_exp$CG, explanation = "CG", names = c("1e4", "1e5", "1e6", "1e7"), density = F)

#density_MSR_correlation(CG_exp$CG, "CG", windows = c("500"), fix = T)
#box_comparison(CG_exp$CG, explanation = "CG", names = c("1e4", "1e5", "1e6", "1e7"), density = T)
#box_comparison(CG_exp$CG, explanation = "CG", names = c("1e4", "1e5", "1e6", "1e7"), density = F)

show_zone(total_exp$`stomach_inverted:_TRUE`[[1]],1000)
show_zone(total_exp$`H1_inverted:_TRUE`[[1]],1000)

######################################################

# macroscopic correlation
macroscopic_correlation <- function(exp1, exp2, names)
{
  windows = c("1e3", "1e4", "1e5", "1e6")
  
  for(i in 1:length(windows))
  {
    density1= exp1[[i]]$data$fragments_infos_array[,2]
    msr1 = exp1[[i]]$data$fragments_infos_array[,3]
    
    density2= exp2[[i]]$data$fragments_infos_array[,2]
    msr2 = exp2[[i]]$data$fragments_infos_array[,3]
    
    title = paste("densities comaparison", windows[i], "  corr:", round(cor.test(density1, density2)$estimate,2))
    plot(density1, density2, main = title, xlab = names[1], ylab = names[2])
    print(cor.test(density1,density2))
    
    par(ask=TRUE)
    
    title = paste("msr comaparison", windows[i], "  corr:", round(cor.test(msr1, msr2)$estimate,2))
    plot(msr1, msr2, main = title, xlab = names[1], ylab = names[2])
    print(cor.test(msr1, msr2))
  
  }
  par(ask=FALSE)
  
}


macroscopic_correlation(total_exp$`H1_inverted:_FALSE`, total_exp$`stomach_inverted:_FALSE`, names)
###################################################################

############################################
# correlation at base level methylation

data_H1 = sum_strands(data_H1)
data_stomach = sum_strands(data_stomach)
min_reads = 5
mask = data_H1$reads>=min_reads & data_stomach$reads>=min_reads
p1 = data_H1[mask, prop]
p2 = data_stomach[mask, prop]
cor.test(p1,p2)
############################################


############################################
CpG_positions <- nucleotides_pattern_positions("chr1", "CG", BSgenome.Hsapiens.UCSC.hg38)
CpG_densities <- get_CpG_densities(dinucleotides_neighborhood_ranges = c(50,100,300,500,1000,5000),data = data_stomach)
############################################



































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





