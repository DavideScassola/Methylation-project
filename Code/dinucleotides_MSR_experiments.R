
################################################# CpG sites analysis
source("MSR_analysis_functions.R", chdir = T)
setwd("Scrivania/Tesi/Methylation-project/Code/")

MSR_fixed_M_exp <- function(positions, M, offset=0, names)
{
  rr_list <- lapply(positions, function(p)
  {
    genome_MSR(p[(1+offset):(M+offset)], verbose = T, minimum_bin_size = 2)
  })
  
  plotter <- function()
  {
    rr_plots(rr_list, legend = names, title = paste("dinucleotides comparison, fixed M:",M))
  }
  
  M_densities <- lapply(positions, function(p)
  {
    M/(max(p[(1+offset):(M+offset)])-min(p[(1+offset):(M+offset)]))
  })
  
  
  return(List(rr_list=rr_list, plotter=plotter, M_densities = M_densities))
}

result <- mouse_dinucleotides_experiment(c("CG", "TA", "CA", "CC"), chromosome = "chr1")
rr_plots(result$rr, title = "dinucleotides MSR comparison", legend_names = c("CG", "TA", "CA", "TG"))


a <- nucleotides_pattern_positions("chr1", "CG", Genome = BSgenome.Hsapiens.UCSC.hg38)

b <- nucleotides_pattern_positions("chr1", "CG", Genome = BSgenome.Mmusculus.UCSC.mm10)

ranges <- matchPattern(pattern,(BSgenome.Hsapiens.UCSC.hg38))


patterns = c("CG", "TA", "CA", "CC")
chromosome <- "chr1"
positions = lapply(patterns, function(p) nucleotides_pattern_positions(chromosome, p, Genome = BSgenome.Hsapiens.UCSC.hg38))

rr_fixed_e4 <- MSR_fixed_M_exp(positions, 1E4, offset = 0, names = patterns)
rr_fixed_e6 <- MSR_fixed_M_exp(positions, 1E6, offset = 0, names = patterns)

rr_fixed_M$plotter()
