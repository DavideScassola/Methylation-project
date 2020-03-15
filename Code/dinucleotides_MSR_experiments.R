
################################################# CpG sites analysis
source("MSR_analysis_functions.R", chdir = T)

result <- mouse_dinucleotides_experiment(c("CG", "TA", "CA", "CC"), chromosome = "chr1")
rr_plots(result$rr, title = "dinucleotides MSR comparison", legend_names = c("CG", "TA", "CA", "TG"))

a <- nucleotides_pattern_positions("chr1", "CG", Genome = BSgenome.Hsapiens.UCSC.hg38)
b <- nucleotides_pattern_positions("chr1", "CG", Genome = BSgenome.Mmusculus.UCSC.mm10)

ranges <- matchPattern(pattern,(BSgenome.Hsapiens.UCSC.hg38))

