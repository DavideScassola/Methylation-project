#################################################################################################
setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
source("WGBS_analysis_functions.R", chdir = T)

data_H1 <- sum_strands(readRDS("../../MethylationCode/MethylationData/wgbs/H1.Rda"))
data_stomach <- sum_strands(readRDS("../../MethylationCode/MethylationData/wgbs/stomach.Rda"))
data_K562 <- sum_strands(readRDS("../../MethylationCode/MethylationData/wgbs/K562.Rda"))
data_HeLa <- sum_strands(readRDS("../../MethylationCode/MethylationData/wgbs/HeLa_S3.Rda"))
gc(full=T)

names = c("H1", "stomach", "HeLa-S3", "K562")
data_list = List(data_H1, data_stomach, data_HeLa, data_K562)
####################################################

###################################################################
H1_islands      = readRDS("../../Rexperiments/H1_island_data.Rda")
astomach_islands = readRDS("../../Rexperiments/stomach_island_data.Rda")
K562_islands    = readRDS("../../Rexperiments/K562_island_data.Rda")
HeLa_islands    = readRDS("../../Rexperiments/HeLa_S3_island_data.Rda")
#####################################################################

###################################################################
H1_enhancers      = readRDS("../../Rexperiments/H1_enhancers_data.Rda")
stomach_enhancers = readRDS("../../Rexperiments/stomach_enhancers_data.Rda")
K562_enhancers    = readRDS("../../Rexperiments/K562_enhancers_data.Rda")
HeLa_enhancers    = readRDS("../../Rexperiments/HeLa_S3_enhancers_data.Rda")
#####################################################################

# BASE LEVEL

# Histograms
compare_meth_sites_histograms(names, min_reads = 5, data_H1, data_stomach, data_HeLa, data_K562)

# Correlation
min_reads = 10
base_level_meth_correlation(data_H1, data_stomach, min_reads)
base_level_meth_correlation(data_H1, data_Hela,    min_reads)
base_level_meth_correlation(data_H1, data_K562,    min_reads)

##### base level correlation
#             Corr       coherence
# H1-stomach: 0.6083763, 0.9021882
# H1-Hela_S3: 0.2146171, 0.6566614
# H1-K562:    0.1398305, 0.3369727




# CG ISLAND LEVEL

# Histograms
compare_meth_annotation_histograms(names, reads_name = "valid sites", title = "CG islands methylation", min_reads = 10, H1_islands, stomach_islands, HeLa_islands, K562_islands)

# Correlation
min_reads = 10
annotation_level_meth_correlation(H1_islands, stomach_islands, reads_name = "valid sites", min_reads=100, names = c("H1", "stomach"), main = "CG islands: ")
annotation_level_meth_correlation(H1_islands, HeLa_islands,    reads_name = "valid sites", min_reads=100, names = c("H1", "HeLa-S3"), main = "CG islands: ")
annotation_level_meth_correlation(H1_islands, K562_islands,    reads_name = "valid sites", min_reads=100, names = c("H1", "K562"), main = "CG islands: ")

##### island level correlation
#             Corr       coherence
# H1-stomach: 0.9245507, 0.9495616
# H1-Hela_S3: 0.6001392, 0.7220958
# H1-K562:    0.5477659, 0.7831474




# ENHANCERS LEVEL

# Histograms
par(mfrow=c(2,2))
compare_meth_annotation_histograms(names, reads_name = "valid sites", title = "Enhancers methylation", min_reads = 5, y_max = 0.08, H1_enhancers, stomach_enhancers, HeLa_enhancers, K562_enhancers)

# Correlation
min_reads = 5
annotation_level_meth_correlation(H1_enhancers, stomach_enhancers, reads_name = "valid sites", min_reads, names = c("H1", "stomach"), main = "Enhancers: ")
annotation_level_meth_correlation(H1_enhancers, HeLa_enhancers,    reads_name = "valid sites", min_reads, names = c("H1", "HeLa-S3"), main = "Enhancers: ")
annotation_level_meth_correlation(H1_enhancers, K562_enhancers,    reads_name = "valid sites", min_reads, names = c("H1", "K562"), main = "Enhancers: ")


##### island level correlation
#             Corr       coherence
# H1-stomach: 0.527691   0.8009343
# H1-Hela_S3: 0.2983441  0.7712936
# H1-K562:    0.1993023  0.3655242





#############################################################

cgi_anno <- readRDS("../../Rexperiments/cgi_improved.Rda")
cgi_anno = cgi_anno[!is.infinite(i_start)]


d1 = data_H1
d2 = data_K562
name1 = "H1"
name2 = "K562"

min_reads = 10






