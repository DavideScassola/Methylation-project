

setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
source("CpG_islands_functions.R", chdir = T)

############## ULISSE ################
setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
source("CpG_islands_functions.R", chdir = T)
file_stomach = "../../MethylationCode/MethylationData/wgbs/ENCFF844EFX_stomach_man_51.bed.gz"
data_stomach <- read_ENCODE_bed(file_stomach, verbose = T)

a = islands_meth_data(data_stomach, island_meth_counter1, cores = 18)
saveRDS(a, "stomach_island.Rda")
####################################### 

p = data_stomach[data_stomach$strand=="+"]$prop
m = data_stomach[data_stomach$strand=="-"]$prop
p_reads = data_stomach[data_stomach$strand=="+"]$reads
m_reads = data_stomach[data_stomach$strand=="-"]$reads

rm(data_stomach)
gc()
c = data.frame(plus=p, minus=m, p_reads=p_reads, m_reads=m_reads)

min_reads = 10
ps = data_stomach[reads>min_reads, ]$prop
ph = data_h1[reads>min_reads, ]$prop
pl = data_lung[reads>min_reads, ]$prop




hist(ph, col=rgb(0,0,1,0.5),probability = T, main = "methylation proportion distribution", xlab = "proportion", breaks = 20)
hist(ps, col=rgb(1,0,0,0.5), probability = T, add=T, breaks = 20)
box()

s = sum_strands(data_stomach)
mr = 1
s_sig = s[reads>mr]


a = sapply(1:100, function(lag){
  cat(lag, " ")
  cor.test(s$prop[1:(length(s$prop)-lag)], s$prop[(1+lag):length(s$prop)])$estimate
})