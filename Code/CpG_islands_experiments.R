
setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
load("../../MethylationCode/MethylationData/CpGislands.Rdata")
source("CpG_islands_functions.R", chdir = T)

###################################################################
H1_islands = readRDS("../../Rexperiments/H1_island_data.Rda")
stomach_islands = readRDS("../../Rexperiments/stomach_island_data.Rda")

#####################################################################
# correlation study
cor.test(stomach_islands$prop, H1_islands$prop)
cor.test(round(stomach_islands$prop/100), round(H1_islands$prop/100))
island_table = table(round(stomach_islands$prop/100), round(H1_islands$prop/100))

linear_model = lm(stomach_islands$prop ~ H1_islands$prop)
#glm_model = glm(stomach_islands$prop/10 ~ H1_islands$prop/100, family = "binomial")
summary(linear_model)

island_table
prop.table(island_table)
prop.table(island_table,1)
prop.table(island_table,2)
coherent_sites_prop = prop.table(island_table)[1] + prop.table(island_table)[4]
cat("coherent_sites_prop: ", coherent_sites_prop)
different_sites_num = (island_table)[2] + (island_table)[3]
cat("different_sites_num: ", different_sites_num)

######################################################################
# Methylation heterogeneity
names = c("Stomach", "H1")
colors = c(alpha(10,0.5), alpha(5,0.5))

breaks = 40
min_valid_sites = 10
hist(H1_islands$prop[H1_islands$`valid sites`>min_valid_sites], col = colors[2], probability = T, breaks = breaks, xlab = "methylation level", main = ("Methylation level on CpG islands"))
hist(stomach_islands$prop[stomach_islands$`valid sites`>min_valid_sites], col = colors[1], probability = T, add = T, breaks = breaks)
legend("top", legend=names, col=colors, fill = colors)

breaks = 4
hist(H1_islands$prop[H1_islands$`valid sites`>min_valid_sites], col = colors[2], probability = T, breaks = breaks, xlab = "methylation level", main = ("Methylation level on CpG islands"))
hist(stomach_islands$prop[stomach_islands$`valid sites`>min_valid_sites], col = colors[1], probability = T, breaks = breaks, add = T)
legend("top", legend=names, col=colors, fill = colors)

######################################################################




file_stomach = "../../MethylationCode/MethylationData/wgbs/ENCFF844EFX_stomach_man_51.bed.gz"
data_stomach <- read_ENCODE_bed(file_stomach, verbose = T)

#####################################################################
# sanity checks
n = 300
show_island(n,data_stomach, 3)
stomach_islands[n,]

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