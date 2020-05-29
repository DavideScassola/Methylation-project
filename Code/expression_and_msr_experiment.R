setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")

#source("WGBS_analysis_functions.R", chdir = T)
source("Methylation-project/Code/WGBS_analysis_functions.R", chdir = T)
##################################
#stomach_fragments_rna_tables = lapply(c(1e3,1e4,1e5,1e6), function(window)
#{
#  make_rna_window_data_frame(wgbs_data, rna_data, genebody_annotation, window)
#})
##################################

library("PerformanceAnalytics")
##########################
file_rna_stomach = "../../MethylationCode/MethylationData/rna-seq/ENCFF918KPC_stomach.tsv"
stomach_rna <- read_rna_file(file_rna_stomach)
load("../../Rexperiments/stomach_fragments_table.Rdata")
load("../../Rexperiments/stomach_fragments_table_discrete.Rdata")
load("../../Rexperiments/stomach_fragments_table_adaptive.Rdata")
load("../../Rexperiments/stomach_rna_fragment_tables.Rdata")
data_stomach <- readRDS("../../MethylationCode/MethylationData/wgbs/stomach.Rda")
data_stomach <- sum_strands(data_stomach)
genebody = readRDS("../../Rexperiments/genebody_improved.Rda")
#########################


#------------------------
y <- c("gene_count", "genes_nucleotides_count", "total_TPM")
basic_predictors <- c("nucleotides", "CpG_density", "meth rate")
essentials <- c(y, basic_predictors)
essential_msr_predictors <- c("msr", "inverted_msr")
fancy_msr_predictors <- c("ecdf(msr, density)", "inverted ecdf(msr, density)", "residual", "inverted_residual")
msr_predictors <- c(essential_msr_predictors, fancy_msr_predictors)
#------------------------


table <- join_rna_and_msr_tables(stomach_rna_fragment_tables, stomach_fragments_table_adaptive, 2)
table <- exclude_outliers(table, lim = 2e6)

gene_mask =table$gene_count>0
boxplot(table$msr[!gene_mask], table$msr[gene_mask])
boxplot(table$inverted_msr[!gene_mask], table$inverted_msr[gene_mask])


chart.Correlation(table[,y])
chart.Correlation(table[,essentials])
chart.Correlation(table[,c(essentials, essential_msr_predictors)])
chart.Correlation(table[,c(essentials, msr_predictors)])


binary_predictivity_hist(table$inverted_msr, table$gene_count>0, xlab = "inverted MSR", sep_names = c("genes", "no genes"), colors = c(5,7), breaks = 10, main = "Can inverted msr predict the presence of a gene?", max_y = 30)

#gene_mask = table$gene_count>0
#feature = table$inverted_msr
#feature = table$`ecdf(msr, density)`
#feature = table$`inverted ecdf(msr, density)`


gene_0 = data_table$gene_count==0
gene_1 = data_table$gene_count==1
gene_2 = data_table$gene_count==2
gene_many = data_table$gene_count>2

boxplot(data_table$inverted_msr[gene_0], data_table$inverted_msr[gene_1], data_table$inverted_msr[gene_2], data_table$inverted_msr[gene_many], names = c("No genes", "one gene", "two genes", "more than 2 genes"), ylab = "inverted MSR")
title("Inverted MSR and number of genes")


modello = gam( data_table[ data_table$gene_count>0,]$log_tpm ~
                 s(data_table[ data_table$gene_count>0, 4]) +
                 s(data_table[ data_table$gene_count>0, 9]) +
                 s(data_table[ data_table$gene_count>0, 10]) +
                 s(data_table[ data_table$gene_count>0, 11]) +
                 s(data_table[ data_table$gene_count>0, 12]) +
                 s(data_table[ data_table$gene_count>0, 14]) +
                 s(data_table[ data_table$gene_count>0, 18])
               )




stomach_tables = produce_fragments_rna_tables(data_stomach, stomach_rna, genebody, sizes = c(1e3,1e4,1e5))
H1_tables = produce_fragments_rna_tables(data_stomach, stomach_rna, genebody, sizes = c(1e3,1e4,1e5))

