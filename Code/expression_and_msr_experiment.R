setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
source("WGBS_analysis_functions.R", chdir = T)

remove_version_from_gene <- Vectorize(function(s)
{
  if(substr(s, 1, 4)=="ENSG")
    return(substr(s, 1, 15))
  else(return(s))
})

read_rna_file <- function(file_rna, reduced = T, correct_gene_id = T)
{
  stomach_rna <- fread(file = file_rna,verbose=F, showProgress=T, stringsAsFactors = T)
  if(reduced)
    stomach_rna = stomach_rna[,c("gene_id", "TPM")]
  
  if(correct_gene_id)
  {
    stomach_rna$gene_id = remove_version_from_gene(as.character(stomach_rna$gene_id))
  }
  
  return(stomach_rna)
}

get_genes_by_region <- function(start_chr, start_position, end_position, genebody_annotation)
{
  # fully included
  #genebody_annotation[chr==start_chr & start>=start_position & end<=end_position]$id
  
  # start included
  genebody_annotation[chr==start_chr & start>=start_position & start<=end_position]$id
}

get_TPM <- function(ids, rna_data)
{
  rna_data[gene_id %in% (ids), TPM]
}

make_msr_rna_data_frame <- function(wgbs_data, rna_data, genebody_annotation, msr_experiment_data_frame)
{
  
  wgbs_data = sum_strands(wgbs_data, verbose = F)
  i_starting_points = msr_experiment_data_frame$start
  window = i_starting_points[2]-i_starting_points[1]
  start_positions = wgbs_data[i_starting_points]$Cpos
  start_chr = wgbs_data[i_starting_points]$chr
  end_positions = wgbs_data[i_starting_points+window]$Cpos
  #end_chr = wgbs_data[i_starting_points+window]$chr
  l = length(start_positions)
  
  rna_data = rna_data[gene_id %in% genebody_annotation$id]
  
  gene_info = sapply(1:l, function(i)
  {
    cat(i, " ")
    genes = get_genes_by_region(start_chr[i], start_positions[i], end_positions[i], genebody_annotation)
    tpm = get_TPM(genes, rna_data)
    total_TPM = sum(tpm)
    c(length(genes), sum(tpm))
  })
  
  nucleotides = end_positions-start_positions
  nucleotides[nucleotides<=0] = NA
  
  out = cbind(data.frame(start_chr, start_positions, end_positions,
                         nucleotides, gene_count = gene_info[1,], total_TPM = gene_info[2,]
  ), msr_experiment_data_frame)
  
  
  return(out)
}

library("PerformanceAnalytics")
##########################
file_rna_stomach = "../../MethylationCode/MethylationData/rna-seq/ENCFF918KPC_stomach.tsv"
stomach_rna <- read_rna_file(file_rna_stomach)
load("../../Rexperiments/stomach_fragments_table.Rdata")
load("../../Rexperiments/stomach_fragments_table_discrete.Rdata")
data_stomach <- sum_strands(readRDS("../../MethylationCode/MethylationData/wgbs/stomach.Rda"))
genebody = readRDS("../../Rexperiments/genebody_improved.Rda")

##########################


a = make_msr_rna_data_frame(data_stomach, stomach_rna, genebody, stomach_fragments_table[[1]])
b = make_msr_rna_data_frame(data_stomach, stomach_rna, genebody, stomach_fragments_table_discrete[[1]])

chart.Correlation(a[, -c(1,2,3,7)])

