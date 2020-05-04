setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
source("WGBS_analysis_functions.R", chdir = T)

# HOW TO GENERATE CPG ISL DATA
############################################################
#file = "../../MethylationCode/MethylationData/cpgIslandExt.txt.gz"
#CpGislands <- fread(cmd=sprintf("zcat < %s",file), verbose=F, showProgress=T, stringsAsFactors = T)
#colnames(CpGislands) <- c("bin","chr", "start", "end", "name", "length", "cpgNum", "gcNum", "perCpG", "perGc", "obsExp")
#save(file = "../../MethylationCode/MethylationData/CpGislands.Rdata", CpGislands)
############################################################

load("../../MethylationCode/MethylationData/CpGislands.Rdata")

filter_island <- function(chromosome, start, end, data)
{
  data[chr==chromosome & Cpos>=start & Cpos<=end]
}

get_island <- function(n, data)
{
  i = CpGislands[n]
  filter_island(as.character(i$chr), i$start, i$end, data)
}

show_island <- function(n, data, min_reads = 2)
{
  d = get_island(n, data)
  d[reads<min_reads]$prop = NA
  minus <- d$strand=="-"
  plus <- d$strand=="+"
  xlabel = paste("nucleotide position on ", as.character(CpGislands[n]$chr))
  plot(d[plus]$Cpos, d[plus]$prop, ylim = c(0,100), col = 2, xlab = xlabel, ylab = "observed meth proportion")
  points(d[minus]$Cpos, d[minus]$prop, col = 4)
  points(d[reads<min_reads]$Cpos, rep_len(50, length(d[reads<min_reads]$Cpos)), col = "gray65", pch = "|")
  lines(x = c(0,9999999999), y = c(50,50), lty = 2)
  
  cat("island infos:\n ")
  print(CpGislands[n])
}

compare_island <- function(n, min_reads = 5, ...)
{
  par(mfrow=c(length(list(...)),1)) 
  
  for(d in list(...))
  {
    show_island(n, d, min_reads = min_reads)
  }

}


island_meth_counter1 <- function(island) 
{
  isl = island[reads>2, prop]
  meth_count = sum(isl>50)
  c(mean(isl), meth_count, length(isl))
}


islands_meth_data <- function(data, island_meth_counter, cores = 1)
{
  load("../../MethylationCode/MethylationData/CpGislands.Rdata")
  
    d = sum_strands(data)
    l = length(CpGislands$start)
    l = 10
    cat("\nprocessing . . .\n")
    methylation_prop = mcmapply(1:l, mc.preschedule = T, mc.cores = cores, FUN =  function(n)
    {
      cat(n," ")
      island_meth_counter(get_island(n, d))
    })
    
    r = data.frame(t(methylation_prop))
    colnames(r) <- c("prop", "meth count", "valid sites")
    gc()
    return(r)
}


aggregate <- function(d, bin_size, chromosome, limiter = 99999999)
{
  df = filter_chromosome(d,chromosome)
  l = max(df$Cpos)
  f = min(df$Cpos)
  len = floor((l-f)/bin_size)
  len = min(len,limiter)
  cat("len: ", len)
  starting_points = (((1:len)-1)*bin_size)+f
  sapply(starting_points, function(s)
    {
    mean(df[df$Cpos>=s & df$Cpos<(s+bin_size)]$prop, na.rm = T)
  })
}

autocor <- function(v, lag)
{
  l = length(v)
  plot(v[1:(l-lag)],v[(1+lag):l])
  cor.test(v[1:(l-lag)],v[(1+lag):l])
}

get_island_ranges <- function(data, min_reads = 2)
{
  cat("\npreprocessing . . .\n")
  d = sum_strands(data)
  d = d[reads>=min_reads, chr, Cpos]
  l_isl = length(CpGislands$start)
  l_data = length(data$prop)
  ranges = array(dim = c(l_isl,2))
  
  cat("computing indexes\n")
  
  indexes = mcmapply(1:l_isl, mc.silent = F, mc.cores = 1, FUN =  function(n)
  {
    cat(n," ")
    i_start = CpGislands[n, start]
    i_end = CpGislands[n, end]
    i_chr = CpGislands[n, chr]
    s = match()
    c(d[Cpos = ])
  })
  
}
