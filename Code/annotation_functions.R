setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
source("WGBS_analysis_functions.R", chdir = T)


# plot utils

to_chr_factor <- Vectorize(function(chr_number)
{
  paste("chr", as.character(chr_number), sep="")
})

filter_annotation_region <- function(chromosome, start, end, data)
{
  data[chr==chromosome & Cpos>=start & Cpos<=end]
}

get_annotation_region <- function(n, data, Annotations_dataframe)
{
  i = Annotations_dataframe[n]
  filter_annotation_region(as.character(i$chr), i$start, i$end, data)
}

show_annotation_region <- function(n, data, Annotations_dataframe,  min_reads = 2)
{
  d = get_annotation_region(n, data, Annotations_dataframe)
  d[reads<min_reads]$prop = NA
  minus <- d$strand=="-"
  plus <- d$strand=="+"
  xlabel = paste("nucleotide position on ", as.character(Annotations_dataframe[n]$chr))
  plot(d[plus]$Cpos, d[plus]$prop, ylim = c(0,100), col = 2, xlab = xlabel, ylab = "observed meth proportion")
  points(d[minus]$Cpos, d[minus]$prop, col = 4)
  points(d[reads<min_reads]$Cpos, rep_len(50, length(d[reads<min_reads]$Cpos)), col = "gray65", pch = "|")
  lines(x = c(0,9999999999), y = c(50,50), lty = 2)
  
  cat("annotation region infos:\n ")
  print(Annotations_dataframe[n])
}

compare_annotation_region <- function(n, Annotations_dataframe, min_reads = 5, ...)
{
  par(mfrow=c(length(list(...)),1)) 
  
  for(d in list(...))
  {
    show_annotation_region(n, d, Annotations_dataframe, min_reads = min_reads)
  }

}

####################################


# annotation utils

add_wgbs_indexes <- function(Annotations_dataframe, wgbs_data, cores = 1)
{
  wgbs_data = sum_strands(wgbs_data)
  data_scheme <- (wgbs_data)[, c("chr","Cpos")]
  l = length(Annotations_dataframe$start)
  
  chr_set <- levels(data_scheme$chr[1])
  
  anno = Annotations_dataframe[Annotations_dataframe$chr %in% chr_set]
  anno$chr = factor(anno$chr)

  wgbs_chr = (data_scheme$chr)
  anno_chr = anno$chr
  anno_start = anno$start
  anno_end = anno$end
  ranges = mcmapply(1:l, mc.preschedule = F, mc.cores = cores, FUN =  function(n)
  {
    cat(n," ")
    chr_mask = wgbs_chr==as.character(anno_chr[n])
    indexes = which(chr_mask & data_scheme$Cpos>=anno_start[n]  & data_scheme$Cpos<=anno_end[n])
    i_start = min(indexes)
    i_end = max(indexes )
    return(c(i_start, i_end))
    })
  
  r = data.frame(t(ranges))
  colnames(r) <- c("i_start", "i_end")
  #print(r)
  
  out = cbind(anno, i_start = r$i_start, i_end = r$i_end)
  gc()
  return(out[order(out$i_start)])
  
}

add_wgbs_indexes2 <- function(Annotations_dataframe, wgbs_data, cores = 1)
{
  wgbs_data = sum_strands(wgbs_data)
  data_scheme <- (wgbs_data)[, c("chr","Cpos")]
  l = length(Annotations_dataframe$start)
  
  anno = Annotations_dataframe[Annotations_dataframe$chr %in% levels(data_scheme$chr[1])]
  anno$chr = factor(anno$chr)
  
  wgbs_chr = (data_scheme$chr)
  anno_chr = anno$chr
  anno_start = anno$start
  anno_end = anno$end
  
  anno$i_start <- anno$start
  anno$i_start <- NA
  anno$i_end <-anno$end
  anno$i_end <- NA
  
  s = 1
  inside_gene <- F
  for(i in 1:length(wgbs_data$Cpos))
  {
    if(!inside_gene & (wgbs_data$Cpos[i]>=anno_start[s]) & (as.character(anno_chr[s])==wgbs_chr[i]))
    {
      anno$i_start[s] <- i
      inside_gene <- T
    }
    
    if(inside_gene & (((wgbs_data$Cpos[i]>anno_end[s]) & (as.character(anno_chr[s])==wgbs_chr[i])) | (as.character(anno_chr[s])!=wgbs_chr[i])))
    {
      anno$i_end[s] <- i
      inside_gene <- F
      s <- s + 1
    }
    if(s>l)
      break
  }

  gc()
  return(anno)
  
}

add_wgbs_indexes3 <- function(Annotations_dataframe, wgbs_data, cores = 1)
{
  wgbs_data = sum_strands(wgbs_data)
  data_scheme <- (wgbs_data)[, c("chr","Cpos")]
  data_scheme$index <- 1:length(data_scheme$Cpos)
  
  
  chr_set <- levels(data_scheme$chr[1])
  anno = Annotations_dataframe[Annotations_dataframe$chr %in% chr_set, ]
  l = length(anno$start)
  anno$chr = factor(anno$chr)
  
  wgbs_chr = (data_scheme$chr)
  anno_chr = anno$chr
  anno_start = anno$start
  anno_end = anno$end
  
  anno$i_start <- anno$start
  anno$i_start <- NA
  anno$i_end <-anno$end
  anno$i_end <- NA
  
  actual_chr <- "ciao"
  
  for(g in 1:l)
  {
    if(anno_chr[g]!=actual_chr)
    {
      actual_chr <- anno_chr[g]
      actual_data_scheme <- data_scheme[chr==as.character(actual_chr), ]
    }
    
    cat(" ", g)
    indexes <- actual_data_scheme$index[actual_data_scheme$Cpos>=anno_start[g]  & actual_data_scheme$Cpos<=anno_end[g]]
    if(length(indexes)>0)
    {
      anno$i_start[g] = indexes[1]
      #actual_data_scheme <- actual_data_scheme[indexes[1]:length(actual_data_scheme$Cpos)]
      anno$i_end[g] = indexes[length(indexes)]
    }

  }
  
  return(anno)
}

from_bed_to_annotation_with_wgbs_indexes <- function(bed_file, wgbs_data_file, cores = 1)
{
  wgbs_data <- NULL
  if(file_ext(wgbs_data_file)=="gz")
    wgbs_data <- read_ENCODE_bed(wgbs_data_file)
  else
    wgbs_data <- readRDS(wgbs_data_file)

  anno <- fread(file = bed_file,verbose=F, showProgress=T, stringsAsFactors = T)
  anno$chr <- to_chr_factor(anno$chr)

  anno_improved = add_wgbs_indexes(anno, wgbs_data, cores)
  new_name = paste(substring(bed_file, 1, nchar(bed_file)-4), "_improved.Rda", sep = "")
  saveRDS(anno_improved, file = new_name)
}
