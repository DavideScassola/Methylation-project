#setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
source("WGBS_analysis_functions.R", chdir = T)

# HOW TO GENERATE CPG ISL DATA
############################################################
#file = "../../MethylationCode/MethylationData/cpgIslandExt.txt.gz"
#CpGislands <- fread(cmd=sprintf("zcat < %s",file), verbose=F, showProgress=T, stringsAsFactors = T)
#colnames(CpGislands) <- c("bin","chr", "start", "end", "name", "length", "cpgNum", "gcNum", "perCpG", "perGc", "obsExp")
#save(file = "../../MethylationCode/MethylationData/CpGislands.Rdata", CpGislands)
############################################################

############################################################
#file = "../../MethylationCode/MethylationData/enhancer.bed"
#Enhancers <- fread(file = file,verbose=F, showProgress=T, stringsAsFactors = T)
#colnames(CpGislands) <- c("bin","chr", "start", "end", "name", "length", "cpgNum", "gcNum", "perCpG", "perGc", "obsExp")
#save(file = "../../MethylationCode/MethylationData/Enhancers.Rdata", Enhancers)
############################################################


to_chr_factor <- Vectorize(function(chr_number)
{
  paste("chr", as.character(chr_number), sep="")
})


#load("../../MethylationCode/MethylationData/CpGislands.Rdata")

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


annotation_region_meth_counter1 <- function(annotation_region) 
{
  isl = annotation_region[reads>2, prop]
  meth_count = sum(isl>50)
  c(mean(isl), meth_count, length(isl))
}


annotation_regions_meth_data <- function(data, annotation_region_meth_counter, Annotations_dataframe, cores = 1)
{
  #load("../../MethylationCode/MethylationData/CpGislands.Rdata")
  
    d = sum_strands(data)
    l = length(Annotations_dataframe$start)
    cat("\nprocessing . . .\n")
    methylation_prop = mcmapply(1:l, mc.preschedule = T, mc.cores = cores, FUN =  function(n)
    {
      cat(n," ")
      annotation_region_meth_counter(get_annotation_region(n, d, Annotations_dataframe))
    })
    
    r = data.frame(t(methylation_prop))
    colnames(r) <- c("prop", "meth count", "valid sites")
    #r = cbind(r, id=Annotations_dataframe$id)
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

add_wgbs_indexes <- function(Annotations_dataframe, wgbs_data, cores = 1)
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

#########################################


from_bed_to_annotation_with_wgbs_indexes <- function(bed_file, wgbs_data_file, cores = 1)
{
  wgbs_data = read_ENCODE_bed(wgbs_data_file)
  anno <- fread(file = bed_file,verbose=F, showProgress=T, stringsAsFactors = T)
  anno$chr <- to_chr_factor(anno$chr)
  
  anno_improved = add_wgbs_indexes(anno, wgbs_data, cores)
  new_name = paste(substring(bed_file, 1, nchar(bed_file)-4), "_improved.Rda", sep = "")
  saveRDS(anno_improved, file = new_name)
}

