source("MSR_analysis_functions.R", chdir = T)

filter_strand     <- function(data, sign) { return(data[strand==sign, -"strand"]) }
pick_plus_strand  <- function(data)       { return(filter_strand(data,"+"))}
pick_minus_strand <- function(data)       { return(filter_strand(data,"-"))}
ignore_strands    <- function(data)       { return(data[, -"strand"]) }
chromosomes <- function(Genome = BSgenome.Hsapiens.UCSC.hg38) { return(BSgenome.Hsapiens.UCSC.hg38@seqinfo@seqnames[1:24])}
meth_proportion <- function(data, minimum_reads = 1) { mean(data[reads>=minimum_reads]$prop, rm.na = T) }
filter_chromosome <- function(data, chromosome) { return(data[chr == chromosome, -"chr"])}


meth_proportion_chromosome <- function(data, minimum_reads = 1, chromosome)
{
  meth_proportion(data[chr==chromosome, ], minimum_reads)
}

sum_strands <- function(data, verbose = T)
{
  unstranded_data <- pick_plus_strand(data)
  
  minus_reads <- pick_minus_strand(data)[, reads]
  minus_props <- pick_minus_strand(data)[, prop]
  plus_reads <- unstranded_data[, reads]
  plus_props <- unstranded_data[, prop]
  
  if(verbose)
  {
    #covered_positions = as.logical((plus_reads!=0) + (minus_reads!=0))
    double_strand_reads = as.logical((plus_reads!=0) * (minus_reads!=0))
    hemymethylated_sites = as.logical(double_strand_reads * (abs(plus_props-minus_props)>50))
    hemymethylation_prop <- sum(hemymethylated_sites)/sum(double_strand_reads)
    cat("Double reads proportion: ", sum(double_strand_reads)/length(double_strand_reads))
    cat("\nHemimethylation proportion on double reads: ", hemymethylation_prop, "\n")
  }
  
  new_reads <- minus_reads + plus_reads
  unstranded_data$reads = new_reads
  
  #new_reads[new_reads==0] <- 1 # avoid following division by zero
  new_prop <- ((minus_reads*minus_props) + (plus_reads*plus_props))/(new_reads)
  
  unstranded_data$prop = new_prop
  
  gc()
  #print(unstranded_data)
  #return(unstranded_data[,-"strand"])
  return(unstranded_data)
}

read_ENCODE_bed <- function(file, chromosome = "all", verbose = T)
{
  data <- fread(cmd=sprintf("zcat < %s",file), verbose=F, showProgress=verbose, select = c(1,2,6,10,11), stringsAsFactors = T)
  colnames(data) <- c("chr","Cpos", "strand", "reads", "prop")
  if(chromosome=="all")
  {
    return(droplevels(data[(chr %in% chromosomes()),]))
  }
  else return(data[chr==chromosome, -c(1)])
}

threshold_binaryzer <- function(data, threshold = 50)
{
  binarized_rates <- data
  binarized_rates$prop = (data$prop >= threshold)
  return(binarized_rates)
}

find_rate_preserving_threshold <- function(props, verbose = T)
{
  rate <- mean(props, na.rm = T)/100
  best_error <- abs(mean(props >= 50, na.rm = T)-rate)
  best_t <- 50
  
  for(t in 51:80)
  {
    new_error <- abs(mean(props >= t, na.rm = T)-rate)
    if(new_error>best_error) return(t)
    else 
      {
        best_error <- new_error
        best_t <- t
      }
  }
  
  return(best_t)
}

rate_preserving_threshold_binaryzer <- function(data, verbose = T)
{
  t <- find_rate_preserving_threshold(data$prop)
  if(verbose) cat("Threshold: ", t, "\n")
  return(threshold_binaryzer(data, t))
}


stochastic_binaryzer <- function(data)
{
  binarized_rates <- data
  not_nas <- !is.na(data$prop)
  binarized_rates$prop[not_nas] <- as.logical(rbinom(sum(not_nas), 1, data$prop[not_nas]/100))
  return(binarized_rates)
}

standard_binaryzer <- function(data) { threshold_binaryzer(data, 50)}


discard_few_reads_entries <- function(data, threshold = 1, verbose = T) 
{
  if(verbose) cat("Missing data: ", sum(data$reads<threshold)/length(data$reads), "\n")
  return(data[(reads)>=threshold,])
}

discard_no_reads_entries <- function(data) { return(discard_few_reads_entries(data,1)) }

replace_few_reads_entries <- function(data, threshold = 1, verbose = T)
{
  if(verbose) cat("Missing data: ", sum(data$reads<threshold)/length(data$reads), "\n")
  mean_prop <- mean(data$prop,na.rm = T)/max(data$prop,na.rm = T)
  data_with_replacements <- data
  selector <- (data$reads<threshold)
  data_with_replacements$prop[selector] = as.logical(rbinom(sum(selector), 1, mean_prop))
  data_with_replacements$prop = as.logical(data_with_replacements$prop)
  gc()
  return(data_with_replacements[, -"reads"])
}

replace_no_reads_entries <- function(data) {replace_few_reads_entries(data, 1)}

keep_nas <- function(data)
{
  data[(data$reads == 0),]$prop <- NA
  cat("Keeping missing data proportion: ", sum((data$reads == 0)/length(data$reads)), "\n")
  return(data[,-"reads"])
}

clean_bed_file <- function(data, strands_handler, methylation_assigner, missing_read_handler)
{
  return(missing_read_handler(methylation_assigner(strands_handler(data))))
}

ouget_methylation_positions <- function(data, chromosome, strands_handler, methylation_assigner, missing_read_handler)
{
  #filter a chromosome
  out <- filter_chromosome(data, chromosome)
  
  # solving various problems
  out <- clean_bed_file(out, strands_handler, methylation_assigner, missing_read_handler)
  
  cat("Methylation proportion: ", sum(out$prop)/length(out$prop),"\n")
  
  positions <- out[prop==T]$Cpos
  
  gc()
  return(positions)
}

get_methylation_CpG_binary_vector <- function(data, chromosome = "all", strands_handler, methylation_assigner, missing_read_handler)
{
  #filter a chromosome
  if(chromosome == "all")
    out <- data[,-"chr"]
  else
    out <- filter_chromosome(data, chromosome)
  
  #print((out))
  # solving various problems
  out <- clean_bed_file(out, strands_handler, methylation_assigner, missing_read_handler)
  
  cat("Methylation proportion: ", mean(out$prop, na.rm = T),"\n")
  
  v <- out$prop
  
  gc()
  return(v)
}


methylation_experiment_by_chromosome <- function(data_list, names, strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = keep_nas, invert = F)
{
  rr_list = lapply(data_list, function(d) 
  {
    pos <- get_methylation_positions(d, chromosome = chromosome, strands_handler = strands_handler, methylation_assigner = methylation_assigner, missing_read_handler = missing_read_handler)
    return(genome_MSR(pos,minimum_bin_size,T,invert = invert))
  })
  
  plotter <- function()
  {
    rr_plots(rr_list, legend = names, title = paste("MSR curve on", chromosome))
  }
  
  return(List(rr_list=rr_list, plotter=plotter))
}

methylation_experiment_CpGlist <- function(data_list, names, na_tolerance, strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = keep_nas, invert = F)
{
  rr_list = lapply(data_list, function(d) 
  {
    binary <- get_methylation_CpG_binary_vector(d, chromosome = "all", strands_handler = strands_handler, methylation_assigner = methylation_assigner, missing_read_handler = missing_read_handler)
    calculate_relevance_resolution_vector(binary, na_tolerance = na_tolerance, na_values_handler = replace_nas_hybrid_stochastic, invert = invert)
  })
  
  plotter <- function()
  {
    rr_plots(rr_list, legend = names, title = paste("MSR curve on CpG list of whole genome"))
  }
  
  return(List(rr_list=rr_list, plotter=plotter))
}

find_few_na_window <- function(binary_vector, size, rate = 0.1)
{
  start <- 1
  l <- length(binary_vector)
  max_nas <- rate*size
  while(start<(l-size))
  {
    nas <- sum(is.na(binary_vector[start:(start+size)]))
    #cat("start:", start, " nas: ", nas)
    if(nas>max_nas) start = start + (nas-max_nas)
    else return(start)
  }
  
  cat("\nnot found")
  return(0)
}

all_chromosome_methylation_experiment <- function(data_list, names, strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = replace_no_reads_entries, minimum_bin_size = 20)
{
  pdf("genome_wide_experiment.pdf")
  
  mclapply(chromosomes(), function(c)
    {
      rr_list = lapply(data_list, function(d) 
      {
        pos <- get_methylation_positions(d, chromosome = c, strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = replace_no_reads_entries)
        return(genome_MSR(pos,minimum_bin_size,T))
      })
      
      compare_resolution_relevance_plot(rr_list, names, title = sprintf("%s: MSR on genome"))
    })
  
  dev.off()
}

different_scales_experiment <- function(data_list, size_list, names, na_tolerance, strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = keep_nas, invert = F)
{
  
  rr_list = lapply(data_list, function(d) 
  {
    binary <- get_methylation_CpG_binary_vector(d, chromosome = "all", strands_handler = strands_handler, methylation_assigner = methylation_assigner, missing_read_handler = missing_read_handler)
    rr_dim_list = lapply(size_list, function(s)
    {
      start <- find_few_na_window(binary, s, rate = 0.05)
      M = sum(binary[(1+start):(s+start)], na.rm = T)
      cat("M:", M, " M/size", M/s, "\n")
      calculate_relevance_resolution_vector(binary[(1+start):(s+start)], na_tolerance = na_tolerance, na_values_handler = replace_nas_hybrid_stochastic, invert = invert)
    })
    return(rr_dim_list)
  })
  
  return(rr_list)
}


different_scales_experiment_genomewide <- function(data_list, size_list, offset=0, chromosome, names, na_tolerance, strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = replace_no_reads_entries, invert = F)
{
  
  rr_list = lapply(data_list, function(d) 
  {
    pos <- get_methylation_positions(d, chromosome, strands_handler, methylation_assigner, missing_read_handler)
    print(max(pos))
    rr_dim_list = mclapply(size_list, mc.cores = 1, function(s)
    {
      new_pos = (pos[pos>offset & pos<(offset+s)])
      cat("M: ", length(new_pos), "M/size:", length(new_pos)/s, "\n")
      genome_MSR(new_pos, minimum_bin_size = 20, invert = F, verbose = F)
    })
    
  })
  
  return(rr_list)
}

different_positions_scale_CpG_list_experiment <- function(data_list, size, na_tolerance, strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = keep_nas, invert = F, undersample = 0)
{
  rr_list = lapply(data_list, function(d) 
  {
    binary <- get_methylation_CpG_binary_vector(d, chromosome = "all", strands_handler = strands_handler, methylation_assigner = methylation_assigner, missing_read_handler = missing_read_handler)
    
    l <- length(binary)
    fragments <- floor(l/size)
    start_list <- ((0:(fragments-1))*size)+1
    cat("fragments: ", fragments, "\n")
    
    if(undersample!=0)
      start_list <- sample(start_list, undersample, replace=F)
    
    rr_fragments_list = lapply(start_list, function(s)
    {
      calculate_relevance_resolution_vector(binary[(s):(size+s)], na_tolerance = na_tolerance, na_values_handler = replace_nas_hybrid_stochastic, invert = invert)
    })
    return(rr_fragments_list)
  })
  
  return(rr_list)
}

different_positions_scale_by_chromosome_experiment <- function(d, size, chromosome, minimum_bin_size = 30, strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = replace_no_reads_entries, invert = F, undersample = 0)
{
  
  pos <- get_methylation_positions(d, chromosome,  strands_handler, methylation_assigner, missing_read_handler)
  
  l <- max(pos)
  fragments <- floor(l/size)
  start_list <- ((0:(fragments-1))*size)+1
  cat("fragments: ", fragments, "\n")
  
  if(undersample!=0)
    start_list <- sample(start_list, undersample, replace=F)
  
  rr_fragments_list = mclapply(start_list, mc.cores = 1, function(s)
  {
    cat("start: ", s, "end: ", size+s, "max: ", l, "\n")
    if(length(pos[pos>s & pos<=(s+size)])<100) return(NA)
    genome_MSR(pos[pos>s & pos<=(s+size)], invert = invert,minimum_bin_size = minimum_bin_size )
  })
  return(rr_fragments_list[!is.na(rr_fragments_list)])
}



