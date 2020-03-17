source("MSR_analysis_functions.R", chdir = T)

filter_strand     <- function(data, sign) { return(data[strand==sign, -c(2)]) }
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
  return(unstranded_data[,-"strand"])
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

get_methylation_positions <- function(data, chromosome, strands_handler, methylation_assigner, missing_read_handler)
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
  
  # solving various problems
  out <- clean_bed_file(data, strands_handler, methylation_assigner, missing_read_handler)
  
  cat("Methylation proportion: ", mean(out$prop, na.rm = T),"\n")
  
  v <- out$prop
  
  gc()
  return(v)
}


methylation_experiment_by_chromosome <- function(data_list, names, chromosome, strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = replace_no_reads_entries, minimum_bin_size = 20, invert = F)
{
  rr_list = lapply(data_list, function(d) 
  {
    pos <- get_methylation_positions(d, chromosome = chromosome, strands_handler = sum_strands, methylation_assigner = standard_binaryzer, missing_read_handler = replace_no_reads_entries, invert = invert)
    return(genome_MSR(pos,minimum_bin_size,T))
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

