source("MSR_analysis_functions.R", chdir = T)
suppressMessages(library("tools"))

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
  if(length(data$strand)==0)
  {
    if(verbose)
      cat("no strands")
    return(data)
  }
  
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

ordered_data  <- function(data)
{
  separated = lapply(chromosomes(), function(c)
  {
    data[chr == c]
  })
  
  ordered = separated[[1]]
  for(i in 2:length(separated))
  {
    ordered = rbind(ordered, separated[[i]])
  }
  
  return(ordered)
  
}


read_ENCODE_bed <- function(file, chromosome = "all", verbose = T)
{
  data <- fread(cmd=sprintf("zcat < %s",file), verbose=F, showProgress=verbose, select = c(1,2,6,10,11), stringsAsFactors = T)
  colnames(data) <- c("chr","Cpos", "strand", "reads", "prop")
  if(chromosome=="all")
  {
    return(ordered_data(droplevels(data[(chr %in% chromosomes()),])))
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

keep_nas <- function(data, minimum_reads = 1)
{
  data[(data$reads < minimum_reads),]$prop <- NA
  data[(data$prop == 50),]$prop <- NA
  cat("Keeping missing data proportion: ", sum((is.na(data$prop))/length(data$reads)), "\n")
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

spatial_MSR_experiment_CpG_list <- function(binary, invert, window_size, fake_data, na_tolerance, na_values_handler = replace_nas_hybrid_stochastic)
{
  l <- length(binary)
  
  if(fake_data)
  {
    prop <- mean(binary, na.rm=T)
    binary <- rbinom(l, 1, prop)
  }
  
  
  fragments <- floor(l/window_size)
  start_list <- ((0:(fragments-1))*window_size)+1
  cat("fragments: ", fragments, "\n")
  
  fragments_infos_array = array(dim = c(fragments, 3))
  
  rr_fragments_list = lapply(1:fragments, function(i)
  {
    cat(i, "...\n")
    sub <- binary[start_list[i]:(start_list[i]+window_size)]
    #fragments_infos_array[i,1] <- start_list[i]
    #fragments_infos_array[i,2] <- mean(sub,na.rm=T)
    rr <- calculate_relevance_resolution_vector(sub, verbose=F, na_tolerance = na_tolerance, na_values_handler = na_values_handler, invert = invert)
    #fragments_infos_array[i,3] <- MSR_area(rr)
    return(rr)
  })
  
  fragments_infos_array[,1] <- start_list
  fragments_infos_array[,2] <-  sapply(1:fragments, function(i) 
  {
    mean(binary[start_list[i]:(start_list[i]+window_size)],na.rm=T)
  })
  if(invert) fragments_infos_array[,2]= 1 - fragments_infos_array[,2]
  
  fragments_infos_array[,3] <-  sapply(rr_fragments_list, function(rr) MSR_area(rr))
  
  return(List(fragments_infos_array=fragments_infos_array, rr_list=rr_fragments_list))
}

bernoulli_positions <- function(max, prop)
{
  unique(sort(round(runif(min = 1, max=max, n = prop*(max-1)))))
}


spatial_MSR_experiment_by_chromosome <- function(pos, window_size, fake_data, minimum_bin_size)
{
  l = max(pos)-min(pos)+1
  
  if(fake_data)
  {
    prop <- length(pos)/l
    pos  <- bernoulli_positions(l, prop)
  }
  
  fragments <- floor(l/window_size)
  start_list <- ((0:(fragments-1))*window_size)+1
  cat("fragments: ", fragments, "\n")
  
  fragments_infos_array = array(dim = c(fragments, 3))
  
  rr_fragments_list = lapply(1:fragments, function(i)
  {
    cat(i, "...\n")
    
    sub <- subset_positions(pos, start_list[i], window_size)
    
    if(length(sub)<5)
      rr <- NA
    else
      rr <- genome_MSR(sub,minimum_bin_size,verbose=F)
    
  })
  
  fragments_infos_array[,1] <- start_list
  fragments_infos_array[,2] <-  sapply(1:fragments, function(i) 
  {
    sub <- subset_positions(pos, start_list[i], window_size)
    length(sub)/window_size
  })
  
  fragments_infos_array[,3] <-  sapply(rr_fragments_list, function(rr) 
  {
    if(length(rr)<=1) return(NA)
    else MSR_area(rr)})
  
  valids <- (!is.na(fragments_infos_array[,3]))
  get_valid_rrs <- function() {return(rr_fragments_list[valids])}
  return(List(fragments_infos_array=fragments_infos_array, rr_list=rr_fragments_list, get_valid_rrs=get_valid_rrs))
}

total_spatial_experiment <- function(files, sizes, inversion, names, methylation_assigner, na_tolerance, fake_data = F, minimum_reads = 1, na_values_handler = replace_nas_with_bin_prop)
{
  result = List()
  
  for(i in 1:length(files))
  {
    data <- NULL
    if(file_ext(files[i])=="gz")
      data <- read_ENCODE_bed(files[i], verbose = T)
    else
      data <- readRDS(files[i])
    
    mrh <- function(x) keep_nas(x, minimum_reads)
    binary <- get_methylation_CpG_binary_vector(data,strands_handler = sum_strands, methylation_assigner = methylation_assigner, missing_read_handler = mrh)
    remove(data)
    gc()
    
    for(inv in inversion)
    {
      result_si = lapply(sizes, function(s)
      {
        gc()
        rrs <- spatial_MSR_experiment_CpG_list(binary, inv, s, fake_data, na_tolerance, na_values_handler)
        return(List(name=names[i], file= files[i],inverted=inv, window_size=s, data=rrs))
      })
      result[[paste(names[i], "inverted:", inv, sep = "_")]] <- result_si
    }
  }
  return(result)
}

total_spatial_experiment_by_chromosome <- function(files, sizes, chromosome, names, methylation_assigner, fake_data, minimum_bin_size = 20, mc = F)
{
  result = List()
  
  for(i in 1:length(files))
  {
    data <- read_ENCODE_bed(files[i], verbose = T)
    pos <- get_methylation_positions(data, chromosome, sum_strands, methylation_assigner = methylation_assigner, missing_read_handler = replace_no_reads_entries)
    remove(data)
    gc()
    
    cores = 1
    if(mc) cores = length(sizes)
    result_si = mclapply(sizes, mc.cores = cores, function(s)
    {
      gc()
      rrs <- spatial_MSR_experiment_by_chromosome(pos, s, fake_data, minimum_bin_size)
      return(List(name=names[i], window_size=s, data=rrs))
    })
    result[[(names[i])]] <- result_si
    
  }
  return(result)
}



total_spatial_experiment_for_dinucleotides_locations_by_chromosome <- function(pattern_list = c("CG"), sizes, chromosome, fake_data, minimum_bin_size = 20, mc = F)
{
  result = List()
  
  for(i in 1:length(pattern_list))
  {
    pos <- nucleotides_pattern_positions(chromosome, pattern_list[i], Genome = BSgenome.Hsapiens.UCSC.hg38)
    gc()
    
    cores = 1
    if(mc) cores = length(sizes)
    result_si = mclapply(sizes, mc.cores = cores, function(s)
    {
      gc()
      rrs <- spatial_MSR_experiment_by_chromosome(pos, s, fake_data, minimum_bin_size)
      return(List(name=pattern_list[i], window_size=s, data=rrs))
    })
    result[[(pattern_list[i])]] <- result_si
    
  }
  return(result)
}

rda_convert <- function(file, new_name)
{
  data <- read_ENCODE_bed(file, verbose = T)
  saveRDS(data, file = new_name)
  system(paste("rm", file))
  gc()
}

############################################
# correlation at base level methylation
base_level_meth_correlation <- function(data1, data2, min_reads = 10)
{
  mask = data1$reads>=min_reads & data2$reads>=min_reads
  p1 = data1[mask, prop]
  p2 = data2[mask, prop]
  
  if(length(p1)<1e5)
    plot(p1,p2, col=alpha(1, 0.5))
  print(cor.test(p1,p2))
  
  t = table(round(p1/100),round(p2/100))
  cat("Table:\n")
  print(t)
  
  cat("\n\nproportion table:\n")
  print(prop.table(t))
  
  coherent_sites_prop = prop.table(t)[1] + prop.table(t)[4]
  cat("\ncoherent_sites_prop: ", coherent_sites_prop)
  different_sites_num = (t)[2] + (t)[3]
  cat("\ndifferent_sites_num: ", different_sites_num)
  
}


compare_meth_sites_histograms <- function(names, min_reads = 1, ...)
{
  i = 1
  for( d in list(...))
  {
    p = d$prop[(d[,reads])>=min_reads]
    hist(p, main = paste(names[i], ": methylation proportion on single sites,", "mean: ", round(mean(p, na.rm = T), 2)), probability = T, ylim = c(0, 0.12), xlab = "methylation level")
    i = i + 1
    par(ask=TRUE)
  }
  
  par(ask=FALSE)
  
}

compare_meth_annotation_histograms <- function(names, reads_name = "reads", min_reads = 1, y_max = 0.12, title = "", ...)
{
  i = 1
  for( d in list(...))
  {
    p = d$prop[d[,reads_name]>=min_reads]
    main = paste(names[i], title, " mean: ", round(mean(p, na.rm = T), 2))
    hist(p, main = main, probability = T, ylim = c(0, y_max), xlab = "methylation level")
    i = i + 1
    par(ask=TRUE)
  }
  
  par(ask=FALSE)
  
}

annotation_level_meth_correlation <- function(data1, data2, reads_name, min_reads = 10, names, main)
{
  mask = data1[, reads_name]>=min_reads & data2[, reads_name]>=min_reads
  p1 = data1[mask, "prop"]
  p2 = data2[mask, "prop"]
  
  print(cor.test(p1,p2))
  
  t = table(round(p1/100),round(p2/100))
  cat("Table:\n")
  print(t)
  
  cat("\n\nproportion table:\n")
  print(prop.table(t))
  
  coherent_sites_prop = prop.table(t)[1] + prop.table(t)[4]
  cat("\ncoherent_sites_prop: ", coherent_sites_prop)
  different_sites_num = (t)[2] + (t)[3]
  cat("\ndifferent_sites_num: ", different_sites_num)
  
  if(length(p1)<1e5)
  {
    plot(p1,p2, col=alpha(1, 0.5), xlab = names[1], ylab = names[2], main = sprintf("%s correlation: %s, coherence: %s", main, round(cor.test(p1,p2)$estimate,2), round(coherent_sites_prop,2)))
    lines(x = c(-100,200), y = c(50,50), lty = 2, col = 2)
    lines(x = c(50,50), y = c(-100,200), lty = 2, col = 2)
  }
  
  
  get_transition_matrix <- function(wgbs_data, min_reads = 10)
  {
    wgbs_data = sum_strands(wgbs_data)
    p = wgbs_data$prop[data$reads>=min_reads]/100
    bin <- rbinom(n = length(p), size=1, prob = p)
    matrix(prop.table(autotable(bin,1),1), nrow = 2)
  }
  
  
}


#####################################################

remove_version_from_gene <- Vectorize(function(s)
{
  if(substr(s, 1, 4)=="ENSG")
    return(substr(s, 1, 15))
  if(substr(s, 1, 7)=="ENSMUSG")
    return(substr(s, 1, 18))
  if(substr(s, 1, 4)=="ENST")
    return(substr(s, 1, 15))
  return(s)
})

read_rna_file <- function(file_rna, reduced = T, correct_gene_id = T)
{
  # name stomach is arbitrary here
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
  #indexes <- which(genebody_annotation$chr==start_chr & genebody_annotation$start>=start_position & genebody_annotation$start<=end_position)
}

get_genes_indexes_by_region <- function(start_chr, start_position, end_position, genebody_annotation)
{
  # start included
  which(genebody_annotation$chr==start_chr & genebody_annotation$start>=start_position & genebody_annotation$start<=end_position)
}

get_genes_nucleotides_intersection <- function(start_chr, start_position, end_position, genebody_annotation)
{
  # fully included
  
  # start included
  i = which(genebody_annotation$chr==start_chr & genebody_annotation$start>=start_position & genebody_annotation$start<=end_position)
  sum(pmin(genebody_annotation$end[i], end_position)-genebody_annotation$start[i])
  
}

get_TPM <- function(ids, rna_data)
{
  rna_data[gene_id %in% (ids), TPM]
}

get_expected_count <- function(ids, rna_data)
{
  rna_data[gene_id %in% (ids), expected_count]
}

make_rna_window_data_frame <- function(wgbs_data, rna_data, genebody_annotation, window)
{
  wgbs_data = sum_strands(wgbs_data, verbose = F); gc()
  
  l = length(wgbs_data$prop)
  fragments <- floor(l/window)
  i_starting_points <- ((0:(fragments-1))*window)+1
  
  window = i_starting_points[2]-i_starting_points[1]
  start_position = wgbs_data[i_starting_points]$Cpos
  start_chr = wgbs_data[i_starting_points]$chr
  end_position = wgbs_data[i_starting_points+window]$Cpos
  #end_chr = wgbs_data[i_starting_points+window]$chr
  l = length(start_position)
  
  cat("  not found genes:", sum(! (rna_data$gene_id %in% genebody_annotation$id)), "\n")
  
  rna_data = rna_data[gene_id %in% genebody_annotation$id]
  
  colnames(genebody_annotation)[5] = "gene_id"
  merged = merge(x = genebody_annotation, y = rna, by = "gene_id", all = F)
  cat("total gene set TPM:", sum(merged$TPM), "\n")
  merged$TPM[is.na(merged$TPM)] <- 0
  
  gc()
  gene_info = sapply(1:l, function(i)
  {
    if((i %% floor(l/100)) == 0)
      cat(floor(100*i/l), "%  ")
    
    genes_indexes = get_genes_indexes_by_region(start_chr[i], start_position[i], end_position[i], merged)
    genes = merged[genes_indexes, gene_id]
    genes_nucleotides_count <- get_genes_nucleotides_intersection(start_chr[i], start_position[i], end_position[i], genebody_annotation)
    tpm = merged[genes_indexes, TPM] #get_TPM(genes, rna_data)
    total_TPM = sum(tpm)
    std_TPM = sd(tpm)
    
    c(length(genes), total_TPM, genes_nucleotides_count, std_TPM)
  })
  
  nucleotides = end_position-start_position
  nucleotides[nucleotides<=0] = NA
  
  data.frame(start_chr, start_position, end_position,
             nucleotides, gene_count = gene_info[1,], genes_nucleotides_count = gene_info[3,], total_TPM = gene_info[2,], std_TPM = gene_info[4,])
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
    genes_nucleotides_count <- get_genes_nucleotides_intersection(start_chr[i], start_positions[i], end_positions[i], genebody_annotation)
    tpm = get_TPM(genes, rna_data)
    expected_count <- get_expected_count(genes, rna_data)
    total_TPM = sum(tpm)
    total_expected_count = sum(expected_count)
    c(length(genes), total_TPM, genes_nucleotides_count, total_expected_count)
  })
  
  nucleotides = end_positions-start_positions
  nucleotides[nucleotides<=0] = NA
  
  out = cbind(data.frame(start_chr, start_positions, end_positions,
                         nucleotides, gene_count = gene_info[1,], genes_nucleotides_count = gene_info[3,], total_TPM = gene_info[2,], total_expected_count = gene_info[3,]
  ), msr_experiment_data_frame)
  
  
  return(out)
}

join_rna_and_msr_tables <- function(rna_tables, msr_tables, i)
{
  df = cbind(rna_tables[[i]], msr_tables[[i]])
  l = length(colnames(rna_tables[[i]]))
  
  df$log_tpm = log(df$total_TPM+ 1e-3)
  df$log_std_tpm = log(df$std_TPM)
  
  colnames(df)[l+1] <- "i_start"
  colnames(df)[l+3] <- "meth rate"
  colnames(df)[l+6] <- "ecdf"
  colnames(df)[l+7] <- "inverted ecdf"
  
  df$CpG_density = (df$i_start[2]-df$i_start[1])/df$nucleotides
  
  df
}

join_rna_and_msr_table <- function(rna_table, msr_table)
{
  df = cbind(rna_table, msr_table)
  l = length(colnames(rna_table))
  
  df$log_tpm = log(df$total_TPM+ 1e-3)
  #df$log_std_tpm = log(df$std_TPM)
  
  colnames(df)[l+1] <- "i_start"
  colnames(df)[l+3] <- "meth rate"
  colnames(df)[l+6] <- "ecdf"
  colnames(df)[l+7] <- "inverted ecdf"
  
  df$CpG_density = (df$i_start[2]-df$i_start[1])/df$nucleotides
  
  df
}

produce_fragments_rna_tables <- function(wgbs_data, rna_data, genebody_annotation, sizes = c(1e3,1e4,1e5,1e6), file = "")
{
  tables = lapply(sizes, function(window)
  {
    make_rna_window_data_frame(wgbs_data, rna_data, genebody_annotation, window)
  })
  if(nchar(file>1))
    saveRDS(tables, file = file)
  return(tables)
}

exclude_outliers <- function(table, lim = 1e6, verbose = T)
{
  too_much_nucleotides <- (table$nucleotides > lim)
  if(verbose) cat(sprintf("%d rows had too many nucleotides", sum(too_much_nucleotides, na.rm=T)))
  return(table[!too_much_nucleotides, ])
}
#####################################################

binary_predictivity_hist <- function(feature, mask, feature_name, sep_names, colors = c(5,7), breaks = 10, main = "", max_y = 10, xlab = "")
{
  yes = feature[mask]
  no = feature[!mask]
  
  colors = c(alpha(colors[1],0.5), alpha(colors[2],0.5))
  
  hist(yes, col = colors[1], probability = T, breaks = breaks, xlab = xlab, main = main, ylim = c(0,max_y))
  hist(no, col = colors[2], probability = T, add = T, breaks = breaks)
  
  legend("top", legend=c(sep_names[1], sep_names[2]), col=colors, fill = colors, cex = 0.75)
}


difference_meth_plot <- function(d1, d2, name1, name2, min_reads, cgi_anno, island_mask, meth_island_mask)
{
  
  reads_mask = (d1$reads >= min_reads) & (d2$reads >= min_reads)
  
  d1p_islands = d1[reads_mask & island_mask, ]$prop
  d2p_islands = d2[reads_mask & island_mask, ]$prop
  diff_islands = d1p_islands-d2p_islands
  
  d1p_meth_islands = d1[reads_mask & meth_island_mask, ]$prop
  d2p_meth_islands = d2[reads_mask & meth_island_mask, ]$prop
  diff_meth_islands = d1p_meth_islands-d2p_meth_islands
  
  d1p_sea = d1[reads_mask & !island_mask, ]$prop
  d2p_sea = d2[reads_mask & !island_mask, ]$prop
  diff_sea = d1p_sea-d2p_sea
  
  par(ask=TRUE)
  hist(diff_sea, breaks = 50, main = sprintf("Methylation difference in CG sites outside CpG island, mean: %s", round(mean(diff_sea, na.rm=T),2)), xlab = sprintf("%s - %s", name1, name2), col = 5)
  hist(diff_islands, breaks = 50, main = sprintf("Methylation difference in CG sites inside CpG island, mean: %s", round(mean(diff_islands, na.rm=T),2)), xlab = sprintf("%s - %s", name1, name2), col = "yellow")
  
  par(mfrow = c(1,1))
  hist(diff_meth_islands, breaks = 50, main = sprintf("Methylation difference in CG sites belonging to methylated CpG island, mean: %s", round(mean(diff_meth_islands, na.rm=T),2)), xlab = sprintf("%s - %s", name1, name2), col = "red")
  
  #par(ask=TRUE)
  #par(mfrow = c(2,1))
  
  #plot(diff_sea[0:1e3], type = "l", main = "outside CpG island", xlab = sprintf("%s - %s", name1, name2), col = 1)
  #lines(x = c(-100, 1e7), y = c(0,0), lty = 2, col = alpha("red", 0.7))
  #lines(x = c(-100, 1e7), y = c(50,50), lty = 2, col = alpha("orange", 0.7))
  #lines(x = c(-100, 1e7), y = c(-50,-50), lty = 2, col = alpha("orange", 0.7))
  
  #plot(diff_islands[0:1e3], type = "l", main = "inside CpG island", xlab = sprintf("%s - %s", name1, name2), col = 1)
  #lines(x = c(-100, 1e7), y = c(0,0), lty = 2, col = alpha("red", 0.7))
  #lines(x = c(-100, 1e7), y = c(50,50), lty = 2, col = alpha("orange", 0.7))
  #lines(x = c(-100, 1e7), y = c(-50,-50), lty = 2, col = alpha("orange", 0.7))
  
  par(ask=FALSE)
  
  
  
}


noNa <- function(d, verbose = T)
{
  nn <- d[complete.cases(d),]
  if(verbose)
    cat("NAs: ", length(d[,1])/length(nn[,1]),"%\n")
  return(nn)
}


suppressMessages(library(infotheo))

fast.mutual.information <- function(a, b, method = "mm", normalized = T)
{
  d <- data.frame(a,b)
  d <- discretize(d)
  
  mi <-  mutinformation(d[,1],d[,2], method)
  if(!normalized)
    return(mi)
  else
    return(mi/sqrt(entropy(d[,1],method)*entropy(d[,2],method)))
  
}


suppressMessages(library(corrplot))

mi.plot <- function(data, method = "mm")
{
  names = colnames(data)
  size = length(names)
  M = array(dim = c(size,size),0)
  
  for(i in 1:size)
  {
    for(j in 1:i)
    {
      M[i,j] = fast.mutual.information(data[,i], data[,j], method)
      M[j,i] = M[i,j]
    }
  }
  
  colnames(M) <- names
  rownames(M) <- names
  corrplot(M, "number")
  return(M)
}

pairwise.plot <- function(data, binary_function)
{
  names = colnames(data)
  size = length(names)
  M = array(dim = c(size,size),0)
  
  for(i in 1:size)
  {
    for(j in 1:i)
    {
      M[i,j] = binary_function(data[,i], data[,j])
      M[j,i] = M[i,j]
    }
  }
  
  colnames(M) <- names
  rownames(M) <- names
  corrplot(M, "number")
  return(M)
}

discretized_density <- function(wgbs_data, indexes)
{
  discretized_prop <- wgbs_data$prop>=50
  cs <- corrected_cumsum(discretized_prop)
  b <- indexes[2:length(indexes)]
  a <- indexes[1:(length(indexes)-1)]
  (cs[b]-cs[a])/(b-a)
  
}

discretized_density <- function(wgbs_data, indexes)
{
  discretized_prop <- wgbs_data$prop>=50
  step <- indexes[2]-indexes[1]
  
  sapply(1:length(indexes), function(i)
  {
    mean(discretized_prop[(indexes[i]):(indexes[i]+step)], na.rm = T)
  })
}


significance_measure <- function(msr, density, msr_ecdf, inverted)
{
  if(inverted)
    density = 1-density
  
  l = length(density)
  zero = function(x) {0}
  
  if(msr_ecdf[[1]]$prop==0) msr_ecdf[[1]]$cdf = zero
  if(msr_ecdf[[length(msr_ecdf)]]$prop==1) msr_ecdf[[length(msr_ecdf)]]$cdf = zero
  
  sapply( 1:l ,function(n)
  {
    #cat(n, " ")
    if(is.na(density[n]) || is.na(msr[n])) return(NA)
    general_msr_cdf(msr_ecdf, density[n], msr[n])
  })
  
}

get_residuals <- function(msr, density, inverted)
{
  if(inverted)
    density = 1-density
  res = lm(msr~density)$residuals
  mask = !is.na(density) & !is.na(msr)
  out = array(dim = length(mask))
  out[mask] = res
  return(out)
}

show_region <- function(i, data_table, wgbs_data, size = 1000, prop = F)
{
  data_table_row <- data_table[i,] 
  start = data_table_row$i_start
  end = start+size
  if(prop)
  {
    plot((wgbs_data$prop[start:end]/100))
  }
  else
  {
    d <- data_table_row[,c("msr", "inverted_msr", "ecdf", "msr_density")]
    m <- paste("MSR(1):", round(d[1], 3) ," MSR(0):", round(d[2], 3), " ecdf:", round(d[3], 4), " rate:", round(d[4], 3))
    plot(round(wgbs_data$prop[start:end]/100), pch = "|", ylab = " ", main = m)
  }
  
  print(data_table_row[,c("msr", "inverted_msr", "residual", "total_TPM", "inverted ecdf", "msr_density")])
}


get_file_names <- function(dir, patterns, full.names = F)
{
  a <- list.files(path = dir, full.names =  full.names)
  for(p in patterns)
  {
    a <- a[grepl(p, a, ignore.case = T)]
  }
  
  return(a)
}

# mean distance from 1 or 0
drift <- function(x) {mean(0.5-abs(x-0.5), na.rm=T)}

tmse <- function(model, test_model_data, y_name)
{
  1-(var(predict(model, test_model_data,  type="response")-test_model_data[,y_name]))/var(test_model_data[,y_name])
}

tmse.glmnet <- function(model, test_model_data, y_name, x_names)
{
  newx = as.matrix.data.frame(test_model_data[, x_names])
  1-(var(predict.glmnet(model, newx = newx,  type="response")-test_model_data[,y_name]))/var(test_model_data[,y_name])
}

get_genomewide_meth_positions <- function(wgbs, start, end, chr)
{
  wgbs$Cpos[wgbs$chr==chr & wgbs$reads>0 & wgbs$prop>=50 & wgbs$Cpos>=start & wgbs$Cpos<=end]
}

autocorrelation <- function(v, method = "pearson")
{
  x <- 2:length(v)
  # nas <- sum(is.na(v))
  # if(nas>((1-na_tolerance)*length(v)) | length(v)-nas<5)
  #   return(NA)
  # cor.test(v[x-1], v[x], method = method)$estimate
  
  tryCatch(
    expr = {
      return(cor.test(v[x-1], v[x], method = method)$estimate)
    },
    error = function(e){
      return(NA)
    },
    warning = function(w){
      return(NA)
    },
    finally = {
      
    }
  )    
}  

meth_autoc <- function(wgbs, data_table, na_tolerance = 0.2, minimum_reads = 1, method = "pearson")
{
  l <- data_table$i_start[2]-data_table$i_start[1]
  
  prop <- wgbs$prop
  prop[wgbs$reads<minimum_reads] <- NA
  autoc <- sapply(data_table$i_start, function(i)
  {
    autocorrelation(prop[i:(i+l)], na_tolerance, method)
  })
  
  return(autoc)
}


# doesn't work
auto.mi <- function(wgbs, data_table, na_tolerance = 0.2, minimum_reads = 1, binaryze = T)
{
  l <- data_table$i_start[2]-data_table$i_start[1]
  
  prop <- wgbs$prop
  prop[wgbs$reads<minimum_reads] <- NA
  
  if(binaryze)
    prop <- round((prop/100)+1e-4)
  
  sapply(data_table$i_start[1:(length(data_table$i_start))], function(i)
  {
    x <- i:(i+l-1)
    #cat(i, i+l-1, "\n")
    if(sum(is.na(prop[x]))>((1-na_tolerance)*l))
      return(NA)
    fast.mutual.information((prop[x]), (prop[x+1]), method = "emp")
  })
  
}

show_region_meth <- function(wgbs, chromosome, start, end, prop = T, genomewide = F)
{
  w <- wgbs[wgbs$chr == chromosome & wgbs$Cpos>=start & wgbs$Cpos<=end, ]
  p <- w$prop/100
  m <- length(p)
  
  if(m<1 | all(is.na(p)))
    return(F)
  
  m <- paste("CG: ", m, "  nucleotides: ", end-start, "\n", chromosome, start,"-", end,
             "\nautocorrelation: ", round(autocorrelation(p), 2), "  drift", round(drift(p), 2))
  #cat(p, m)
  
  x <- 1:length(p)
  if(genomewide)
    x <- w$Cpos
  
  if(prop)
  {
    plot(x, p, main=m)
  }
  else
  {
    plot(x, round(p), pch = "|", ylab = " ", main=m)
  }
}

show_gene_meth <- function(wgbs, genebody_annotation, i, prop = T, genomewide = F, margins = 0)
{
  g <- genebody_annotation[i,]
  print(g)
  show_region_meth(wgbs, g$chr, g$start - margins, g$end + margins, prop, genomewide)
}

# chr start end strand gene_id gene_type i_start i_end CG_count CG_density meth_rate meth_rate_binary CGsites_msr CGmeth_msr CGunmeth_msr msr inverted_msr

make_genes_table <- function(wgbs, genebody_annotation, gene_type_filter = NA, na_tolerance = 0.3, no_msr = F)
{
  data_table <- genebody_annotation
  
  if(!is.na(gene_type_filter))
    data_table <- data_table[data_table$gene_type==gene_type_filter, ]
  
  l <- length(data_table$start)
  data_table$CG_count <- data_table$i_end - data_table$i_start + 1
  data_table$nucleotides <- data_table$end - data_table$start
  data_table$CG_density <- data_table$CG_count/data_table$nucleotides
  
  data_table$missing_prop <- NA
  data_table$meth_rate <- NA
  data_table$meth_rate_binary <- NA
  data_table$CGsites_msr <- NA
  data_table$meth_msr <- NA
  data_table$unmeth_msr <- NA
  data_table$CG_list_msr <- NA
  data_table$CG_list_inverted_msr <- NA
  data_table$meth_autocorrelation <- NA
  data_table$drift <- NA

  for(i in 1:l)
  {
    if(!is.na(data_table$i_start[i]))
    {
      cat(i, "\n")
      w <- wgbs[data_table$i_start[i]:data_table$i_end[i], ]
      CG_pos <- w$Cpos
      meth_CG_pos <- w[w$prop>=50,]$Cpos
      unmeth_CG_pos <- w[w$prop<50,]$Cpos
      CG_lits <- (w$prop>=50)
      
      data_table$missing_prop[i] <- mean(w$reads==0)
      data_table$meth_rate[i] <- mean(w$prop, na.rm = T)/100
      data_table$meth_rate_binary[i] <- mean(CG_lits, na.rm = T)
      if(data_table$CG_count[i]>=20)
        data_table$meth_autocorrelation[i] <- autocorrelation(w$prop)
      data_table$drift[i] <- drift(w$prop/100)
      
      # cat(w$prop)
      # cat("\nCG_pos", CG_pos)
      # cat("\nmeth_CG_pos", meth_CG_pos)
      # cat("\nCG_lits", CG_lits)
      if(data_table$CG_count[i]>=100 & data_table$missing_prop[i]<=na_tolerance & !no_msr)
      {
        data_table$CGsites_msr[i] <- genome_MSR(CG_pos, minimum_bin_size = 1, max_bins = 40)
        data_table$meth_msr[i] <- genome_MSR(meth_CG_pos, minimum_bin_size = 1, max_bins = 40)
        data_table$unmeth_msr[i] <- genome_MSR(unmeth_CG_pos, minimum_bin_size = 1, max_bins = 40)
        data_table$CG_list_msr[i] <- MSR_area(calculate_relevance_resolution_vector(CG_lits, na_tolerance = na_tolerance, na_values_handler = replace_nas_with_bin_prop, invert = F, verbose = F))
        data_table$CG_list_inverted_msr[i] <- MSR_area(calculate_relevance_resolution_vector(CG_lits, na_tolerance = na_tolerance, na_values_handler = replace_nas_with_bin_prop, invert = T, verbose = F))
      }
    }
  }
  
  return(data_table)
}

make_genes_table2 <- function(wgbs, genebody_annotation, gene_type_filter = NA, na_tolerance = 0.3, no_msr = F)
{
  data_table <- genebody_annotation
  
  if(!is.na(gene_type_filter))
    data_table <- data_table[data_table$gene_type==gene_type_filter, ]
  
  l <- length(data_table$start)
  data_table$CG_count <- data_table$i_end - data_table$i_start + 1
  data_table$nucleotides <- data_table$end - data_table$start
  data_table$CG_density <- data_table$CG_count/data_table$nucleotides
  
  data_table$missing_prop <- NA
  data_table$meth_rate <- NA
  data_table$meth_rate_binary <- NA
  data_table$CGsites_msr <- NA
  data_table$meth_msr <- NA
  data_table$unmeth_msr <- NA
  data_table$CG_list_msr <- NA
  data_table$CG_list_inverted_msr <- NA
  data_table$meth_autocorrelation <- NA
  data_table$drift <- NA
  
  mask <- !is.na(data_table$i_start + data_table$i_end)
  valid_rows <- (1:length(data_table$i_start))[mask]
  cat("\n precomputing")
  wList <- lapply(valid_rows, function(i){wgbs[data_table$i_start[i]:data_table$i_end[i], -"chr"]})
  
  cat("\ncomputing basic features")
  data_table$missing_prop[valid_rows] <- sapply(wList, function(w){mean(w$reads==0)}) 
  data_table$meth_rate[valid_rows] <- sapply(wList, function(w){mean(w$prop, na.rm = T)/100})
  data_table$meth_rate_binary[valid_rows] <- sapply(wList, function(w){mean(w$prop>=50, na.rm = T)})
  cat("\ncomputing autocorrelation")
  data_table$meth_autocorrelation[valid_rows] <- sapply(wList, function(w){autocorrelation(w$prop)})
  cat("\ncomputing drift")
  data_table$drift[valid_rows] <- sapply(wList, function(w){drift(w$prop/100)})
  
  if(!no_msr)
  {
    max_bins <- 40
    mask <- data_table$CG_count>=100 & data_table$missing_prop<=na_tolerance & !is.na(data_table$i_start + data_table$i_end)
    valid_rows <- (1:length(data_table$i_start))[mask]
    cat("\n precomputing")
    wList <- lapply(valid_rows, function(i){wgbs[data_table$i_start[i]:data_table$i_end[i], -"chr"]})
    r <- 1:length(wList)
    
    cat("\ncomputing CGsites_msr")
    data_table$CGsites_msr[valid_rows] <- sapply(r, function(i){cat(i," "); w<-wList[[i]]; genome_MSR(w$Cpos, minimum_bin_size = 1, max_bins = max_bins)})
    cat("\ncomputing meth_msr")
    data_table$meth_msr[valid_rows] <- sapply(r, function(i){cat(i," "); w<-wList[[i]]; genome_MSR(w[w$prop>=50]$Cpos, minimum_bin_size = 1, max_bins = max_bins)})
    cat("\ncomputing unmeth_msr")
    data_table$unmeth_msr[valid_rows] <- sapply(r, function(i){cat(i," "); w<-wList[[i]]; genome_MSR(w[w$prop<50]$Cpos, minimum_bin_size = 1, max_bins = max_bins)})
    cat("\ncomputing CG_list_msr")
    data_table$CG_list_msr[valid_rows] <- sapply(r, function(i){cat(i," "); w<-wList[[i]]; MSR_area(calculate_relevance_resolution_vector((w$prop>=50), na_tolerance = na_tolerance, na_values_handler = replace_nas_with_bin_prop, invert = F, verbose = F))})
    cat("\ncomputing CG_list_inverted_msr")
    data_table$CG_list_inverted_msr[valid_rows] <- sapply(r, function(i){cat(i," "); w<-wList[[i]]; MSR_area(calculate_relevance_resolution_vector((w$prop>=50), na_tolerance = na_tolerance, na_values_handler = replace_nas_with_bin_prop, invert = T, verbose = F))})
    
  }
  
  
  return(data_table)
}

model_this <- function(response_variable_name)
{
  as.formula(paste(response_variable_name, "~ ."))
}


suppressMessages(library(glmnet))

lasso <- function(response_variable_name, df, lambda = 0.01, alpha = 1)
{
  x_vars <- (as.matrix.data.frame(df[, names(df) != response_variable_name]))
  y_var <- (df[, response_variable_name])
  lambda_seq <- 10^seq(2, log(lambda,10), by = -.1)
  cv_output <- cv.glmnet(x_vars, y_var, alpha = alpha, lambda = lambda_seq, nfolds = 5, standardize = T)
  glmnet(x_vars, y_var, alpha = alpha, lambda = cv_output$lambda.min)
}

produce_and_save_genes_msr_table <- function(wgbs_file, short_name, genebody_annotation_file, gene_type_filter = NA, na_tolerance = 0.3, no_msr = F, dir = "../../Rexperiments/")
{
  if(is.character(wgbs_file)) { wgbs <- sum_strands(readRDS(wgbs_file)); gc() }
  else{wgbs <- sum_strands(wgbs_file)}
  genebody_annotation <- readRDS(genebody_annotation_file)
  gbt <- make_genes_table(wgbs, genebody_annotation, gene_type_filter = gene_type_filter, na_tolerance = na_tolerance, no_msr = no_msr)
  saveRDS(gbt, paste(dir,short_name, "_genes_msr_table.Rda", sep = ""))
}

produce_and_save_fragments_msr_table <- function(wgbs_file, short_name, size, msr_ecdf_file, na_tolerance = 0.4, minimum_reads=1, methylation_assigner = standard_binaryzer, bed = NA, dir = "../../Rexperiments/")
{
  msr_ecdf_ref <- readRDS(msr_ecdf_file)
  
  # convert bed file into RDA
  if(!is.na(bed))
    rda_convert(bed, wgbs_file)
  
  # read rda file
  if(is.character(wgbs_file))
  {
    wgbs <- sum_strands(readRDS(new_name)); gc()
  }
  else{wgbs <- sum_strands(wgbs_file)}
  
  # produce msr table
  rr_table <- total_spatial_experiment(c(new_name), c(size), c(F,T), c(short_name), methylation_assigner, na_tolerance, F, minimum_reads)
  start = rr_table[[1]][[1]]$data$fragments_infos_array[,1]
  msr_density = rr_table[[1]][[1]]$data$fragments_infos_array[,2]
  true_density = sapply(start, function(x){mean(wgbs$prop[x:(x+size)], na.rm = T)})/100
  msr = rr_table[[1]][[1]]$data$fragments_infos_array[,3]
  inverted_msr = rr_table[[2]][[1]]$data$fragments_infos_array[,3]
  inverted_msr[is.na(inverted_msr)] <- msr[is.na(inverted_msr)]
  msr[is.na(msr)] <- inverted_msr[is.na(msr)]
  sig <- significance_measure(msr, msr_density, msr_ecdf_ref, inverted = F)
  inverted_sig <- significance_measure(inverted_msr, msr_density, msr_ecdf_ref, inverted = T)
  median_function <- extract_ecdf_function(msr_ecdf_ref, 0.5)
  residual <- msr-median_function(msr_density)
  inverted_residual <- inverted_msr-median_function(1-msr_density)
  meth_autocorrelation <- sapply(start, function(x){autocorrelation(wgbs$prop[x:(x+size)])})
  mean_drift <- sapply(start, function(x){drift(wgbs$prop[x:(x+size)])})
  
  # save msr tables
  msr_table = data.frame(start,msr_density,true_density,msr,inverted_msr, sig, inverted_sig, residual, inverted_residual, meth_autocorrelation, drift = mean_drift)
  saveRDS(msr_table, file = paste("../../Rexperiments/",short_name, "_msr_table_", size, ".Rda", sep = ""))
}

produce_and_save_fragments_expression_table <- function(wgbs_file, short_name, rna_file, size, genebody_annotation_file = "../../Rexperiments/detailed_genebody_improved.Rda", correct_gene_id = T)
{
  genebody_annotation <- readRDS(genebody_annotation_file); gc()
  rna <- read_rna_file(rna_file, reduced = F, correct_gene_id = correct_gene_id)
  
  # read rda file
  if(is.character(wgbs_file)) { wgbs <- sum_strands(readRDS(wgbs_file)); gc() }
  else{wgbs <- sum_strands(wgbs_file)}
  
  if(correct_gene_id)
    genebody_annotation$id <- remove_version_from_gene(genebody_annotation$id)
  
  rna_table <- make_rna_window_data_frame(wgbs, rna, genebody_annotation, size)
  saveRDS(rna_table, file = paste("../../Rexperiments/",short_name, "_rna_table_", size, extension_tag, ".Rda", sep = ""))
}