source("MSR_analysis_functions.R", chdir = T)
suppressMessages(library("tools"))
suppressMessages(library(infotheo))
suppressMessages(library(corrplot))
suppressMessages(library(glmnet))

# WGBS MANIPULATION UTILS

filter_strand     <- function(data, sign) { return(data[strand==sign, -"strand"]) }
pick_plus_strand  <- function(data)       { return(filter_strand(data,"+"))}
pick_minus_strand <- function(data)       { return(filter_strand(data,"-"))}
chromosomes <- function(Genome = BSgenome.Hsapiens.UCSC.hg38) { return(BSgenome.Hsapiens.UCSC.hg38@seqinfo@seqnames[1:24])}
meth_proportion <- function(data, minimum_reads = 1) { mean(data[reads>=minimum_reads]$prop, rm.na = T) }
filter_chromosome <- function(data, chromosome) { return(data[chr == chromosome, -"chr"])}

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

standard_binaryzer2 <- function(data){
  v <- data$prop
  v[v==50] <- NA
  v>50
  }

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

clean_bed_file <- function(data, strands_handler, methylation_assigner, missing_read_handler)
{
  return(missing_read_handler(methylation_assigner(strands_handler(data))))
}

rda_convert <- function(file, new_name)
{
  data <- read_ENCODE_bed(file, verbose = T)
  saveRDS(data, file = new_name)
  system(paste("rm", file))
  gc()
}

#####################################################


# DATASET GENERATION UTILS

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
  # start included
  genebody_annotation[chr==start_chr & start>=start_position & start<=end_position]$id
}

get_genes_indexes_by_region <- function(start_chr, start_position, end_position, genebody_annotation)
{
  # start included
    which(genebody_annotation$chr==start_chr & genebody_annotation$start>=start_position & genebody_annotation$start<=end_position)
}

get_genes_nucleotides_intersection <- function(start_chr, start_position, end_position, genebody_annotation)
{
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
  
  #window = i_starting_points[2]-i_starting_points[1]
  start_position = wgbs_data[i_starting_points]$Cpos
  start_chr = wgbs_data[i_starting_points]$chr
  end_position = wgbs_data[i_starting_points+window]$Cpos
  l = length(start_position)
  
  cat("  not found genes:", sum(! (rna_data$gene_id %in% genebody_annotation$gene_id)), "\n")
  
  #rna_data = rna_data[gene_id %in% genebody_annotation$gene_id]
  
  merged = merge(x = genebody_annotation, y = rna_data, by = "gene_id", all = F)
  cat("total gene set TPM:", sum(merged$TPM), "\n")
  if(any(is.na(merged$TPM))) cat("some nas")
  #merged$TPM[is.na(merged$TPM)] <- 0
  merged$start[merged$strand == "-"] <- merged$end[merged$strand == "-"]
  gc()
  gene_info = sapply(1:l, function(i)
  {
    show_perc_progress(i,l)
    
    genes_indexes = get_genes_indexes_by_region(start_chr[i], start_position[i], end_position[i], merged)
    genes = merged[genes_indexes, gene_id]
    genes_nucleotides_count <- get_genes_nucleotides_intersection(start_chr[i], start_position[i], end_position[i], genebody_annotation)
    total_TPM = sum(merged[genes_indexes, TPM])
    total_pme_TPM <- sum(merged[genes_indexes, pme_TPM])
    total_FPKM <- sum(merged[genes_indexes, FPKM])
    total_pme_FPKM <- sum(merged[genes_indexes, pme_FPKM])
    total_expected_count <- sum(merged[genes_indexes, expected_count])
    
    c(length(genes), total_TPM, genes_nucleotides_count, total_pme_TPM, total_FPKM, total_pme_FPKM, total_expected_count)
  })
  
  nucleotides = end_position-start_position
  nucleotides[nucleotides<=0] = NA
  
  data.frame(start_chr, start_position, end_position,
             nucleotides, gene_count = gene_info[1,], genes_nucleotides_count = gene_info[3,], total_TPM = gene_info[2,], total_pme_TPM = gene_info[4,],
             total_FPKM = gene_info[5,], total_pme_FPKM = gene_info[6,], total_expected_count = gene_info[7,])
}

join_rna_and_msr_table <- function(rna_table, msr_table)
{
  rna_table <- rna_table[(rna_table$end_position-rna_table$start_position)>0,]
  msr_table <- msr_table[(msr_table$end-msr_table$start)>0,]
  df = cbind(msr_table,rna_table)
  l = length(colnames(rna_table))
  
  df$log_tpm = log(df$total_TPM+ 1e-3)
  #df$log_std_tpm = log(df$std_TPM)
  
  df$CpG_density = (df$i_start[2]-df$i_start[1])/df$nucleotides
  
  df
}

exclude_outliers <- function(table, lim = 1e6, verbose = T)
{
  too_many_nucleotides <- (table$nucleotides > lim)
  if(verbose) cat(sprintf("%d rows had too many nucleotides", sum(too_many_nucleotides, na.rm=T)))
  return(table[!too_many_nucleotides, ])
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

get_file_names <- function(dir, patterns, full.names = F)
{
  a <- list.files(path = dir, full.names =  full.names)
  for(p in patterns)
  {
    a <- a[grepl(p, a, ignore.case = T)]
  }
  
  return(a)
}

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
  data_table$meth_sd <- NA

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
      data_table$meth_sd[i] <- sd(w$prop/100, na.rm=T)
      
      # cat(w$prop)
      # cat("\nCG_pos", CG_pos)
      # cat("\nmeth_CG_pos", meth_CG_pos)
      # cat("\nCG_lits", CG_lits)
      if(data_table$CG_count[i]>=100 & data_table$missing_prop[i]<=na_tolerance & !no_msr)
      {
        data_table$CGsites_msr[i] <- genome_MSR(CG_pos, minimum_bin_size = 2, max_gbins = 40)
        data_table$meth_msr[i] <- genome_MSR(meth_CG_pos, minimum_bin_size = 2, max_gbins = 40)
        data_table$unmeth_msr[i] <- genome_MSR(unmeth_CG_pos, minimum_bin_size = 2, max_gbins = 40)
        data_table$CG_list_msr[i] <- MSR_area(calculate_relevance_resolution_vector(CG_lits, na_tolerance = na_tolerance, na_values_handler = replace_nas_with_bin_prop, invert = F, verbose = F))
        data_table$CG_list_inverted_msr[i] <- MSR_area(calculate_relevance_resolution_vector(CG_lits, na_tolerance = na_tolerance, na_values_handler = replace_nas_with_bin_prop, invert = T, verbose = F))
      }
    }
  }
  
  return(data_table)
}

produce_and_save_genes_msr_table <- function(wgbs_file, short_name, genebody_annotation_file, gene_type_filter = NA, na_tolerance = 0.3, no_msr = F, dir = "../../Rexperiments/")
{
  if(is.character(wgbs_file)) { wgbs <- sum_strands(readRDS_if_needed(wgbs_file)); gc() }
  else{wgbs <- sum_strands(wgbs_file)}
  genebody_annotation <- readRDS(genebody_annotation_file)
  gbt <- make_genes_table(wgbs, genebody_annotation, gene_type_filter = gene_type_filter, na_tolerance = na_tolerance, no_msr = no_msr)
  saveRDS(gbt, paste(dir,short_name, "_genes_msr_table.Rda", sep = ""))
}

produce_and_save_fragments_msr_table <- function(wgbs_file, short_name, size, msr_ecdf_file, na_tolerance = 0.4, minimum_reads=1, dir = "../../Rexperiments/")
{
  msr_ecdf_ref <- readRDS(msr_ecdf_file)
  
  wgbs <- sum_strands(readRDS_if_needed(wgbs_file)); gc(verbose = F)

  l = length(wgbs$prop)
  fragments <- floor(l/size)
  i_start <- ((0:(fragments-1))*size)+1
  
  binary <- wgbs$prop>=50
  binary[wgbs$reads<minimum_reads] <- NA
  binary[wgbs$prop==50] <- NA
  
  cat("\nmsr_density...")
  msr_density = sapply(i_start, function(x){mean(binary[x:(x+size)], na.rm = T)})
  cat("\nmeth_rate...")
  meth_rate = sapply(i_start, function(x){mean(wgbs$prop[x:(x+size)], na.rm = T)})/100
  cat("\nmsr...")
  msr = sapply(i_start, function(x){ cat(((x-1)/size)+1, "\n"); MSR_area(calculate_relevance_resolution_vector(binary[x:(x+size)], verbose=F, na_tolerance = na_tolerance, na_values_handler = replace_nas_with_bin_prop, invert = F))})
  cat("\ninverted_msr...")
  inverted_msr = sapply(i_start, function(x){ cat(((x-1)/size)+1, "\n"); MSR_area(calculate_relevance_resolution_vector(binary[x:(x+size)], verbose=F, na_tolerance = na_tolerance, na_values_handler = replace_nas_with_bin_prop, invert = T))})

  inverted_msr[is.na(inverted_msr)] <- msr[is.na(inverted_msr)]
  msr[is.na(msr)] <- inverted_msr[is.na(msr)]
  
  cat("\necdf...")
  ecdf <- significance_measure(msr, msr_density, msr_ecdf_ref, inverted = F)
  cat("\ninverted_ecdf...")
  inverted_ecdf <- significance_measure(inverted_msr, msr_density, msr_ecdf_ref, inverted = T)
  median_function <- extract_ecdf_function(msr_ecdf_ref, 0.5)
  cat("\nresidual...")
  residual <- msr-median_function(msr_density)
  cat("\ninverted_residual...")
  inverted_residual <- inverted_msr-median_function(1-msr_density)
  cat("\nmeth_autocorrelation...")
  meth_autocorrelation <- sapply(i_start, function(x){autocorrelation(wgbs$prop[x:(x+size)])})
  cat("\nmean_drift...")
  mean_drift <- sapply(i_start, function(x){drift(wgbs$prop[x:(x+size)]/100)})
  cat("\nmeth_sd...")
  meth_sd <- sapply(i_start, function(x){sd(wgbs$prop[x:(x+size)]/100, na.rm=T)})
  cat("\nmissing_prop...")
  missing_prop <- sapply(i_start, function(x){mean(is.na(wgbs$prop[x:(x+size)]))})
  
  # save msr tables
  msr_table = data.frame(i_start,msr_density,meth_rate,msr,inverted_msr, ecdf, inverted_ecdf, residual, inverted_residual, meth_autocorrelation, drift = mean_drift, meth_sd, missing_prop)
  saveRDS(msr_table, file = paste(dir,short_name, "_msr_table_", size, ".Rda", sep = ""))
}

expression_file_type <- function(file_name)
{
  polyA <- grepl("polyA", file_name, fixed = T) | grepl("poly-A", file_name, fixed = T) | grepl("poly_A", file_name, fixed = T)
  total <- grepl("total", file_name, ignore.case = T)
  type <- ""
  if(polyA) type = "polyA"
  if(total) type = "total"
  if(total & polyA) type = ""
  return(type)
}

produce_and_save_fragments_expression_table <- function(wgbs_file, expression_file, genebody_annotation_file, size, ignore_gene_version, short_name, dir = "../../Rexperiments/", tag = "", filter_gene_type = NA)
{
  wgbs <- sum_strands(readRDS_if_needed(wgbs_file))
  genebody_annotation <- readRDS_if_needed(genebody_annotation_file)
  expression_table <- read_rna_file(expression_file, reduced = F, correct_gene_id = F)
  
  if(!is.na(filter_gene_type)) genebody_annotation <- genebody_annotation[genebody_annotation$gene_type==filter_gene_type, ]
    
  if(ignore_gene_version)
  {
    expression_table$gene_id <- remove_version_from_gene(expression_table$gene_id)
    genebody_annotation$gene_id <- remove_version_from_gene(genebody_annotation$gene_id)
  }
  
  
  cat("genes in anno but not in expr file: ", sum(!(genebody_annotation$gene_id %in% expression_table$gene_id)), "\n")
  
  type_tag <- expression_file_type(expression_file)
  fragments_expression_table <- make_rna_window_data_frame(wgbs, expression_table, genebody_annotation, size)
  saveRDS(fragments_expression_table, file = paste(dir,short_name, "_expression_table_", size, "_",type_tag,tag, ".Rda", sep = ""))
}

exclude_sex_chromosomes <- function(data_table)
{
  data_table[data_table$chr!="chrX" & data_table$chr!="chrY", ]
}

msr_exp_complete <- function(wgbs, short_name, CG_window_size, msr_ecdf_file, na_tolerance = 0.4, minimum_reads=1, dir = "../../Rexperiments/")
{
  max_gbins <- 100
  msr_ecdf_ref <- readRDS_if_needed(msr_ecdf_file)
  
  wgbs <- sum_strands(readRDS_if_needed(wgbs)); gc(verbose = F)
  
  fragments <- floor(length(wgbs$prop)/CG_window_size)
  i_start <- ((0:(fragments-1))*CG_window_size)+1

  binary <- wgbs$prop>=50
  binary[wgbs$reads<minimum_reads] <- NA
  binary[wgbs$prop==50] <- NA
  
  i_start <- i_start[wgbs$chr[i_start]==wgbs$chr[i_start+CG_window_size]]
  msr_table = data.frame(i_start = i_start)
  l <- length(i_start)
  i_end <- i_start+CG_window_size-1
  
  
  msr_table$i_end <- i_end
  msr_table$start <- wgbs$Cpos[i_start]
  msr_table$end <- wgbs$Cpos[i_start+CG_window_size-1]
  msr_table$chr <- wgbs$chr[i_start]
  msr_table$nucleotides <- msr_table$end - msr_table$start

  cat("\nmsr_density...")
  msr_table$msr_density = sapply(i_start, function(x){mean(binary[x:(x+CG_window_size-1)], na.rm = T)})
  cat("\nmeth_rate...")
  msr_table$meth_rate = sapply(i_start, function(x){mean(wgbs$prop[x:(x+CG_window_size-1)], na.rm = T)})/100
  cat("\nmsr...")
  msr_table$msr = sapply(i_start, function(x){ cat(((x-1)/CG_window_size)+1, "\n"); MSR_area(calculate_relevance_resolution_vector(binary[x:(x+CG_window_size-1)], verbose=F, na_tolerance = na_tolerance, na_values_handler = replace_nas_with_bin_prop, invert = F))})
  cat("\ninverted_msr...")
  msr_table$inverted_msr = sapply(i_start, function(x){ cat(((x-1)/CG_window_size)+1, "\n"); MSR_area(calculate_relevance_resolution_vector(binary[x:(x+CG_window_size-1)], verbose=F, na_tolerance = na_tolerance, na_values_handler = replace_nas_with_bin_prop, invert = T))})
  
  cat("\ncomputing CGsites_msr")
  msr_table$CGsites_msr <- sapply(i_start, function(x){ cat(((x-1)/CG_window_size)+1, "\n"); w<-wgbs$Cpos[x:(x+CG_window_size-1)]; genome_MSR(w, minimum_bin_size = 2, max_gbins = max_gbins)})
  cat("\ncomputing meth_msr")
  msr_table$meth_msr <- sapply(i_start, function(x){ cat(((x-1)/CG_window_size)+1, "\n"); genome_MSR(wgbs$Cpos[x:(x+CG_window_size-1)][binary[x:(x+CG_window_size)]], minimum_bin_size = 2, max_gbins = max_gbins)})
  cat("\ncomputing unmeth_msr")
  msr_table$unmeth_msr <- sapply(i_start, function(x){ cat(((x-1)/CG_window_size)+1, "\n"); genome_MSR(wgbs$Cpos[x:(x+CG_window_size-1)][!binary[x:(x+CG_window_size)]], minimum_bin_size = 2, max_gbins = max_gbins)})

  msr_table$inverted_msr[is.na(msr_table$inverted_msr)] <- msr_table$msr[is.na(msr_table$inverted_msr)]
  msr_table$msr[is.na(msr_table$msr)] <- msr_table$inverted_msr[is.na(msr_table$msr)]
  
  cat("\necdf...")
  msr_table$ecdf <- significance_measure(msr_table$msr, msr_table$msr_density, msr_ecdf_ref, inverted = F)
  cat("\ninverted_ecdf...")
  msr_table$inverted_ecdf <- significance_measure(msr_table$inverted_msr, msr_table$msr_density, msr_ecdf_ref, inverted = T)
  median_function <- extract_ecdf_function(msr_ecdf_ref, 0.5)
  cat("\nresidual...")
  msr_table$residual <- msr_table$msr-median_function(msr_table$msr_density)
  cat("\ninverted_residual...")
  msr_table$inverted_residual <- msr_table$inverted_msr-median_function(1-msr_table$msr_density)
  
  
  #other_stats_to_calculate <- c(autocorrelation, mean_bernoulli_entropy, meth_sd, missing_prop)
  cat("\nmeth_autocorrelation...")
  msr_table$meth_autocorrelation <- sapply(i_start, function(x){autocorrelation(wgbs$prop[x:(x+CG_window_size-1)])})
  #cat("\nmean_drift...")
  #msr_table$mean_drift <- sapply(i_start, function(x){drift(wgbs$prop[x:(x+CG_window_size-1)]/100)})
  cat("\nmeth_sd...")
  msr_table$meth_sd <- sapply(i_start, function(x){sd(wgbs$prop[x:(x+CG_window_size-1)]/100, na.rm=T)})
  cat("\nmissing_prop...")
  msr_table$missing_prop <- sapply(i_start, function(x){mean(is.na(binary[x:(x+CG_window_size-1)]))})
  cat("\nmean_bernoulli_entropy ...")
  msr_table$mean_bernoulli_entropy <- sapply(i_start, function(x){mean_be(wgbs$prop[x:(x+CG_window_size-1)]/100)})
  
  # save msr tables
  #msr_table = data.frame(i_start,msr_density,meth_rate,msr,inverted_msr, ecdf, inverted_ecdf, residual, inverted_residual, meth_autocorrelation, drift = mean_drift, meth_sd, missing_prop)
  saveRDS(msr_table, file = paste(dir,short_name, "_msr_complete_experimental_table_", CG_window_size, ".Rda", sep = ""))
}

methylation_summary_feature <- function(stat, data_table, wgbs)
{
  wgbs <- sum_strands(wgbs)
  wgbs$prop <- wgbs$prop/100
  if(! "i_end" %in% colnames(data_table)) # for fragments when
    data_table$i_end <- data_table$i_start + (data_table$i_start[2]-data_table$i_start[1])-1
  sapply(1:length(data_table$i_start), function(i) {
    if(is.na(data_table$i_start[i]+data_table$i_end[i])) return(NA)
    stat( wgbs$prop[data_table$i_start[i]:data_table$i_end[i]] )
    })
}

add_feature_to_table <- function(data_file, relative_wgbs, stat, stat_name)
{
  data_table <- readRDS(data_file)
  wgbs <- sum_strands(readRDS_if_needed(relative_wgbs))
  data_table[,stat_name] <- methylation_summary_feature(stat, data_table, wgbs)
  saveRDS(data_table,file=data_file)
}

#####################################################


# OTHER

noNa <- function(d, verbose = T)
{
  nn <- d[complete.cases(d),]
  if(verbose)
    cat("NAs: ", length(d[,1])/length(nn[,1]),"%\n")
  return(nn)
}

readRDS_if_needed <- function(file_or_table)
{
  if(is.character(file_or_table))
  {
    file_or_table <- readRDS(file_or_table); gc()
  }
  return(file_or_table)
}

show_perc_progress <- function(i,l)
{
  #cat('\014')
  #cat(paste0(round(i/n*100), '% completed'))
  
  if((i %% floor(l/100)) == 0)
    cat(floor(100*i/l), "%  ")
}

validation <- function(n,train_prop,data,lambda_model,evaluation_metric=tmse)
{
  sapply(1:n, function(i)
  {
    df <- train_and_test_split(data,train_prop,random = T)
    model <- lambda_model(df$train)
    evaluation_metric(model, df$test)
  })
}

from_table_to_data_frame <- function(M, label=NA, label_names=NA)
{
  l <- prod(dim(M))
  
  rn <- rownames(M)
  dffs <- lapply(rn, function(r)
  {
    dff <- data.frame(result=M[r,])
    dff$name <- r
    dff
  })
  
  df <- dffs[[1]]
  for(i in 2:length(dffs))
  {
    df <- rbind(df,dffs[[i]])
  }
  
  if(!is.na(label[1]))
  {
    for(i in 1:length(label))
      df[, label_names[i]] = as.factor(label[i])
  }
  
  
  df$name = as.factor(df$name)
  df
}

validation_on_subsets <- function(n,train_prop,data,lambda_model,evaluation_metric=tmse, id, merged)
{
  L <- length(data[,1])
  id_vector <- unique(data[,id])
  l <- length(id_vector)
  train_size <- l*train_prop
  cell_names <- unique(data$cell_name)
  M <- sapply(1:n, function(i)
  {
    id_train <- sample(id_vector, train_size, replace = F)
    train_rows <- which(data[,id] %in% id_train)
    test_rows = (1:L)[-train_rows]
    df_train <- data[train_rows,]
    df_test  <- data[test_rows,]
    
    
    sapply(cell_names, function(cn)
    {
      if(merged)
        model <- lambda_model(df_train)
      else model <- lambda_model(df_train[df_train$cell_name==cn,])
      evaluation_metric(model, df_test[df_test$cell_name==cn,])
    })
  })
  
  rownames(M) <- cell_names
  M
  
}

validation_on_subsets_data_frame <- function(n,train_prop,data,lambda_model,evaluation_metric=tmse, id, label=NA, label_name="model")
{
  M_single <- validation_on_subsets(n,train_prop,data,lambda_model,evaluation_metric, id, merged=F)
  M_merged <- validation_on_subsets(n,train_prop,data,lambda_model,evaluation_metric, id, merged=T)
  
  d1 <- from_table_to_data_frame(M_single, label="single", label_names="type")
  d2 <- from_table_to_data_frame(M_merged, label="shared", label_names="type")
  df <- rbind(d1,d2)
  
  if(!is.na(label))
    df[,label_name] <- as.factor(label)
  
  return(df)
  
}

table_to_frame <- function(M,cnames=NA)
{
  l <- prod(dim(M))
  r <- array(dim = l)
  c <- array(dim = l)
  v <- array(dim = l, data = 0.1)
  df <- data.frame(v,r,c)
  
  i = 1
  
  for(cl in colnames(M))
  {
    for(rw in rownames(M))
    {
      df[i,"v"] <- M[rw,cl]
      df[i,"r"] <- rw
      df[i,"c"] <- cl
      i <- i+1
    }
  }
  
  df$r = as.factor(df$r)
  df$c = as.factor(df$c)
  if(!is.na(cnames[1]))
    colnames(df)<-cnames
  
  return(df)
}

add_log_features <- function(data_table, names, epsilon)
{
  new <- log(data_table[, names]+epsilon,2)
  colnames(new) <- paste0("log_", colnames(new))
  cbind(data_table, new)
}

#####################################################


# STAT UTILS

drift <- function(x) {mean(0.5-abs(x-0.5), na.rm=T)}

Bernoulli_entropy <- function(p) {(-p*log(p+1e-20,2)-(1-p)*log((1-p)+1e-20,2))}

mean_be <- function(p) {mean(Bernoulli_entropy(p),na.rm=T)}

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

tmse <- function(model, test_model_data)
{
  y_name <- as.character(formula(model)[2])
  y <- test_model_data[,y_name]
  y_mean <- mean(y, na.rm=T)
  y_pred <- predict(model, test_model_data,  type="response")
  RSS <- sum((y-y_pred)^2)
  TSS <- sum((y-y_mean)^2)
  1-RSS/TSS
}

prediction_correlation_score <- function(model, test_data)
{
  y_name <- as.character(formula(model)[2])
  c <- cor.test(predict(model, test_data,  type="response"), test_data[,y_name])$estimate
  names(c) <- NULL
  c
}

tmse.glmnet <- function(model, test_model_data, y_name, x_names)
{
  newx = as.matrix.data.frame(test_model_data[, x_names])
  1-(var(predict.glmnet(model, newx = newx,  type="response")-test_model_data[,y_name]))/var(test_model_data[,y_name])
}

prediction_correlation_score.glmnet <- function(model, test_data, y_name, x_names)
{
  newx = as.matrix.data.frame(test_data[, x_names])
  c <- cor.test(predict.glmnet(model, newx = newx,  type="response"), test_data[,y_name])$estimate
  names(c) <- NULL
  c
}

autocorrelation <- function(v, method = "pearson")
{
  x <- 2:length(v)
  
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

model_this <- function(response_variable_name, predictors = NA)
{
  p <- "."
  if(!any(is.na(predictors)))
    p <- paste(predictors, collapse = " + ")
  as.formula(paste(response_variable_name, "~ ", p))
}

lasso <- function(response_variable_name, df, lambda = 0.01, alpha = 1)
{
  x_vars <- (as.matrix.data.frame(df[, names(df) != response_variable_name]))
  y_var <- (df[, response_variable_name])
  lambda_seq <- 10^seq(2, log(lambda,10), by = -.1)
  cv_output <- cv.glmnet(x_vars, y_var, alpha = alpha, lambda = lambda_seq, nfolds = 5, standardize = T)
  glmnet(x_vars, y_var, alpha = alpha, lambda = cv_output$lambda.min)
}

lasso_normalized_coefficients <- function(response_variable_name, df, lambda = 0.01, alpha = 1)
{
  x_vars <- scale(as.matrix.data.frame(df[, names(df) != response_variable_name]))
  y_var <- scale(df[, response_variable_name])
  lambda_seq <- 10^seq(2, log(lambda,10), by = -.1)
  cv_output <- cv.glmnet(x_vars, y_var, alpha = alpha, lambda = lambda_seq, nfolds = 5, standardize = T)
  m <- glmnet(x_vars, y_var, alpha = alpha, lambda = cv_output$lambda.min)
  m$beta
}

single_correlations <- function(data_table, response, predictors)
{
  sapply(predictors, function(p){cor.test(data_table[,response], data_table[,p])$estimate})
}

train_and_test_split <- function(df, train_prop, random, test_chr = NA)
{
  l = length(df[,1])
  train_length = train_prop*l
  
  if(random) 
  {
    train_rows = sample(1:l, size = train_length, replace = F)
  } else train_rows = 1:train_length
  
  if(!is.na(test_chr[1]))
  {
    train_rows = which(!(df$chr %in% test_chr))
    cat("\ntrain_prop:", length(train_rows)/l)
  }
  
  test_rows = (1:l)[-train_rows]
  
  List(train= df[train_rows,], test=df[test_rows,])
}

undersample_data_frame <- function(df, n)
{
  df[sample(1:length(df[,1]),replace = F, size = n),]
}

check_linear_model <- function(train_model_data, test_model_data, response_variable, predictors, summary = T)
{
  formula <- as.formula(paste(response_variable, "~ ."))
  model = lm(formula, data = train_model_data[, c(response_variable, predictors)])
  if(summary)
    print(summary(model))
  else
  {
    cat("predictors: ", predictors, "\n")
  }
  cat("Test data R squared: ", tmse(model, test_model_data))
  #1-(var(predict(model, test_model_data, type="response")-test_model_data$log_tpm))/var(test_model_data$log_tpm) )
}

#####################################################


# PLOTS

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
    plot(x, round(p), pch = "|", ylab = " ", main=m, xlab="")
  }
}

show_region_meth2 <- function(wgbs, chromosome, start, end, prop = T, genomewide = F, title = NA, title_pos = -6)
{
  w <- wgbs[wgbs$chr == chromosome & wgbs$Cpos>=start & wgbs$Cpos<=end, ]
  p <- w$prop/100
  m <- length(p)
  
  if(m<1 | all(is.na(p)))
    return(F)
  
  if(is.na(title))
    m <- paste("CG: ", m, "  nucleotides: ", end-start, "\n", chromosome, start,"-", end)
  else{
    m <- title
  }
  #cat(p, m)
  
  x <- 1:length(p)
  if(genomewide)
    x <- w$Cpos
  
  if(prop)
  {
    plot(x, p)
  }
  else
  {
    p <- standard_binaryzer2(w)*1
    plot(x, p, pch = "|", ylab = m, xlab="", yaxt="none")
    axis(2, c(0,1), las=2, font=2)
  }
  title(m, line = title_pos)
  
  
  nas <- which(is.na(p))
  y <- rep(1, length(nas))
  points(nas,y,pch = "|", col="gray")
  y <- rep(0, length(nas))
  points(nas,y,pch = "|", col="gray")
  
}

show_gene_meth <- function(wgbs, genebody_annotation, i, prop = T, genomewide = F, margins = 0)
{
  g <- genebody_annotation[i,]
  print(g)
  show_region_meth(wgbs, g$chr, g$start - margins, g$end + margins, prop, genomewide)
}

show_fragment_meth <- function(wgbs, data_table, i = 1, prop = T, genomewide = F, title = NA, title_pos = -6)
{
  g <- data_table[i,]
  print(g)
  show_region_meth2(wgbs, data_table$chr[i], data_table$start[i], data_table$end[i], prop, genomewide, title, title_pos = title_pos)
}

ggplot_meth_fragment <- function(wgbs, i_start, size, prop=T,genomewide=F, minimum_reads=1, loess_span = NA)
{
  logit <- function(p){log(p/(1-p))}
  epsilon <- 1e-3
  
  fragment <- wgbs[i_start:(i_start+size-1),]
  fragment$x <- 1:size
  xscale <- "CpG index"
  yscale <- "meth prop"
  if(genomewide)
  {
    xscale <- "nucleotide position"
    fragment$x <- fragment$Cpos
  }
  
  fragment$prop <- fragment$prop/100; fragment$prop[fragment$prop==0] <- epsilon; fragment$prop[fragment$prop==1] <- 1-epsilon
  fragment <- fragment[fragment$reads>=minimum_reads, ]
  
  if(!prop)
  {
    fragment<-fragment[fragment$prop!=0.5,]
    fragment$prop<-(fragment$prop>0.5)+0
    yscale <- "meth state"
  }
  
  
  if(!is.na(loess_span))
    model <- loess(logit(prop) ~ x, data = fragment, span = loess_span)
  else{
    model <- gam(logit(prop) ~ s(x), data = fragment, gamma = 0)
  }
  
  
  fragment$fit <- sigmoid(predict(model, type="response"))
  
  pl <- ggplot(data=fragment, aes(x = x, y = fit), alpha = 0.9) +
    geom_point(data=fragment, aes(x = x, y = prop), alpha=0.5)+
    geom_line(data=fragment, aes(x = x, y = fit), alpha = 0.9, color = "blue")+
    ylim(0,1)+
    ylab(yscale)+
    xlab(xscale)
  
  return(pl)
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

ggscatter <- function(data_table, xn, yn, alpha = 0.1)
{
  df <- data_table[, c(xn, yn)]
  colnames(df) <- c("x","y")
  ((ggplot(df,aes(x=x,y=y)) + geom_point(alpha = alpha))) + ylab(yn) + xlab(xn)
}

gg_expression_scatterplot <- function(data_table, xname, yname, response_variable, undersample = NA, points_alpha = 1, title="")
{
  c1 <- '#666699'
  c2 <- '#00b300'
  
  if(!is.na(undersample))
    data_table <- undersample_data_frame(data_table, undersample)
  
  df <- data.frame(x=data_table[,xname], y=data_table[,yname], ltpm=data_table[,response_variable])
  df <- df[complete.cases(df),]
  #data_table <- data_table[,c(xname, yname, response_variable)]
  ggplot(df, aes(x, y, color=ltpm)) +
    geom_point(alpha = points_alpha) +
    scale_colour_gradient(low = c1, high = "green", na.value = NA, name = response_variable) +
    theme(legend.position=c(0,1), legend.justification=c(0,1)) +
    labs(x = xname, y=yname, title = title)
}

ggplot_meth_fragment_comparison <- function(wgbs1, wgbs2, table, i,prop=T,genomewide=F, minimum_reads=1, loess_span = NA)
{
  table <- table[i,]
  i_start <- table$i_start
  size <- table$i_end - table$i_start
  
  p1 <- ggplot_meth_fragment(wgbs1, i_start, size, prop,genomewide, minimum_reads, loess_span)
  p2 <- ggplot_meth_fragment(wgbs2, i_start, size, prop,genomewide, minimum_reads, loess_span)
  
  grid.arrange(p1,p2,nrow=2)
  print(table)
}

ggplot_meth_fragment2 <- function(wgbs, table_row, prop=T,genomewide=F, minimum_reads=1, loess_span = NA)
{
  i_start <- table_row$i_start
  size <- table_row$i_end - table_row$i_start - 1
  ggplot_meth_fragment(wgbs, i_start, size, prop,genomewide, minimum_reads, loess_span)
    
}

gg_CpG_density <- function(wgbs, table_row)
{
  data <- data.frame(x=wgbs$Cpos[table_row$i_start:(table_row$i_start+999)])
  ggplot(data, aes(x=x)) +
    geom_density(adjust=1, alpha=.4, fill="gray")+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
}

MSR_CpG_visualization <- function(wgbs, table_row, size=1000)
{
  x = wgbs$Cpos[table_row$i_start:(table_row$i_start+size-1)]
  MSR_visualization(x)
}

#####################################################