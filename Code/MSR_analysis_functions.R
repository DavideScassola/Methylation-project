suppressMessages(library(data.table))
suppressMessages(library(pracma))
suppressMessages(library(parallel))
suppressMessages(library(scales))
#suppressMessages(library("markovchain"))

calculate_relevance_from_counts <- function(vector, M)
{
  M <- sum(vector) #efficiency ?
  t <- table(vector)
  m_k <- as.integer(t)
  k <- as.integer(names(t))
  
  # probability that a methyle belongs to a bin with k methyles
  p <- (k*m_k)/M
  return(-sum(p*log(p+1E-10, M)))
}


calculate_resolution_from_counts <- function(counts, M)
{
  M <- sum(counts)
  p <- counts/M
  return(-sum(p*log(p+1E-10,base = M)))
}

meth_prop <- function(meth_vector)
{
  valid_meth <- meth_vector[!is.na(meth_vector)]
  sum(valid_meth)/length(valid_meth)
}

stochastic_round <- function(v)
{
  v_floor <- floor(v)
  return(v_floor + rbinom(length(v), 1, v-v_floor))
}

replace_nas_with_zeros <- function(effective_counts_and_na, bin_size, methylation_vector)
{
  return(effective_counts_and_na[1,])
}

replace_nas_with_meth_prop <- function(effective_counts_and_na, bin_size, methylation_vector)
{
  return(stochastic_round(effective_counts_and_na[1,] + effective_counts_and_na[2,]*meth_prop(methylation_vector)))
}

replace_nas_with_bin_prop <- function(effective_counts_and_na, bin_size, methylation_vector)
{
  estimated_prop <- effective_counts_and_na[1,]/(bin_size-effective_counts_and_na[2,])
  estimated_prop[is.na(estimated_prop)] <- 0
  return(stochastic_round(effective_counts_and_na[1,] + effective_counts_and_na[2,]*estimated_prop))
}

replace_nas_hybrid <- function(effective_counts_and_na, bin_size, methylation_vector)
{
  bin_estimated_prop <- effective_counts_and_na[1,]/(bin_size-effective_counts_and_na[2,])
  bin_estimated_prop[is.na(bin_estimated_prop)] <- 0
  weights <- effective_counts_and_na[2,]/bin_size
  weighted_estimated_prop <- weights*meth_prop(methylation_vector) + (1-weights)*bin_estimated_prop
  return(stochastic_round(effective_counts_and_na[1,] + effective_counts_and_na[2,]*weighted_estimated_prop))
}

replace_nas_hybrid_stochastic <- function(effective_counts_and_na, bin_size, methylation_vector)
{
  bin_estimated_prop <- effective_counts_and_na[1,]/(bin_size-effective_counts_and_na[2,])
  bin_estimated_prop[is.na(bin_estimated_prop)] <- 0
  weights <- (effective_counts_and_na[2,]/bin_size)
  weighted_estimated_prop <- weights*meth_prop(methylation_vector) + (1-weights)*bin_estimated_prop
  add <- rbinom(length(effective_counts_and_na[2,]), size = effective_counts_and_na[2,], weighted_estimated_prop)
  return(effective_counts_and_na[1,] + add)
}

corrected_cumsum <- function(vector)
{
  l <- length(vector)+1
  c = rep_len(0, l)
  c[2:l] <- cumsum(vector)
  return(c)
}

# PLEASE USE A SPARSE MATRIX IF POSSIBLE
calculate_relevance_and_resolution_ignoring_nas <- function(methylation_vector, cumulative_sum_vector, bin_size) 
{

  l <- length(methylation_vector)
  #if(bin_size==1) return(c(0,1,1,1))
  #print(bin_size)
  if(bin_size==l) return(c(0,0,l,NA))
  
  starting_points <- ((1:ceiling(l/bin_size))*bin_size)-bin_size+1
  counts  <-  cumulative_sum_vector[pmin(starting_points+bin_size, l+1)]-cumulative_sum_vector[starting_points]
  
  if(l/bin_size > 1E8)
  {
    # it's useful to save memory, but for large bin size this would just slow down
    remove(starting_points)
    gc(verbose = F)
  }

  M <- cumulative_sum_vector[l+1]
  relevance <- calculate_relevance_from_counts(counts, M)
  resolution <- calculate_resolution_from_counts(counts, M)

  return(c(relevance, resolution, bin_size, 1))
}


calculate_relevance_and_resolution <- function(methylation_vector, cumulative_sum_vector, cumulative_na_vector, bin_size, na_tolerance, no_data_out = F, na_values_handler, verbose = T, throw_away_threshold = 0) 
{
  # cumulative_sum_vector[end]-cumulative_sum_vector[start]=sum(methylation_vector[start:(end-1)])
  
  l <- length(methylation_vector)
  data_size <- l-cumulative_na_vector[l+1]
  cumulative_vector_length <- length(cumulative_sum_vector)
  
  if(bin_size==1) return(c(0,1,1,1))
  if(bin_size==l) return(c(0,0,l,NA))
  
  counts_and_na <- array(dim = c(2,(l/bin_size)+1))
  valid_bins_found <- 0
  start <- 1
  max_n_na <- floor(bin_size*na_tolerance)
  
  while(start <= l)
  {
    end <- min(start+bin_size, cumulative_vector_length)
    n_of_na <- cumulative_na_vector[end]-cumulative_na_vector[start]
      
    if(n_of_na>max_n_na)
    {
      start = start + n_of_na - max_n_na
    }
    else
    {
      valid_bins_found <- valid_bins_found + 1
      end <- min(start+bin_size, cumulative_vector_length)
      counts_and_na[,valid_bins_found] <- c(cumulative_sum_vector[end]-cumulative_sum_vector[start], n_of_na)
      start <- start + bin_size
    }
    
  }
  
  if(valid_bins_found==0) return(c(NA,NA,bin_size, 0))
  effective_data <- (bin_size*valid_bins_found)/data_size
  if(effective_data<throw_away_threshold) return(c(NA,NA,bin_size, effective_data))

  #cat("ciao\n")
  effective_counts_and_na <- counts_and_na[, 1:(valid_bins_found),drop=F]

  # na values handling
  effective_counts <- na_values_handler(effective_counts_and_na, bin_size, methylation_vector)

  
  if(no_data_out)
  {
    effective_counts <- counts_and_na[1, 1:(valid_bins_found+1)]
    # adding to counts the "bin" including all excluded observation
    effective_counts[valid_bins_found+1] = cumulative_sum_vector[l+1]-sum(effective_counts[1:valid_bins_found])
  }
  

  
  M <- sum(effective_counts)

  relevance <- calculate_relevance_from_counts(effective_counts, M)
  resolution <- calculate_resolution_from_counts(effective_counts, M)

  if(verbose) cat("bin_size:", bin_size, " valid_bins_found:", valid_bins_found, "effective data:", (bin_size*valid_bins_found)/data_size, "\n")

  return(c(relevance, resolution, bin_size, effective_data))
}


good_bin_sizes <- function(n, max_bins = 100, losing_data_prop = 0.1)
{
  start <- 1
  end <- n
  
  log_spaced_bins <- unique(round((exp(linspace(log(start), log(end), n = max_bins)))))
  log_spaced_bins <- sort(unique(floor(end/log_spaced_bins)))
  b <- log_spaced_bins[((n%%log_spaced_bins)/n)<=losing_data_prop]
  #cat(b, "\n")
  return(b)
}

good_bin_sizes <- function(n, max_bins = 100, losing_data_prop = NA)
{
  return(unique(round((exp(linspace(log(1), log(n), n = max_bins))))))
}


calculate_relevance_resolution_vector <- function(methylation_vector, na_tolerance = 0, max_bins = 100, no_data_out = F, na_values_handler = replace_nas_hybrid, invert = F, verbose = T)
{
  if(invert) methylation_vector <- !methylation_vector
  
  l <- length(methylation_vector)
  bin_sizes <- good_bin_sizes(l, max_bins)
  
  nas <- is.na(methylation_vector)
  cumulative_nas_vector <- rep_len(0,l+1)
  cumulative_nas_vector[2:(l+1)] <- cumsum(nas)
  
  replaced_nas_vector <- methylation_vector
  replaced_nas_vector[nas] <- 0
  
  cumulative_sum_vector <- rep_len(0,l+1)
  cumulative_sum_vector[2:(l+1)] <- cumsum(replaced_nas_vector)
  
  return(sapply(bin_sizes, function(x) 
  { 
    calculate_relevance_and_resolution(methylation_vector, cumulative_sum_vector, cumulative_nas_vector, x, na_tolerance =na_tolerance , no_data_out = no_data_out, na_values_handler = na_values_handler, verbose = verbose)
  }))
}

add_max_resolution_point_if_needed <- function(rr)
{
  if(rr[1, 1]==(0))
    return(rr)
  l <- length(rr[1,])
  x <- array(NA, dim = c(4,l+1))
  x[, 2:(l+1)] <- rr
  x[, 1] <- c(0,1,NA, NA)
  x
}

calculate_relevance_resolution_vector_ignoring_nas <- function(methylation_vector, max_bins = 100, verbose = T, minimum_bin_size = 1, invert = F)
{
  if(invert) methylation_vector <- !methylation_vector
  
  start_time <- Sys.time()
  l <- length(methylation_vector)
  bin_sizes <- good_bin_sizes(l, max_bins)
  bin_sizes <- bin_sizes[(bin_sizes>minimum_bin_size) | (bin_sizes==1) ]
  
  
  nas <- is.na(methylation_vector)
  replaced_nas_vector <- methylation_vector
  replaced_nas_vector[nas] <- 0
  
  if(verbose) cat("Calculating cumulative sum vector \n")
  cumulative_sum_vector <- corrected_cumsum(replaced_nas_vector)

  #################################
  # I don't want my RAM to explode
  remove(replaced_nas_vector, nas)
  if(l>1e7) gc()
  #################################
    
  out <- (sapply(bin_sizes, function(x) 
  { 
    if(verbose) cat(x, "...\n")
    calculate_relevance_and_resolution_ignoring_nas(methylation_vector, cumulative_sum_vector, bin_size = x)
  }))
  
  out <- add_max_resolution_point_if_needed(out)
  
  end_time <- Sys.time()
  if(verbose) print(end_time-start_time)
  
  return(out)
}

calculate_relevance_resolution_vector_with_positions <- function(positions, max_bins = 100, verbose = T, pre_binning = 1)
{
  start_time <- Sys.time()
  l <- length(positions)
  pre_binned <- hist(positions,breaks = l/(pre_binning), plot = F)$counts
  bin_sizes <- good_bin_sizes(length(pre_binned), max_bins)

  if(verbose) cat("Calculating cumulative sum vector \n")
  cumulative_sum_vector <- corrected_cumsum(pre_binned)
  
  out <- (sapply(bin_sizes, function(x) 
  { 
    if(verbose) cat(x, "...\n")
    calculate_relevance_and_resolution_ignoring_nas(pre_binned, cumulative_sum_vector, bin_size = x)
  }))
  
  end_time <- Sys.time()
  if(verbose) print(end_time-start_time)
  
  return(out)
}

genome_MSR <- function(methylation_positions, minimum_bin_size = 10, verbose = F, invert = F, msr = T, max_bins = 100)
{
  if(length(methylation_positions)<2) return(NA)
  v = methylation_positions - min(methylation_positions) + 1
  v = sparseVector(i = v, x = T, length = max(v))
  rr <- calculate_relevance_resolution_vector_ignoring_nas(v, minimum_bin_size = minimum_bin_size, invert = invert, verbose = verbose, max_bins = max_bins)
  if(msr)
    return(MSR_area(rr))
  return(rr)
}

positions_MSR <- function(positions, discretization_bin_size = NA, verbose = F, msr = T, max_bins = 1e6, discrete = F,  max_gbins = 100)
{
  if(length(positions)<2) return(NA)
  positions <- positions[!is.na(positions)]
  positions <- positions-min(positions)
  
  if(is.na(discretization_bin_size))
  {
    s <- sort(positions)
    discretization_bin_size <- min(s[2:length(s)]-s[1:(length(s)-1)])
  }
  #print(discretization_bin_size)
  breaks <- floor(min(max_bins, (max(positions)-min(positions))/discretization_bin_size)*1)
  #print(breaks)
  v <- hist(positions, plot = F, breaks = breaks)$counts
  #print(v)
  
  if(discrete)
  {
    gc()
    positions <- sort(positions+1)
    v <- sparseVector(i = positions, x = T, length = max(positions))
  }    

  rr <- calculate_relevance_resolution_vector_ignoring_nas(v, minimum_bin_size = 1, verbose = verbose,  max_bins = max_gbins)
  if(msr)
    return(MSR_area(rr))
  cat(MSR_area(rr))
  return(rr)
}


# supposing resolution is in decreasing order
MSR_area<- function(rr_vector, M_correction = 10)
{
  return(trapz(x = rev(rr_vector[2,]), y = rev(rr_vector[1,]))*log10(M_correction))
}

MSR_area<- function(rr_vector)
{
  ord = order(rr_vector[2,])
  return(trapz(x = rr_vector[2,][ord], y = rr_vector[1,][ord] ))
}

calculate_MSR_area<- function(v)
{
  rr <- calculate_relevance_resolution_vector_ignoring_nas(v, verbose = F)
  MSR_area(rr)
}


############################### PLOT FUNCTIONS ###############################

resolution_relevance_plot <- function(relevance_resolution_vector, color = rgb(0, 200, 50, maxColorValue = 255), plotter = plot, area_alpha = 0.5)
{
  ord = order(relevance_resolution_vector[2,])
  xx = relevance_resolution_vector[2,][ord]
  yy <- relevance_resolution_vector[1,][ord]
  max_rel <- max(0.5, max(yy))
  plotter(x = xx, y = yy, type="o" , lwd=3, xlab = "resolution", ylab = "relevance", ylim = c(0,max_rel), xlim = c(0,1), col = alpha(color, 0.9))
  polygon( 
    c(min(xx), xx , max(xx)) , 
    c( min(yy) , yy , min(yy)) , 
    col=alpha(color, area_alpha) , border=F)
}

bin_size_relevance_plot <- function(relevance_resolution_vector)
{
  plot(x = relevance_resolution_vector[3,], y = relevance_resolution_vector[1,], type = "l", log = "x")
}

bin_size_resolution_plot <- function(relevance_resolution_vector)
{
  plot(x = relevance_resolution_vector[3,], y = relevance_resolution_vector[2,], type = "l", log = "x")
}

compare_resolution_relevance_plot <- function(relevance_resolution_vector_list, legend_names, title)
{
  plot(x = (relevance_resolution_vector_list[[1]])[2,], y = (relevance_resolution_vector_list[[1]])[1,], type = "l", col = 1, xlab = "resolution", ylab = "relevance", ylim = c(0,0.5), xlim = c(0,1))
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  n <- length(relevance_resolution_vector_list)
  if(n>1)
  {
    for(i in 2:n)
    {
      lines(x = (relevance_resolution_vector_list[[i]])[2,], y = (relevance_resolution_vector_list[[i]])[1,], type = "l", col = i)
    }
  }
  title(title)
  legend("topleft", legend=legend_names, col=1:n, lty = 1, cex = 0.7, y.intersp = 0.5)

}

compare_bin_size_relevance_plot <- function(bsrvl, legend_names, title)
{
  l <- length(bsrvl[[1]][1,])
  n <- length(bsrvl)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  plot(x = (bsrvl[[1]])[3,], y = (bsrvl[[1]])[1,], type = "p", col = 1, xlab = "bin_size", ylab = "relevance", ylim = c(0,0.6), log = "x")
  lines(x = (bsrvl[[1]])[3,], y = (bsrvl[[1]])[4,], type = "l", col = 1)

  if(n>1)
  {
    for(i in 2:n)
    {
      lines(x = (bsrvl[[i]])[3,], y = (bsrvl[[i]])[1,], type = "p", col = i)
      lines(x = (bsrvl[[i]])[3,], y = (bsrvl[[i]])[4,], type = "l", col = i)
    }
  }
  title(title)
  legend("topright", legend=legend_names, col=1:n, lty = 1, cex = 0.8, y.intersp = 0.8)
  
}

compare_bin_size_resolution_plot <- function(bsrvl, legend_names, title)
{
  l <- length(bsrvl[[1]][2,])
  n <- length(bsrvl)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  
  plot(x = (bsrvl[[1]])[3,], y = (bsrvl[[1]])[2,], type = "p", col = 1, xlab = "bin_size", ylab = "resolution", ylim = c(0,1), xlim = c(1, l*100000), log = "x")
  lines(x = (bsrvl[[1]])[3,], y = (bsrvl[[1]])[4,], type = "l", col = 1)
  
  if(n>1)
  {
    for(i in 2:n)
    {
      lines(x = (bsrvl[[i]])[3,], y = (bsrvl[[i]])[2,], type = "p", col = i)
      lines(x = (bsrvl[[i]])[3,], y = (bsrvl[[i]])[4,], type = "l", col = i)
    }
  }
  title(title)
  legend("topright", legend=legend_names, col=1:n, lty = 1, cex = 0.8, y.intersp = 0.8)
  
}

compare_resolution_relevance_plot_confidence <- function(relevance_resolution_vector_list, legend_names, title)
{
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  plot(x = (relevance_resolution_vector_list[[1]])[2,], y = (relevance_resolution_vector_list[[1]])[1,], type = "p", col = 1, xlab = "resolution", ylab = "relevance", ylim = c(0,0.8))
  lines(x = (relevance_resolution_vector_list[[1]])[2,], y = (relevance_resolution_vector_list[[1]])[4,], type = "h", col = 1)
  n <- length(relevance_resolution_vector_list)
  if(n>1)
  {
    for(i in 2:n)
    {
      lines(x = (relevance_resolution_vector_list[[i]])[2,], y = (relevance_resolution_vector_list[[i]])[1,], type = "p", col = i)
      lines(x = (relevance_resolution_vector_list[[i]])[2,], y = (relevance_resolution_vector_list[[i]])[4,], type = "h", col = i)
      
    }
  }
  title(title)
  legend("topleft", legend=legend_names, col=1:n, lty = 1, cex = 0.8, y.intersp = 0.8)
  
}

remove_low_confidence_points <- function(rr_vector, threshold = 0.01)
{
  return(rr_vector[,rr_vector[4,]>=threshold])
}

remove_low_confidence_points_list <- function(rr_vector_list, threshold = 0.01)
{
  return(lapply(rr_vector_list, function(x) { return(remove_low_confidence_points(x, threshold))}))
}

rr_plots <- function(rr_list, title = "", legend_names = "", threshold = 0.05)
{
  compare_resolution_relevance_plot(rr_list, legend_names, title)
  par(ask=TRUE)
  compare_resolution_relevance_plot_confidence(rr_list, legend_names, title)
  compare_resolution_relevance_plot(remove_low_confidence_points_list(rr_list, threshold = threshold), legend_names = legend_names, title = paste(title, "\n removing points using less than ", threshold*100, "% effective data", sep = ""))
  compare_bin_size_relevance_plot(rr_list, legend_names, title)
  compare_bin_size_resolution_plot(rr_list, legend_names, title)
  par(ask=FALSE)
}




############################### AUX FUNCTIONS ###############################

get_max_information_vector <- function(size)
{
  probs = rep_len(0, size)
  n = floor(sqrt(size))
  for(i in 1:n)
  {
    probs[(1+((i-1)*n)):(i*n)] = 1/i
  }

  return(rbinom(n = size, size = 1, prob = probs))
}

mutilated <- function(v, missing_proportion)
{
  l <- length(v)
  mask <- as.logical(rbinom(l, size = 1, prob = missing_proportion))
  m <- v
  m[mask] = NA
  return(m)
}

nameof <- function(v1) {
  deparse(substitute(v1))
}

missing_data_experiment <- function(full_vector, missing_data_rates = c(0,0.1,0.2,0.5,0.7,0.9), max_bins = 100, no_data_out = F, na_values_handler, info = "", na_tolerance = 0)
{
  cat("missing_data_rates: ", missing_data_rates, "\n")
  vectors <- lapply(missing_data_rates, function(r) mutilated(full_vector, r))
  rr_missings = mclapply(vectors, function(v) calculate_relevance_resolution_vector(v,na_tolerance = na_tolerance, max_bins = max_bins, no_data_out = no_data_out, na_values_handler = na_values_handler), mc.cores = 1)
  legend_names <- paste(missing_data_rates*100, "% miss.", sep = "")
  rr_plots(rr_missings, title = info, legend_names = legend_names)
  return(list(rr_data = rr_missings, missing_rates = missing_data_rates, na_tolerance = na_tolerance, max_bins = max_bins, info = info))
}



apply_same_missing_data_pattern <- function(vector, mask_vector)
{
  v <- vector
  v[is.na(mask_vector)] = NA
  return(v)
}


##############################################################################
suppressMessages(library(Matrix))
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(require(Biostrings))

#load("DataCleaning/CpG_sites_dataframe.RData")

nucleotides_pattern_positions <- function(chromosome, pattern, Genome = BSgenome.Hsapiens.UCSC.hg38)
{
  ranges <- matchPattern(pattern,(Genome[[chromosome]]))
  return(IRanges(ranges)@start)
}

binary_nucleotides_pattern_positions <- function(chromosome, pattern, Genome = BSgenome.Hsapiens.UCSC.hg38)
{
  positions <- IRanges(matchPattern(pattern,(Genome[[chromosome]])))@start
  chromosome_size <- length(Genome[[chromosome]])
  return(sparseVector(i = positions, length = chromosome_size, x = T))
}

mouse_dinucleotides_experiment <- function(pattern_list, chromosome)
{
  
  rr_list <- lapply(pattern_list, function(p) {
    chr_p <- binary_nucleotides_pattern_positions(chromosome,p)
    rr_chr_p <- calculate_relevance_resolution_vector_ignoring_nas(chr_p)
    remove(chr_p)
    gc()
    return(rr_chr_p)
  })
  
  gc()
  
  areas <- array(length(rr_list))
  for(i in 1:length(rr_list))
  {
    areas[i] <- MSR_area(rr_list[[i]], M_correction = sum(binary_nucleotides_pattern_positions(chromosome,pattern_list[i])))
  }
  
  cat("\nMSR corrected areas: ", areas)
  
  rr_plots(rr_list, title = "dinucleotides_MSR_comparison", legend_names = pattern_list)
  
  return((List(rr = rr_list, areas = areas)))
}


#############################################################################

subset_positions <- function(pos, start, size) { pos[pos>=start & pos<=(start+size)] }

pos_to_neighbor_counts  <- function(pos, dinucleotides_neighborhood_ranges, centered = T)
{
  max_half = max(dinucleotides_neighborhood_ranges)/2
  v = (pos - min(pos) + 1)+max_half
  v2 = sparseVector(i = v, x = T, length = max(v)+max_half)
  cumulative_sum_vector <- corrected_cumsum(v2)
  #cat("\ncumulative_sum_vector: ", cumulative_sum_vector)
  #cat("\npos: ", pos, "\n")
  #cat("v: ", v, "\n")
  result = sapply(dinucleotides_neighborhood_ranges, function(r)
  {
    cumulative_sum_vector[v+r/2]-cumulative_sum_vector[v-r/2]
  })
  remove(cumulative_sum_vector)
  colnames(result) <- (dinucleotides_neighborhood_ranges)
  gc()
  return(as.data.frame(result))
  #sum(pos<pos+half & pos>pos-half)
}

filter_odd_position_values <- function(v)
{
  v[2*(0:(ceiling(length(v)/2)-1))+1]
}

get_CpG_densities <- function(dinucleotides_neighborhood_ranges, Genome = BSgenome.Hsapiens.UCSC.hg38, data = NA)
{
  chr_list = chromosomes(Genome = BSgenome.Hsapiens.UCSC.hg38)
  n = length(chr_list)

  CpGlists <- lapply(chr_list, function(chr)
    {
     cat("processing ", chr, "\n")
    if(type(data)=="logical")
      pos <- nucleotides_pattern_positions(chr, "CG", Genome = Genome)
    else
      pos <- filter_odd_position_values(filter_chromosome(data,chr)$Cpos)
     gc()
     pos_to_neighbor_counts(pos, dinucleotides_neighborhood_ranges)
    })
     gc()
  
  result <- data.frame()
  
  for(p in CpGlists)
  {
    result = rbind(result,p)
  }

  return(result)
}

autocor <- function(v, lag)
{
  l = length(v)
  plot(v[1:(l-lag)],v[(1+lag):l])
  cor.test(v[1:(l-lag)],v[(1+lag):l])
}

autotable <- function(v, lag)
{
  l = length(v)
  table(v[1:(l-lag)],v[(1+lag):l])
}

random_binary_vector <- function(l, prop)
{
  M <- round(l*prop)
  bin <- array(dim = l, 0)
  bin[sample(1:l, M)] <- 1
  return((bin))
}

sample_msr <- function(l, prop)
{
  v <- random_binary_vector(l, prop)
  rr_v <- calculate_relevance_resolution_vector_ignoring_nas(v, verbose = F)
  MSR_area(rr_v)
}

msr_samples <- function(l, prop, sample_size, cores = 1, verbose = F)
{
  start_time <- Sys.time()
  msr_samples = mcmapply(1:sample_size, mc.preschedule = T, mc.cores = cores, FUN =  function(n)
  {
    #if(verbose) cat(n, " ")
    sample_msr(l, prop)
  })
  if(verbose) cat("time: ", Sys.time()-start_time)
  
  return(msr_samples)
}

msr_significance_experiment <- function(l,prop_list, samples_size, cores = 1, verbose = F)
{
  start_time <- Sys.time()
  exp = mcmapply(prop_list, mc.preschedule = T, mc.cores = cores, FUN =  function(p)
  {
    if(verbose) cat(p, " ")
    msr_samples(l, p, samples_size, cores = 1, verbose = F)
  })
  cat("time: ", Sys.time()-start_time)
  
  exp = t(exp)
  rownames(exp) <- prop_list
  return(List(data=exp, prop_list=prop_list))
}

get_msr_ecdf <- function(l, prop, sample_size, cores = 1)
{
  return(ecdf(msr_samples(l, prop, sample_size, verbose = F, cores = cores)))
}

get_msr_ecdfs <- function(l, prop_list, sample_size, cores = 1)
{
  sapply(prop_list, function(prop)
    {
      if(prop==0 || prop==1) e = NA
      else e = get_msr_ecdf(l, prop, sample_size, cores = cores)
      return(List(cdf = e, prop=prop))
    })
}

inverse = function(f, lower = -100, upper = 100) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]$root
}

msr_confidence_line <- function(ecdfs, confidence = 0.99)
{
  sapply(ecdfs, function(f)
    {
    if( f$prop==0 || f$prop==1) return(c(f$prop,NA))
    c(f$prop, inverse(f$cdf, 0, 0.4)(confidence))
  })
}



add_msr_confidence_line <- function(ecdfs, confidence = 0.99, col = 1, lty = 2)
{
  conf = msr_confidence_line(ecdfs, confidence)
  lines(conf[1,], conf[2,], col=col, lty=lty)
}

general_msr_cdf <- function(ecdfs, density, msr)
{
  for( i in 1:length(ecdfs))
  {
    if(density>= ecdfs[[i]]$prop & density<= ecdfs[[i+1]]$prop)
    {
      dist1 = density - ecdfs[[i]]$prop
      dist2 = ecdfs[[i+1]]$prop - density
      w = dist1/(dist1+dist2)
      return(ecdfs[[i]]$cdf(msr)*(1-w)+(w)*ecdfs[[i+1]]$cdf(msr))
    }
  }
}

extract_ecdf_function <- function(ecdfs, density)
{
  d = msr_confidence_line(ecdfs, density)
  d[2] <- 0
  d[length(d)] <- 0
  approxfun(x = d[1,], y = d[2,])
}

general_msr_cdf <- function(ecdfs, density, msr)
{
  
  for( i in 1:length(ecdfs))
  {
    if(density>= ecdfs[[i]]$prop & density<= ecdfs[[i+1]]$prop)
    {
      dist1 = density - ecdfs[[i]]$prop
      dist2 = ecdfs[[i+1]]$prop - density
      w = dist1/(dist1+dist2)
      return(ecdfs[[i]]$cdf(msr)*(1-w)+(w)*ecdfs[[i+1]]$cdf(msr))
    }
  }
}



#markovchain_fake_data <- function(size, transition_matrix)
#{
#  model <- new("markovchain", states = c("0","1"), transitionMatrix = transition_matrix)
#  as.integer(rmarkovchain(size, object = meth_model, t0 = "0"))
#}


prop_msr_samples <- function(l, props, sample_size, cores = 1, verbose = F)
{
  start_time <- Sys.time()
  msr_samples = mcmapply(1:sample_size, mc.preschedule = T, mc.cores = cores, FUN =  function(n)
  {
    v = rbinom(n = l, size = 1, prob = probs)
    rr_v <- calculate_relevance_resolution_vector(methylation_vector = v, verbose = F, na_tolerance = 0.1)
    MSR_area(rr_v)
  })
  if(verbose) cat("time: ", Sys.time()-start_time)
  
  return(msr_samples)
}

plot_positions <- function(positions)
{ 
  positions <- positions-min(positions)
  plot(x = c(0, 1), y = c(0,0), col = alpha(1, 0.4), ylab = "", xlab = "Positions", type = "l", yaxt='n')
  points(x = positions/max(positions), y=array(data = 0, dim = length(positions)), pch = "|", ylab = "")
}

MSR_example <- function(positions, perturb = 0)
{
  if(perturb)
  {
    s <- sort(positions)
    sd <- mean(s[2:length(s)]-s[1:(length(s)-1)])*perturb
    positions <- s + rnorm(length(s), sd = sd)
  }
  
  rr <- positions_MSR(positions, msr = F)
  msr <- MSR_area(rr)
  plot_positions(positions - min(positions))
  resolution_relevance_plot(rr)
  title(paste("MSR:", round(msr, digits = 3)))
}

MSR_visualization <- function(positions, discretization_bin_size = NA, perturb = 0)
{
  if(perturb)
  {
    s <- sort(positions)
    sd <- mean(s[2:length(s)]-s[1:(length(s)-1)])*perturb
    positions <- s + rnorm(length(s), sd = sd)
  }
    
  rr <- positions_MSR(positions, discretization_bin_size, F, F, max_bins = 1e6)
  msr <- MSR_area(rr)
  resolution_relevance_plot(rr)
  title(paste("MSR:", round(msr, digits = 4)))
}

rr_curve_comparison <- function(v_list, legend_names, add_msr = T, discrete = F)
{
  rr_list <- lapply(v_list, function(x)
    {
    positions_MSR(x, 1, F, F, max_bins = 1e6, discrete = discrete)
  })
  
  
  msr_list <- sapply(rr_list, MSR_area)
  
  if(add_msr)
    legend_names <- as.character(sapply(1:length(legend_names), function(i){paste(legend_names[i], round(msr_list[i],3), sep = ": ")}))
  
  ord <- order(msr_list, decreasing = T)
  
  plotter = plot
  for(i in ord)
  {
    resolution_relevance_plot(rr_list[[i]], alpha(i+1, 0.5), plotter, area_alpha = 0.1)
    plotter = points
  }
  
  legend("topleft", legend=legend_names[ord], col=ord+1, lty = 1, cex = 0.8, y.intersp = 0.5)
}

filter_by_range <- function(v, min, max)
{
  v <- v[v>=min]
  v[v<=max]
}

basis_pattern_rr_curve_comparison <- function(patterns, chromosome, min, max, Genome = BSgenome.Hsapiens.UCSC.hg38, same_M = F, discrete = F)
{
  v_list <- lapply(patterns, function(p) {
    x <- nucleotides_pattern_positions(chromosome, p, Genome)
    if(same_M)
      return(x[min:max])
    filter_by_range(x, min, max)
  })
  
  rr_curve_comparison(v_list, legend_names = patterns, discrete = discrete)
  title(paste(chromosome, "from basis", min, "to", max))
}

#setwd("./Scrivania/Tesi/MethylationCode/")
#directory <- "MethylationData/binary_rate/"

#cell_files <- list.files(directory,pattern="(_converted.Rda)$")
#cell_names <- sub("_converted.Rda","", cell_files)
#methylation_vector <- as.logical(readRDS("MethylationData/binary_rate/GSM3436261_O1_TA_Hi_10_converted.Rda"))


make_random_msr_data_frame <- function(length, samples, lw = 0, hp = 1)
{
  p <- runif(samples, lw, hp)
  msr <- array(dim = samples)
  for(i in 1:samples)
  {
    v <- rbinom(length, 1, p[i])
    msr[i] <- calculate_MSR_area(v)
    p[i] <- mean(v)
  }

  return(data.frame(p=p,msr=msr))
}

make_random_msr_data_frame2 <- function(samples, lw = 100, hp = 10000, verbose = F)
{
  M <- round(runif(samples, lw, hp))
  msr <- array(dim = samples)
  for(i in 1:samples)
  {
    v <- runif(M[i])
    msr[i] <- positions_MSR(v, max_gbins = 20)
    if(verbose) cat(i, "")
  }
  
  return(data.frame(M=M,msr=msr))
}
