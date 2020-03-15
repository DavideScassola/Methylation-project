suppressMessages(library(data.table))
suppressMessages(library(hashmap))
suppressMessages(library(doParallel))
#suppressMessages(library(argparse))

# returns a dataframe with the name of the file
cov.gz_to_data <- function(folder_directory, sample_name, input_format, rate = T, position_refers_to_C = F) {
  
  # Load data
  data <- fread(cmd=sprintf("zcat < %s/%s.cov.gz",folder_directory,sample_name), verbose=F, showProgress=F)
  
  # Selecting relevant columns
  # Input format 1 (chr,pos, pos2, rate, m, nm) 
  if (input_format == 1) {
    colnames(data) <- c("chr","Cpos", "Cpos2", "rate", "m", "nm")
    if(rate){ 
      data = subset(data, select = c(chr,Cpos,rate) )
      data$rate <- as.integer(data$rate)}
  
    
    else {
      data = subset(data, select = c(chr,Cpos, m, nm) )
      data$m = as.integer(data$m)
      data$nm = as.integer(data$nm)
    }
    
  }
  
  # Selecting relevant columns
  # Input format 2 (chr,pos, rate) 
  if (input_format == 2) {
    colnames(data) <- c("chr","Cpos","rate")
    data = subset(data, select = c(chr,Cpos,rate) )
  }
  
  if(!position_refers_to_C) {
    # Transforming actual Gpos into Cpos
    data$Cpos <- as.integer(data$Cpos - 1)
    #data$rate <- as.integer(data$rate)
  }
  
  return((data))
}

#cov.gz_to_data <- Vectorize(cov.gz_to_data, USE.NAMES = F)


convert_to_met_rate_vector <- function(data, verbose = T, tolerance = 1) {
  
  l <- CpG_kit$length
  bmv <- array(dim = l)
  actual_data_length <- length(data$Cpos)
  not_found_sites <- 0
  
  indexes <- CpG_kit$tolerant_genetic_pos_to_index(data$chr, data$Cpos, tolerance = tolerance)
  wich_are_not_NA <- !is.na(indexes)
  valid_indexes <- indexes[wich_are_not_NA]
  valid_positions <- (1:actual_data_length)[wich_are_not_NA]
  bmv[valid_indexes] <- data[valid_positions]$rate
  
  not_found_sites <- length(indexes) - length(valid_indexes)

  if(not_found_sites!=0) {
    print(sprintf("%d sites were not found", not_found_sites))
  }
  
  return(bmv)
  
}

convert_to_met_counts_vector <- function(data, verbose = T, tolerance = 1) {
  
  l <- CpG_kit$length
  bmv <- list(m = array(dim = l), nm = array(dim = l))
  actual_data_length <- length(counts_data$Cpos)
  not_found_sites <- 0
  
  indexes <- CpG_kit$tolerant_genetic_pos_to_index(data$chr, data$Cpos, tolerance = tolerance)
  wich_are_not_NA <- !is.na(indexes)
  valid_indexes <- indexes[wich_are_not_NA]
  valid_positions <- (1:actual_data_length)[wich_are_not_NA]
  bmv$m[valid_indexes] <- data[valid_positions]$m
  bmv$nm[valid_indexes] <- data[valid_positions]$nm
  
  not_found_sites <- length(indexes) - length(valid_indexes)
  
  if(not_found_sites!=0) {
    print(sprintf("%d sites were not found", not_found_sites))
  }
  
  return(bmv)
  
}

#convert_to_met_rate_vector <- Vectorize(convert_to_met_rate_vector)

get_met_rate_matrix <- function(folder_directory, format = 1, verbose = T) {
  samples <- sub(".cov.gz","",list.files(indir, pattern="(.cov.gz)$"))
  dataset_list <- cov.gz_to_data(folder_directory, samples, format)  
  rate_vector_list <- convert_to_met_rate_vector(dataset_list)
  remove(dataset_list)
  gc(full = T, verbose = F)
  return(rate_vector_list)
}


get_CpG_kit <- function() {
  
  print("loading CpG sites file . . .")
  load(file="CpG_sites_dataframe.RData")
  
  print("loading CpG_hashmap.keys file . . .")
  load(file="CpG_hashmap.keys.Rdata")
  
  l <- length(CpG_hashmap.keys)
  values <- 1:l
  
  print("builiding map . . .")
  CpG_hashmap <- hashmap(CpG_hashmap.keys, values)
  
  chromosome_to_int <- function(chromosome) {
    l <- nchar(chromosome)
    chr_n <- chromosome
    if(l>2) {
      chr_n <- substr(chromosome,4,l)
    }
    
    if(chr_n=="X") return(20)
    if(chr_n=="Y") return(21)
    else return(strtoi(chr_n))
  }
  
  chromosome_to_int <- Vectorize(chromosome_to_int, USE.NAMES = F)
  
  genetic_pos_to_index <- function(chromosome, position) {
    chr_n <- chromosome
    if(typeof(chromosome)=="character") chr_n <- chromosome_to_int(chr_n)
    return(CpG_hashmap[[position + chr_n*1000000000]])
  }
  
  genetic_pos_to_index <- Vectorize(genetic_pos_to_index, USE.NAMES = F)
  
  index_to_genetic_pos <- function(index) {
    return(list(sites$chr_name[index], sites$position[index]))
  }
  
  index_to_genetic_pos <- Vectorize(index_to_genetic_pos, USE.NAMES = F)
  
  tolerant_genetic_pos_to_index  <- function(chromosome, position, tolerance = 1)
  {
    index <- genetic_pos_to_index(chromosome, position)
    if(!is.na(index)) return(index)
    
    for(i in 1:tolerance) {
      index <- genetic_pos_to_index(chromosome, position+i)
      if(!is.na(index)) return(index)
      
      index <- genetic_pos_to_index(chromosome, position-i)
      if(!is.na(index)) return(index)
    }
    
    return(NA)
  }
  
  tolerant_genetic_pos_to_index <- Vectorize(tolerant_genetic_pos_to_index, USE.NAMES = F)
  
  CpG_kit <- list(tolerant_genetic_pos_to_index = tolerant_genetic_pos_to_index, genetic_pos_to_index = genetic_pos_to_index, keys = CpG_hashmap.keys, sites = sites, index_to_genetic_pos = index_to_genetic_pos, chromosome_to_int = chromosome_to_int, values = values, CpG_hashmap = CpG_hashmap, length = l)
  
  print("done")
  return(CpG_kit)
  
}

cov.gz_to_met_rate_vector <- function(folder_directory, sample_name, input_format = 2, tolerance=1, position_refers_to_C = F, verbose = T)
{
  data <- cov.gz_to_data(folder_directory, sample_name, input_format, position_refers_to_C)
  print(sample_name)
  return(convert_to_met_rate_vector(data, tolerance = tolerance))
}

cov.gz_to_met_counts_vector <- function(folder_directory, sample_name, tolerance=1, position_refers_to_C = F, verbose = T)
{
  data <- cov.gz_to_data(folder_directory, sample_name, input_format = 1, position_refers_to_C, rate = F)
  print(sample_name)
  return(convert_to_met_counts_vector(data, tolerance = tolerance))
}

convert_cov.gz <- function(directory, input_format = 2, tolerance = 1, position_refers_to_C = F) {
  samples <- sub(".cov.gz","",list.files(directory,pattern="(.cov.gz)$"))
  lapply(samples, function(x) {
    rv <- cov.gz_to_met_rate_vector(directory, x, input_format, tolerance, position_refers_to_C)
    saveRDS(rv, file = sprintf("%s_converted.Rda", x))
  } )
}

convert_counts_cov.gz <- function(directory, tolerance = 1, position_refers_to_C = F) {
  samples <- sub(".cov.gz","",list.files(directory,pattern="(.cov.gz)$"))
  lapply(samples, function(x) {
    cv <- cov.gz_to_met_counts_vector(directory, x, tolerance, position_refers_to_C)
    saveRDS(cv, file = sprintf("directory/met_counts/%s_count_converted.Rda", x))
  } )
}


cov.gz_to_met_rate_vector_vectorized <- function(folder_directory, samples, input_format = 2, tolerance=1, position_refers_to_C = F)
{
  return(lapply(samples, function(x) cov.gz_to_met_rate_vector(folder_directory, x, input_format, tolerance, position_refers_to_C)))
}

#save(get_CpG_kit, file = "get_CpG_kit.Rdata")

build_met_matrix_from_rda <- function(directory, samples_names = NA) {
  if(sum(is.na(samples_names))>0) {
    samples_names <- sub("_converted.Rda","",list.files(indir,pattern="(_converted.Rda)$"))
  }
  
  l <- 21867550
  n <- length(samples_names)
  x = array(dim = c(n,l))
  
  for(i in 1:n) {
    sample_name = samples_names[i]
    x[i,] <- readRDS(file = sprintf("%s/%s_converted.Rda", directory, samples_names[i]))
    gc(verbose = F)
  }
  
  #x = sapply(1:n , function(i) {
  #  sample_name = samples_names[i]
  #  return(readRDS(file = sprintf("%s/%s_converted.Rda", directory, sample_name)))
  #})
  
  #gc(verbose = F)
  #return(t(x))
  return(x)
  
}

########################################################################

#l <- length(sites$position)
#values <- 1:l
#chr_numbers <- chromosome_to_int(as.character(sites$chr_name))
#CpG_hashmap.keys <- sites$position + chr_numbers*1000000000
#remove(chr_numbers)
#CpG_hashmap <- hashmap(CpG_hashmap.keys, values)


setwd("./Scrivania/Tesi/MethylationCode/")
#load(file="CpG_sites_dataframe.RData")

indir <- "./MethylationData"
cores <- 2
input_format <- 2 #(chr,pos, pos2, rate, m, nm) 
position_refers_to_C <- FALSE # the position of the CpG sites in the files refers to the C position
sample_index <- 1
samples <- sub(".cov.gz","",list.files(indir,pattern="(.cov.gz)$"))
#sample_name <- samples[sample_index]

# Parallelise processing
registerDoParallel(cores=3)

CpG_kit <- get_CpG_kit()
convert_cov.gz(indir, input_format = 1, tolerance = 1, position_refers_to_C = F)
#convert_counts_cov.gz(indir, tolerance = 2)

#data<- cov.gz_to_data(indir, samples[1], 2, position_refers_to_C)
#datalist <- cov.gz_to_data(indir, samples, 2, position_refers_to_C)
#matrix = cov.gz_to_met_rate_vector_vectorized(folder_directory = indir, samples = samples)
#bmw <- convert_to_met_rate_vector(data)
#get_met_rate_matrix(indir, format = 2)
#l <- CpG_kit$length

#samples_names <- sub("_converted.Rda","",list.files(indir,pattern="(_converted.Rda)$"))

#x = build_met_matrix_from_rda(indir, samples_names[1:5])
