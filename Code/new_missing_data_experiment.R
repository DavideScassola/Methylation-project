setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
source("WGBS_analysis_functions.R", chdir = T)

reconstruction_test <- function(v, missing_data_prop, n, na_tolerance = 0, na_values_handler = replace_nas_hybrid)
{
  true_msr = MSR_area(calculate_relevance_resolution_vector(v, verbose = F))
  samples <- sapply(1:n, function(i)
  {
    MSR_area(calculate_relevance_resolution_vector(mutilated(v, missing_data_prop), na_tolerance = na_tolerance, verbose = F, na_values_handler = na_values_handler))
  })
  
  hist(samples, xlim = c(min(samples-0.02), max(samples+0.02)))
  lines(x = c(true_msr, true_msr), y=c(0,1e4), col = "red")
  
  mn =  mean(samples, na.rm=T)
  sdev = sqrt(var(samples, na.rm=T))
  bias = abs(true_msr-mn)
  cat("true: ", true_msr)
  cat("  error: ", sqrt(bias^2+sdev^2))
  cat("\n mean: ", mn)
  cat("  std: ", sdev)
  cat("  bias: ", bias)
  
}

l <- 1e5
v = rbinom(1e3, 1, 0.05)
v = get_max_information_vector(1e3)
v <- rbinom(l, 1, (1:l)/l)

v <- rbinom(1e3, 1, 0.1)
v2 <- v
v2[4:30] <- 1
v2[400:500] <- 1
v2[600:800] <- 0


reconstruction_test(v2, 0.06, 100, 1, replace_nas_hybrid_stochastic)
reconstruction_test(v2, 0.1, 100, 1, replace_nas_with_bin_prop)

reconstruction_test(v, 0.1, 200, 1, replace_nas_hybrid_stochastic)
reconstruction_test(v, 0.1, 200, 1, replace_nas_with_bin_prop)
