setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
source("WGBS_analysis_functions.R", chdir = T)

reconstruction_test <- function(v, missing_data_prop, n, na_tolerance = 0, na_values_handler = replace_nas_hybrid)
{
  true_msr = MSR_area(calculate_relevance_resolution_vector(v, verbose = F))
  samples <- sapply(1:n, function(i)
  {
    MSR_area(calculate_relevance_resolution_vector(mutilated(v, missing_data_prop), na_tolerance = na_tolerance, verbose = F, na_values_handler = na_values_handler))
  })
  
  hist(samples, xlim = c(min(samples-0.02), max(samples+0.02)), breaks = linspace(min(samples), max(samples), sqrt(n/1.5)))
  lines(x = c(true_msr, true_msr), y=c(0,1e4), col = "red", lwd = 5)
  
  mn =  mean(samples, na.rm=T)
  sdev = sqrt(var(samples, na.rm=T))
  bias = abs(true_msr-mn)
  cat("true: ", true_msr)
  #cat("  error: ", sqrt(bias^2+sdev^2))
  cat("\n mean: ", mn)
  #cat("  std: ", sdev)
  #cat("  bias: ", bias)
  cat("\nmean absolute error: ", mean(abs(samples-true_msr)))
  
}


make_reconstruction_df <- function(samples, vector_generator, missing_data_prop, na_tolerance, na_values_handler, verbose = T)
{
  msr <- array(dim=samples)
  rec_msr <- msr
  p <- msr
  
  for(i in 1:samples)
  {
    if(verbose)
      cat(i, " ")
    v <- vector_generator()
    msr[i] <- calculate_MSR_area(v)    
    rec_msr[i] <- MSR_area(calculate_relevance_resolution_vector(mutilated(v, missing_data_prop), na_tolerance = na_tolerance, verbose = F, na_values_handler = na_values_handler))
    p[i] <- mean(v)
  }
  
  par(pty="s")
  lim = c(min(c(msr,rec_msr)), max(c(msr,rec_msr)))
  plot(msr, rec_msr, xlim=lim, ylim=lim)
  abline(coef = c(0,1))
  cat("\ncorrelation:",(cor(msr,rec_msr)))
  cat("\nmean absolute error:", mean(abs(msr-rec_msr)))
  cat("\n rmse:", sqrt(mean((msr-rec_msr)^2)))
  return(data.frame(msr,rec_msr,p))
}

make_sensibiity_df <- function(samples, vector_generator, modified_prop, verbose = T)
{
  msr <- array(dim=samples)
  mod_msr <- msr
  p <- msr
  
  for(i in 1:samples)
  {
    if(verbose)
      cat(i, " ")
    v <- vector_generator()
    msr[i] <- calculate_MSR_area(v)    
    mod_msr[i] <- calculate_MSR_area(xor(v, rbinom(length(v), 1, modified_prop)))
    p[i] <- mean(v)
  }
  
  par(pty="s")
  plot(msr, mod_msr, xlim = c(0, 0.3), ylim = c())
  abline(coef = c(0,1))
  cat("\ncorrelation:",(cor(msr,mod_msr)))
  return(data.frame(msr,mod_msr,p))
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


reconstruction_test(v2, 0.1, 500, 1, replace_nas_hybrid_stochastic)
reconstruction_test(v2, 0.1, 200, 1, replace_nas_with_bin_prop)

reconstruction_test(v, 0.1, 500, 1, replace_nas_hybrid_stochastic)
reconstruction_test(v, 0.1, 500, 1, replace_nas_with_bin_prop)

# bernoulli 1e3 low p
l <- 1e3
v <- rbinom(l,1,0.1)
reconstruction_test(v, 0.1, 300, 0.2, replace_nas_with_bin_prop)
reconstruction_test(v, 0.1, 300, 0.2, replace_nas_with_bin_prop_deterministic)

# bernoulli 1e3 medium p
l <- 1e3
v <- rbinom(l,1,0.5)
reconstruction_test(v, 0.1, 500, 0.33, replace_nas_with_bin_prop)
reconstruction_test(v, 0.1, 300, 0.2, replace_nas_with_bin_prop_deterministic)

# bernoulli 1e3 high p
l <- 1e3
v <- rbinom(l,1,0.9)
reconstruction_test(v, 0.1, 100, 0.2, replace_nas_with_bin_prop)
reconstruction_test(v, 0.1, 300, 0.2, replace_nas_with_bin_prop_deterministic)

# sinusoidal
l <- 1e3
p <- ((sin(linspace(0,2*pi, l)))^2)/2
v <- rbinom(l,1,p)
cat(mean(p))
reconstruction_test(v, 0.1, 300, 0.2, replace_nas_with_bin_prop)
df <- make_reconstruction_df(300, function(){rbinom(l,1,p)}, 0.1, 0.2, replace_nas_with_bin_prop)
reconstruction_test(v, 0.1, 300, 0.2, replace_nas_with_bin_prop_deterministic)

# linear
l <- 1e3
p <- linspace(0,1, l)
v <- rbinom(l,1,p)
cat(mean(p))
reconstruction_test(v, 0.1, 300, 0.2, replace_nas_with_bin_prop)
df <- make_reconstruction_df(300, function(){rbinom(l,1,p)}, 0.1, 0.2, replace_nas_with_bin_prop)
reconstruction_test(v, 0.1, 300, 0.2, replace_nas_with_bin_prop_deterministic)


######
l <- 1e3
f <- function() {rbinom(l,1,0.1)}
df <- make_reconstruction_df(300, f, 0.1, 0.2, replace_nas_with_bin_prop)

######
l <- 1e3
f <- function() {rbinom(l,1,0.5)}
df <- make_reconstruction_df(200, f, 0.1, 0.2, replace_nas_with_bin_prop)

######
l <- 1e3
f <- function() {rbinom(l,1,0.9)}
df <- make_reconstruction_df(200, f, 0.1, 0.2, replace_nas_with_bin_prop)



######
l <- 1e3
f <- function() {rbinom(l,1,0.1)}
df <- make_reconstruction_df(500, f, 0.1, 0.2, replace_nas_with_bin_prop)

######
l <- 1e3
f <- function() {rbinom(l,1,linspace(0,1, l))}
f <- function() {
  if(rbinom(1,1,0.5))
  return(rbinom(l,1,linspace(0,1, l)))
  
    rbinom(l,1,0.5) }

df <- make_reconstruction_df(200, f, 0.1, 0.2, replace_nas_with_bin_prop)
df <- make_sensibiity_df(100,f,0.1)
######


