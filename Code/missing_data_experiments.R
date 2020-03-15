source("MSR_analysis_functions.R", chdir = T)
setwd("./Scrivania/Tesi/MethylationCode/")
directory <- "MethylationData/binary_rate/"

cell_files <- list.files(directory,pattern="(_converted.Rda)$")
cell_names <- sub("_converted.Rda","", cell_files)
methylation_vector <- as.logical(readRDS("MethylationData/binary_rate/GSM3436261_O1_TA_Hi_10_converted.Rda"))
l <- length(methylation_vector)

meth_rr_vectors <- lapply(cell_names, function(x)
  {
  cat("reading", x, "\n")
  methylation_vector <- as.logical(readRDS(sprintf("%s%s_converted.Rda", directory, x)))
  rr_vec <- calculate_relevance_resolution_vector(methylation_vector, na_tolerance = 0.05)
  remove(methylation_vector)
  gc()
  return(rr_vec)
})


# fake methylation vectors
l = 120000
f1 <- get_max_information_vector(l)
f2 <- rbinom(l, 1, (1:l)/l) # the mountain
f3 <- rbinom(l, 1, 0.5)
f4 <- rbinom(l, 1, (1:50)/200)

fake_methylation_vectors <- list(f1, !f1, f2, f3, f4, !f4)

info = "Picking only bins without NAs"
info = "Vertical lines are proportion of available data actually used, picking only bins without NAs"
missing_data_experiment(f1, missing_data_rates = c(0, 0.1, 0.1, 0.1, 0.1), na_tolerance = 0, na_values_handler = replace_nas_hybrid, info = info)
missing_data_experiment(!f1, missing_data_rates = c(0, 0.1, 0.1, 0.1, 0.1), na_tolerance = 0, na_values_handler = replace_nas_hybrid, info = info)

info = "Picking only bins with at maximum 10% of NAs"
missing_data_experiment(f1, missing_data_rates = c(0, 0.1, 0.1, 0.1, 0.1), na_tolerance = 0.1, na_values_handler = replace_nas_hybrid_stochastic, info = info)
missing_data_experiment(!f1, missing_data_rates = c(0, 0.1, 0.1, 0.1, 0.1), na_tolerance = 0.1, na_values_handler = replace_nas_hybrid_stochastic, info = info)
missing_data_experiment(f1, missing_data_rates = c(0, 0.2, 0.5, 0.9), na_tolerance = 0.1, na_values_handler = replace_nas_hybrid_stochastic, info = info)

info = "Picking only bins with at maximum 50% of NAs"
missing_data_experiment(f1, missing_data_rates = c(0, 0.1, 0.5, 0.7, 0.9), na_tolerance = 0.5, na_values_handler = replace_nas_hybrid_stochastic, info = info)

info = "Picking every possible bin and replacing all NAs"
missing_data_experiment(f1, missing_data_rates = c(0, 0.7, 0.9), na_tolerance = 1, na_values_handler = replace_nas_hybrid_stochastic, info = info)
missing_data_experiment(f2, missing_data_rates = c(0, 0.7, 0.9), na_tolerance = 1, na_values_handler = replace_nas_hybrid_stochastic, info = info)


info = "Picking every possible bin, but replacing NAs with zeros"
missing_data_experiment( f1, missing_data_rates = c(0, 0.2, 0.5, 0.9), na_tolerance = 1, na_values_handler = replace_nas_with_zeros, info = info)
missing_data_experiment(!f1, missing_data_rates = c(0, 0.2, 0.5, 0.9), na_tolerance = 1, na_values_handler = replace_nas_with_zeros, info = info)


# Now experiments with meth vectors
methylation_vector1 <- as.logical(readRDS("MethylationData/binary_rate/GSM3436261_O1_TA_Hi_10_converted.Rda"))
methylation_vector2 <- as.logical(readRDS("MethylationData/binary_rate/GSM3436264_O1_TA_Hi_13_converted.Rda"))
methylation_vector3 <- as.logical(readRDS("MethylationData/binary_rate/GSM3436271_O1_TA_Hi_1_converted.Rda"))
methylation_vector4 <- as.logical(readRDS("MethylationData/binary_rate/GSM3436313_O8_TA_Hi_19_converted.Rda"))
methylation_vector5 <- as.logical(readRDS("MethylationData/binary_rate/GSM3436360_Y2_TA_Hi_1_converted.Rda"))

mv <- list(methylation_vector1, methylation_vector2, methylation_vector3, methylation_vector4, methylation_vector5)
l <- length(methylation_vector1)
fm1 <- get_max_information_vector(l)
fm1_masked <- lapply(mv, function(v) {apply_same_missing_data_pattern(fm1,v)})

rr_fm1 <- calculate_relevance_resolution_vector(fm1, na_tolerance = 0, na_values_handler = replace_nas_hybrid)


# picking only bins without NAs
rr_fm1_masked <- lapply(fm1_masked, function(v) calculate_relevance_resolution_vector(v,na_tolerance = 0, na_values_handler = replace_nas_hybrid))
rr_plots(list(rr_fm1, rr_fm1_masked[[1]], rr_fm1_masked[[2]], rr_fm1_masked[[3]], rr_fm1_masked[[4]]), threshold = 0.05, title = "Applying cells missing data pattern, picking only bins without NAs", legend_names = c("original data","cell1 pattern", "cell2 pattern", "cell3 pattern", "cell4 pattern"))

# picking every possible bin
rr_fm1_masked_tolerance <- lapply(fm1_masked, function(v) calculate_relevance_resolution_vector(v,na_tolerance = 1, na_values_handler = replace_nas_hybrid_stochastic))
rr_plots(list(rr_fm1, rr_fm1_masked_tolerance[[1]], rr_fm1_masked_tolerance[[2]], rr_fm1_masked_tolerance[[3]], rr_fm1_masked_tolerance[[4]]), title = "Applying cells missing data pattern, replacing every NAs", legend_names = c("original data","cell1 pattern", "cell2 pattern", "cell3 pattern", "cell4 pattern"))

# picking 0.8 bins
rr_fm1_masked_tolerance08 <- lapply(fm1_masked, function(v) calculate_relevance_resolution_vector(v,na_tolerance = 0.8, na_values_handler = replace_nas_hybrid_stochastic))
rr_plots(list(rr_fm1, rr_fm1_masked_tolerance08[[1]], rr_fm1_masked_tolerance08[[2]], rr_fm1_masked_tolerance08[[3]], rr_fm1_masked_tolerance08[[4]]), title = "Applying cells missing data pattern, picking only bins with at maximum 80% of NAs", legend_names = c("original data","cell1 pattern", "cell2 pattern", "cell3 pattern", "cell4 pattern"))


# analysing our data
rr_mv_0t <- lapply(mv, function(v) {calculate_relevance_resolution_vector(v,na_tolerance = 0, na_values_handler = replace_nas_hybrid)})
rr_mv_0.8t <- lapply(mv, function(v) {calculate_relevance_resolution_vector(v,na_tolerance = 0.8, na_values_handler = replace_nas_hybrid)})
rr_mv_0.2t <- lapply(mv, function(v) {calculate_relevance_resolution_vector(v,na_tolerance = 0.2, na_values_handler = replace_nas_hybrid)})
legend_names <- c("cell1", "cell2", "cell3", "cell4", "cell5")
rr_plots(rr_mv_0t, title = "Cell data, picking only bins without NAs", legend_names = legend_names)
rr_plots(rr_mv_0.2t, title = "Cell data, picking only bins with at maximum 20% of NAs", legend_names = legend_names)
rr_plots(rr_mv_0.8t, title = "Cell data, picking only bins with at maximum 80% of NAs", legend_names = legend_names)

rr_plots(list(rr_mv_0t[[1]], rr_mv_0.2t[[1]],rr_mv_0.8t[[1]]), title = "Cell1 data", legend_names = c("no replacing", "20% replacing", "80% replacing"))

#################################################################################
# is our data different from a repetition of bernoulli extraction?
simulated <- rbinom(length(methylation_vector1), 1, meth_prop(methylation_vector1))
simulated_masked <- apply_same_missing_data_pattern(simulated, methylation_vector1)
exp_vectors <- list(simulated, simulated_masked, methylation_vector1)
rr_experiment0 <- lapply(exp_vectors, function(v) {calculate_relevance_resolution_vector(v,na_tolerance = 0, na_values_handler = replace_nas_hybrid_stochastic)})
rr_experiment02 <- lapply(exp_vectors, function(v) {calculate_relevance_resolution_vector(v,na_tolerance = 0.2, na_values_handler = replace_nas_hybrid_stochastic)})
rr_experiment08 <- lapply(exp_vectors, function(v) {calculate_relevance_resolution_vector(v,na_tolerance = 0.8, na_values_handler = replace_nas_hybrid_stochastic)})

legend_names <- c("binomial data", "binomial data with cell miss pattern", "cell data")
rr_plots(rr_experiment0, title = "Binomial data with cell missing data pattern, no NAs", legend_names = legend_names)
rr_plots(rr_experiment02, title = "Binomial data with cell missing data pattern, 20% accepted NAs", legend_names = legend_names)
rr_plots(rr_experiment08, title = "Binomial data with cell missing data pattern, 80% accepted NAs", legend_names = legend_names)
#################################################################################
