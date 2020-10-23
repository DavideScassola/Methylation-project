setwd(dir = "Scrivania/Tesi/Methylation-project/Code/")
source("WGBS_analysis_functions.R", chdir = T)

# single_correlations_matrix(train_datas, response, predictors)
# lasso_normalized_coefficients_matrix(train_datas, response, predictors, lambda)
# generic_model(train_datas, response, predictors) -> 
# self_test_r2(models, test_datas)
# cross_matrix_test_r2(models, test_datas)
# mixed_test_r2(model, test_datas)

names <- c("H1", "K562", "GM12878", "GM23248", "Hela", "endodermal", "lung", "stomach")

# general paramenters
train_prop <- 0.6
random_train_test_split <- T
epsilon <- 1e-3

# load data
# prepare data
  # join expression data
  # add log features
  # remove sex chromosomes

response_variable <- "log_total_pme_TPM"

# predictors selection
basic_predictors <- c("log_nucleotides", "CpG_density", "meth_rate")
advanced_predictors <- c("meth_autocorrelation", "drift", "meth_sd")
msr_related_predictors <- c("msr", "inverted_msr", "ecdf", "inverted_ecdf", "residual", "inverted_residual")
basic_and_advanced_predictors <- c(basic_predictors, advanced_predictors)
all_predictors <- c(basic_predictors, advanced_predictors, msr_related_predictors)


# show single_correlations_matrix(train_datas, response, predictors)

# show lasso_normalized_coefficients_matrix(train_datas, response, predictors, lambda)
# for basic_and_advanced_predictors, all_predictors


# LINEAR MODEL BASIC_PREDICTORS
# self_test_r2(models, test_datas)
# cross_matrix_test_r2(models, test_datas)
# mixed_test_r2(model, test_datas)

# LINEAR MODEL BASIC AND ADVANCED PREDICTORS
# self_test_r2(models, test_datas)
# cross_matrix_test_r2(models, test_datas)
# mixed_test_r2(model, test_datas)

# LINEAR MODEL ALL PREDICTORS
# self_test_r2(models, test_datas)
# cross_matrix_test_r2(models, test_datas)
# mixed_test_r2(model, test_datas)

# LASSO MODEL BASIC AND ADVANCED PREDICTORS
# self_test_r2(models, test_datas)
# cross_matrix_test_r2(models, test_datas)
# mixed_test_r2(model, test_datas)

# LASSO MODEL ALL PREDICTORS
# self_test_r2(models, test_datas)
# cross_matrix_test_r2(models, test_datas)
# mixed_test_r2(model, test_datas)

# RANDOM FOREST MODEL BASIC AND ADVANCED PREDICTORS
# self_test_r2(models, test_datas)
# cross_matrix_test_r2(models, test_datas)
# mixed_test_r2(model, test_datas)

# RANDOM FOREST MODEL ALL PREDICTORS
# self_test_r2(models, test_datas)
# cross_matrix_test_r2(models, test_datas)
# mixed_test_r2(model, test_datas)
