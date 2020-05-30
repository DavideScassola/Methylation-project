MSR and expression, for stomach cells and window size of 10,000
================

Here I will investigate if there is a relationship between the presence of genes and genes expression in a certain genomic region and the MSR (with some of its derivate statistics).

I chose stomach cells data in order to do this, and CpG windows of size 10,000, that corresponds to a variable window size in term of nucleotides (on average about 1,000,000).

This is an example of total-rna-seq file, that shows for each "gene" its transcripts and some measures of expression. In this case I just kept two colums. The first one indicates the "gene"", the second one is the Transcript Per Million that is a relative measure of how much a gene is expressed.

    ##                    gene_id    TPM
    ##     1:     ENSG00000000003   0.49
    ##     2:     ENSG00000000005   0.00
    ##     3:     ENSG00000000419   0.24
    ##     4:     ENSG00000000457   0.28
    ##     5:     ENSG00000000460   2.18
    ##    ---                           
    ## 60818: gSpikein_ERCC-00165   8.98
    ## 60819: gSpikein_ERCC-00168   0.04
    ## 60820: gSpikein_ERCC-00170   0.56
    ## 60821: gSpikein_ERCC-00171 530.85
    ## 60822:    gSpikein_phiX174  98.55

This is the annotation file that store the position occupied by each human gene.

    ##         chr     start       end strand              id                    anno
    ##     1: chr1     65419     71585      + ENSG00000186092 genebody_protein_coding
    ##     2: chr1    450703    451697      - ENSG00000284733 genebody_protein_coding
    ##     3: chr1    685679    686673      - ENSG00000284662 genebody_protein_coding
    ##     4: chr1    923928    944581      + ENSG00000187634 genebody_protein_coding
    ##     5: chr1    944204    959309      - ENSG00000188976 genebody_protein_coding
    ##    ---                                                                        
    ## 19801: chrY  24763069  24813492      - ENSG00000187191 genebody_protein_coding
    ## 19802: chrY  24833843  24907040      + ENSG00000205916 genebody_protein_coding
    ## 19803: chrY  25030901  25062548      - ENSG00000185894 genebody_protein_coding
    ## 19804: chrY  25622162  25624902      + ENSG00000172288 genebody_protein_coding
    ## 19805: chrX 135309480 135309659      + ENSG00000283644 genebody_protein_coding

The number of genes is much less than the ones in the total-rna-seq file, since the first one also contains so called pseudogenes and other stuff.

So the final dataFrame is the following (excluding some columns for readability):

    ##   start_chr start_position end_position gene_count total_TPM meth rate
    ## 2      chr1         921648      1151546         10     12.65 0.5677261
    ## 3      chr1        1151546      1390325         17     32.56 0.6524979
    ## 4      chr1        1390325      1647121         13     16.74 0.6697962
    ## 5      chr1        1647121      1981023          8     10.67 0.7296424
    ## 6      chr1        1981023      2270621          4      5.49 0.7102574
    ## 7      chr1        2270621      2559502          7      9.64 0.6749857

The full scheme includes:

**nucleotides**: number of nucleotides in the window

**CpG density**: fraction of nucleotides that is a C of a CpG site (= 10000/nucleotides)

**meth rate**: ratio of methylated CpG sites

**gene\_count**: number of genes included (even partially) inside the interval

**total\_TPM**: sum of the TPMs of the genes in the interval

then the MSR and some related statistics: **msr**, **inverted msr**, **msr ecdf **, **inverted msr ecdf**, **residual** (residual of the linear regression between msr and meth rate), **inverted residual**.

First let's see if there are pairwise correlations between the features.

###### Basic features:

![](MSR_and_expression_stomach_1e4_files/figure-markdown_github/unnamed-chunk-7-1.png)

log(TPM) is considered only for fragments with at least a gene.

###### Comparison with simple MSR statistics:

![](MSR_and_expression_stomach_1e4_files/figure-markdown_github/unnamed-chunk-8-1.png)

###### Comparison with other MSR statistics:

![](MSR_and_expression_stomach_1e4_files/figure-markdown_github/unnamed-chunk-9-1.png)

![](MSR_and_expression_stomach_1e4_files/figure-markdown_github/unnamed-chunk-10-1.png)

The correlation between the total TPM with the standard deviation of the TPM is:

    ##      cor 
    ## 0.892561

![](MSR_and_expression_stomach_1e4_files/figure-markdown_github/unnamed-chunk-12-1.png)

![](MSR_and_expression_stomach_1e4_files/figure-markdown_github/unnamed-chunk-13-1.png)

#### Predicting gene number

![](MSR_and_expression_stomach_1e4_files/figure-markdown_github/unnamed-chunk-14-1.png)

Negative binomial regression for gene number with basic predictors (nucleotides, CpG\_density, meth rate):

    ## 
    ## Call:
    ## glm.nb(formula = model_data$gene_count ~ nucleotides + CpG_density + 
    ##     `meth rate`, data = model_data, init.theta = 2.360023795, 
    ##     link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.4647  -0.8359  -0.1575   0.4520   5.0829  
    ## 
    ## Coefficients:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  8.034e+00  2.745e-01  29.267   <2e-16 ***
    ## nucleotides -1.027e-06  7.245e-08 -14.170   <2e-16 ***
    ## CpG_density -4.230e+01  4.712e+00  -8.977   <2e-16 ***
    ## `meth rate` -6.212e+00  3.168e-01 -19.610   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.36) family taken to be 1)
    ## 
    ##     Null deviance: 4187.9  on 2883  degrees of freedom
    ## Residual deviance: 3286.1  on 2880  degrees of freedom
    ## AIC: 16268
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.3600 
    ##           Std. Err.:  0.0910 
    ## 
    ##  2 x log-likelihood:  -16257.6530

Negative binomial regression with inverted\_msr as predictor

    ## 
    ## Call:
    ## glm.nb(formula = model_data$gene_count ~ model_data$inverted_msr, 
    ##     data = model_data, init.theta = 2.753414566, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.8162  -0.8493  -0.2670   0.3648   7.0595  
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)               7.4171     0.1709   43.39   <2e-16 ***
    ## model_data$inverted_msr -22.0763     0.6825  -32.34   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.7534) family taken to be 1)
    ## 
    ##     Null deviance: 4090.9  on 2766  degrees of freedom
    ## Residual deviance: 3078.7  on 2765  degrees of freedom
    ## AIC: 15533
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.753 
    ##           Std. Err.:  0.110 
    ## 
    ##  2 x log-likelihood:  -15526.564

Negative Binomial Regression with several predictors

    ## 
    ## Call:
    ## glm.nb(formula = gene_count ~ ., data = model_data, init.theta = 4.191147656, 
    ##     link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -4.3272  -0.7845  -0.2033   0.3758   8.4455  
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                     Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)       -6.864e+01  3.590e+00 -19.119  < 2e-16 ***
    ## nucleotides        1.224e-07  7.476e-08   1.637   0.1017    
    ## CpG_density       -1.855e+01  4.136e+00  -4.485 7.30e-06 ***
    ## `meth rate`        3.273e+01  1.597e+00  20.499  < 2e-16 ***
    ## msr                2.716e+02  1.213e+01  22.384  < 2e-16 ***
    ## inverted_msr      -1.706e+01  2.216e+00  -7.700 1.36e-14 ***
    ## ecdf               1.882e-01  9.802e-02   1.920   0.0548 .  
    ## `inverted ecdf`    9.660e-02  7.458e-02   1.295   0.1952    
    ## residual          -2.444e+02  1.290e+01 -18.950  < 2e-16 ***
    ## inverted_residual         NA         NA      NA       NA    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(4.1911) family taken to be 1)
    ## 
    ##     Null deviance: 5214.5  on 2766  degrees of freedom
    ## Residual deviance: 2965.3  on 2758  degrees of freedom
    ## AIC: 14730
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  4.191 
    ##           Std. Err.:  0.190 
    ## 
    ##  2 x log-likelihood:  -14709.734

![](MSR_and_expression_stomach_1e4_files/figure-markdown_github/unnamed-chunk-17-1.png)

#### Predicting log(TPM)

Distribution of TPM values (only for regions that contains some genes). ![](MSR_and_expression_stomach_1e4_files/figure-markdown_github/unnamed-chunk-18-1.png)

Linear model for log(TPM) with standard predictors:

    ## 
    ## Call:
    ## lm(formula = log_tpm ~ ., data = model_data)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -12.6703  -0.4418   0.3675   1.0995   5.2812 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -2.343e+00  8.573e-01  -2.734  0.00631 ** 
    ## nucleotides -1.974e-06  1.973e-07 -10.001  < 2e-16 ***
    ## CpG_density -2.577e+01  1.272e+01  -2.026  0.04284 *  
    ## `meth rate`  5.510e+00  9.615e-01   5.731 1.11e-08 ***
    ## gene_count   1.803e-01  7.440e-03  24.233  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.975 on 2650 degrees of freedom
    ##   (23 observations deleted due to missingness)
    ## Multiple R-squared:  0.3085, Adjusted R-squared:  0.3074 
    ## F-statistic: 295.6 on 4 and 2650 DF,  p-value: < 2.2e-16

Linear model for TPM with all features and MSR statistics:

    ## 
    ## Call:
    ## lm(formula = log_tpm ~ ., data = model_data)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -8.9883 -0.5826  0.2021  0.9444  5.7676 
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)             -1.322e+02  1.036e+01 -12.766  < 2e-16 ***
    ## nucleotides              1.459e-07  2.080e-07   0.701  0.48312    
    ## CpG_density              2.303e+01  1.177e+01   1.957  0.05045 .  
    ## `meth rate`              7.530e+01  4.567e+00  16.486  < 2e-16 ***
    ## msr                      4.819e+02  3.551e+01  13.572  < 2e-16 ***
    ## inverted_msr            -5.691e+01  6.212e+00  -9.161  < 2e-16 ***
    ## ecdf                    -1.328e-01  2.747e-01  -0.483  0.62900    
    ## `inverted ecdf`          6.355e-01  2.047e-01   3.104  0.00193 ** 
    ## residual                -3.930e+02  3.770e+01 -10.424  < 2e-16 ***
    ## inverted_residual               NA         NA      NA       NA    
    ## gene_count               1.075e-01  7.212e-03  14.904  < 2e-16 ***
    ## genes_nucleotides_count  9.412e-07  1.436e-07   6.553 6.77e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.712 on 2593 degrees of freedom
    ##   (74 observations deleted due to missingness)
    ## Multiple R-squared:  0.4486, Adjusted R-squared:  0.4464 
    ## F-statistic: 210.9 on 10 and 2593 DF,  p-value: < 2.2e-16

![](MSR_and_expression_stomach_1e4_files/figure-markdown_github/unnamed-chunk-20-1.png)

Linear model for TPM with all features and MSR statistics, without information about genes:

    ## 
    ## Call:
    ## lm(formula = log_tpm ~ ., data = model_data[, c(to_predict, basic_predictors, 
    ##     msr_predictors)])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -9.6612 -0.6083  0.2365  1.0532  5.8986 
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       -1.796e+02  1.064e+01 -16.881   <2e-16 ***
    ## nucleotides        3.592e-07  2.185e-07   1.644   0.1003    
    ## CpG_density        1.603e+00  1.241e+01   0.129   0.8973    
    ## `meth rate`        9.572e+01  4.688e+00  20.417   <2e-16 ***
    ## msr                6.699e+02  3.605e+01  18.583   <2e-16 ***
    ## inverted_msr      -6.296e+01  6.563e+00  -9.593   <2e-16 ***
    ## ecdf               2.068e-01  2.908e-01   0.711   0.4771    
    ## `inverted ecdf`    5.155e-01  2.170e-01   2.376   0.0176 *  
    ## residual          -5.808e+02  3.847e+01 -15.097   <2e-16 ***
    ## inverted_residual         NA         NA      NA       NA    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.817 on 2595 degrees of freedom
    ##   (74 observations deleted due to missingness)
    ## Multiple R-squared:  0.3784, Adjusted R-squared:  0.3764 
    ## F-statistic: 197.4 on 8 and 2595 DF,  p-value: < 2.2e-16

Linear model for TPM with some features:

    ## 
    ## Call:
    ## lm(formula = log_tpm ~ (model_data$inverted_msr) + (model_data$CpG_density) + 
    ##     (meth_rate) + (model_data$gene_count) + (model_data$msr), 
    ##     data = model_data)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -10.2501  -0.5352   0.2486   0.9994   4.6921 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)             -32.543413   2.616437 -12.438  < 2e-16 ***
    ## model_data$inverted_msr -85.147389   3.883254 -21.927  < 2e-16 ***
    ## model_data$CpG_density   44.821984   6.489619   6.907 6.21e-12 ***
    ## meth_rate                38.451478   2.285999  16.820  < 2e-16 ***
    ## model_data$gene_count     0.140415   0.006967  20.155  < 2e-16 ***
    ## model_data$msr          132.336954   9.362353  14.135  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.773 on 2598 degrees of freedom
    ##   (74 observations deleted due to missingness)
    ## Multiple R-squared:  0.4074, Adjusted R-squared:  0.4062 
    ## F-statistic: 357.2 on 5 and 2598 DF,  p-value: < 2.2e-16
