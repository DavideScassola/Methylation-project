MSR and expression for K562, windows of 1000 sites
================

Here I will investigate if there is a relationship between the presence of genes and genes expression in a certain genomic region and the MSR (with some of its derivate statistics).

I chose K562 cells data in order to do this, and CpG windows of size 1000, that corresponds to a variable window size in term of nucleotides (on average about 100,000).

This is an example of total-rna-seq file, that shows for each "gene" its transcripts and some measures of expression. In this case I just kept two colums. The first one indicates the "gene"", the second one is the Transcript Per Million that is a relative measure of how much a gene is expressed.

    ##                    gene_id   TPM
    ##     1:     ENSG00000000003  0.49
    ##     2:     ENSG00000000005  0.01
    ##     3:     ENSG00000000419 78.32
    ##     4:     ENSG00000000457  8.23
    ##     5:     ENSG00000000460 43.52
    ##    ---                          
    ## 60818: gSpikein_ERCC-00165  0.00
    ## 60819: gSpikein_ERCC-00168  0.09
    ## 60820: gSpikein_ERCC-00170 44.26
    ## 60821: gSpikein_ERCC-00171 29.73
    ## 60822:    gSpikein_phiX174  0.00

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

    ## 21 rows had too many nucleotides

    ##    start_chr start_position end_position gene_count total_TPM meth rate
    ## 12      chr1         940826       961902          2     75.90 0.4715903
    ## 13      chr1         961902       982731          2      0.10 0.7012655
    ## 14      chr1         982731      1007283          2     13.36 0.4854497
    ## 23      chr1        1206432      1228380          2     31.60 0.8951453
    ## 24      chr1        1228380      1246900          2      5.97 0.6799055
    ## 26      chr1        1274217      1295503          2      9.00 0.8296808

The full scheme includes:

**nucleotides**: number of nucleotides in the window

**CpG density**: fraction of nucleotides that is a C of a CpG site (= 1000/nucleotides)

**meth rate**: ratio of methylated CpG sites

**gene\_count**: number of genes included (even partially) inside the interval

**total\_TPM**: sum of the TPMs of the genes in the interval

then the MSR and some related statistics: **msr**, **inverted msr**, **msr ecdf **, **inverted msr ecdf**, **residual** (residual of the linear regression between msr and meth rate), **inverted residual**.

First let's see if there are pairwise correlations between the features.

###### Basic features:

![](MSR_and_expression_K562_1e3_files/figure-markdown_github/unnamed-chunk-6-1.png)

log(TPM) is considered only for fragments with at least a gene.

###### Comparison with simple MSR statistics:

![](MSR_and_expression_K562_1e3_files/figure-markdown_github/unnamed-chunk-7-1.png)

###### Comparison with other MSR statistics:

![](MSR_and_expression_K562_1e3_files/figure-markdown_github/unnamed-chunk-8-1.png)

inverted msr vs log(tpm): ![](MSR_and_expression_K562_1e3_files/figure-markdown_github/unnamed-chunk-9-1.png)

#### Predicting gene presence

Check if features can predict gene presence:

![](MSR_and_expression_K562_1e3_files/figure-markdown_github/unnamed-chunk-10-1.png)

![](MSR_and_expression_K562_1e3_files/figure-markdown_github/unnamed-chunk-11-1.png)

![](MSR_and_expression_K562_1e3_files/figure-markdown_github/unnamed-chunk-12-1.png)

![](MSR_and_expression_K562_1e3_files/figure-markdown_github/unnamed-chunk-13-1.png)

![](MSR_and_expression_K562_1e3_files/figure-markdown_github/unnamed-chunk-14-1.png)

The fraction of fragments that have at least one gene inside is

    ## [1] 0.4420778

Logistic Regression Model for gene presence with basic predictors (nucleotides, CpG\_density, meth rate):

    ##        prediction
    ## actual      FALSE      TRUE
    ##   FALSE 0.4524840 0.1054078
    ##   TRUE  0.1903731 0.2517350

    ## 
    ## accuracy:  0.7042191

Logistic Regression Model with inverted\_msr as predictor

    ##        prediction
    ## actual      FALSE      TRUE
    ##   FALSE 0.3780359 0.1660416
    ##   TRUE  0.1862506 0.2696719

    ## 
    ## accuracy:  0.6477078

Adding other predictors doesn't significantly improve the accuracy.

#### Predicting log(TPM)

Distribution of TPM values (only for regions that contains some genes) ![](MSR_and_expression_K562_1e3_files/figure-markdown_github/unnamed-chunk-18-1.png)

Linear model for TPM with standard predictors:

    ## 
    ## Call:
    ## lm(formula = log_tpm ~ ., data = model_data[, standard_predictors])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -12.7421  -1.8543   0.3118   2.0978  11.8434 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)             -4.034e-01  3.275e-01  -1.232   0.2181    
    ## gene_count               8.269e-01  3.334e-02  24.799  < 2e-16 ***
    ## nucleotides             -2.786e-06  1.578e-06  -1.766   0.0775 .  
    ## CpG_density             -5.129e+01  7.033e+00  -7.293 3.21e-13 ***
    ## `meth rate`             -4.151e+00  2.944e-01 -14.098  < 2e-16 ***
    ## genes_nucleotides_count  4.770e-06  1.077e-06   4.428 9.58e-06 ***
    ## msr_density              9.063e+00  1.194e-01  75.931  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 3.188 on 11806 degrees of freedom
    ## Multiple R-squared:  0.4051, Adjusted R-squared:  0.4048 
    ## F-statistic:  1340 on 6 and 11806 DF,  p-value: < 2.2e-16

Linear model for TPM with all features and MSR statistics:

    ## 
    ## Call:
    ## lm(formula = log_tpm ~ ., data = model_data[, c(standard_predictors, 
    ##     msr_predictors)])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -11.855  -1.887   0.291   2.055  11.617 
    ## 
    ## Coefficients: (2 not defined because of singularities)
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)             -1.785e+00  4.460e-01  -4.002 6.31e-05 ***
    ## gene_count               7.649e-01  3.271e-02  23.383  < 2e-16 ***
    ## nucleotides             -4.679e-08  1.566e-06  -0.030    0.976    
    ## CpG_density             -3.900e+01  6.968e+00  -5.597 2.23e-08 ***
    ## `meth rate`             -2.630e+00  2.980e-01  -8.825  < 2e-16 ***
    ## genes_nucleotides_count  5.291e-06  1.054e-06   5.018 5.29e-07 ***
    ## msr_density              4.377e+00  3.298e-01  13.271  < 2e-16 ***
    ## msr                     -6.565e+00  1.150e+00  -5.707 1.18e-08 ***
    ## inverted_msr             2.275e+01  1.205e+00  18.883  < 2e-16 ***
    ## ecdf                    -7.969e-01  1.446e-01  -5.512 3.63e-08 ***
    ## `inverted ecdf`         -1.893e+00  1.424e-01 -13.299  < 2e-16 ***
    ## residual                        NA         NA      NA       NA    
    ## inverted_residual               NA         NA      NA       NA    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 3.115 on 11802 degrees of freedom
    ## Multiple R-squared:  0.432,  Adjusted R-squared:  0.4315 
    ## F-statistic: 897.5 on 10 and 11802 DF,  p-value: < 2.2e-16

Linear model for TPM with some features:

    ## 
    ## Call:
    ## lm(formula = log_tpm ~ (model_data$inverted_msr) + (model_data$CpG_density) + 
    ##     (meth_rate) + (model_data$gene_count), data = model_data)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -13.604  -2.351   0.621   2.488  11.119 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)             -4.31986    0.29114 -14.838  < 2e-16 ***
    ## model_data$inverted_msr 29.46298    0.51787  56.893  < 2e-16 ***
    ## model_data$CpG_density  -3.84026    4.65028  -0.826    0.409    
    ## meth_rate               -2.42098    0.32025  -7.560 4.34e-14 ***
    ## model_data$gene_count    1.00531    0.03563  28.213  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 3.515 on 11808 degrees of freedom
    ## Multiple R-squared:  0.2766, Adjusted R-squared:  0.2764 
    ## F-statistic:  1129 on 4 and 11808 DF,  p-value: < 2.2e-16
