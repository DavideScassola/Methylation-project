MSR and expression for H1 genes
================

    ## TPM fraction:  0.7970407

    ## protein_coding TPM fraction:  0.6840699

###### Basic features:

![](MSR_and_expression_genebodies_H1_files/figure-markdown_github/unnamed-chunk-5-1.png)

###### Comparison with MSR statistics:

![](MSR_and_expression_genebodies_H1_files/figure-markdown_github/unnamed-chunk-6-1.png)

meth\_autocorrelation vs log(tpm): ![](MSR_and_expression_genebodies_H1_files/figure-markdown_github/unnamed-chunk-8-1.png)

CG\_list\_inverted\_msr vs log(tpm): ![](MSR_and_expression_genebodies_H1_files/figure-markdown_github/unnamed-chunk-9-1.png)

meth\_rate\_binary vs log(tpm): ![](MSR_and_expression_genebodies_H1_files/figure-markdown_github/unnamed-chunk-10-1.png)

drift vs log(tpm): ![](MSR_and_expression_genebodies_H1_files/figure-markdown_github/unnamed-chunk-11-1.png)

    ## missing data:  14.58641 %

    ## 
    ## train_data_proportion:  0.6

    ## 
    ## 
    ## basic missing data:  0 %

    ## 
    ## train_data_proportion:  0.6

Linear model for log\_tpm with basic features:

    ## 
    ## Call:
    ## lm(formula = formula, data = train_model_data[, c(response_variable, 
    ##     predictors)])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -10.5007  -1.4012   0.3873   1.5646   7.7899 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       13.7439     0.5047   27.23   <2e-16 ***
    ## log_nucleotides   -6.1614     0.2269  -27.15   <2e-16 ***
    ## CG_density      -109.8607     4.5037  -24.39   <2e-16 ***
    ## log_CG_count       7.1666     0.2408   29.77   <2e-16 ***
    ## meth_rate         -2.8489     0.1692  -16.84   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.275 on 10075 degrees of freedom
    ## Multiple R-squared:  0.1043, Adjusted R-squared:  0.1039 
    ## F-statistic: 293.2 on 4 and 10075 DF,  p-value: < 2.2e-16
    ## 
    ## Test data R squared:  0.1077674

Linear model for log\_tpm with basic features with meth\_autocorrelation and drift:

    ## 
    ## Call:
    ## lm(formula = formula, data = train_model_data[, c(response_variable, 
    ##     predictors)])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -8.7069 -1.0301  0.1415  1.1676  8.1023 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            0.47119    0.49322   0.955    0.339    
    ## log_nucleotides       -1.84034    0.19861  -9.266  < 2e-16 ***
    ## CG_density           -20.72545    4.00866  -5.170 2.38e-07 ***
    ## log_CG_count           2.25243    0.21248  10.601  < 2e-16 ***
    ## meth_rate              2.10195    0.15527  13.538  < 2e-16 ***
    ## meth_autocorrelation   5.30003    0.09633  55.021  < 2e-16 ***
    ## drift                -23.68195    0.91807 -25.795  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.861 on 10073 degrees of freedom
    ## Multiple R-squared:  0.4008, Adjusted R-squared:  0.4004 
    ## F-statistic:  1123 on 6 and 10073 DF,  p-value: < 2.2e-16
    ## 
    ## Test data R squared:  0.419495

    ## 
    ## keeping also data with NA msr features:

    ## Test data R squared:  0.4812646

Linear model for TPM with all predictors:

    ## 
    ## Call:
    ## lm(formula = formula, data = train_model_data[, c(response_variable, 
    ##     predictors)])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -8.6982 -1.0285  0.1362  1.1676  8.1017 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            6.8933     1.2605   5.469 4.64e-08 ***
    ## log_nucleotides       -2.0079     0.2017  -9.955  < 2e-16 ***
    ## CG_density           -26.2950     4.1345  -6.360 2.11e-10 ***
    ## log_CG_count           2.2802     0.2173  10.494  < 2e-16 ***
    ## meth_rate              2.1718     0.2384   9.111  < 2e-16 ***
    ## meth_autocorrelation   5.3317     0.1387  38.436  < 2e-16 ***
    ## drift                -23.1976     1.1156 -20.793  < 2e-16 ***
    ## CGsites_msr          -17.3206     3.6138  -4.793 1.67e-06 ***
    ## meth_msr              -2.7155     1.2366  -2.196   0.0281 *  
    ## unmeth_msr             0.2860     0.7745   0.369   0.7120    
    ## CG_list_msr            0.2733     0.8741   0.313   0.7545    
    ## CG_list_inverted_msr  -0.6933     0.6520  -1.063   0.2877    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.858 on 10068 degrees of freedom
    ## Multiple R-squared:  0.403,  Adjusted R-squared:  0.4024 
    ## F-statistic: 617.9 on 11 and 10068 DF,  p-value: < 2.2e-16
    ## 
    ## Test data R squared:  0.4239681

Linear model with few predictors

    ## 
    ## Call:
    ## lm(formula = formula, data = train_model_data[, c(response_variable, 
    ##     predictors)])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -9.1625 -1.0552  0.1969  1.2300  7.3576 
    ## 
    ## Coefficients:
    ##                        Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          -9.978e-01  1.190e-01  -8.385   <2e-16 ***
    ## meth_autocorrelation  5.236e+00  8.707e-02  60.137   <2e-16 ***
    ## drift                -1.983e+01  8.937e-01 -22.188   <2e-16 ***
    ## nucleotides           1.166e-06  1.340e-07   8.699   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.9 on 10076 degrees of freedom
    ## Multiple R-squared:  0.3748, Adjusted R-squared:  0.3747 
    ## F-statistic:  2014 on 3 and 10076 DF,  p-value: < 2.2e-16
    ## 
    ## Test data R squared:  0.398728

Lasso:

    ## lambda: 0.1

    ## 11 x 1 sparse Matrix of class "dgCMatrix"
    ##                               s0
    ## log_nucleotides        .        
    ## CG_density             .        
    ## log_CG_count           0.2509730
    ## meth_rate              0.6150052
    ## meth_autocorrelation   4.9620363
    ## drift                -17.9304300
    ## CGsites_msr           -1.4757285
    ## meth_msr               .        
    ## unmeth_msr             .        
    ## CG_list_msr           -0.8794890
    ## CG_list_inverted_msr   .

    ## 
    ## Test data R squared:  0.4057738
