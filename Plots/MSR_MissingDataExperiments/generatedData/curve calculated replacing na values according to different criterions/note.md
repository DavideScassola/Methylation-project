Where not specified, the criterion for sobstituting na values is the following: 
If a certain bin has x ones and y missing values, then the total count is given by summing to x a  binomial random variable with n=y and p= an estimated proportion.
The estimated proportion is given by a weighted average of the proportion of ones in the bin, and the proportion of ones in the full binary string, with the weight depending on the proportion of missing values in the bin.

Anyway, I think this method give biased results when replacing too much missing values.





