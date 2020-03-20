### Methylation state assignment methods
50% threshold: a base is methylated if the percentage of methylated reads is greater than 50%

Adaptive (rate preserving) threshold:  first the overall methylation rate is estimated by averaging the rates of all the reads, then a base is methylated if the percentage of methylaed reads is above a certain threshold such that the rate is the same as the estimated one (For example in our case the new threshold was between 60% and 80%).

Stochastic assignment: the mehylation state of a base is a Bernoulli random variable with p equal to the rate of methylated reads for that base.


