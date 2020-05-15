### Density - MSR Plots
In these plots each point represents a window of a specified size in the genome.

The confidence interval is built simulating for each density value (proportion of methylated sites) 
several windows in a random way (reshuffling a fixed number of methylated sites).
Since the amount of data that can be simulated for each density value is limited (10^5 - 10^6)
the given confidence interval is almost equivalent to the min and max values obtained from the simulations.

The "inverted" version is the same plot but MSR and density are calculated for the negation of the given binary string.

### Fragments comparison
In these plots are shown some of the windows of size 1000 that have a different "significance",
the cdf is the empirical cumulative distribution function evaluated for a certain (density-msr) pair; a value that's almost 0 or 1 is almost impossible to be obtained with random data.