
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


### Some considerations
- The cdfs are correlated with the residuals of a linear regression.
- The cdfs of the (msr, density) pairs of the windows are correlated with the inverted variants of the windows.


#### There is a biological meaning to this significance measure ?

The data does not fit the distribution of a random generated data, this is because the generative process is different. The methylation of the genome is characterized by a methylation probability that is around 80%, but this regularity is sometimes interrupted by demethylated CG islands, that are regions of contiguous demethylated sites (<10%) (and also methylated islands, but it's less frequent).

__About the windows of size 1000:__
- Windows that have a density that is >80% probably does not contain any island, in this case the values are in the expected bounds.
- Windows that have a density lower than 0.2 are probably inside large CG islands, even in this case the values are in the expected bounds with few exceptions.
- Windows that have a density between 0.2 and 0.8 instead are probably regions that contains one or more island of various size. In this case the atypical pattern caused by the islands may be the cause of unexpected low values of MSR.

__About the windows of size 10^5 and 10^6:__
For bigger windows (10^5, 10^6) the MSR is bigger than the expected, in this case all windows contains some islands so there is more homogeneity in the data. In those case the presence of various sized islands may be the cause of a bigger MSR.

__About the windows of size 10^4:__
The situation is similar to the previous one, but in this case probably there is more variability in the number of islands that are contained in each window.
