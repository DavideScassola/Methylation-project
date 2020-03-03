### Curves calculated avoiding bins with na values
In these case the MSR curves were calculate selecting only bins containing no missing data.
This implied that not all of the available data were used to calculate resolution and relevance, especially for large bins (bins containing at least one missing value were not used).

### Curve calculated replacing na values according to different criterions
In these cases instead bins with a percentage of missing values below a certain threshold were accepted. 
The total count of ones in a certain bin was then adjusted according to a weighted average of the proportion of ones in the bin and the global proportion of ones in the entire string.


