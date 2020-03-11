## Dinucleotides experiments

### Dinucleotides MSR comparison
This data is about chromosome 1 of a mice.

The strings of bits that are used to calculate the MSR are such that at the position i there is a 1 if at the position i of the genome it starts a certain dinucleotide, 0 otherwise.

The total number of nucleotides in the chromosome 1 of a mice (mm10) is about 195,471,971.

M is the total count of sites found.

chr1, CG, M=1477561	MSR area = 0.2736858

chr1, CC, M=9741343	MSR area = 0.2566359

chr1, TA, M=12527811	MSR area = 0.251933

chr1, CA, M=14258405	MSR area = 0.2415733

I think MSR areas are not significantly different, moreover a low M correspond often to a higher MSR area.

for CG (but it's similar for the other dinucleotides)
maximum sum of resolution and relevance is at bin size of 394,
maximum relevance is at bin size of 128430

### CpG sites vs bernoulli
This data is about chromosome 1 of a human.

The "independent Bernoulli" MSR curve is generated from data with a length equal to the CpG sites length, and every bit is a Bernoulli trial with p = CpG sites proportion. In this way the M of the two binary strings is almost the same.
The two MSR are significantly different.
