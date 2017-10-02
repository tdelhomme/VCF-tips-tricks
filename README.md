# VCF-tips-tricks

## Piece of codes to analyse VCF files

### 1. Control the distribution of substitutions

The distribution of substitutions is expected to keep the same whatever the range of a particular variable you are dealing with.  
A good test on the quality of the variant calling is to have a look on these distributions given a specific variable (allelic fraction, samples, ...). This gives a measure of quality of the variant calling for the particular variable

#### 1.1 In function of the allelic fraction

[Here]() is a R code to plot the distribution of substitutions for multiple classes of allelic fractions.  
The input parameters are :
  * __input_VCF__ (no default)
  * __min_AF__ (default=0)
  * __max_AF__ (default=1)
  * __nCPU__ (default=1)
