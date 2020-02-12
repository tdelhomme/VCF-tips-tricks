# VCF-tips-tricks

## Piece of codes to analyse VCF files

### 1. Control the distribution of substitutions

The distribution of substitutions is expected to keep the same whatever the range of a particular variable you are dealing with.  
A good test on the quality of the variant calling is to have a look on these distributions given a specific variable (allelic fraction, samples, ...). This gives a measure of quality of the variant calling for the particular variable

#### 1.1 In function of the allelic fraction

[Here](https://github.com/tdelhomme/VCF-tips-tricks/blob/master/code/cumulative_plot_vcf.r) is a R code to plot the distribution of substitutions for multiple classes of allelic fractions.  
The input parameters are :
  * __input_vcf__ (no default)
  * __min_AF__ (default=0)
  * __max_AF__ (default=1)
  * __nCPU__ (default=1)
  * __nb_classes__ (default=100)
  * __AF_field__ (default=AF)

Command line example:
```
cumulative_plot_vcf.r --input_vcf=path_input_vcf_file.vcf
```

[Here](https://github.com/tdelhomme/VCF-tips-tricks/blob/master/plots/substitutions_proportion_by_AF.png) is an example of output.  


## Piece of codes to analyse annotated VCF files (with annovar to keep only positive samples)

Arguments are the same than previously described.  

Command line example:
```
cumulative_plot_annotation.r --input_txt=path_input_vcf_file.vcf
```

Note: for the moment this only output graph per sample, need to add an option for AF scale.

## Piece of code to analyze raw VCF files in R
