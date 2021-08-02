# VCF-tips-tricks

## Piece of code to analyze VCF files in R

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

### 2. Tips for strelka2 VCFs analyzes

[Here]() is a piece of R code that computes the variant allelic fraction of each individual mutation found in a VCF file output of Strelka2 variant caller. Note that the VCF should be compressed with bgzip and indexed with tabix.

## Piece of code to analyze annotated VCF files (with annovar to keep only positive samples)

Arguments are the same than previously described.  

Command line example:
```
cumulative_plot_annotation.r --input_txt=path_input_vcf_file.vcf
```

Note: for the moment this only output graph per sample, need to add an option for AF scale.

## 3. Scripts to analyse/modify VCFs

### Correct the REF/ALT in a VCF file
using [This script](https://github.com/tdelhomme/VCF-tips-tricks/blob/master/code/correct_refalt_vcf.R). Some tools like plink does not output (in the vcf) the REF/ALT based on the reference genome but based on the SNP chip fluorescence. This can be problematic when comparing 2 VCF files that should have the same chr-pos-ref (for imputation as an example) or for annovar that requires REF to match the databases. This script takes as input a VCF and an annotation file in order to correct the REF and ALT in the VCF. This annotation file is basically from convert2annovar (or done by hand), for each ID field in the VCF, the informations about the SNP, e.g.:
```
chr1	564477	564477	A	G	rs6650104	SNP_A-8575395
chr1	564621	564621	C	T	rs10458597	SNP_A-8575125
chr1	565400	565400	C	T	rs8179414	SNP_A-8575389
```

### Extract a list of SNPs (Id field in input VCF) from a VCF (output a new VCF)
The script [extract_snps.R](https://github.com/tdelhomme/VCF-tips-tricks/blob/master/code/extract_snps.R) takes 3 input parameters:
- --input_vcf
- --output_vcf
- --snps

It creates a new vcf 'output_vcf' by retaining the snps present in the txt file 'snps' from the 'input_vcf' (in the ID field).

## 4. BCFtools tricks

### 4.1 Annotate a VCF with BAD_REGION from positions in a bed file

2 files are needed to start:
* a VCF file (let's assume it is named `variants.vcf`)
* a BED file (let's assume it is named `bad_regions.bed`)
The BED file looks like:
```
10	60680	60690	1
```
First compress and index the bed if it is not already done:
```
bgzip -c bad_regions.bed > bad_regions.bed.gz
tabix -p bed bad_regions.bed.gz
```
Then create the INFO line that would be added in the header:
```
echo "##INFO=<ID=BAD_REGION,Number=0,Type=Flag,Description="My bad region for some reason">" > bad_regions.bed.hdr
```
Finally, annotate the VCF by adding a new field `BAD_REGION` if the variant is in the bed regions:
```
bcftools annotate -a bad_regions.bed.gz -h bad_regions.bed.hdr -c CHROM,FROM,TO,BAD_REGION variants.vcf -Oz -o variants_annotated.vcf.gz
```
One can also want to filter out the bad regions, for that you can do this:
```
bcftools filter variants_annotated.vcf.gz -s BAD_REGION -m+ -e 'INFO/BAD_REGION=1' | bcftools view -f PASS -Oz > variants_annotated_filter.vcf.gz
```
