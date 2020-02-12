vcf = "path/to/file.vcf.bgz"

genome = "hg38"

calls = readVcf( open(VcfFile(vcf,  yieldSize=100000)), genome)

# if you need to filter on PASS mutations
calls = calls[which(rowRanges(calls)$FILTER == "PASS"),]

# get alternative bases to search in the VCF columns
raw_alt = unlist(lapply(alt(calls), as.character))

alt = paste(raw_alt,"U",sep="")
vaf = unlist(lapply(1:length(alt2), function(i){
  a = alt2[i]
  ao = geno(calls[i,],a)[[2]]
  dp =geno(calls[i,],"DP")[,"TUMOR"]
  ao / dp
}))

# if you want to plot the distribution of variant allelic fractions
hist(vaf, xlab="VAF", ylab="counts", main="VAF distribution", col="grey")
