library(VariantAnnotation)

input_vcf = "input.vcf.gz"

allelic_fractions_field = "AF"
af_length = 100

# contruct a prop barplot for substitutions stratified by AF

vcf <- open(VcfFile(input_vcf,  yieldSize=1000000))
vcf_chunk = readVcf(vcf, "hg19")
AF_df = as.data.frame(geno(vcf_chunk, allelic_fractions_field))
AF_df$ref = as.character(ref(vcf_chunk))
AF_df$alt = as.character(unlist(alt(vcf_chunk)))
AF_df_snv = AF_df[which(AF_df$ref %in% c("A","G","C","T") & AF_df$alt %in% c("A","G","C","T")),]
AF_df_snv$subst = paste(AF_df_snv$ref, AF_df_snv$alt, sep="")

# compute only the 6 non-equivalent substitutions
AF_df_snv[which(AF_df_snv$subst == "TC"),"subst"] = "AG"
AF_df_snv[which(AF_df_snv$subst == "TG"),"subst"] = "AC"
AF_df_snv[which(AF_df_snv$subst == "TA"),"subst"] = "AT"
AF_df_snv[which(AF_df_snv$subst == "GT"),"subst"] = "CA"
AF_df_snv[which(AF_df_snv$subst == "GC"),"subst"] = "CG"
AF_df_snv[which(AF_df_snv$subst == "GA"),"subst"] = "CT"

AF_classes = seq(0, 1, length = af_length + 1)

l = lapply(1:af_length, function(i){
  nb_af_ok = apply(AF_df_snv[,1:(ncol(AF_df_snv)-3)], 1, function(x) {
    sum(x>AF_classes[i] & x<=AF_classes[i+1], na.rm = T)
  })
  t = table(AF_df_snv[which(nb_af_ok>0),"subst"])[c("AG","AC","AT","CA","CG","CT")]  # compute number of mutations in the range of AF
  t[which(is.na(t))] = 0 # if there is not mutation for a particular subsitution, set the count as 0
  t = t / length(which(nb_af_ok>0)) # transform counts to proportions
  names(t) = c("AG","AC","AT","CA","CG","CT") # correct names, because if NA in the table, the name was not kept
  t # return the result
})

prop_subst = matrix(unlist(l), ncol = af_length, byrow = F)
rownames(prop_subst) = c("AG","AC","AT","CA","CG","CT")

pdf("substitutions_proportion_by_AF.pdf",10,5)
par(lwd = 0.25)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
barplot(prop_subst, col = rainbow(length(rownames(prop_subst))), space=0)
axis(1,labels=AF_classes, at=c(0, seq(1:af_length)))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("topright", legend=c("AG","AC","AT","CA","CG","CT"), col = rainbow(length(rownames(prop_subst))), pch=19)
dev.off()
