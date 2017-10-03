#! /usr/bin/env Rscript

library(VariantAnnotation)

# parse input parameters

args <- commandArgs(TRUE)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$input_vcf)) {stop("no input VCF file")} else {input_vcf = args$input_vcf}
if(is.null(args$nb_classes)) {nb_classes = 100} else {nb_classes = as.numeric(args$nb_classes)}
if(is.null(args$AF_field)) {AF_field = "AF"} else {AF_field = args$AF_field}
if(is.null(args$min_AF)) {min_AF = 0} else {min_AF = as.numeric(args$min_AF)}
if(is.null(args$max_AF)) {max_AF = 1} else {max_AF = as.numeric(args$max_AF)}
if(is.null(args$nCPU)) {nCPU = 1} else {nCPU = as.numeric(args$nCPU)}

# contruct a prop barplot for substitutions stratified by AF

vcf <- open(VcfFile(input_vcf,  yieldSize=1000000))
vcf_chunk = readVcf(vcf, "hg19")
AF_df = as.data.frame(geno(vcf_chunk, AF_field))
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

AF_classes = seq(0, 1, length = nb_classes + 1)

l = lapply(1:nb_classes, function(i){
  nb_af_ok = apply(AF_df_snv[,1:(ncol(AF_df_snv)-3)], 1, function(x) {
    sum(x>AF_classes[i] & x<=AF_classes[i+1], na.rm = T)
  })
  t = table(AF_df_snv[which(nb_af_ok>0),"subst"])[c("AG","AC","AT","CA","CG","CT")]  # compute number of mutations in the range of AF
  t[which(is.na(t))] = 0 # if there is not mutation for a particular subsitution, set the count as 0
  t = t / length(which(nb_af_ok>0)) # transform counts to proportions
  names(t) = c("AG","AC","AT","CA","CG","CT") # correct names, because if NA in the table, the name was not kept
  t # return the result
})

prop_subst = matrix(unlist(l), ncol = nb_classes, byrow = F)
rownames(prop_subst) = c("AG","AC","AT","CA","CG","CT")
prop_subst = prop_subst[names(sort(apply(prop_subst, 1, var), decreasing = T)),]

# output the resulting PDF plot

pdf("substitutions_proportion_by_AF.pdf",10,5)
par(lwd = 0.25)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
cols = c("purple1","red1","chartreuse2","dodgerblue3","yellow1","chocolate2")
barplot(prop_subst, col = cols, space=0, xlab="AF", ylab="Proportion")
axis(1,labels=AF_classes, at=c(0, seq(1:nb_classes)))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("topright", legend=rownames(prop_subst), col = cols, pch=19)
dev.off()
