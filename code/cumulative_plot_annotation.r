#! /usr/bin/env Rscript

source("https://gist.githubusercontent.com/mfoll/a4dfbb92068dc559f130/raw/714dc8c2e97987fd4385dcef2722b3ef986d38d6/get_vcf_data.r")

args <- commandArgs(TRUE)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$input_txt)) {stop("no input annotated VCF file")} else {input_txt = args$input_txt}
if(is.null(args$nb_classes)) {nb_classes = 100} else {nb_classes = as.numeric(args$nb_classes)}
if(is.null(args$AF_field)) {AF_field = "AF"} else {AF_field = args$AF_field}
if(is.null(args$Ref_field)) {Ref_field = "Ref"} else {Ref_field = args$Ref_field}
if(is.null(args$Alt_field)) {Alt_field = "Alt"} else {Alt_field = args$Alt_field}
if(is.null(args$Sample_field)) {Sample_field = "SM"} else {Sample_field = args$Sample_field}
if(is.null(args$min_AF)) {min_AF = 0} else {min_AF = as.numeric(args$min_AF)}
if(is.null(args$max_AF)) {max_AF = 1} else {max_AF = as.numeric(args$max_AF)}
if(is.null(args$nCPU)) {nCPU = 1} else {nCPU = as.numeric(args$nCPU)}

variants = read.table(input_txt,quote="\"",stringsAsFactors=F,sep="\t",header=T,na.strings = c(".","NA"))

variants$AF = get_genotype(variants$GENOTYPE,variants$FORMAT[1], AF_field)

vars = variants[,c(Ref_field, Alt_field, Sample_field)]

vars = vars[which(vars[, Ref_field] %in% c("A","G","C","T") & vars[, Alt_field] %in% c("A","G","C","T")),]
vars$subst = paste(vars[, Ref_field], vars[, Alt_field], sep="")

vars[which(vars$subst == "TC"),"subst"] = "AG"
vars[which(vars$subst == "TG"),"subst"] = "AC"
vars[which(vars$subst == "TA"),"subst"] = "AT"
vars[which(vars$subst == "GT"),"subst"] = "CA"
vars[which(vars$subst == "GC"),"subst"] = "CG"
vars[which(vars$subst == "GA"),"subst"] = "CT"

l = lapply(unique(vars[, Sample_field]), function(sm){
  t = table(vars[which(vars[, Sample_field]==sm),"subst"])[c("AG","AC","AT","CA","CG","CT")]  # compute number of mutations in the range of AF
  t[which(is.na(t))] = 0 # if there is not mutation for a particular subsitution, set the count as 0
  names(t) = c("AG","AC","AT","CA","CG","CT") # correct names, because if NA in the table, the name was not kept
  t / sum(t) # return the result
})

# construct matrix of proportion and order it by variance for substitution, and by increasing order for samples

prop_subst = matrix(unlist(l), ncol = length(unique(vars[, Sample_field])), byrow = F)
rownames(prop_subst) = c("AG","AC","AT","CA","CG","CT")
prop_subst = prop_subst[names(sort(apply(prop_subst, 1, var), decreasing = T)),]
ord = order(prop_subst[1,])
prop_subst = prop_subst[, ord]
colnames(prop_subst) = unique(vars[, Sample_field])[ord]

# output the resulting plot

cols = c("purple1","red1","chartreuse2","dodgerblue3","yellow1","chocolate2")
pdf("substitutions_proportion_by_sample.pdf",10,5)
par(lwd = 0.25)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
barplot(prop_subst, col = adjustcolor(cols, 0.75), space=0, xlab="", ylab="Proportion", las=2)
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("topright", legend=rownames(prop_subst), col = cols, pch=19)
dev.off()

print("DONE")
