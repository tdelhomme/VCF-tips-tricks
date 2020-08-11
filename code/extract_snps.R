args <- commandArgs(TRUE)
parseArgs <- function(x) {
  res = strsplit(sub("^--", "", x), "=")
  if(length(unlist(res))==1) res[[1]][2]=""
  return(res)
}

argsL <- as.list(do.call("cbind", parseArgs(args))[c(F,T)])
names(argsL) <- as.list(do.call("cbind", parseArgs(args))[c(T,F)])
args <- argsL;rm(argsL)

if(! is.null(args$help)) {
  cat("
      Mandatory arguments:
      --input_vcf       Input .gz vcf file  
      --snps            Text file containing the list of RS snps we want to keep in the VCF
      --output_vcf      Output .vcf file
      --help \n\n")
  q(save="no")
}

if(is.null(args$snps)) {stop("Option --snps should be provided")} else{snps=args$snps}
if(is.null(args$input_vcf)) {stop("Option --input_vcf should be provided")} else{input_vcf=args$input_vcf}
if(is.null(args$output_vcf)) {stop("Option --output_vcf should be provided")} else{output_vcf=args$output_vcf}

snps_chip = read.table("/g/strcombio/fsupek_cancer1/TCGA/annovar/GenomeWideSNP_6.na35.annot_IDs.csv", h=F)
snps_tcga = read.table(snps, h=F)
rownames(snps_chip) = snps_chip[,2]
kept_snps = as.character(snps_chip[snps_tcga$V1,1])

suppressMessages(library(VariantAnnotation))

vcf <- open(VcfFile(input_vcf,  yieldSize=10000))
vcf_chunk = readVcf(vcf, "hg19")
nbc = 1

#and continue
while(dim(vcf_chunk)[1] != 0) {
  print(paste(date(), " INFO: working on chunk number ", nbc, sep=""))
  ids_vcf_snps = 1:nrow(vcf_chunk) #to efficiently search for snps: use names() and get the row ids to keep
  names(ids_vcf_snps) = names(vcf_chunk)
  ids_to_keep = ids_vcf_snps[kept_snps]
  ids_to_keep = ids_to_keep[which(!is.na(ids_to_keep))] # NA is when the snp is not in this vcf chunk
  
  #remove snps already seen (we had duplicates in vcf output...)
  kept_snps = setdiff(kept_snps, names(ids_to_keep))
  
  # keep only our snps in the vcf
  vcf_chunk = vcf_chunk[ids_to_keep,]
  
  #write out the selected VCF
  con = file(output_vcf, open = "a")
  writeVcf(vcf_chunk, con)
  vcf_chunk = readVcf(vcf, "hg19")
  close(con)
  nbc = nbc + 1
}

# zip and index the final vcf
system(paste(" bgzip -c ", output_vcf, " > ", output_vcf, ".gz", sep=""))
