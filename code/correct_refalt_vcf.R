#! /usr/bin/env Rscript

# Copyright (C) 2020 BIST/IRB Barcelona
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Script to correct REF and ALT in a VCF file based on the RS id from dbSNP
# WARNING: this version works for the moment only for hg19 -- should add a genome option

args <- commandArgs(TRUE)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$input_vcf)) {stop("no input VCF file")} else {input_vcf = args$input_vcf}
if(is.null(args$input_annot)) {stop("no input from annovar")} else {input_annot = args$input_annot}
annot = read.table(input_annot, h=F, sep="\t", stringsAsFactors=F)
if(is.null(args$chunk_size)) {chunk_size = 10000} else {chunk_size = as.numeric(args$chunk_size)}

# which column of the input_annot does contain the ID to look for in the VCF
if(is.null(args$col_vcf_id)) {col_vcf_id = ncol(annot)} else {col_vcf_id = as.numeric(args$col_vcf_id)}
ids_dup = paste( annot[,2], annot[,col_vcf_id], sep="__") # we can have dup snp (same pos-ref but different alt)
annot = annot[which(!duplicated(ids_dup)),] #take the first, hope it is ok
rownames(annot) = paste( annot[,2], annot[,col_vcf_id], sep="__")

if(is.null(args$out_vcf)) {args$out_vcf = out_vcf = paste(gsub(".vcf.gz","",input_vcf),
                                                          "_corrected.vcf",sep="")}else{out_vcf=args$out_vcf}
if(is.null(args$genome)) {genome = "hg19"} else {genome = as.numeric(args$genome)}

library(VariantAnnotation)

#initiate the first chunk
vcf <- open(VcfFile(input_vcf,  yieldSize=chunk_size))
vcf_chunk = readVcf(vcf, genome)

#and continue
chunkid = 1
while(dim(vcf_chunk)[1] != 0) {
  print(paste(date(), " INFO: reading the chunk ", chunkid, sep=""))
  # in case we have 2 dbnp id in the annot file, we should check the position in the vcf
  ids = paste( start(ranges(rowRanges(vcf_chunk))), names(rowRanges(vcf_chunk)), sep="__")
  true_ref = annot[ids,4] #rarely the snp has not been found in dbsnp
  true_ref[which(is.na(true_ref))] = as.character(ref(vcf_chunk)[which(is.na(true_ref))])
  true_alt = annot[ids,5]
  true_alt[which(is.na(true_alt))] = as.character(unlist(alt(vcf_chunk)[which(is.na(true_alt))]))
  ref(vcf_chunk) = DNAStringSet(x=true_ref)
  alt(vcf_chunk) = DNAStringSetList(as.list(true_alt))
  con = file(out_vcf, open = "a")
  writeVcf(vcf_chunk, con)
  vcf_chunk = readVcf(vcf, genome)
  close(con)
  chunkid = chunkid + 1
}


