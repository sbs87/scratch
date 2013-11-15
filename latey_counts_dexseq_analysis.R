source("http://bioconductor.org/biocLite.R")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(Rsamtools)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
load("/Users/stsmith/Documents/Analysis/Latey/counts.Rda")
prefix<-"/Users/stsmith/Documents/Data/"
suffix<-"_thout/accepted_hits.bam"
samples.file.meta.pre<-read.table("/Users/stsmith/Documents/Data/files",sep="\t")
design_table<-with(samples.file.meta.pre,data.frame(V2,V3,V4,V5,row.names=1))
names(design_table)<-c("infection_3cats","infection","subject")

#Load Ref Genome
txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene

#Define exons
exons<-unlist(exonsBy(txdb,by="gene"))
exon_id<-exons$exon_id
length(exon_id)

exon_gene_mapping<-select(txdb,keys=exon_id,cols="GENEID",keytype="EXONID")

exon_ID<-exon_gene_mapping$EXONID
gene_ID<-exon_gene_mapping$GENEID
exon_range<-select(txdb,keys=exon_id,cols=c("EXONCHROM","EXONSTART","EXONEND","EXONSTRAND"),keytype="EXONID")
exon_range<-with(exon_range,data.frame(EXONCHROM,EXONSTART,EXONEND,EXONSTRAND))
names(exon_range)<-c("chr", "start", "end", "strand")


