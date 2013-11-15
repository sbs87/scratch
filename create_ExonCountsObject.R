#source("http://bioconductor.org/biocLite.R")
library(DEXSeq)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Rsamtools)

#Define Variables
prefix<-"/Volumes/40GIG/"
suffix<-"_thout/accepted_hits.bam"
samples.file.meta.pre<-read.table("/Users/stevensmith/Documents/School/Maryland/Classes/Comp_Bio_Systems/project/files.txt",sep="\t")
design_table<-with(samples.file.meta.pre,data.frame(V2,V4,V5,V7,V8,row.names=1))
names(design_table)<-c("stability","infection","activity","subject")

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

#Load gapped alignment (bam file from thout)
#count.names<-list(c(),rownames(design_table))
#counts<-matrix(nrow=length(exons),ncol=nrow(samples.file.meta),dimnames=count.names)

#for (sample in seq(1:nrow(samples.file.meta))){
#  file.name<-paste(prefix,rownames(design_table)[sample],suffix,sep="")
#  aln<-readGappedAlignments(file.name)
#  new.counts<- countOverlaps(exons, aln)
#  counts[,sample]<-new.counts
#  print(rownames(design_table)[sample])
#}
write.table(counts,file="counts_as_file",quote=F)
head(counts)
#Make ExonCountSet
#exon_count_set<-newExonCountSet(countData=counts,design=design_table,geneIDs=gene_ID,exonIDs=exon_ID)
vaginal.model<-newExonCountSet(countData=counts,design=design_table$stability,geneIDs=gene_ID,exonIDs=exon_ID,exonIntervals=exon_range)

#Normalize the exon counts
vaginal.model<-estimateSizeFactors(vaginal.model)
sizeFactors(vaginal.model)

#Models
#formuladispersion<-count~sample + (stability + health)*exon
#formula0<-count ~ sample + health*exon + stability
#formula1<-count ~ sample + health*exon + stability*I(exon==exonID)
vaginal.model
vaginal.model<-estimateDispersions(vaginal.model)
vaginal.model<-fitDispersionFunction(vaginal.model)
plotDispEsts(vaginal.model)

vaginal.model<-testForDEU(vaginal.model,formula0=formula0,formula1=formula1,nCores=6)
tail(fData(vaginal.model)[,c("pvalue","padjust")])

vaginal.model<-estimatelog2FoldChanges(vaginal.model)
res3<-DEUresultTable(vaginal.model)
res3<-res3[!is.na(res3$padjust),]
write.table(res3[res3$padjust<0.1,]$geneID, "genes.txt",quote=F)
DEUresultTable(exon_count_set)

table(res1$padjust<0.1,res3$padjust<0.1)
bottom=function(x) pmax(x,1e-6)
plot(bottom(res1$padjust),bottom(res3$padjust),log="xy",pch=20)
table(tapply(res1$padjust<0.1,geneIDs(vaginal.model),any))
plotDEXSeq(vaginal.model,"1")
vaginal.model
design(vaginal.model)
vaginal.model@designColumns
plotMA(with(res1,data.frame(baseMean=meanBase,log2FoldChange=`log2fold(Unstable/Mid)`,padj=padjust)))

vaginal.model<-makeCompleteDEUAnalysis(vaginal.model)

plotDispEsts = function( cds, ymin, linecol="#ff000080",
                         xlab = "mean of normalized counts", ylab = "dispersion",
                         log = "xy", cex = 0.45, ... )
{
  px = rowMeans( counts( cds, normalized=TRUE ) )
  sel = (px>0)
  px = px[sel]
  
  py = fData(cds)$dispBeforeSharing[sel]
  if(missing(ymin))
    ymin = 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1)
  
  plot(px, pmax(py, ymin), xlab=xlab, ylab=ylab,
       log=log, pch=ifelse(py<ymin, 6, 16), cex=cex, ... )
  xg = 10^seq( -.5, 5, length.out=100 )
  fun = function(x) { cds@dispFitCoefs[1] + cds@dispFitCoefs[2] / x }
  lines( xg, fun(xg), col=linecol, lwd=4)
}

plotMA = function(x, ylim,
                  col = ifelse(x$padj>=0.1, "gray32", "red3"),
                  linecol = "#ff000080",
                  xlab = "mean of normalized counts", ylab = expression(log[2]~fold~change),
                  log = "x", cex=0.45, ...)
{
  if (!(is.data.frame(x) && all(c("baseMean", "log2FoldChange") %in% colnames(x))))
    stop("'x' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")
  
  x = subset(x, baseMean!=0)
  py = x$log2FoldChange
  if(missing(ylim))
    ylim = c(-1,1) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1
  plot(x$baseMean, pmax(ylim[1], pmin(ylim[2], py)),
       log=log, pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
       cex=cex, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
  abline(h=0, lwd=4, col=linecol)
}


