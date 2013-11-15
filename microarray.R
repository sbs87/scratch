source("http://bioconductor.org/biocLite.R")
biocLite()
library(affy)
library(limma)

setwd("/Users/stsmith/Documents/Data/microarray")
celfile.list<-read.table("cel_file_list",stringsAsFactors=FALSE)$V1
celfiles<-ReadAffy(filenames=celfile.list)
getwd()


source("http://bioconductor.org/workflows.R")
workflowInstall("arrays")
meta<-read.table("Sample_Meta_Table",header=T,sep="\t")
## import 'phenotype' data, describing the experimental design
#phenoData <- read.AnnotatedDataFrame(system.file("extdata", "pdata.txt", package = "arrays"))

meta<-meta[meta$File_Name!="missing"& meta$Group=="Control",]
pd<-data.frame(group=meta$Group,ID=meta$ID,BV_Stat=meta$Nugent,fn=meta$File_Name)

rownames(pd)<-meta$File_Name
vl<-data.frame(labelDescription=c("Group","blah","blah","blah"))
phenoData<-new("AnnotatedDataFrame",data=pd,varMetadata=vl)
## RMA normalization
#celfiles <- system.file("extdata", package = "arrays")
eset <- justRMA(phenoData = phenoData, celfile.path = "/Users/stsmith/Documents/Data/microarray/cel_files")

## differential expression
group<-pData(phenoData)[, 1]
status<-pData(phenoData)[, 3]
combn <- factor(paste(pData(phenoData)[, 1], pData(phenoData)[, 3], sep = "_"))
design <- model.matrix(~status)  # describe model to be fit

fit <- lmFit(y, design)  # fit each probeset to model
efit <- eBayes(fit)  # empirical Bayes adjustment
topTable(efit, coef = 2)  # table of differentially expressed probesets
?topTable
y<-log2(exprs(eset))
boxplot(y,col=as.numeric(pData(eset)$BV_Stat))

Index<-as.numeric(pData(eset)$BV_Stat)
d<-rowMeans(y[,Index==2])-rowMeans(y[,Index==1])
a<-rowMeans(y)
smoothScatter(a,d,ylim=c(1,-1))
abline(h=c(-0.5,0.5),col="red")
colnames(tab)
dimnames(tab)
tab<-topTable(efit,coef=2,adjust="fdr",n=150)
labCol<-c("BV","Intermediate","Normal")[Index]

heatmap(y[dimnames(tab)[[1]],],labCol=labCol)
library("annotate")
library("hgu133plus2.db")