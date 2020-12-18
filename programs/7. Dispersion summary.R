##Demo
source("D:/Course/BIOS7695/Final projects/programs/1. Demographics.R")


##################DESeq#####################################################

#source("D:/Course/BIOS7695/Final projects/programs/2. DESeq2.R")
obj.deseq2.distp<-read.csv('D:/Course/BIOS7695/Final projects/data/obj.deseq2.distp.csv')
obj.deseq2.distp<-obj.deseq2.distp[,2]
load("D:/Course/BIOS7695/Final projects/data/fit.deseq2.rdata")
results.deseq2<-results(fit.deseq2, contrast=c("case","Case","Control"))

df.res.deseq2<-data.frame(results.deseq2$log2FoldChange,results.deseq2$pvalue,results.deseq2$padj,row.names = rownames(results(fit.deseq2)))
colnames(df.res.deseq2)<-c("logFC.deseq2","pvalue.deseq2","padj.deseq2")
sig.deseq2<-df.res.deseq2[c(df.res.deseq2$padj<0.05),]
results.deseq2<-df.res.deseq2

##################edgeR#####################################################

#source("D:/Course/BIOS7695/Final projects/programs/3. edgeR.R")
obj.edgeR.distp<-read.csv('D:/Course/BIOS7695/Final projects/data/obj.edgeR.distp.csv')
obj.edgeR.distp<-obj.edgeR.distp[,2]

load("D:/Course/BIOS7695/Final projects/data/LRT.fit.edgeR.rdata")
results.edgeR<-LRT.fit.edgeR$table

results.edgeR$padj<-p.adjust(LRT.fit.edgeR$table$PValue,"BH")
results.edgeR$logCPM<-NULL
results.edgeR$LR<-NULL
colnames(results.edgeR)<-c("logFC.edgeR","pvalue.edgeR","padj.edgeR")

sig.edgeR<-results.edgeR[c(results.edgeR$padj<0.05),]


source("D:/Course/BIOS7695/Final projects/programs/4. Normalization.R")
source("D:/Course/BIOS7695/Final projects/programs/5. NBQsummary.R")
source("D:/Course/BIOS7695/Final projects/programs/6. Result summary.R")

NB2glm.disp<-read.csv('D:/Course/BIOS7695/Final projects/data/NB2glm.disp.csv')
NB2glm.disp<-NB2glm.disp[,2]
NB2glm.disp[NB2glm.disp=="Inf"]<-NA

NB2glm.disp2<-read.csv('D:/Course/BIOS7695/Final projects/data/NB2glm.disp2.csv')
NB2glm.disp2<-NB2glm.disp2[,2]
NB2glm.disp2[NB2glm.disp2=="Inf"]<-NA

NB2glm.disp[ave.gene.count>2^12]<-NB2glm.disp2

#Plot dispersion vs mean
ave.gene.count<-rowMeans(counts/norm.factors)


# x and y axis are transformed and formatted


#dispersion from NB2 for each gene

p5 <- ggplot(data=NULL, aes(x = ave.gene.count, y = NB2glm.disp)) + geom_point(shape = 1, alpha = 0.5) +
  scale_x_log2(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_y_log2(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x)))+
  labs(x = "Mean of normalized gene count",y = "Estimate of NB2 Dispersion")


#dispersion from DESeq2
p6 <- ggplot(data=NULL, aes(x = ave.gene.count, y = obj.deseq2.distp)) + geom_point(shape = 1, alpha = 0.5) +
  scale_x_log2(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_y_log2(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x)))+
  labs(x = "Mean of normalized gene count",y = "Estimate of NB2 Dispersion") 

p7 <- ggplot(data=NULL, aes(x = ave.gene.count, y = obj.edgeR.distp)) + geom_point(shape = 1, alpha = 0.5) +
  scale_x_log2(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_y_log2(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x)))+
  labs(x = "Mean of normalized gene count",y = "Estimate of NB2 Dispersion")

p8 <- ggplot(data=NULL, aes(x = unlist(obj.NBPSeq.avect[seq(1,3000),]), y = unlist(obj.NBPSeq.distp[seq(1,3000),]))) + geom_point(shape = 1, size=2, alpha = 0.3) +
  scale_x_log2(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_y_log2(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x)))+
  labs(x = "Mean of normalized gene count",y = "Estimate of NB2 Dispersion")

#p9 <- ggplot(data=NULL, aes(x = unlist(obj.NBPSeq.avect[seq(3501,9500),]), y = unlist(obj.NBPSeq.distp[seq(3501,9500),]))) + geom_point(shape = 1, size=2, alpha = 0.3) +
p9 <- ggplot(data=NULL, aes(x = unlist(obj.NBPSeq.avect[seq(3001,9000),]), y = unlist(obj.NBPSeq.distp[seq(3001,9000),]))) + geom_point(shape = 1, size=2, alpha = 0.3) +
  
    scale_x_log2(breaks = trans_breaks("log2", function(x) 2^x),
               labels = trans_format("log2", math_format(2^.x))) +
  scale_y_log2(breaks = trans_breaks("log2", function(x) 2^x),
               labels = trans_format("log2", math_format(2^.x)))+
  labs(x = "Mean of normalized gene count",y = "Estimate of NB2 Dispersion")

#p10 <- ggplot(data=NULL, aes(x = unlist(obj.NBPSeq.avect[seq(10001,13081),]), y = unlist(obj.NBPSeq.distp[seq(10001,13081),]))) + geom_point(shape = 1, size=2, alpha = 0.3) +
p10 <- ggplot(data=NULL, aes(x = unlist(obj.NBPSeq.avect[seq(9001,12081),]), y = unlist(obj.NBPSeq.distp[seq(9001,12081),]))) + geom_point(shape = 1, size=2, alpha = 0.3) +
    scale_x_log2(breaks = trans_breaks("log2", function(x) 2^x),
               labels = trans_format("log2", math_format(2^.x))) +
  scale_y_log2(breaks = trans_breaks("log2", function(x) 2^x),
               labels = trans_format("log2", math_format(2^.x)))+
  labs(x = "Mean of normalized gene count",y = "Estimate of NB2 Dispersion")

#p11 <- ggplot(data=NULL, aes(x = unlist(obj.NBPSeq.avect[c(seq(3001,3500),seq(9501,10000),seq(13082,14081)),]), y = unlist(obj.NBPSeq.distp[c(seq(3001,3500),seq(9501,10000),seq(13082,14081)),]))) + geom_point(shape = 1, size=2, alpha = 0.3) +
p11 <- ggplot(data=NULL, aes(x = unlist(obj.NBPSeq.avect[seq(12082,14081),]), y = unlist(obj.NBPSeq.distp[seq(12082,14081),]))) + geom_point(shape = 1, size=2, alpha = 0.3) +
  scale_x_log2(breaks = trans_breaks("log2", function(x) 2^x),
               labels = trans_format("log2", math_format(2^.x))) +
  scale_y_log2(breaks = trans_breaks("log2", function(x) 2^x),
               labels = trans_format("log2", math_format(2^.x)))+
  labs(x = "Fitted mean of gene count under null model",y = "Estimate of NB2 Dispersion")

p12 <- ggplot(data=NULL, aes(x = unlist(obj.NBPSeq.avect[seq(14082,16081),]), y = unlist(obj.NBPSeq.distp[seq(14082,16081),]))) + geom_point(shape = 1, size=2, alpha = 0.3) +
  scale_x_log2(breaks = trans_breaks("log2", function(x) 2^x),
               labels = trans_format("log2", math_format(2^.x))) +
  scale_y_log2(breaks = trans_breaks("log2", function(x) 2^x),
               labels = trans_format("log2", math_format(2^.x)))+
  labs(x = "Fitted mean of gene count under null model",y = "Estimate of NB2 Dispersion")
