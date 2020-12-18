
#source("D:/Course/BIOS7695/Final projects/programs/NBQ1.R")
#source("D:/Course/BIOS7695/Final projects/programs/NBQ2.R")
#source("D:/Course/BIOS7695/Final projects/programs/NBQ3.R")
#source("D:/Course/BIOS7695/Final projects/programs/NBQ4.R")
#source("D:/Course/BIOS7695/Final projects/programs/NBQ5.R")

load("D:/Course/BIOS7695/Final projects/programs/dispersion1.rdata")
dispersion1<-x
load("D:/Course/BIOS7695/Final projects/programs/dispersion2.RData")
dispersion2<-x
load("D:/Course/BIOS7695/Final projects/programs/dispersion3.RData")
dispersion3<-x
load("D:/Course/BIOS7695/Final projects/programs/dispersion4.RData")
dispersion4<-x
load("D:/Course/BIOS7695/Final projects/programs/dispersion5.RData")
dispersion5<-x

obj.NBPSeq.distp<-rbind(
  data.frame(#id=seq(1,3000),
             dispersion1$estimates),
  data.frame(#id=seq(3501,9500),
             dispersion2$estimates),
  data.frame(#id=seq(10001,13081),
             dispersion3$estimates),
  data.frame(#id=c(seq(3001,3500),seq(9501,10000),seq(13082,14081)),
             dispersion4$estimates),
  data.frame(#id=seq(14082,16081),
             dispersion5$estimates)
)
#obj.NBPSeq.distp <- obj.NBPSeq.distp[order(obj.NBPSeq.distp$id),]

#Summary the results
load("D:/Course/BIOS7695/Final projects/programs/res1.rdata")
res1<-x
gene.id<-seq(1,3000)
df.res.NBQ1<-data.frame(gene.id, rownames(counts[seq(1,3000),]), res1$beta.hat[,3]/log(2),res1$LR$p.values)
colnames(df.res.NBQ1)<-c("Gene id","Gene name","logFC.NBQ","pvalue.NBQ")

load("D:/Course/BIOS7695/Final projects/programs/res2.rdata")
res2<-x
gene.id<-seq(3501,9500)
df.res.NBQ2<-data.frame(gene.id, rownames(counts[seq(3501,9500),]), res2$beta.hat[,3]/log(2),res2$LR$p.values)
colnames(df.res.NBQ2)<-c("Gene id","Gene name","logFC.NBQ","pvalue.NBQ")

load("D:/Course/BIOS7695/Final projects/programs/res3.rdata")
res3<-x
gene.id<-seq(10001,13081)
df.res.NBQ3<-data.frame(gene.id, rownames(counts[seq(10001,13081),]), res3$beta.hat[,3]/log(2),res3$LR$p.values)
colnames(df.res.NBQ3)<-c("Gene id","Gene name","logFC.NBQ","pvalue.NBQ")

load("D:/Course/BIOS7695/Final projects/programs/res4.rdata")
res4<-x
gene.id<-c(seq(3001,3500),seq(9501,10000),seq(13082,14081))
df.res.NBQ4<-data.frame(gene.id, rownames(counts[c(seq(3001,3500),seq(9501,10000),c(seq(13082,14081))),]), res4$beta.hat[,2]/log(2),res4$LR$p.values)
colnames(df.res.NBQ4)<-c("Gene id","Gene name","logFC.NBQ","pvalue.NBQ")

load("D:/Course/BIOS7695/Final projects/programs/res5.rdata")
res5<-x
gene.id<-seq(14082,16081)
df.res.NBQ5<-data.frame(gene.id, rownames(counts[seq(14082,16081),]), res5$beta.hat[,3]/log(2),res5$LR$p.values)
colnames(df.res.NBQ5)<-c("Gene id","Gene name","logFC.NBQ","pvalue.NBQ")

obj.NBPSeq.avect<-rbind(
  data.frame(#id=seq(1,3000),
             res1$mu.tilde),
  data.frame(#id=seq(3501,9500),
             res2$mu.tilde),
  data.frame(#id=seq(10001,13081),
             res3$mu.tilde),
  data.frame(#id=c(seq(3001,3500),seq(9501,10000),seq(13082,14081)),
             res4$mu.tilde),
  data.frame(#id=seq(14082,16081),
             res5$mu.tilde)
)
#obj.NBPSeq.avect <- obj.NBPSeq.avect[order(obj.NBPSeq.avect$id),]

results.NBPSeq<-rbind(df.res.NBQ1,df.res.NBQ2,df.res.NBQ3,df.res.NBQ4,df.res.NBQ5)
#results.NBPSeq <- results.NBPSeq[order(results.NBPSeq$"Gene id"),]

results.NBPSeq$padj.NBPSeq<-p.adjust(results.NBPSeq$pvalue.NBQ, method = "fdr")
rownames(results.NBPSeq)<-results.NBPSeq$`Gene name`
results.NBPSeq$`Gene name`<-NULL
results.NBPSeq$`Gene id`<-NULL
sig.NBPSeq<-results.NBPSeq[c(results.NBPSeq$padj<0.05),]