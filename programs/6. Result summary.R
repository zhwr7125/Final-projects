

#Is gene name duplicated?
table(duplicated(rownames(counts)))
#So we can call them using gene name

#How many significane were detected?
dim(sig.deseq2)[1]
dim(sig.edgeR)[1]
dim(sig.NBPSeq)[1]
print(c("Sig in DESeq2",dim(sig.deseq2)[1]))
print(c("Sig in edgeR",dim(sig.edgeR)[1]))
print(c("Sig in NBPSeq",dim(sig.NBPSeq)[1]))

#Find out only
only.deseq2<-sig.deseq2[!(rownames(sig.deseq2) %in% rownames(sig.edgeR))&!(rownames(sig.deseq2) %in% rownames(sig.NBPSeq)),]
only.edgeR<-sig.edgeR[!(rownames(sig.edgeR) %in% rownames(sig.deseq2))&!(rownames(sig.edgeR) %in% rownames(sig.NBPSeq)),]
only.NBPSeq<-sig.NBPSeq[!(rownames(sig.NBPSeq) %in% rownames(sig.edgeR))&!(rownames(sig.NBPSeq) %in% rownames(sig.deseq2)),]
dim(only.deseq2)[1]
dim(only.edgeR)[1]
dim(only.NBPSeq)[1]
print(c("Only in DESeq2",dim(only.deseq2)[1]))
print(c("Only in edgeR",dim(only.edgeR)[1]))
print(c("Only in NBPSeq",dim(only.NBPSeq)[1]))

#Find out all
all.methods.name<-intersect(intersect(rownames(sig.deseq2),rownames(sig.edgeR)),rownames(sig.NBPSeq))
length(all.methods.name)
print(c("All three",length(all.methods.name)))

#Find out both in but one not in
both.deseq2.edgeR<-cbind(sig.deseq2[intersect(rownames(sig.deseq2),rownames(sig.edgeR)),],sig.edgeR[intersect(rownames(sig.deseq2),rownames(sig.edgeR)),])
both.edgeR.NBPSeq<-sig.edgeR[intersect(rownames(sig.NBPSeq),rownames(sig.edgeR)),]
both.NBPSeq.deseq2<-sig.NBPSeq[intersect(rownames(sig.NBPSeq),rownames(sig.deseq2)),]
dim(both.deseq2.edgeR)[1]-length(all.methods.name)
dim(both.edgeR.NBPSeq)[1]-length(all.methods.name)
dim(both.NBPSeq.deseq2)[1]-length(all.methods.name)
print(c("both.deseq2.edgeR but not NBP",dim(both.deseq2.edgeR)[1]-length(all.methods.name)))
print(c("both.edgeR.NBPSeq but not deseq2",dim(both.edgeR.NBPSeq)[1]-length(all.methods.name)))
print(c("both.NBPSeq.deseq2 but not edgeR",dim(both.NBPSeq.deseq2)[1]-length(all.methods.name)))

#test if correct
dim(only.deseq2)[1]+dim(both.deseq2.edgeR)[1]+dim(both.NBPSeq.deseq2)[1]-length(all.methods.name)==dim(sig.deseq2)[1]



#Estimated log fold change of gene counts
#Merge average gene count to sig.deseq2, sig.edgeR, and sig.NBPSeq
norm.factors = estimate.norm.factors(counts)
ave.gene.count<-rowMeans(counts/norm.factors)


sig.deseq2<-merge(sig.deseq2, ave.gene.count, by="row.names")   
sig.edgeR<-merge(sig.edgeR, ave.gene.count, by="row.names")   
sig.NBPSeq<-merge(sig.NBPSeq, ave.gene.count, by="row.names")   


background<-merge(results.edgeR, ave.gene.count, by="row.names")   
p1 <- ggplot(data=background, aes(x = y, y = logFC.edgeR)) + geom_point(color="grey",shape = 4) + 
  geom_point(data=sig.deseq2,aes( x = y, y = logFC.deseq2),color="magenta",shape = 1, alpha = 0.9)+
  geom_point(data=sig.edgeR,aes( x = y, y = logFC.edgeR),color="green",shape = 2, alpha = 0.5)+
  geom_point(data=sig.NBPSeq,aes( x = y, y = logFC.NBQ),color="cyan",shape = 3, alpha = 0.3)+
  scale_x_log2(breaks = trans_breaks("log2", function(x) 2^x),
               labels = trans_format("log2", math_format(2^.x))) +
  labs(y = "Estimated Log Fold Change of Gene Counts",x = "Average Gene Counts") 

background<-merge(results.deseq2, ave.gene.count, by="row.names") 
p2 <- ggplot(data=background, aes(x = y, y = logFC.deseq2)) + geom_point(color="grey",shape = 4) + 
  geom_point(data=sig.deseq2,aes( x = y, y = logFC.deseq2),color="magenta",shape = 1)+
#  geom_point(data=sig.edgeR,aes( x = y, y = logFC.edgeR),color="green",shape = 2, alpha = 0.1)+
#  geom_point(data=sig.NBPSeq,aes( x = y, y = logFC.NBQ),color="cyan",shape = 3, alpha = 0.1)+
  scale_x_log2(breaks = trans_breaks("log2", function(x) 2^x),
               labels = trans_format("log2", math_format(2^.x))) +
  labs(y = "Estimated Log Fold Change of Gene Counts",x = "Average Gene Counts") 


background<-merge(results.edgeR, ave.gene.count, by="row.names") 
p3 <- ggplot(data=background, aes(x = y, y = logFC.edgeR)) + geom_point(color="grey",shape = 4) + 
  # geom_point(data=sig.deseq2,aes( x = y, y = logFC.deseq2),color="magenta",shape = 1, alpha = 0.8)+
  geom_point(data=sig.edgeR,aes( x = y, y = logFC.edgeR),color="green",shape = 2)+
  #  geom_point(data=sig.NBPSeq,aes( x = y, y = logFC.NBQ),color="cyan",shape = 3, alpha = 0.1)+
  scale_x_log2(breaks = trans_breaks("log2", function(x) 2^x),
               labels = trans_format("log2", math_format(2^.x))) +
  labs(y = "Estimated Log Fold Change of Gene Counts",x = "Average Gene Counts") 


background<-merge(results.NBPSeq, ave.gene.count, by="row.names") 
p4 <- ggplot(data=background, aes(x = y, y = logFC.NBQ)) + geom_point(color="grey",shape = 4) + 
  # geom_point(data=sig.deseq2,aes( x = y, y = logFC.deseq2),color="magenta",shape = 1, alpha = 0.8)+
  # geom_point(data=sig.edgeR,aes( x = y, y = logFC.edgeR),color="green",shape = 2, alpha = 0.1)+
  geom_point(data=sig.NBPSeq,aes( x = y, y = logFC.NBQ),color="cyan",shape = 3)+
  scale_x_log2(breaks = trans_breaks("log2", function(x) 2^x),
               labels = trans_format("log2", math_format(2^.x))) +
  labs(y = "Estimated Log Fold Change of Gene Counts",x = "Average Gene Counts") 

