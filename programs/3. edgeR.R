
#No need to exclude because this is from the exclude result from DESeq2
obj.edgeR<- as.DGEList(obj.deseq2)

#Calculate the normalized factor
obj.edgeR <- calcNormFactors(obj.edgeR,method="RLE")
obj.edgeR.sizefac<-obj.edgeR$samples$norm.factors
names(obj.edgeR.sizefac) <- colnames(counts)

#Calculate the dispersion
design<-model.matrix(~female+case+PC1_MEGA+PC2_MEGA+PC3_MEGA+Ever.Smoked,data=demo.dat)
obj.edgeR<-estimateDisp(obj.edgeR,design,robust=TRUE)
obj.edgeR.distp<-obj.edgeR$tagwise.dispersion
write.csv(obj.edgeR.distp,'D:/Course/BIOS7695/Final projects/data/obj.edgeR.distp.csv')

#Perform test
fit.edgeR<-glmFit(obj.edgeR,design)
LRT.fit.edgeR<-glmLRT(fit.edgeR, coef=3)

save(LRT.fit.edgeR, file='D:/Course/BIOS7695/Final projects/data/LRT.fit.edgeR.rdata')



