

obj.deseq2 <- ipf

#Exclude gene rows with counts>=10
detect.zero<-rowSums(matrix(counts(obj.deseq2) %in% seq(from=0,to=10,by=1),nrow=nrow(counts(obj.deseq2))))
table(detect.zero==0)
#KEEP detect.zero==0
filter10<-counts(obj.deseq2)[detect.zero==0,]
test<-rowSums(matrix(filter10 %in% seq(from=0,to=10,by=1),nrow=nrow(filter10)))
table(test==0)

#Exclude participants with missing ever smoke
keeplist<-rownames(demo.dat)[is.na(demo.dat$Ever.Smoked)==FALSE]
demo.dat<-demo.dat[is.na(demo.dat$Ever.Smoked)==FALSE,]
filter10<-filter10[,colnames(filter10) %in% keeplist]

#Convert to factors
obj.deseq2@colData@listData$female<-as.factor(obj.deseq2@colData@listData$female)
obj.deseq2@colData@listData$case<-as.factor(obj.deseq2@colData@listData$case)
obj.deseq2@colData@listData$Ever.Smoked<-as.factor(obj.deseq2@colData@listData$Ever.Smoked)

obj.deseq2 <- DESeqDataSetFromMatrix(countData = filter10,
                                     colData = DataFrame(demo.dat),
                                     design= ~ case)




#Calculate the normalized factor
obj.deseq2  <- estimateSizeFactors(obj.deseq2)
obj.deseq2.sizefac<-sizeFactors(obj.deseq2)


#Calculate the dispersion and base mean
obj.deseq2<-estimateDispersions(obj.deseq2,fitType = c("local"))
obj.deseq2.distp<-dispersions(obj.deseq2)
write.csv(obj.deseq2.distp,'D:/Course/BIOS7695/Final projects/data/obj.deseq2.distp.csv')

ave.deseq2.count<-rowMeans(counts/obj.deseq2.sizefac)
basemean.deseq2<-obj.deseq2@rowRanges@elementMetadata@listData$baseMean
  
#Perform test

design(obj.deseq2)<-formula(~female+case+PC1_MEGA+PC2_MEGA+PC3_MEGA+Ever.Smoked)
#design(obj.deseq2)<-formula(~female+case+PC1_MEGA+PC2_MEGA+PC3_MEGA)
#design(obj.deseq2)<-formula(~female+case)

#fit.deseq2 <- nbinomLRT(obj.deseq2,reduced=~female+PC1_MEGA+PC2_MEGA+PC3_MEGA+Ever.Smoked,maxit = 10000000)
fit.deseq2 <- nbinomWaldTest(obj.deseq2)


save(fit.deseq2, file='D:/Course/BIOS7695/Final projects/data/fit.deseq2.rdata')
#69 rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT


