library(EDASeq)


seq.counts = newSeqExpressionSet(counts=counts)


#Diagnose
#meanVarPlot(seq.counts,log=TRUE)

#Box plot and normalization
counts.between.median<-betweenLaneNormalization(seq.counts,which=c("median"))
counts.between.upper<-betweenLaneNormalization(seq.counts,which=c("upper"))
counts.between.full<-betweenLaneNormalization(seq.counts,which=c("full"))
#The result is not too bad

