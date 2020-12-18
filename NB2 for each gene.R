source("D:/Course/BIOS7695/Final projects/programs/1. Demographics.R")
library(Hmisc)
library(rlist)
library(DESeq2)
library(DEFormats)
library(EDASeq)
library(edgeR)
library(MASS)
library(gnlm)
library(NBPSeq)
library(dplyr)
library(jrnoldmisc)
library(scales)
library(gnlm)

norm.factors = estimate.norm.factors(counts)
dat<-demo.dat
dat$female<-as.character(dat$female)
dat$case<-as.character(dat$case)
dat$case[dat$case=="Case"]<-1
dat$case[dat$case=="Control"]<-0
dat$case<-as.numeric(dat$case)
dat$female<-as.numeric(dat$female)
dat$Ever.Smoked<-as.numeric(dat$Ever.Smoked)-1
dat$PC1_MEGA<-scale(dat$PC1_MEGA)
dat$PC2_MEGA<-scale(dat$PC2_MEGA)
dat$PC3_MEGA<-scale(dat$PC3_MEGA)
log2norm<-log2(norm.factors)

attach(dat)

haha<-function(x){
  y<-x
  nbmodel<-gnlr(y=y,
                distribution = "negative binomial",
                mu = ~ 2^(b0 + b1*female+b2*case+b3*PC1_MEGA+b4*PC2_MEGA+b5*PC3_MEGA+b6*Ever.Smoked + log2norm),
                pmu=list(b0=100,b1=1,b2=1,b3=1,b4=1,b5=1,b6=1), pshape=2
  )
#  if (sum(is.na(nbmodel$corr))!=0){
#    NB2glm.disp<-NA
#  }else{
    NB2glm.disp<-exp(nbmodel$coefficients[8])*(sum(is.na(nbmodel$corr))==0)
#  }
  NB2glm.disp<-1/NB2glm.disp
  NB2glm.disp
}
NB2glm.disp<-apply(counts,1,FUN=haha)
NB2glm.disp2<-apply(counts[rowMeans(counts/norm.factors)>2^12,],1,FUN=haha)

write.csv(NB2glm.disp,'NB2glm.disp.csv')
write.csv(NB2glm.disp2,'NB2glm.disp.csv')
