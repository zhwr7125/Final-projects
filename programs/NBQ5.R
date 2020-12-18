
## Estimate normalization factors (we want to use the entire data set)


#counts.less<-counts[13082:16081,]
counts.less<-counts[14082:16081,]


norm.factors = estimate.norm.factors(counts.less);
## Prepare the data
nb.data = prepare.nb.data(counts.less, lib.sizes = colSums(counts.less), norm.factors = norm.factors);


## Specify the model matrix (experimental design)
#design<-model.matrix(~case,data=demo.dat)
design<-model.matrix(~female+case+PC1_MEGA+PC2_MEGA+PC3_MEGA+Ever.Smoked,data=demo.dat)
## Estimate dispersion model. 
dispersion5 = estimate.dispersion(nb.data, design);

## Specify the null hypothesis
## The null hypothesis is beta[2]=0 (beta[2] is the log fold change).
#beta0 = c(NA, 0);
beta0 = c(NA, NA, 0, NA, NA,NA, NA);
## Test regression coefficient
res5 = test.coefficient(nb.data, dispersion5, design, beta0);


list.save(dispersion5, 'D:/Course/BIOS7695/Final projects/data/dispersion5.rdata')
list.save(res5, 'D:/Course/BIOS7695/Final projects/data/res5.rdata')