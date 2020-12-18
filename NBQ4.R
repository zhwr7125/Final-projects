
## Estimate normalization factors (we want to use the entire data set)


counts.less<-counts[c(seq(3001,3500),seq(9501,10000),c(seq(13082,14081))),]
#counts.less<-counts[c(seq(13082,14081)),]
norm.factors = estimate.norm.factors(counts.less);
## Prepare the data
nb.data = prepare.nb.data(counts.less, lib.sizes = colSums(counts.less), norm.factors = norm.factors);


## Specify the model matrix (experimental design)
design<-model.matrix(~case,data=demo.dat)
#design<-model.matrix(~female+case+PC1_MEGA+PC2_MEGA+PC3_MEGA+Ever.Smoked,data=demo.dat)
## Estimate dispersion model. 
dispersion4 = estimate.dispersion(nb.data, design);


print(dispersion4);
plot(dispersion4);
## Specify the null hypothesis
## The null hypothesis is beta[2]=0 (beta[2] is the log fold change).
beta0 = c(NA, 0);
## Test regression coefficient
res4 = test.coefficient(nb.data, dispersion4, design, beta0);


list.save(dispersion4, 'dispersion4.rdata')
list.save(res4, 'res4.rdata')