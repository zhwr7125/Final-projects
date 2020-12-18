

## Estimate normalization factors (we want to use the entire data set)


#counts.less<-counts[10001:16081,]
counts.less<-counts[10001:13081,]
norm.factors = estimate.norm.factors(counts.less);
## Prepare the data
nb.data = prepare.nb.data(counts.less, lib.sizes = colSums(counts.less), norm.factors = norm.factors);


## Specify the model matrix (experimental design)
#design<-model.matrix(~female+case,data=demo.dat)
design<-model.matrix(~female+case+PC1_MEGA+PC2_MEGA+PC3_MEGA+Ever.Smoked,data=demo.dat)
## Estimate dispersion model. 
dispersion3 = estimate.dispersion(nb.data, design);


print(dispersion3);
plot(dispersion3);
## Specify the null hypothesis
## The null hypothesis is beta[2]=0 (beta[2] is the log fold change).
beta0 = c(NA, NA, 0, NA, NA,NA, NA);
## Test regression coefficient
res3 = test.coefficient(nb.data, dispersion3, design, beta0);



list.save(dispersion3, 'dispersion3.rdata')
list.save(res3, 'res3.rdata')