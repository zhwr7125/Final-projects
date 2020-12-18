library(NBPSeq)

## Estimate normalization factors (we want to use the entire data set)


counts.less<-counts[1:3000,]
norm.factors = estimate.norm.factors(counts.less);
## Prepare the data
nb.data = prepare.nb.data(counts.less, lib.sizes = colSums(counts.less), norm.factors = norm.factors);


## Specify the model matrix (experimental design)

design<-model.matrix(~female+case+PC1_MEGA+PC2_MEGA+PC3_MEGA+Ever.Smoked,data=demo.dat)
## Estimate dispersion model. 
dispersion1 = estimate.dispersion(nb.data, design);



## Specify the null hypothesis
## The null hypothesis is beta[2]=0 (beta[2] is the log fold change).
beta0 = c(NA, NA, 0, NA, NA,NA, NA);
## Test regression coefficient
res1 = test.coefficient(nb.data, dispersion1, design, beta0);

list.save(dispersion1, 'D:/Course/BIOS7695/Final projects/data/dispersion1.rdata')
list.save(res1, 'D:/Course/BIOS7695/Final projects/data/res1.rdata')

