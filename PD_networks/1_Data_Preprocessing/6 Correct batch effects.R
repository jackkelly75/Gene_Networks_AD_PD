library(sva)

############################
#correcting all data       #
############################
pheno = pData
edata = exprs(data.rma.norm)
mod = model.matrix(~as.factor(`disease label:ch1`), data=pheno)

# reference-batch version, with covariates
batch = pheno$`batch:ch1`
data.rma.ComBat = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
batch = pheno$final_gender
data.rma.ComBat = ComBat(dat=data.rma.ComBat, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
