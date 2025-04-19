args=(commandArgs(TRUE))
i <- as.numeric(args[1]) ### index for data split
l <- as.numeric(args[2]) ### index for parallel sketching
p <- as.numeric(args[3]) ### sketching dimension for FADI
library(gdsfmt)
library(here)
gfile <- openfn.gds("/users/shutingshen/fast_pca/real_data/1KG_pruned_forPCA.gds") #### 1000 Genome data, publicly available at https://www.internationalgenome.org/data/
geno <- index.gdsn(gfile, "genotype")

d <- 2504
N <- 168047
n <- round(N/100)
ni <- if (i < 99) n else (N - 99*n)
xi <- read.gdsn(geno)
ind <- apply(xi,2, function(x){any(x>2 & x <0)})
xi<-xi[,!ind]
xi <- scale(xi)
ts <- Sys.time()
xi <- xi[,(((i-1)*n + 1):( (i-1)*n + ni) )]
####Step 1: compute parallel sketches corresponding to each data split
set.seed(l)
omega <- matrix(rnorm(d*p),d , p)
Y <- t(xi)%*%omega
Y <- xi %*% Y/N
te <- Sys.time()
runtime <- as.numeric(te - ts, unit = "secs")
##### Saving parallel sketches corresponding to each data split
dir.create(here("1000g_est_Results"))
fname<-paste(c(here("1000g_est_Results","results_"),args,".RData"),collapse = '_')
save(Y, runtime,file=fname)