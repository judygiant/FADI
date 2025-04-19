args=(commandArgs(TRUE))
l <- as.numeric(args[1]) ### index for parallel sketching
p <- as.numeric(args[2]) ### sketching dimension for FADI
d <- 2504
N <- 168047
mu= (d/(N*p)^{.5}*log(d))^{0.75}/12 #### treshold parameter for estimating the rank K (see Section 3.5 of paper)
library(here)
###### Step 2 aggregating sketches across distributed splits
Yt <- 0 
rt <- c()
for (i in 1:100) {
  fname <- paste(c(here("1000g_est_Results","results_"),c(i,l,50),".RData"),collapse = '_')
  load(fname)
  Yt <- Yt+Y
  rt <- c(rt, runtime)
}

rt <- max(rt)
ts<- Sys.time()
set.seed(l)
omega <- matrix(rnorm(d*p),d , p)

Yt <- Yt - 0.730152*omega

##### Compute parallel PCA
vk_hat <- prcomp(Yt)$x[,1:25]
te <- Sys.time()

rt <- rt + te - ts

##### Local estimation of K
dd <- svd(Yt/sqrt(p))$d
diff <- dd-dd[p]
k <- 1
flag <- FALSE
while (!flag) {
  flag <- all(diff[k:p] <= mu)
  k<-k+1
}
kl<-k-2

fname<-paste(c(here("1000g_est_Results","results_collect1"),args,".RData"),collapse = '_')
### vk_hat: top K eigenspace estimate by parallel sketching l
### rt: total runtime
### kl: rank estimate by parallel sketching l
save(vk_hat,rt,kl,file=fname)
