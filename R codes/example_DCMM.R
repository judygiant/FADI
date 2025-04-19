args=(commandArgs(TRUE))
d <- as.numeric(args[1]) #### Data dimension/number of variables
mc <- as.numeric(args[2]) #### Monte Carlo index
rt <- as.numeric(args[3]) #### rt = L*p / d
library(here)
#####Define fast power PCA method for aggregation of parallel sketching results
fast_pca_final<-function(y_ls,p0,qq){
  ####y_ls is the list of parallel SVD results, i.e., hat{V}^l in FADI algorithm
  ####p0 is the dimension of the power sketching
  ####qq is the number of iteration for the power method
  d <- dim(y_ls[[1]])[1]
  omega<-matrix(rnorm(d*p0),d,p0) ####Gassian sketching matrix
  
  for (q in 1:qq){
    Yt <- 0
    L <- length(y_ls)
    for (l in 1:L){
      Y <- t(y_ls[[l]])%*%omega
      Y <- y_ls[[l]] %*% Y
      Yt <- Yt + Y
    }
    Yt <- Yt / L
    omega <- Yt
  }
  
  
  Q<-svd(Yt)$u ##### output: final SVD estimator through powered sketching
  return(list('u'=Q))
}

##### Parameters setup for FADI
K =3
p = 4*K
p0 = 4*K ####p0 is the dimension of the power sketching
L = round(rt*d/p)
lmode = 1  ##### Determine if rt = L*p / d < 1 or rt = L*p / d >= 1. Different covariance estimation methods for the FADI estimator under two scenarios
q = 7  #### the number of iteration for the power method in the aggregation step


######### Generating the stacked community assignment probability vectors
set.seed(153)
Pi <- matrix(runif(K*K),K,K)
Pi <- Pi + diag(rep(10,K))
Pi <- Pi / rowSums(Pi)
Pi1 <- matrix(0,d,K)
nk <- round(d/K) 
for (k in 1:(K-1)){
  Pi1[(((k-1)*nk+1):((k-1)*nk+round(nk/3))),] <- rep(Pi[k,],each = round(nk/3))
  Pi1[(((k-1)*nk+round(nk/3)+ 1):(k*nk)),k] <- 1
  
}
Pi1[(((K-1)*nk+1):((K-1)*nk+round(nk/3))),] <- rep(Pi[K,],each = round(nk/3))
Pi1[(((K-1)*nk+round(nk/3)+1):d),K] <- 1

##### Generating connection probability matrix among clusters
P <- diag(rep(0.3,3)) + 0.1*c(2,1,1)%*%t(c(2,1,1)) + 0.1 * c(1,-1,-1)%*%t(c(1,-1,-1))
svd1 <- svd(Pi1)
mid <- diag(svd1$d) %*% t(svd1$v) %*% P %*% svd1$v %*% diag(svd1$d)
svd2 <- svd(mid)
###### True top K PCs
Vk <- svd1$u %*% svd2$u


M <- Pi1 %*% P %*% t(Pi1) ##### True d by d matrix of interest
set.seed(mc)
Mh <- matrix(rbinom(rep(1,d^2),1,c(M)),d,d) #### Observed adjacency matrix of the graph
Mh[lower.tri(Mh)] <- t(Mh)[lower.tri(Mh)]




runtime <- 0

y_ls <- list()
t_ls <- c()

###########distributed fast sketches
for (l in 1:L){
  ts <- Sys.time()
  ####### generate parallel gaussian sketchings
  omega <- matrix(rnorm(d*p),d,p)
  Y <- Mh%*%omega
  vl <- svd(Y)$u[,1:K]
  y_ls <- append(y_ls,list(vl))
  te <- Sys.time()
  t_ls <- c(t_ls,  difftime(te,ts,units="secs"))
}

runtime <- runtime + max(t_ls)
ts <- Sys.time()
###### Aggregating parallel PCA results
vtild <- fast_pca_final(y_ls, p0,q)$u[,1:K]

if (lmode == 1){
  ##### When rt = L*p / d < 1
  
  lmtild <- (t(vtild)%*%Mh%*%vtild)
  Mtild <- vtild%*%lmtild%*%t(vtild)
  
  Sig_hat1 <- solve(lmtild)%*%t(vtild)%*%( (    (Mtild[1,])*(1-Mtild[1,])   ) *vtild  )%*%solve(lmtild)
  Sig_hat2 <- solve(lmtild)%*%t(vtild)%*%( (   (Mtild[2,])*(1-Mtild[2,])   ) *vtild  )%*%solve(lmtild)
  Sig_hat <- Sig_hat1 + Sig_hat2
  }
te <- Sys.time()

#### Total runtime for FADI
runtime <- runtime +  difftime(te,ts,units="secs")
##### Final FADI estimator 
vkh <- svd(Mh)$u[,1:K]


dir.create(here("DCMM_Results"))
fname<-paste(c(here("DCMM_Results","results_"),args,".RData"),collapse = '_')
#### Output: 
### vtild: FADI estimator for top K PCs
### Sig_hat: covariance estimate for difference of row 1 and row 2 of FADI estimator
### Sig_hat1: covariance estimate for row 1 of FADI estimator
### Sig_hat2: covariance estimate for row 2 of FADI estimator
### vkh: traditional PCA estimator
### runtime: total runtime for each Monte Carlo

save(vtild,Sig_hat,Sig_hat1,Sig_hat2, vkh,runtime,file=fname)
