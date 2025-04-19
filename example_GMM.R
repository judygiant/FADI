args=(commandArgs(TRUE))

d <- as.numeric(args[1])  #### Data dimension/number of variables
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
K = 3
m = 20
ni = 1000 ### sample size for each Gaussian distribution
n = ni*m
Delta2 = n^{2/3} ### Scaling parameter of mean matrix


p = 4*K
p0 = 4*K ####p0 is the dimension of the power sketching
L = round(rt*d/p)
q = 7  #### the number of iteration for the power method in the aggregation step
lmode <- rt < 1 ##### Determine if rt = L*p / d < 1 or rt = L*p / d >= 1. Different covariance estimation methods for the FADI estimator under two scenarios

##### Generating means for K Gaussian distributions
set.seed(152)
Theta = matrix(rnorm(n*K), n, K)
Theta <- Theta * (Delta2/2/n)^{0.5}

###### Generating membership labels
F_st <- matrix(0,d,K)
nk <- round(d/K) #### sample size of each distribution
for (k in 1:(K-1)){
  F_st[(((k-1)*nk+1):(k*nk)),k] <- 1
}
F_st[(((K-1)*nk+1):d),K] <- 1

svd1 <- svd(F_st)

mid <- t(Theta) %*% Theta 
svd2 <- svd(mid)

Vk <- svd1$u %*% svd2$u ##### True top K eigenspace

EX <- Theta %*% t(F_st) #### stacking of Gaussian means

#### Generate Gaussian data from K distributions
set.seed(mc+10086)

X <- matrix(rnorm(n*d),n,d) + EX

runtime <- 0




omega_ls <- list()
yl_ls <- list()
y_ls <- list()

tL_ls <- c()

###########distributed fast sketches
for (l in 1:L){
  ####### generate parallel gaussian sketchings
  omega <- matrix(rnorm(d*p),d,p)
  omega_ls <- append(omega_ls,list(omega))
  Y<- 0 
  tm_ls <- c()
  ####### Step 1: distributed fast sketching on local machines
  for(i in 1:m){
    ts <- Sys.time()
    xi <- X[(((i-1)*ni+1):(i*ni)),]
    y <- t(xi)%*%(xi%*%omega)
    Y <- Y+y
    te <- Sys.time()
    tm_ls <- c(tm_ls,difftime(te,ts,units="secs"))
  }
  tm <- max(tm_ls)
  ts <- Sys.time()
  ###### Step 2: parallel PCA for each gaussian sketching
  Y <- Y-n*omega
  yl_ls <- append(yl_ls,list(Y))
  vl <- svd(Y)$u[,1:K]
  y_ls <- append(y_ls,list(vl))
  te <- Sys.time()
  tL_ls <- c(tL_ls, tm + difftime(te,ts,units="secs"))
}

tL <- max(tL_ls)
runtime <- runtime + tL
ts <- Sys.time()
###### Aggregating parallel PCA results
vtild <- fast_pca_final(y_ls, p0,q)$u[,1:K]
te <- Sys.time()
runtime <- runtime + difftime(te,ts,units="secs")



########estimation of covariance

if(lmode == 1){
  ##### When rt = L*p / d < 1
  Bo <- c()
  Ysig <- c()
  Omega <- c()
  ts <- Sys.time()
  for (l in 1:L){
    Bl <- t(yl_ls[[l]])%*%vtild/sqrt(p)/n
    svdl <- svd(Bl)
    Bl <- svdl$u %*% diag((svdl$d)^{-1}) %*% t(svdl$v)/n
    Bo <- rbind(Bo, Bl)
    Ysig <- cbind(Ysig, (n*omega_ls[[l]]+yl_ls[[l]])/sqrt(p))
    Omega <- cbind(Omega, omega_ls[[l]]/sqrt(p))
  }
  Sig_hat <- t(Omega)%*%Ysig
  Sig_hat <- t(Bo)%*%Sig_hat %*% Bo *(d*n^{1/3}/L)
  sc <- (L*d*n^{1/3})^{0.5}
  te <- Sys.time()
  tsig <- difftime(te,ts,units="secs")
}else{
  ##### When rt = L*p / d >= 1
  lmtilde <- 0
  tsig <- c()
  for (i in 1:m){
    ts <- Sys.time()
    xi <- X[(((i-1)*ni+1):(i*ni)),]
    lmtilde <- lmtilde + t(vtild)%*% t(xi) %*% (xi %*% vtild)
    te <- Sys.time()
    tsig <- c(tsig, difftime(te,ts,units="secs"))
  }  
  tsig <- max(tsig)
  
  ts <- Sys.time()
  lmtilde <- lmtilde - n*diag(rep(1,K))
  lminv <- solve(lmtilde)
  Sig_hat <- (lminv+n*lminv %*% lminv)*(d^2*n^{1/3})
  sc <- d*n^{1/6}
  te <- Sys.time()
  tsig <- tsig + difftime(te,ts,units="secs")
}
runtime <- runtime + tsig

dir.create(here("GMM_Results"))
fname<-paste(c(here("GMM_Results","results_"),args,".RData"),collapse = '_')
#### Output: 
### vtild: FADI estimator for top K PCs
### Sig_hat: covariance estimate for FADI estimator
### sc: scaling of estimation
### runtime: total runtime for each Monte Carlo
save(vtild,Sig_hat, sc, runtime,file=fname)


