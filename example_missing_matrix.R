args=(commandArgs(TRUE))
d <- as.numeric(args[1]) #### Data dimension/number of variables
mc <- as.numeric(args[2]) #### Monte Carlo index
rt <- as.numeric(args[3]) #### rt = L*p / d
library(here)
m <- 10 ### number of distributed splits
dj <- d/m ### dimension of each split

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
det =2 ### Eigen gap
lambda = det*(K:1) ### top K eigenvalues
p = 4*K
p0 = 4*K  ####p0 is the dimension of the power sketching
L = round(rt*d/p) 
lmode = (rt >= 1) ##### Determine if rt = L*p / d < 1 or rt = L*p / d >= 1. Different covariance estimation methods for the FADI estimator under two scenarios
q = 7 #### the number of iteration for the power method in the aggregation step
theta = 0.4 #### missing probability
sig = d^{-1}*det*4 ##### standard deviation of noise


##### Generating true matrix of interest
set.seed(154)
u <- matrix(rnorm(d*K), d,K)
u <- svd(u)$u
M <- u%*%diag(lambda)%*%t(u)


set.seed(mc)
#### adding entry-wise noise
E <- matrix(rnorm(d^2)*sig,d,d)
Mh <- M+E
#### only observe entries not missing
obv_loc <- rbinom(d^2, 1, 1-theta) == 1
Mh[obv_loc] <- 0
Mh[lower.tri(Mh)] <- t(Mh)[lower.tri(Mh)]

runtime <- 0

yl_ls <- list()
y_ls <- list()
t_ls <- c()
omega_ls <- list()

###########distributed fast sketches
for (l in 1:L){
  tt_ls <- c()
  Y <- 0
  omega <- c()
  for (j in 1:m){
    ts <- Sys.time()
    omegaj <- matrix(rnorm(d*p/m),d/m,p)   
    loc <- (dj*(j-1)+1):(dj*j)
    Yj <- Mh[,loc]%*%omegaj
    Y <- Y + Yj
    te <- Sys.time()
    tt_ls <- c(tt_ls,  difftime(te,ts,units="secs"))
    omega <- rbind(omega,omegaj) 
  }
  
  ts <- Sys.time()
  vl <- svd(Y)$u[,1:K]
  te <- Sys.time()
  y_ls <- append(y_ls,list(vl))
  omega_ls <- append(omega_ls,list(omega))
  yl_ls <- append(yl_ls,list(Y))
  t_ls <- c(t_ls,  difftime(te,ts,units="secs") + max(tt_ls))
  
}

runtime <- runtime + max(t_ls)
ts <- Sys.time()
###### Aggregating parallel PCA results
vtild <- fast_pca_final(y_ls, p0,q)$u[,1:K]
Ss <- sum(obv_loc & !lower.tri(Mh))
that <- 1-2*Ss/d/(d+1) #### missing probability estimate
lmtild <- (t(vtild)%*%Mh%*%vtild)/that ##### Estimating top K eigenvalues
Mtild <- vtild%*%lmtild%*%t(vtild) #### Estimating true matrix of interest by FADI estimator
sig2hat <- mean(((Mh - Mtild)[!obv_loc & !lower.tri(Mh)])^2) #### noise variance estimate


########estimation of covariance
if (lmode == 1){
  ##### When rt = L*p / d >= 1
  Sig_hat1 <- solve(lmtild)%*%t(vtild)%*%(((Mtild[1,]^2*(1-that) + sig2hat    )/that)*vtild)%*%solve(lmtild)
  Sig_hatd <- solve(lmtild)%*%t(vtild)%*%((( Mtild[d,]^2*(1-that) + sig2hat    )/that)*vtild)%*%solve(lmtild)
  Sig_hat2 <- solve(lmtild)%*%t(vtild)%*%((( Mtild[2,]^2*(1-that) + sig2hat    )/that)*vtild)%*%solve(lmtild)
  Sig_hat <- Sig_hat1 + Sig_hat2
}else{
  ##### When rt = L*p / d < 1
  Bo <- c()
  Omega <- c()
  for (l in 1:L){
    Bl <- t(yl_ls[[l]])%*%vtild/sqrt(p)
    svdl <- svd(Bl)
    Bl <- svdl$u %*% diag((svdl$d)^{-1}) %*% t(svdl$v) ### Estimating the matrix B^{(l)} defined below Assumption 2 of paper 
    Bo <- rbind(Bo, Bl) ### Estimating the matrix B_{Omega}= (B^{(1)}, ..., B^{(L)}) defined below Assumption 2 of paper 
    Omega <- cbind(Omega, omega_ls[[l]]/sqrt(p))
  }
  
  Sig_hat1 <-  t(Bo)%*%t(Omega)%*%(((Mtild[1,]^2*(1-that) + sig2hat   )/that)*Omega)%*%Bo/L^2
  Sig_hatd <- t(Bo)%*%t(Omega)%*%((( Mtild[d,]^2*(1-that) + sig2hat    )/that)*Omega)%*%Bo/L^2
  Sig_hat2 <- t(Bo)%*%t(Omega)%*%((( Mtild[2,]^2*(1-that) + sig2hat    )/that)*Omega)%*%Bo/L^2
  Sig_hat <- Sig_hat1 + Sig_hat2
  
}
te <- Sys.time()

runtime <- runtime +  difftime(te,ts,units="secs")

dir.create(here("missing_mat_Results"))
fname<-paste(c(here("missing_mat_Results","results_"),args,".RData"),collapse = '_')
#### Output: 
### vtild: FADI estimator for top K PCs
### lmtild: estimate of top eigenvalues
### Mtild: estimate of true matrix of interest
### Sig_hat1: covariance estimate for difference row 1 of FADI estimator
### Sig_hat2: covariance estimate for difference row 2 of FADI estimator
### Sig_hatd: covariance estimate for difference row d of FADI estimator
### Sig_hat: covariance estimate for difference between row 1 and row 2 of FADI estimator
### sig2hat: estimate of noise variance
### that: estimate of missing probability
### runtime: total runtime for each Monte Carlo
### Estimating the matrix B_{Omega}= (B^{(1)}, ..., B^{(L)}) defined below Assumption 2 of paper 
### Omega: stacking of scaled Gaussian sketching matrices
if (lmode == 1){
  save(vtild,lmtild,Mtild,Sig_hat, Sig_hat1, Sig_hatd, sig2hat,that, runtime,file=fname)
}else{
  save(vtild,Bo,Omega,Mtild,Sig_hat, Sig_hat1, Sig_hatd, sig2hat,that, runtime,file=fname)
}
