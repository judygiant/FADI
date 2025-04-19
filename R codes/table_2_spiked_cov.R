args=(commandArgs(TRUE))
d <- as.numeric(args[1]) ## dimension / number of variables
m <- as.numeric(args[2]) ## number of distributed data splits
mci <- as.numeric(args[3]) ## monte carlo index
p <- 12 ## sketching dimension

#####Define fast power PCA method for aggregation of parallel sketching results
fast_pca_final<-function(y_ls,l,qq){
  d < -dim(y_ls[[1]])[1]
  omega<-matrix(rnorm(d*l),d,l)
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

#### Comparing error rate, runtime among FADI, Fan et al's distributed PCA and the traditional PCA method
err_fast<-function(d,m,n,det,q,L,mc,p,p_0){
  #### d: data dimension
  #### m: number of distributed data splits 
  #### n: sample size on each data split
  #### det: eigengap parameter
  #### q: number of iteration of power method in final PC aggregation
  #### L: number of repeated parallel sketching
  #### mc: number of monte carlos
  #### p: sketching dimension on local machines
  #### p_0: sketching dimension of the power method in final PC aggregation 
  
  N<-m*n ### total sample size
  #### Generating true eigenvalues and top eigenvectors of the covariance
  lm<-4*(det+1)
  sig<-diag(c(lm,lm/2,lm/4,rep(1,d-3)))
  vk<-matrix(0,d,3) ### true top K PCs
  diag(vk)<-1
  err_eigsp<-c()#distributed, fast pca eigenspace error
  err_eigsp_all<-c()
  err_eigsp_fan <- c()
  running_time <- c()
  running_time_all <-c()
  running_time_fan <-c()
  for (rp in 1:mc){
    tt<-0
    #### Generate the distributed Gaussian data
    x<-matrix(rnorm(d*N),N,d)%*%sqrt(sig)
    
    #######this is for full sample traditional PCA########
    ts<-Sys.time()
    sig_hat <- t(x)%*%x/N
    
    vk_all <- svd(sig_hat)$u[,1:3]
    te<-Sys.time()
    running_time_all <- c(running_time_all, difftime(te,ts,units="secs"))
    #################end######################
    
    ##############calculate sig_K for FADI#########
    ts<-Sys.time()
    sig_k <- t(x[,1:4])%*%x[,1:4]/N
    s_hat <- svd(sig_k)$d[4]
    sig_tilde <- list()
    te<-Sys.time()
    tt<-tt+difftime(te,ts,units="secs")
    
    ###########this is for Fan et al's distributed PCA################
    tl<-c()
    t_ch <- c()
    vs <- 0
    for (ch in 1:m){
      ts <- Sys.time()
      xc <- x[((ch-1)*n+1):(ch*n),]
      sig_ch <- t(xc) %*% xc / n
      vc <- svd(sig_ch)$u[,1:3]
      vs <- vs + vc %*% t(vc)
      te <- Sys.time()
      t_ch <- c(t_ch, difftime(te,ts,units="secs"))
    }
    tt_fan <- max(t_ch)
    ts <- Sys.time()
    vs <- vs /m
    vk_fan <- svd(vs)$u[,1:3]
    te <- Sys.time()
    tt_fan <- tt_fan + difftime(te,ts,units="secs")
    running_time_fan <- c(running_time_fan, tt_fan)
    ##############end#################
    
    ########This is FADI#################
    
    ###########distributed fast sketches
    for (j in 1:L){
      ts<-Sys.time()
      set.seed(j)
      omega <-matrix(rnorm(d*p),d,p)
      Y<-0
      te<-Sys.time()
      dt<- difftime(te,ts,units="secs")
      tm <-c()
      ####### Step 1: distributed fast sketching on local machines
      for (i in 1:m){
        ts<-Sys.time()
        xi <- x[((i-1)*n+1):(i*n),]
        Yi<- xi%*% omega
        Yi<- t(xi) %*% Yi
        Y<-Y+Yi
        te<-Sys.time()
        tm<-c(tm,difftime(te,ts,units="secs"))
      }
      dt<-dt+max(tm)
      ts<-Sys.time()
      ###### Step 2: parallel PCA for each gaussian sketching
      Y <- Y/N - s_hat*omega
      vj<- svd(Y)$u[,1:3]
      sig_tilde <- append(sig_tilde,list(vj))
      te<-Sys.time()
      dt <- dt+difftime(te,ts,units="secs")
      tl <- c(tl,dt)
    }
    tt<-tt+max(tl)
    ts<-Sys.time()
    vk_t<- fast_pca_final(sig_tilde, p_0, q)$u[,1:3]
    te<-Sys.time()
    tt<-tt+difftime(te,ts,units="secs")
    running_time<-c(running_time,tt)
    
    er <- norm(vk_t %*% t(vk_t) - vk%*%t(vk),type = "F")
    era <- norm(vk_all %*% t(vk_all) - vk%*%t(vk),type = "F")
    erf <- norm(vk_fan %*% t(vk_fan) - vk%*%t(vk),type = "F")
    err_eigsp <- c(err_eigsp,er)
    err_eigsp_all<-c(err_eigsp_all,era)
    err_eigsp_fan <- c(err_eigsp_fan , erf)
    
    
  }
  ##### error_rate: Error rate of FADI
  ##### running_time: Runtime of FADI
  ##### error_rate_all: Error rate of traditional PCA
  ##### running_time_all: Runtime of traditional PCA
  ##### error_rate_fan: Error rate of Fan et al's distributed PCA
  ##### running_time_fan: Runtime of Fan et al's distributed PCA
  return(list('error_rate'= mean(err_eigsp), 'running_time'=mean(running_time),'error_rate_all'=mean(err_eigsp_all),'running_time_all'=mean(running_time_all), 'error_rate_fan' = mean(err_eigsp_fan), 'running_time_fan' = mean(running_time_fan)))
}


n = 2000
det = 11.5
q = 7
L <- round(d/p*1.2)
mc = 1
#p = 12
p_0 = p


re <- err_fast(d,m,n,det,q,L,mc,p,p_0)

dir.create(here("Spiked_Cov_Results"))
fname<-paste(c(here("Spiked_Cov_Results","results_compare_"),args,".RData"),collapse = '_')
save(re,file=fname)
