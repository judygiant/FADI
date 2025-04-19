######inputting arguments from the shell script
args=(commandArgs(TRUE))
mc <- as.numeric(args[1])
K = 4
p = 50
m = 10
library(here)
#####Loading graph generated from the 1000 Genome Data
load(here("Data", "1000g_sbm95.RData"))
####Loading the population labels
info <- read.delim(here("Data", "1KG_TRACE_pca.txt"), sep = " ")
ll <- order(info$Population.2)
pop_vec <- info$Population.2[ll]
loc_pure <- pop_vec %in% c("AFR", "EAS", "EUR", "SAS")
Mh <- sbm95[loc_pure, loc_pure] #Obtaining hat{M}, including only four super populations: AFR, EAS, EUR and SAS

d = dim(Mh)[1] #dimension of hat{M}
dj <- ceiling(d/m) #the dimension of data split on each machine

#####Estimating parameters from observed data matrix
thetah = sum(Mh)/d^2/2
mu = d*log(d)*sqrt(thetah/p)/12

#####Computing the local sketchings on each machine
t_ls <- c()
Yt <- 0
for (j in 1:m){
  loc <- (dj*(j-1)+1):(min(dj*j,d)) ###data index corresponding to each machine
  djj <- length(loc)
  ts <- Sys.time()
  omega <- matrix(rnorm(djj*p),djj,p)
  Ytj <- Mh[,loc]%*%omega
  Yt <- Yt + Ytj
  te <- Sys.time()
  t_ls <- c(t_ls, difftime(te,ts,units = "secs"))
}

ts <- Sys.time()
svdY <- svd(Yt/sqrt(p)) ##### SVD for each parallel sketching
te <- Sys.time()
rt <- difftime(te,ts,units = "secs") + max(t_ls) ####computing the runtime for step 1 and step 2 of FADI
vkl <- svdY$u[,1:K] ##### Extract the top K left singular vectors for each parallel sketching

#####estimate K
dd <- svdY$d
diff <- dd-dd[p]
k <- 1
flag <- FALSE
while (!flag) {
  flag <- all(diff[k:p] <= mu)
  k<-k+1
}
kl<-k-2

#####Saving data
dir.create(here("1000g_inf_Results"))
fname<-paste(c(here("1000g_inf_Results","results_"),args,".RData"),collapse = '_')

save(vkl,kl,rt, file = fname)