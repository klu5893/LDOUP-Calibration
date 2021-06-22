save.image("BDLP Spec Sample.Rdata")

tt <- proc.time()
ncol <- 13
par.est.all.frame <- array(NA,c(size2,ncol))
#size2 <- 500

registerDoParallel(cores=8)
par.est.all.frame <- foreach(i=1:size2,.combine=rbind) %dopar%
  #for(i in 1:size2)
  {
    library(pracma)
    Psi.Y.V.1D <- Vectorize(Psi.Y.1D,vectorize.args="z1")
    Psi.Y2.V.1D <- Vectorize(Psi.Y2.1D,vectorize.args="z1")
    Psi.Y.V <- Vectorize(Psi.Y,vectorize.args=c("z1","z2"))
    Psi.Y2.V <- Vectorize(Psi.Y2,vectorize.args=c("z1","z2"))
  
  X.sample <- X.sample.all[[i]]
  X.sample1 <- X.sample[,1]
  X.sample2 <- X.sample[,2]
  
  acf1 <- acf(X.sample[,1],plot=FALSE)$acf
  acf2 <- acf(X.sample[,2],plot=FALSE)$acf
  lambda.est <-  optim(NA,lambda.objf,m=1,method="Brent",lower=0.01,upper=10,
                       control=list(maxit=1e3,reltol=1e-8))$par
  Z.Delta1 <- exp(lambda.est*Delta.t)*X.sample1[2:(size+1)]-X.sample1[1:size]
  Z.Delta2 <- exp(lambda.est*Delta.t)*X.sample2[2:(size+1)]-X.sample2[1:size]
  Z.sample1 <- c(X.sample[1,1],Z.Delta1)
  Z.sample2 <- c(X.sample[1,2],Z.Delta2)
  Z.sample <- cbind(Z.sample1,Z.sample2)

  
  adjust.test <- 1
  power.N <- 6

  # Component 1
  parscale <- 1
  par <- c(alpha1,mu1,S11,eta1)*adjust.test
  par.optim <- optim(par*parscale,llf.1D,par.est=lambda.est,X.sample=Z.sample1,method="Nelder-Mead",
                     control=list(maxit=1e3,trace=FALSE,reltol=1e-5),hessian =FALSE)
  par.optim1 <- par.optim
  par.est1 <- par.optim$par/parscale
  
  
  
  # Component 2
  parscale <- 1
  par <- c(alpha2,mu2,S22,eta2)*adjust.test
  par.optim <- optim(par*parscale,llf.1D,par.est=lambda.est,X.sample=Z.sample2,method="Nelder-Mead",
                     control=list(maxit=1e3,trace=FALSE,reltol=1e-5),hessian=FALSE)
  par.optim2 <- par.optim
  par.est2 <- par.optim$par/parscale

  

  
  # Get a, S12 par
  power.N <- 3
  parscale <- 1
  
  (alpha1.est <- par.est1[1])
  (alpha2.est <- par.est2[1])
  (mu1.est <- par.est1[2])
  (mu2.est <- par.est2[2])
  (S11.est <- par.est1[3])
  (S22.est <- par.est2[3])
  (eta1.est <- par.est1[4])
  (eta2.est <- par.est2[4])
  par.est <- c(alpha1.est,alpha2.est,mu1.est,mu2.est,S11.est,S22.est,eta1.est,eta2.est,lambda.est)
  
  vf <- .5*(exp(2*lambda.est*Delta.t)-1)
  cov.Z <- cov(Z.Delta1,Z.Delta2)
  S12.bound <- ((cov.Z/vf)/min(1/alpha1.est,1/alpha2.est)-alpha1.est*alpha2.est*mu1.est*mu2.est)/min(alpha1.est,alpha2.est)
  if(cov.Z>eps){
    S12.lb <- max(S12.bound,-sqrt(S11.est*S22.est))+eps
    S12.ub <- sqrt(S11.est*S22.est)-eps
  }else if(cov.Z<eps){
    S12.lb <- -sqrt(S11.est*S22.est)+eps
    S12.ub <- min(S12.bound,sqrt(S11.est*S22.est))-eps
  }
  
  if(abs(cov.Z)<eps){
    # cov.Z = 0 case (technically, search over "a" along the line S12=S12.bound)
    S12.est <- S12.bound
    a.est <- 0
    par.optim.fn <- NA
  }else if(S12.lb>S12.ub){
    # No exact solution: a takes its upper bound, take closest S12
    S12.est <- (cov.Z>0)*sqrt(S11.est*S22.est)
    a.est <- min(1/alpha1.est,1/alpha2.est)
    par.optim.fn <- NA
  }else{
    # The usual case: optimize along the cov constraint
    par.optim <- optim(NA,llf,par.est=par.est,X.sample=Z.sample,method="Brent",
                       lower=S12.lb,upper=S12.ub,
                       control=list(maxit=1e3,trace=TRUE,reltol=1e-4))
    S12.est <- par.optim$par
    a.est <- cov.Z/(vf*(min(alpha1.est,alpha2.est)*S12.est+alpha1.est*alpha2.est*mu1.est*mu2.est))
    par.optim.fn <- par.optim$value
  }

  par.est.all <- c(a.est,alpha1.est,alpha2.est,mu1.est,mu2.est,S11.est,S22.est,S12.est,
                   eta1.est,eta2.est,lambda.est)
  
  
  # Serial
  # write.val <- par.est.all.frame[i,] <- c(i,par.est.all,par.optim.fn)
  # write.val <- array(write.val,c(1,ncol))
  # write(write.val, file="Results BDLP.csv", append=TRUE,sep=",",ncol=ncol)
  
  # Parallel
  write.val <- returnrow <- c(i,par.est.all,par.optim.fn)
  write.val <- array(write.val,c(1,ncol))
  write(write.val, file="Results BDLP.csv", append=TRUE,sep=",",ncol=ncol)
  returnrow
  
}
proc.time()-tt

par.est.all.frame2 <- par.est.all.frame[,-c(1,ncol)]
#True vale
c1 <- c(a,alpha1,alpha2,mu1,mu2,S11,S22,S12,eta1,eta2,lambda)
#Mean estimate
c2 <- apply(par.est.all.frame2,2,mean,na.rm=TRUE)
c3 <- abs(c1-apply(par.est.all.frame2,2,mean,na.rm=TRUE))
#RMSE
c4 <- sqrt(apply((t(array(c1,c(dim(par.est.all.frame2)[2],size2)))-par.est.all.frame2)^2,2,mean))
cbind(c1,c2,c3,c4)
