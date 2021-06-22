save.image("SD Spec Sample.Rdata")

tt <- proc.time()
ncol <- 8
par.est.all.frame <- array(NA,c(size2,ncol))
#size2 <- 500

registerDoParallel(cores=8)
par.est.all.frame <- foreach(i=1:size2,.combine=rbind) %dopar%
#for(i in 1:size2)
{
  library(pracma)

  X.sample <- X.sample.all[[i]]
  X.sample1 <- X.sample[,1]
  X.sample2 <- X.sample[,2]
  Z.Delta1 <- exp(h2)*X.sample1[2:(size+1)]-X.sample1[1:size]
  Z.Delta2 <- exp(h2)*X.sample2[2:(size+1)]-X.sample2[1:size]
  
  Z.sample1 <- c(X.sample[1,1],Z.Delta1)
  Z.sample2 <- c(X.sample[1,2],Z.Delta2)
  Z.sample <- cbind(Z.sample1,Z.sample2)
  
  
  adjust.test <- 1
  power.N <- 6
  
  # Component 1
  parscale <- 1
  par <- c(alpha1,S11)*adjust.test
  par.optim <- optim(par*parscale,llf.1D,par.est=eta1,X.sample=Z.sample1,method="Nelder-Mead",
                     control=list(maxit=1e3,trace=TRUE,reltol=1e-5),hessian =FALSE)
  par.optim1 <- par.optim
  par.est1 <- par.optim$par/parscale
  
  # Component 2
  parscale <- 1
  par <- c(alpha2,S22)*adjust.test
  par.optim <- optim(par*parscale,llf.1D,par.est=eta2,X.sample=Z.sample2,method="Nelder-Mead",
                     control=list(maxit=1e3,trace=TRUE,reltol=1e-5),hessian=FALSE)
  par.optim2 <- par.optim
  par.est2 <- par.optim$par/parscale
  
  
  
  adjust.test <- 1
  power.N <- 3
  
  # Get a, S12 par
  parscale <- 1
  par <- c(a,S12)*adjust.test
  (alpha1.est <- par.est1[1])
  (alpha2.est <- par.est2[1])
  (S11.est <- par.est1[2])
  (S22.est <- par.est2[2])
  # Ensure initial values are within bounds
  if(!(eps<par[1])){
    par[1] <- 2*eps
  }else if(!(par[1]<min(1/alpha1.est,1/alpha2.est)-eps)){
    par[1] <- min(1/alpha1.est,1/alpha2.est)-2*eps
  }else if(!(par[2]/sqrt(S11.est*S22.est)<1-eps)){
    par[2] <- sqrt(S11.est*S22.est)*(1-2*eps)
  }else if(!(par[2]/sqrt(S11.est*S22.est)>-1+eps)){
    par[2] <- -sqrt(S11.est*S22.est)*(1-2*eps)
  }
  par.est <- c(alpha1.est,alpha2.est,S11.est,S22.est,eta1,eta2)
  
  
  par.optim <- optim(par*parscale,llf,par.est=par.est,X.sample=Z.sample,method="Nelder-Mead",
                     control=list(maxit=1e3,trace=TRUE,reltol=1e-5),hessian =FALSE)
  par.optim0 <- par.optim
  par.est0 <- par.optim$par/parscale
  
  
  
  
  # Optimize over all paramters
  parscale <- 1
  (a.est <- par.est0[1])
  (S12.est <- par.est0[2])
  par <- c(a.est,alpha1.est,alpha2.est,S11.est,S22.est,S12.est)
  par.est <- c(eta1,eta2)
  
  
  par.optim <- optim(par*parscale,llf,par.est=par.est,X.sample=Z.sample,method="Nelder-Mead",
                     control=list(maxit=1e3,trace=TRUE,reltol=1e-5),hessian =FALSE)
  par.est.all <- par.optim$par/parscale
  
  
  # Serial
  # write.val <- par.est.all.frame[i,] <- c(i,par.est.all,par.optim$value)
  # write.val <- array(write.val,c(1,ncol))
  # write(write.val, file="Results SD.csv", append=TRUE,sep=",",ncol=ncol)
  
  # Parallel
  write.val <- returnrow <- c(i,par.est.all,par.optim$value)
  write.val <- array(write.val,c(1,ncol))
  write(write.val, file="Results SD.csv", append=TRUE,sep=",",ncol=ncol)
  returnrow
  
}
proc.time()-tt

par.est.all.frame2 <- par.est.all.frame[,-c(1,ncol)]
#True vale
c1 <- c(a,alpha1,alpha2,S11,S22,S12)
#Mean estimate
c2 <- apply(par.est.all.frame2,2,mean,na.rm=TRUE)
c3 <- abs(c1-apply(par.est.all.frame2,2,mean,na.rm=TRUE))
#RMSE
c4 <- sqrt(apply((t(array(c1,c(dim(par.est.all.frame2)[2],size2)))-par.est.all.frame2)^2,2,mean))
cbind(c1,c2,c3,c4)
