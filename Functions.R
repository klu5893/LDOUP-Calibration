library(pracma)
library(MASS)
library(doParallel)

sim.vag <- function(a,alpha1,alpha2,mu1,mu2,S11,S22,rho,t.max=1000){
  S12 <- rho*sqrt(S11*S22)
  if(a >= min(1/c(alpha1,alpha2))){
    return(NA)
  }
  adjustedmu <- c(mu1,mu2)
  Sigma <- array(c(S11,S12,S12,S22),c(2,2))
  adjustedmu.list <- Sigma.list <- list()
  #V0
  S0 <- rgamma(n=t.max,shape=a,rate=a)
  adjustedmu0 <- a*c(alpha1,alpha2)*adjustedmu
  Sigma0 <- a*outer(c(alpha1,alpha2),c(alpha1,alpha2),pmin)*Sigma
  N0 <- t(array(mvrnorm(n=t.max,mu=c(0,0),Sigma=Sigma0),c(t.max,2)))
  V0 <- rbind(sqrt(S0)*N0[1,]+S0*adjustedmu0[1],sqrt(S0)*N0[2,]+S0*adjustedmu0[2])
  #V1
  S1 <- rgamma(n=t.max,shape=(1-a*alpha1)/alpha1,rate=(1-a*alpha1)/alpha1)
  adjustedmu1 <- S1*(1-a*alpha1)*mu1
  Sigma1 <- S1*(1-a*alpha1)*S11
  V1 <- rnorm(n=t.max,mean=adjustedmu1,sd=sqrt(Sigma1))
  #V2
  S2 <- rgamma(n=t.max,shape=(1-a*alpha2)/alpha2,rate=(1-a*alpha2)/alpha2)
  adjustedmu2 <- S2*(1-a*alpha2)*mu2
  Sigma2 <- S2*(1-a*alpha2)*S22
  V2 <- rnorm(n=t.max,mean=adjustedmu2,sd=sqrt(Sigma2))
  returns <- t(V0)+cbind(V1,V2)
  colnames(returns) <- c("R1","R2")
  return(returns)
}

eval.pdf <- function(R.obs,pdf){
  x.ind.le <- which(pdf$xy[1,]<R.obs[1])
  if(length(x.ind.le)!=0){
    x.ind <- max(x.ind.le)
  }else{
    x.ind <- 1
  }
  y.ind.le <- which(pdf$xy[2,]<R.obs[2])
  if(length(y.ind.le)!=0){
    y.ind <- max(y.ind.le)
  }else{
    y.ind <- 1
  }
  if(x.ind==1 | x.ind==dim(pdf$xy)[2] | y.ind==1 | y.ind==dim(pdf$xy)[2]){
    return(0)
  }
  rv <- pdf$pdf[x.ind,y.ind]
  return(rv)
}

# 1D versions

eval.pdf.1d <- function(R.obs,pdf){
  x.ind.le <- which(pdf$xy<R.obs)
  if(length(x.ind.le)!=0){
    x.ind <- max(x.ind.le)
  }else{
    x.ind <- 1
  }
  if(x.ind==1 | x.ind==length(pdf$xy)){
    return(0)
  }
  rv <- pdf$pdf[x.ind]
  return(rv)
}
eval.pdf.1d <- Vectorize(eval.pdf.1d,vectorize.args="R.obs")

# Lambda objective function

lambda.objf <- function(lambda,m){
  norm1 <- (acf1[2:(m+1)]-exp(-(1:m)*Delta.t*lambda))^2
  norm2 <- (acf2[2:(m+1)]-exp(-(1:m)*Delta.t*lambda))^2
  return(sum(c(norm1,norm2)))
}
