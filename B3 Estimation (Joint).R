########### 2 = Psi Z*(Delta)
Psi.Z2.int <- function(t,z1,z2,par,type){
  a <- par[1]
  alpha1 <- par[2]
  alpha2 <- par[3]
  mu1 <- par[4]
  mu2 <- par[5]
  S11 <- par[6]
  S22 <- par[7]
  S12 <- par[8]
  eta1 <- par[9]
  eta2 <- par[10]
  lambda <- par[11]
  
  Z1 <- exp(t)*z1
  Z2 <- exp(t)*z2
  wvag.cf <- (1+alpha1*(-1i*mu1*Z1+S11*Z1^2/2))^((a-1/alpha1))*
    (1+alpha2*(-1i*mu2*Z2+S22*Z2^2/2))^((a-1/alpha2))*
    (1+alpha1*(-1i*mu1*Z1+S11*Z1^2/2)+alpha2*(-1i*mu2*Z2+S22*Z2^2/2)+min(alpha1,alpha2)*S12*Z1*Z2)^(-a)
  wvag.ce <- log(wvag.cf)+1i*Z1*eta1+1i*Z2*eta2
  if(type=="Re"){
    return(Re(wvag.ce))
  }else if(type=="Im"){
    return(Im(wvag.ce))
  }else{
    return(wvag.ce)
  }
}

Psi.Y2 <- function(z1,z2,par){
  lambda <- par[11]
  inte.Re <- integrate(Psi.Z2.int,lower=0,upper=Delta.t*lambda,z1=z1,z2=z2,par=par,type="Re")
  inte.Im <- integrate(Psi.Z2.int,lower=0,upper=Delta.t*lambda,z1=z1,z2=z2,par=par,type="Im")
  return((inte.Re$val+1i*inte.Im$val))
}
Psi.Y2.V <- Vectorize(Psi.Y2,vectorize.args=c("z1","z2"))





llf <- function(par,par.est,X.sample){
  par <- par/parscale
  
  S12 <- par
  
  alpha1 <- par.est[1]
  alpha2 <- par.est[2]
  mu1 <- par.est[3]
  mu2 <- par.est[4]
  S11 <- par.est[5]
  S22 <- par.est[6]
  eta1 <- par.est[7]
  eta2 <- par.est[8]
  lambda <- par.est[9]
  
  Y.obs <- t(X.sample[1,])
  Z.obs <- t(X.sample[2:(size+1),])
  vf <- .5*(exp(2*lambda*Delta.t)-1)
  a <- cov(Z.obs[1,],Z.obs[2,])/(vf*(min(alpha1,alpha2)*S12+alpha1*alpha2*mu1*mu2))
  par <- c(a,alpha1,alpha2,mu1,mu2,S11,S22,S12,eta1,eta2,lambda)
  
  cat(c("par > ",round(c(a,S12/sqrt(S11*S22)),2),sep="\n"))

  if(!(eps <a & a<min(1/alpha1,1/alpha2)-eps &
       abs(S12/sqrt(abs(S11*S22)))<1-eps & lambda > eps)){
    return(Inf)
  }
  
  v1 <- vf*(alpha1*mu1^2+S11)
  v2 <- vf*(alpha2*mu2^2+S22)
  v.size <- 0.25/(2^power.N)
  Grid.n <- 2^7*(2^power.N)
  Grid.v <- v.size*sqrt(c(v1,v2))
  Grid.s <- 1/(Grid.n*Grid.v)
  Grid.x <- rbind(Grid.v[1]*((-Grid.n/2):(Grid.n/2-1)),Grid.v[2]*((-Grid.n/2):(Grid.n/2-1)))
  Grid.z <- 2*pi*rbind(Grid.s[1]*((-Grid.n/2):(Grid.n/2-1)),Grid.s[2]*((-Grid.n/2):(Grid.n/2-1)))
  Grid.x <-  Grid.x
  Z1 <- array(Grid.z[1,],c(Grid.n,Grid.n))
  Z2 <- t(array(Grid.z[2,],c(Grid.n,Grid.n)))
  sgn <- 2*repmat(diag(2),Grid.n/2)-1
  
  # Test for out of range
  all.o1 <- c(Y.obs[1],Z.obs[1,])
  all.o2 <- c(Y.obs[2],Z.obs[2,])
  if(sum(!(min(Grid.x[1,]) <= all.o1 & all.o1 <= max(Grid.x[1,])) | 
         !(min(Grid.x[2,]) <= all.o2 & all.o2 <= max(Grid.x[2,])))!=0){
    return(Inf)
  }
  

  # density Y
  Y.ce <- try(array(Psi.Y.V(z1=Z1,z2=Z2,par=par),dim(Z1)))
  if(is.character(Y.ce)==TRUE){
    return(Inf)
  }
  Y.cf <- exp(Y.ce)
  charfn <- Y.cf*sgn
  P <- Re(fft(charfn))
  P <- P*sgn*Grid.s[1]*Grid.s[2]
  unit.area <- median(diff(Grid.x[1,]))*median(diff(Grid.x[2,]))
  P[P<=0] <- min(P[P>0])/2
  P.Y <- P/(sum(P)*unit.area)
  pdf.Y <- list(xy=Grid.x,pdf=P.Y)
  
  # density Z
  Y.ce <- try(array(Psi.Y2.V(z1=Z1,z2=Z2,par=par),dim(Z1)))
  if(is.character(Y.ce)==TRUE){
    return(Inf)
  }
  Y.cf <- exp(Y.ce)
  charfn <- Y.cf*sgn
  P <- Re(fft(charfn))
  P <- P*sgn*Grid.s[1]*Grid.s[2]
  unit.area <- median(diff(Grid.x[1,]))*median(diff(Grid.x[2,]))
  P[P<=0] <- min(P[P>0])/2
  P.Z <- P/(sum(P)*unit.area)
  pdf.Z <- list(xy=Grid.x,pdf=P.Z)
  
  # Evaulating density Y and Z
  fY <- -log(eval.pdf(Y.obs,pdf.Y))
  R.obs <- Z.obs
  R.obs.list <- split(R.obs,rep(1:ncol(R.obs), each = nrow(R.obs)))
  evaluated.pdf.val <- unlist(lapply(R.obs.list,eval.pdf,pdf=pdf.Z))
  fZ <- -sum(log(evaluated.pdf.val))
  const <- -2*size*lambda*Delta.t
  
  
  return(fY+const+fZ)
}
