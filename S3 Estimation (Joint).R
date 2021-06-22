Psi.ZD.1D <- function(z1,par){
  a <- par[1]
  alpha1 <- par[2]
  S11 <- par[3]
  
  ce <-
    -(1/alpha1-a)*log((1+alpha1*S11*z1^2/2*exp(2*h2))/(1+alpha1*S11*z1^2/2))
  
  return(ce)
}

Psi.ZD.2D <- function(z1,z2,par){
  a <- par[1]
  alpha1 <- par[2]
  alpha2 <- par[3]
  S11 <- par[4]
  S22 <- par[5]
  S12 <- par[6]
  
  at1 <- alpha1*S11*z1^2/2
  at2 <- alpha2*S22*z2^2/2
  at12 <- min(alpha1,alpha2)*S12*z1*z2
  a0 <- at1+at2+at12
  ex <- exp(2*h2)
  ce <- -a*log((1+a0*ex)/(1+a0))-
          (1/alpha1-a)*log((1+at1*ex)/(1+at1))-(1/alpha2-a)*log((1+at2*ex)/(1+at2))
  
  return(ce)
}


Psi.Y <- function(z1,z2,par){
  a <- par[1]
  alpha1 <- par[2]
  alpha2 <- par[3]
  S11 <- par[4]
  S22 <- par[5]
  S12 <- par[6]
  eta1 <- par[7]
  eta2 <- par[8]
  
  # adjusted drift for the shifted Grid.x
  # Note: h2 depends on lambda which is known
  ce <- 1i*z1*(eta1-eta1*(exp(h2)-1))+1i*z2*(eta2-eta2*(exp(h2)-1))+
        (a-1/alpha1)*log(1+alpha1*S11*z1^2/2)+
        (a-1/alpha2)*log(1+alpha2*S22*z2^2/2)-
        a*log(1+alpha1*S11*z1^2/2+alpha2*S22*z2^2/2+min(alpha1,alpha2)*S12*z1*z2)
  return(ce)
}



llf <- function(par,par.est,X.sample){
  par <- par/parscale
  if(length(par)==2){
    a <- par[1]
    S12 <- par[2]
    
    alpha1 <- par.est[1]
    alpha2 <- par.est[2]
    S11 <- par.est[3]
    S22 <- par.est[4]
    eta1 <- par.est[5]
    eta2 <- par.est[6]
  }else{
    a <- par[1]
    alpha1 <- par[2]
    alpha2 <- par[3]
    S11 <- par[4]
    S22 <- par[5]
    S12 <- par[6]
    
    eta1 <- par.est[1]
    eta2 <- par.est[2]
  }
  

  par <- c(a,alpha1,alpha2,S11,S22,S12,eta1,eta2)

  cat(c("par > ",round(c(a,alpha1,alpha2,S11,S22,S12/sqrt(S11*S22)),2),sep="\n"))

  if(!(eps<a & a<min(1/alpha1,1/alpha2)-eps & alpha1>eps & alpha2>eps & S11>eps & S22>eps &
       abs(S12/sqrt(abs(S11*S22)))<1-eps)){
    return(Inf)
  }
  
  vf <- .5*(exp(2*lambda*Delta.t)-1)
  v1 <- vf*2*S11
  v2 <- vf*2*S22
  v.size <- 0.25/(2^power.N)
  Grid.n <- 2^7*(2^power.N)
  Grid.v <- v.size*sqrt(c(v1,v2))
  Grid.s <- 1/(Grid.n*Grid.v)
  Grid.x <- rbind(Grid.v[1]*((-Grid.n/2):(Grid.n/2-1)),Grid.v[2]*((-Grid.n/2):(Grid.n/2-1)))
  Grid.z <- 2*pi*rbind(Grid.s[1]*((-Grid.n/2):(Grid.n/2-1)),Grid.s[2]*((-Grid.n/2):(Grid.n/2-1)))
  Grid.x <-  Grid.x+array(c(eta1,eta2)*(exp(h2)-1),c(2,Grid.n))
  Z1 <- array(Grid.z[1,],c(Grid.n,Grid.n))
  Z2 <- t(array(Grid.z[2,],c(Grid.n,Grid.n)))
  sgn <- 2*repmat(diag(2),Grid.n/2)-1
  # 1D terms
  Grid.x.1D.1 <- Grid.x[1,]
  Grid.x.1D.2 <- Grid.x[2,]
  Z.1D.1 <- Grid.z[1,]
  Z.1D.2 <- Grid.z[2,]
  sgn.1D <- (-1)^(0:(Grid.n-1))
  
  # Test for out of range
  Y.obs <- t(X.sample[1,])
  Z.obs <- t(X.sample[2:(size+1),])
  all.o1 <- c(Y.obs[1],Z.obs[1,])
  all.o2 <- c(Y.obs[2],Z.obs[2,])
  if(sum(!(min(Grid.x[1,]) <= all.o1 & all.o1 <= max(Grid.x[1,])) | 
         !(min(Grid.x[2,]) <= all.o2 & all.o2 <= max(Grid.x[2,])))!=0){
    return(Inf)
  }
  
  # density Y
  Y.ce <- array(Psi.Y(z1=Z1,z2=Z2,par=par),dim(Z1))
  Y.cf <- exp(Y.ce)
  charfn <- Y.cf*sgn
  P <- Re(fft(charfn))
  P <- P*sgn*Grid.s[1]*Grid.s[2]
  unit.area <- median(diff(Grid.x[1,]))*median(diff(Grid.x[2,]))
  P[P<=0] <- min(P[P>0])/2
  P.Y <- P/(sum(P)*unit.area)
  pdf.Y <- list(xy=Grid.x,pdf=P.Y)


  # density Z1
  Y.cond1.ce <- Psi.ZD.1D(z1=Z.1D.1,par=c(a,alpha1,S11))
  Y.cond1.cf <- exp(Y.cond1.ce)
  pr.nojump1 <- (1-exp(-2*(1/alpha1-a)*h2))*exp(-2*(1/alpha2-a)*h2)*exp(-2*a*h2)
  pr.nojump.1D.1 <- exp(-2*(1/alpha1-a)*h2)
  Y.cf <- (Y.cond1.cf-pr.nojump.1D.1)/(1-pr.nojump.1D.1)
  Y.cf1 <- array(Y.cf,c(Grid.n,Grid.n))
  charfn <- Y.cf*sgn.1D
  P <- Re(fft(charfn))
  P <- P*sgn.1D*Grid.s[1]
  unit.area <- median(diff(Grid.x.1D.1))
  P[P<=0] <- min(P[P>0])/2
  P.Z <- P/(sum(P)*unit.area)
  pdf.Z1 <- list(xy=Grid.x.1D.1,pdf=(pr.nojump1)*P.Z)
  
  # density Z2
  Y.cond2.ce <- Psi.ZD.1D(z1=Z.1D.2,par=c(a,alpha2,S22))
  Y.cond2.cf <- exp(Y.cond2.ce)
  pr.nojump2 <- (1-exp(-2*(1/alpha2-a)*h2))*exp(-2*(1/alpha1-a)*h2)*exp(-2*a*h2)
  pr.nojump.1D.2 <- exp(-2*(1/alpha2-a)*h2)
  Y.cf <- (Y.cond2.cf-pr.nojump.1D.2)/(1-pr.nojump.1D.2)
  Y.cf2 <- t(array(Y.cf,c(Grid.n,Grid.n)))
  charfn <- Y.cf*sgn.1D
  P <- Re(fft(charfn))
  P <- P*sgn.1D*Grid.s[2]
  unit.area <- median(diff(Grid.x.1D.2))
  P[P<=0] <- min(P[P>0])/2
  P.Z <- P/(sum(P)*unit.area)
  pdf.Z2 <- list(xy=Grid.x.1D.2,pdf=(pr.nojump2)*P.Z)
  
  # density Z0
  Y.ce <- array(Psi.ZD.2D(z1=Z1,z2=Z2,par=par),dim(Z1))
  Y.cond.cf <- exp(Y.ce)
  pr.nojump <- exp(-2*(a+(1/alpha1-a)+(1/alpha2-a))*h2)
  Y.cf <- (Y.cond.cf-pr.nojump-pr.nojump1*Y.cf1-pr.nojump2*Y.cf2)/
            (1-pr.nojump-pr.nojump1-pr.nojump2)
  charfn <- Y.cf*sgn
  P <- Re(fft(charfn))
  P <- P*sgn*Grid.s[1]*Grid.s[2]
  unit.area <- median(diff(Grid.x[1,]))*median(diff(Grid.x[2,]))
  P[P<=0] <- min(P[P>0])/2
  P.Z <- P/(sum(P)*unit.area)
  pdf.Z0 <- list(xy=Grid.x,pdf=(1-pr.nojump-pr.nojump1-pr.nojump2)*P.Z)

  # Evaulating density Y and Z
  fY <- -log(eval.pdf(Y.obs,pdf.Y))
  nojump.ind <- abs(Z.obs[1,]-eta1*(exp(h2)-1))<eps & abs(Z.obs[2,]-eta2*(exp(h2)-1))<eps
  jump1.ind <- abs(Z.obs[1,]-eta1*(exp(h2)-1))>=eps & abs(Z.obs[2,]-eta2*(exp(h2)-1))<eps
  jump2.ind <- abs(Z.obs[1,]-eta1*(exp(h2)-1))<eps & abs(Z.obs[2,]-eta2*(exp(h2)-1))>=eps
  jump0.ind <- abs(Z.obs[1,]-eta1*(exp(h2)-1))>=eps & abs(Z.obs[2,]-eta2*(exp(h2)-1))>=eps
  fJ <- -sum(nojump.ind)*log(pr.nojump)
  evaluated.pdf.val <- eval.pdf.1d(Z.obs[1,jump1.ind],pdf=pdf.Z1)
  fZ1 <- -sum(log(evaluated.pdf.val))
  evaluated.pdf.val <- eval.pdf.1d(Z.obs[2,jump2.ind],pdf=pdf.Z2)
  fZ2 <- -sum(log(evaluated.pdf.val))
  Z.obs.list <- split(Z.obs[,jump0.ind],rep(1:sum(jump0.ind), each = 2))
  evaluated.pdf.val <- unlist(lapply(Z.obs.list,eval.pdf,pdf=pdf.Z0))
  fZ0 <- -sum(log(evaluated.pdf.val))

  return(fY+fJ+fZ1+fZ2+fZ0)
}
