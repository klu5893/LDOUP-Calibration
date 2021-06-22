Psi.ZD.1 <- function(z1,par){
  alpha1 <- par[1]
  S11 <- par[2]

  ce <- 
    -(1/alpha1)*log((1+alpha1*S11*z1^2/2*exp(2*h2))/(1+alpha1*S11*z1^2/2))
  
  return(ce)
}

Psi.Y.1 <- function(z1,par){
  alpha1 <- par[1]
  S11 <- par[2]
  eta1 <- par[3]
  
  # adjusted drift for the shifted Grid.x
  # Note: h2 depends on lambda which is known
  ce <- 1i*z1*(eta1-eta1*(exp(h2)-1))-(1/alpha1)*log(1+alpha1*S11*z1^2/2)
  return(ce)
}


llf.1D <- function(par,par.est,X.sample){
  par <- par/parscale
  alpha1 <- par[1]
  S11 <- par[2]
  
  eta1 <- par.est[1]

  par <- c(alpha1,S11,eta1)
  
  cat(c("par >",round(par,2),sep="\n"))
  
  if(!(alpha1>eps & S11>eps)){
    return(Inf)
  }
  
  vf <- .5*(exp(2*lambda*Delta.t)-1)
  v1 <- vf*2*S11
  v.size <- 0.5/(2^power.N)
  Grid.n <- 2^7*(2^power.N)
  Grid.v <- v.size*sqrt(v1)
  Grid.s <- 1/(Grid.n*Grid.v)
  Grid.x <- Grid.v[1]*((-Grid.n/2):(Grid.n/2-1))
  Grid.z <- 2*pi*Grid.s[1]*((-Grid.n/2):(Grid.n/2-1))
  Grid.x <-  Grid.x+eta1*(exp(h2)-1)
  Z1 <- Grid.z
  sgn <- (-1)^(0:(Grid.n-1))
  
  # Test for out of range
  Y.obs <- X.sample[1]
  Z.obs <- X.sample[2:(size+1)]
  all.o1 <- c(Y.obs,Z.obs)
  if(sum(!(min(Grid.x) <= all.o1 & all.o1 <= max(Grid.x)))!=0){
    return(Inf)
  }
  
  # density Y
  Y.ce <- Psi.Y.1(z1=Z1,par=par)
  Y.cf <- exp(Y.ce)
  charfn <- Y.cf*sgn
  P <- Re(fft(charfn))
  P <- P*sgn*Grid.s
  unit.area <- median(diff(Grid.x))
  P[P<=0] <- min(P[P>0])/2
  P.Y <- P/(sum(P)*unit.area)
  pdf.Y <- list(xy=Grid.x,pdf=P.Y)
  
  # density Z
  Y.cond.ce <- Psi.ZD.1(z1=Z1,par=par)
  Y.cond.cf <- exp(Y.cond.ce)
  pr.nojump <- exp(-2/alpha1*h2)
  Y.cf <- (Y.cond.cf-pr.nojump)/(1-pr.nojump)
  charfn <- Y.cf*sgn
  P <- Re(fft(charfn))
  P <- P*sgn*Grid.s
  unit.area <- median(diff(Grid.x))
  P[P<=0] <- min(P[P>0])/2
  P.Z <- P/(sum(P)*unit.area)
  pdf.Z <- list(xy=Grid.x,pdf=(1-pr.nojump)*P.Z)
  
  # Evaulating density Y and Z
  fY <- -log(eval.pdf.1d(Y.obs,pdf.Y))
  nojump.ind <- abs(Z.obs-eta1*(exp(h2)-1))<eps
  fJ <- -sum(nojump.ind)*log(pr.nojump)
  evaluated.pdf.val <- eval.pdf.1d(Z.obs[nojump.ind==FALSE],pdf=pdf.Z)
  fZ <- -sum(log(evaluated.pdf.val))
  
  return(fY+fJ+fZ)
}
