Psi.Z.int.1D <- function(t,z1,par,type){
  alpha1 <- par[1]
  mu1 <- par[2]
  S11 <- par[3]
  eta1 <- par[4]
  lambda <- par[5]
  
  Z1 <- exp(-t)*z1
  wvag.cf <- (1+alpha1*(-1i*mu1*Z1+S11*Z1^2/2))^(-1/alpha1)
  wvag.ce <- log(wvag.cf)+1i*Z1*eta1
  if(type=="Re"){
    return(Re(wvag.ce))
  }else if(type=="Im"){
    return(Im(wvag.ce))
  }else{
    return(wvag.ce)
  }
}

Psi.Y.1D <- function(z1,par){
  inte.Re <- integrate(Psi.Z.int.1D,lower=0,upper=Inf,z1=z1,par=par,type="Re")
  inte.Im <- integrate(Psi.Z.int.1D,lower=0,upper=Inf,z1=z1,par=par,type="Im")
  return((inte.Re$val+1i*inte.Im$val))
}
Psi.Y.V.1D <- Vectorize(Psi.Y.1D,vectorize.args="z1")

########### 2 = Psi Z*(Delta)

Psi.Z2.int.1D <- function(t,z1,par,type){
  alpha1 <- par[1]
  mu1 <- par[2]
  S11 <- par[3]
  eta1 <- par[4]
  lambda <- par[5]
  
  Z1 <- exp(t)*z1
  wvag.cf <- (1+alpha1*(-1i*mu1*Z1+S11*Z1^2/2))^(-1/alpha1)
  wvag.ce <- log(wvag.cf)+1i*Z1*eta1
  if(type=="Re"){
    return(Re(wvag.ce))
  }else if(type=="Im"){
    return(Im(wvag.ce))
  }else{
    return(wvag.ce)
  }
}

Psi.Y2.1D <- function(z1,par){
  lambda <- par[5]
  inte.Re <- integrate(Psi.Z2.int.1D,lower=0,upper=Delta.t*lambda,z1=z1,par=par,type="Re")
  inte.Im <- integrate(Psi.Z2.int.1D,lower=0,upper=Delta.t*lambda,z1=z1,par=par,type="Im")
  return((inte.Re$val+1i*inte.Im$val))
}
Psi.Y2.V.1D <- Vectorize(Psi.Y2.1D,vectorize.args="z1")

llf.1D <- function(par,par.est,X.sample){
  par <- par/parscale
  alpha1 <- par[1]
  mu1 <- par[2]
  S11 <- par[3]
  eta1 <- par[4]
  
  lambda <- par.est[1]
  par <- c(alpha1,mu1,S11,eta1,lambda)
  
  cat(c("par >",round(par,2),sep="\n"))
  
  if(!(alpha1>eps & S11>eps & lambda>eps)){
    return(Inf)
  }
  
  vf <- .5*(exp(2*lambda*Delta.t)-1)
  v1 <- vf*(alpha1*mu1^2+S11)
  v.size <- 0.5/(2^power.N)
  Grid.n <- 2^7*(2^power.N)
  Grid.v <- v.size*sqrt(v1)
  Grid.s <- 1/(Grid.n*Grid.v)
  Grid.x <- Grid.v[1]*((-Grid.n/2):(Grid.n/2-1))
  Grid.z <- 2*pi*Grid.s[1]*((-Grid.n/2):(Grid.n/2-1))
  Grid.x <-  Grid.x
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
  Y.ce <- try(Psi.Y.V.1D(z1=Z1,par=par))
  if(is.character(Y.ce)==TRUE){
    return(Inf)
  }
  Y.cf <- exp(Y.ce)
  charfn <- Y.cf*sgn
  P <- Re(fft(charfn))
  P <- P*sgn*Grid.s
  unit.area <- median(diff(Grid.x))
  P[P<=0] <- min(P[P>0])/2
  P.Y <- P/(sum(P)*unit.area)
  pdf.Y <- list(xy=Grid.x,pdf=P.Y)
  
  # density Z
  Y.ce <- try(Psi.Y2.V.1D(z1=Z1,par=par))
  if(is.character(Y.ce)==TRUE){
    return(Inf)
  }
  Y.cf <- exp(Y.ce)
  charfn <- Y.cf*sgn
  P <- Re(fft(charfn))
  P <- P*sgn*Grid.s
  unit.area <- median(diff(Grid.x))
  P[P<=0] <- min(P[P>0])/2
  P.Z <- P/(sum(P)*unit.area)
  pdf.Z <- list(xy=Grid.x,pdf=P.Z)
  
  # Evaulating density Y and Z
  fY <- -log(eval.pdf.1d(Y.obs,pdf.Y))
  R.obs <- Z.obs
  evaluated.pdf.val <- eval.pdf.1d(R.obs,pdf=pdf.Z)
  fZ <- -sum(log(evaluated.pdf.val))
  const <- -1*size*lambda*Delta.t
  
  return(fY+const+fZ)
}
