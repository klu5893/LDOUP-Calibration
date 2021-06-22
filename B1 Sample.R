set.seed(212084)

space <- 100 # 1/Delta~
size <- 1000 # number of samples: X(0),X(Delta),...,X(size*Delta)
Delta <- 1/space # Delta~
M <- size/Delta # number of points that need to be sampled 
Delta.t <- Delta*space # Delta = 1

a <- 1
alpha1 <- .9
alpha2 <- .5
mu1 <- 0.15
mu2 <- -0.06
S11 <- 0.18
S22 <- 0.08
rho <- 0.75
eta1 <- -.06
eta2 <- 0
S12 <- rho*sqrt(S11*S22)
lambda <- 0.5

eps <- 1e-10




###########################################################################
###########################################################################

############# 1 = Psi Y

Psi.Z.int <- function(t,z1,z2,par,type){
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
  
  Z1 <- exp(-t)*z1
  Z2 <- exp(-t)*z2
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

Psi.Y <- function(z1,z2,par){
  inte.Re <- integrate(Psi.Z.int,lower=0,upper=Inf,z1=z1,z2=z2,par=par,type="Re")
  inte.Im <- integrate(Psi.Z.int,lower=0,upper=Inf,z1=z1,z2=z2,par=par,type="Im")
  return((inte.Re$val+1i*inte.Im$val))
}
Psi.Y.V <- Vectorize(Psi.Y,vectorize.args=c("z1","z2"))


###########################################################################
###########################################################################

par <- c(a,alpha1,alpha2,mu1,mu2,S11,S22,S12,eta1,eta2,lambda)
power.N <- 3

# Fourier inversion for Y
vf <- .5*(exp(2*lambda*Delta.t)-1)
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

Y.ce <- array(Psi.Y.V(z1=Z1,z2=Z2,par=par),dim(Z1))
Y.cf <- exp(Y.ce)
charfn <- Y.cf*sgn
P <- Re(fft(charfn))
P <- P*sgn*Grid.s[1]*Grid.s[2]
unit.area <- median(diff(Grid.x[1,]))*median(diff(Grid.x[2,]))
P[P<=0] <- min(P[P>0])/2
P <- P/(sum(P)*unit.area)

###########################################################################
###########################################################################

# Sample Y
h <- lambda*Delta
size2 <- 500
grid.pts <- expand.grid(Grid.x[1,],Grid.x[2,])
sample.ind <- sample(1:dim(grid.pts)[1],size=size2,prob=c(P),replace=TRUE)
sample.pts <- grid.pts[sample.ind,]

X.sample.all <- list()

for(j in 1:size2){
  # Simualate Z
  return.process <- sim.vag(a*h,alpha1/h,alpha2/h,mu1*h,mu2*h,S11*h,S22*h,rho,t.max=M)
  R1 <- return.process[,1]+eta1*h
  R2 <- return.process[,2]+eta2*h
  
  # Simualate LDOUP
  X <- array(NA,c(M+1,2)); X[1,] <- as.numeric(sample.pts[j,])
  for(i in 1:M){
    X[i+1,] <- exp(-h)*(X[i,]+exp(h)*c(R1[i],R2[i]))
  }
  X.sample <- X[seq(from=1,to=M+1,by=space),]
  X.sample.all[[j]] <- X.sample
}
