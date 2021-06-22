set.seed(958348)

size <- 1000 # number of samples: X(0),X(Delta),...,X(size*Delta)
Delta.t <- 1 # Delta

a <- 1
alpha1 <- .9
alpha2 <- .5
mu1 <- 0 # Note 0 for SD
mu2 <- 0
S11 <- 0.18
S22 <- 0.08
rho <- 0.75
eta1 <- -.06
eta2 <- 0
S12 <- rho*sqrt(S11*S22)
lambda <- 0.5

eps <- 1e-10

h2 <- lambda* Delta.t


############################################################# 

# Sample Z
size2 <- 500
X.sample.all <- list()

# Simulate Y
return.process <- sim.vag(a,alpha1,alpha2,mu1,mu2,S11,S22,rho,t.max=size2)
Y1 <- return.process[,1]+eta1
Y2 <- return.process[,2]+eta2

# Simulate base of the CP process of Z
sim.cpbase <- function(a,alpha1,alpha2,S11,S22,rho,t.max){
  S12 <- rho*sqrt(S11*S22)
  if(a >= min(1/c(alpha1,alpha2))){
    return(NA)
  }
  beta1 <- (1-a*alpha1)/alpha1
  beta2 <- (1-a*alpha2)/alpha2
  Sigma <- array(c(S11,S12,S12,S22),c(2,2))
  total <- a+beta1+beta2
  #V0
  S0 <- rexp(n=t.max,rate=1)
  Sigma0 <- outer(c(alpha1,alpha2),c(alpha1,alpha2),pmin)*Sigma
  N0 <- t(array(mvrnorm(n=t.max,mu=c(0,0),Sigma=Sigma0),c(t.max,2)))
  V0 <- rbind(sqrt(S0)*N0[1,],sqrt(S0)*N0[2,])
  #V1
  S1 <- rexp(n=t.max,rate=1)
  Sigma1 <- S1*Sigma0[1,1]
  V1 <- rnorm(n=t.max,mean=0,sd=sqrt(Sigma1))
  #V2
  S2 <- rexp(n=t.max,rate=1)
  Sigma2 <- S2*Sigma0[2,2]
  V2 <- rnorm(n=t.max,mean=0,sd=sqrt(Sigma2))
  
  sample.ind <- sample(0:2,size=t.max,prob=c(a,beta1,beta2),replace=TRUE)
  returns <- t(rep(sample.ind==0,each=2)*V0 + rep(sample.ind==1,each=2)*rbind(V1,0) +
                 rep(sample.ind==2,each=2)*rbind(0,V2))
  #Add them up for compound Poisson
  return(returns)
}

# Simulate Z*(Delta)
sim.cpint <- function(a,alpha1,alpha2,S11,S22,rho,time,t.max){
  returns <- sim.cpbase(a,alpha1,alpha2,S11,S22,rho,t.max)
  arrivals <- sort(runif(t.max,min=0,max=time))
  integ <- cbind(exp(arrivals)*returns[,1],exp(arrivals)*returns[,2])
  return(apply(integ,2,sum))
}
sim.cpint <- Vectorize(sim.cpint,vectorize.args="t.max")

sim.ZDelta <- function(a,alpha1,alpha2,S11,S22,rho,time,t.max){
  rate <- 2*time*(a+(1-a*alpha1)/alpha1+(1-a*alpha2)/alpha2)
  N <- rpois(n=t.max,lambda=rate)
  nonzero.ind <- which(N!=0)
  cp.nonzero <- t(sim.cpint(a,alpha1,alpha2,S11,S22,rho,time,t.max=N[nonzero.ind]))
  returns <- t(array(c(eta1,eta2)*(exp(time)-1),c(2,t.max)))
  returns[nonzero.ind,] <- returns[nonzero.ind,]+cp.nonzero
  return(returns)
} 

# Simulate X
for(j in 1:size2){
  print(j)
  
  Z.Delta.all <- sim.ZDelta(a,alpha1,alpha2,S11,S22,rho,h2,t.max=size)
  Z.Delta1 <- Z.Delta.all[,1]
  Z.Delta2 <- Z.Delta.all[,2]
  X <- array(NA,c(size+1,2)); X[1,] <- c(Y1[j],Y2[j])
  for(i in 1:size){
    X[i+1,] <- exp(-h2)*(X[i,]+c(Z.Delta1[i],Z.Delta2[i]))
  }
  X.sample <- X
  X.sample.all[[j]] <- X.sample
}
