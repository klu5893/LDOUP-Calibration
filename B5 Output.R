generate.mom.BDLP <- function(par){
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
  
  E1 <- (exp(1*lambda*Delta.t)-1)/1
  E2 <- (exp(2*lambda*Delta.t)-1)/2
  E3 <- (exp(3*lambda*Delta.t)-1)/3
  E4 <- (exp(4*lambda*Delta.t)-1)/4
  
  P11 <- mu1+eta1
  P21 <- S11+alpha1*mu1^2
  P31 <- 3*alpha1*S11*mu1+2*alpha1^2*mu1^3
  P41 <- 3*S11^2+3*alpha1*(S11^2+2*S11*mu1^2)+3*alpha1^2*(4*S11*mu1^2+mu1^4)+6*alpha1^3*mu1^4
  P12 <- mu2+eta2
  P22 <- S22+alpha2*mu2^2
  P32 <- 3*alpha2*S22*mu2+2*alpha2^2*mu2^3
  P42 <- 3*S22^2+3*alpha2*(S22^2+2*S22*mu2^2)+3*alpha2^2*(4*S22*mu2^2+mu2^4)+6*alpha2^3*mu2^4
  C12 <- a*(min(alpha1,alpha2)*S12+alpha1*alpha2*mu1*mu2)
  
  rv <- c(E1*P11,
          E2*P21,
          E3*P31,
          E4*P41-3*E2*P21^2,
          E1*P12,
          E2*P22,
          E3*P32,
          E4*P42-3*E2*P22^2,
          E2*C12)
  return(rv)
}

# space <- 100
# size <- 1000
# Delta <- 1/space
# M <- size/Delta
# Delta.t <- Delta*space
# 
# a <- 1
# alpha1 <- .9
# alpha2 <- .5
# mu1 <- 0.15
# mu2 <- -0.06
# S11 <- 0.18
# S22 <- 0.08
# rho <- 0.75
# eta1 <- -.06
# eta2 <- 0
# S12 <- rho*sqrt(S11*S22)
# lambda <- 0.5

load("BDLP Spec Sample.Rdata")

# BDLP Results
par.est.all.frame <- read.csv("Results BDLP.csv", header=FALSE)
nrow <- dim(par.est.all.frame)[1]
ncol <- dim(par.est.all.frame)[2]

par.est.all.frame2 <- par.est.all.frame[,-c(1,ncol)]
#True vale
c1 <- c(a,alpha1,alpha2,mu1,mu2,S11,S22,S12,eta1,eta2,lambda)
c1 <- c(c1,generate.mom.BDLP(c1))
#Mean estimate
moments <- array(NA,c(nrow,9))
for(i in 1:nrow){
  moments[i,] <- generate.mom.BDLP(as.numeric(par.est.all.frame2[i,]))
}
par.est.all.frame3 <- cbind(par.est.all.frame2,moments)
c2 <- apply(par.est.all.frame3,2,mean,na.rm=TRUE)
c3 <- abs(c1-c2)
#RMSE
c4 <- sqrt(apply((t(array(c1,c(dim(par.est.all.frame3)[2],nrow)))-par.est.all.frame3)^2,2,mean))
cbind(c1,c2,c3,c4)

# Quick summary of avg % error
100*mean((abs((c1-c2)/c1))[is.finite(abs((c1-c2)/c1))])
# Output for paper
round(cbind(c1,c2,c4),4)

#### Simulation
# Here, we set Delta (in the paper) = Delta.t (in the code) = 1/100 to be
# consistent with Fig 1 (note: \tilde{Delta} (in the paper) = Delta (in the
# code) = 1/10000).

set.seed(245780)

TT <- 1000
Delta.t <- 1/100 # Delta
space <- 100 
Delta <- Delta.t/space # Delta~
M <- TT/Delta

h <- lambda*Delta
size2 <- 1
grid.pts <- expand.grid(Grid.x[1,],Grid.x[2,])
sample.ind <- sample(1:dim(grid.pts)[1],size=size2,prob=c(P),replace=TRUE)
sample.pts <- grid.pts[sample.ind,]

X.sample.all <- list()

for(j in 1:size2){
  print(j)
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

# Plot
x1.text = expression(italic("X")[1]*"("*italic("t")*")")
x2.text = expression(italic("X")[2]*"("*italic("t")*")")
t.text = expression(italic("t"))
x.text = expression(bold("X")*"("*italic("t")*")")

pdf("plot2.pdf",width = 5, height = 6) # width=5 is approx text width in LaTeX
layout.matrix <- matrix(c(1,2,3), nrow = 3, ncol = 1)
layout(mat = layout.matrix, heights = c(1,1,2))

# Plot 1
par(mar=c(0,3,3,3)+c(0,.1,0,0),mgp=c(2,.5,0))
tt <- seq(0,TT,by=Delta.t)
n <- TT/Delta.t+1
plot(tt[1:n],X.sample[1:n,1],type="l",xlab="",ylab=x1.text,xlim=c(tt[1],tt[n]),
     xaxt="n",xaxs="i",yaxt="n")
axis(2, at=c(-1,0,1,2))
# Plot 2
par(mar=c(3,3,0,3)+c(0,.1,0,0),mgp=c(2,.5,0))
plot(tt[1:n],X.sample[1:n,2],type="l",xlab=t.text,ylab=x2.text,xlim=c(tt[1],tt[n]),
     xaxs="i",yaxt="n")
axis(2, at=c(-1,0,1))
# Plot 3
par(mar=c(3,3,3,3)+c(0,.1,0,0),mgp=c(2,.5,0))
n <- 10/Delta.t+1
plot(tt[1:n],X.sample[1:n,1],type="l",xlab=t.text,ylab=x.text,xlim=c(tt[1],tt[n]),
     ylim=summary(c(X.sample[1:n,1],X.sample[1:n,2]))[c(1,6)],xaxs="i")
lines(tt[1:n],X.sample[1:n,2],col="blue",lty=5)
legend("topleft", legend=c(x1.text,x2.text),
       col=c("black", "blue"), lty=c(1,5))

dev.off()
