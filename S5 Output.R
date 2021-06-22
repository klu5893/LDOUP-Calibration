generate.mom.SD <- function(par){
  a <- par[1]
  alpha1 <- par[2]
  alpha2 <- par[3]
  S11 <- par[4]
  S22 <- par[5]
  S12 <- par[6]
  
  E1 <- (exp(1*lambda*Delta.t)-1)/1
  E2 <- (exp(2*lambda*Delta.t)-1)/2
  E3 <- (exp(3*lambda*Delta.t)-1)/3
  E4 <- (exp(4*lambda*Delta.t)-1)/4
  
  P11 <- eta1
  P21 <- 2*S11
  P31 <- 0
  P41 <- 12*S11^2*(1+alpha1)
  P12 <- eta2
  P22 <- 2*S22
  P32 <- 0
  P42 <- 12*S22^2*(1+alpha2)
  C12 <- 2*a*min(alpha1,alpha2)*S12
  
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
# mu1 <- 0 # Note 0 for SD
# mu2 <- 0
# S11 <- 0.18
# S22 <- 0.08
# rho <- 0.75
# eta1 <- -.06
# eta2 <- 0
# S12 <- rho*sqrt(S11*S22)
# lambda <- 0.5

load("SD Spec Sample.Rdata")

# BDLP Results
par.est.all.frame <- read.csv("Results SD.csv", header=FALSE)
nrow <- dim(par.est.all.frame)[1]
ncol <- dim(par.est.all.frame)[2]

par.est.all.frame2 <- par.est.all.frame[,-c(1,ncol)]
#True vale
c1 <- c(a,alpha1,alpha2,S11,S22,S12)
c1 <- c(c1,generate.mom.SD(c1))
#Mean estimate
moments <- array(NA,c(nrow,9))
for(i in 1:nrow){
  moments[i,] <- generate.mom.SD(as.numeric(par.est.all.frame2[i,]))
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

# Set sampling interval
TT <- 1000
space <- 100
size <- TT*space # number of points to simulate
Delta.t <- 1/space # Delta
h2 <- lambda* Delta.t
size2 <- 1
set.seed(103824)

# Simulate
return.process <- sim.vag(a,alpha1,alpha2,mu1,mu2,S11,S22,rho,t.max=size2)
Y1 <- return.process[,1]+eta1
Y2 <- return.process[,2]+eta2

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

# Plot
x1.text = expression(italic("X")[1]*"("*italic("t")*")")
x2.text = expression(italic("X")[2]*"("*italic("t")*")")
t.text = expression(italic("t"))
x.text = expression(bold("X")*"("*italic("t")*")")

pdf("plot1.pdf",width = 5, height = 6) # width=5 is approx text width in LaTeX
layout.matrix <- matrix(c(1,2,3), nrow = 3, ncol = 1)
layout(mat = layout.matrix, heights = c(1,1,2))

# Plot 1
par(mar=c(0,3,3,3)+c(0,.1,0,0),mgp=c(2,.5,0))
tt <- seq(0,TT,by=Delta.t)
n <- size+1
plot(tt[1:n],X.sample[1:n,1],type="l",xlab="",ylab=x1.text,xlim=c(tt[1],tt[n]),
     xaxt="n",xaxs="i")
# Plot 2
par(mar=c(3,3,0,3)+c(0,.1,0,0),mgp=c(2,.5,0))
plot(tt[1:n],X.sample[1:n,2],type="l",xlab=t.text,ylab=x2.text,xlim=c(tt[1],tt[n]),
     xaxs="i")
# Plot 3
par(mar=c(3,3,3,3)+c(0,.1,0,0),mgp=c(2,.5,0))
n <- 1000+1
plot(tt[1:n],X.sample[1:n,1],type="l",xlab=t.text,ylab=x.text,xlim=c(tt[1],tt[n]),
     ylim=summary(c(X.sample[1:n,1],X.sample[1:n,2]))[c(1,6)],xaxs="i")
lines(tt[1:n],X.sample[1:n,2],col="blue",lty=5)
legend("topleft", legend=c(x1.text,x2.text),
       col=c("black", "blue"), lty=c(1,5))

dev.off()
