library(VarianceGamma)
library(boot)

###### SD

################################################################################
# Check 1: Marginal distribution

# Get theoretical pdf and cdf and histogram data

check1.sd <- function(a,alpha1,alpha2,mu1,mu2,S11,S22,S12,eta1,eta2,lambda,X1.time1000,X2.time1000,X.sample.all){
  xx <- seq(-5,5,len=1000)
  pdf1 <- dvg(xx,theta=mu1,sigma=sqrt(S11),nu=alpha1,vgC=eta1)
  pdf2 <- dvg(xx,theta=mu2,sigma=sqrt(S22),nu=alpha2,vgC=eta2)
  
  hist1 <- hist(X1.time1000,seq(-10,10,by=0.05),plot=FALSE)
  xlim1.min <- hist1$breaks[min(which(hist1$counts!=0))]
  xlim1.max <- hist1$breaks[max(which(hist1$counts!=0))+1]
  hist2 <- hist(X2.time1000,seq(-10,10,by=0.05),plot=FALSE)
  xlim2.min <- hist2$breaks[min(which(hist2$counts!=0))]
  xlim2.max <- hist2$breaks[max(which(hist2$counts!=0))+1]
  xlim.min <- min(xlim1.min,xlim2.min)
  xlim.max <- max(xlim1.max,xlim2.max)
  ylim.max <- max(c(hist1$density,hist2$density))
  
  # Plot histogram and pdf
  x1.text = expression(italic("X")[1]*"(1000)")
  x2.text = expression(italic("X")[2]*"(1000)")
  
  # X1
  par(mar=c(3,3,3,3)+c(0,.1,0,0),mgp=c(2,.5,0),cex=2/3,lwd=0.5)
  hist(X1.time1000,seq(-10,10,by=0.05),prob=TRUE,xlim=c(xlim.min,xlim.max),ylim=c(0,1.05*ylim.max),
       xaxs="i",xlab = x1.text, ylab = "Density", main="WVAG-OU",yaxs="i")
  box(lwd=1)
  lines(xx,pdf1,type="l",col="blue",lwd=1)
  # legend(x="topright",legend="Theoretical density",lty=1,col="blue")
  # X2
  par(mar=c(3,3,3,3)+c(0,.1,0,0),mgp=c(2,.5,0),cex=2/3,lwd=0.5)
  hist(X2.time1000,seq(-10,10,by=0.05),prob=TRUE,xlim=c(xlim.min,xlim.max),ylim=c(0,1.05*ylim.max),
       xaxs="i",xlab = x2.text, ylab = "Density", main="",yaxs="i")
  box(lwd=1)
  lines(xx,pdf2,type="l",col="blue",lwd=1)
  # legend(x="topright",legend="Theoretical density",lty=1,col="blue")
  
  kstest.table <- rbind(c(
    ks.test(X1.time1000[1:100],pvg,theta=mu1,sigma=sqrt(S11),nu=alpha1,vgC=eta1)$p.val,
    ks.test(X1.time1000[1:1000],pvg,theta=mu1,sigma=sqrt(S11),nu=alpha1,vgC=eta1)$p.val,
    ks.test(X1.time1000[1:10000],pvg,theta=mu1,sigma=sqrt(S11),nu=alpha1,vgC=eta1)$p.val),
    c(ks.test(X2.time1000[1:100],pvg,theta=mu2,sigma=sqrt(S22),nu=alpha2,vgC=eta2)$p.val,
      ks.test(X2.time1000[1:1000],pvg,theta=mu2,sigma=sqrt(S22),nu=alpha2,vgC=eta2)$p.val,
      ks.test(X2.time1000[1:10000],pvg,theta=mu2,sigma=sqrt(S22),nu=alpha2,vgC=eta2)$p.val))
  print(round(kstest.table,4))
}

################################################################################
# Check 2: Autocorrelation function

check2.sd <- function(a,alpha1,alpha2,mu1,mu2,S11,S22,S12,eta1,eta2,lambda,X.sample.all){
  X1 <- X.sample.all[[1]][,1]
  X2 <- X.sample.all[[1]][,2]
  xx <- seq(0,30,len=1000)
  x1.text = expression(paste("ACF of ",italic("X")[1]))
  x2.text = expression(paste("ACF of ",italic("X")[2]))
  
  par(mar=c(3,3,3,3)+c(0,.1,0,0),mgp=c(2,.5,0),cex=2/3)
  acf(X1,ylab=x1.text,xlab="Lag",main="")
  lines(xx,exp(-lambda*xx),col="red")
  title(main="WVAG-OU")
  # legend(x="topright",legend="Theoretical ACF",lty=1,col="red")
  par(mar=c(3,3,3,3)+c(0,.1,0,0),mgp=c(2,.5,0),cex=2/3)
  acf(X2,ylab=x2.text,xlab="Lag",main="")
  lines(xx,exp(-lambda*xx),col="red")
  # legend(x="topright",legend="Theoretical ACF",lty=1,col="red")
}

################################################################################
# Check 3: Cov(X_1(1000),X_2(1000))

check3.sd <- function(a,alpha1,alpha2,mu1,mu2,S11,S22,S12,eta1,eta2,lambda,X1.time1000,X2.time1000){
  # Theoretical cov
  true.cov <- a*(min(alpha1,alpha2)*S12+alpha1*alpha2*mu1*mu2)
  # Sample cov and bootstrap CI
  set.seed(402935)
  sample.cov <- cov(X1.time1000,X2.time1000)
  boot.cov <- boot(cbind(X1.time1000,X2.time1000), cov.boot, R=10000)
  boot.cov.ci<- boot.ci(boot.cov,type="perc")
  
  print(round(c(true.cov,sample.cov,boot.cov.ci$percent[4:5]),4))
}
cov.boot <- function(d, i){
  cov(d[i,1],d[i,2])
}

###### BDLP

################################################################################
# Check 1: Marginal distribution

check1.bdlp <- function(a,alpha1,alpha2,mu1,mu2,S11,S22,S12,eta1,eta2,lambda,X1.time1000,X2.time1000,X.sample.all){
  # Get theoretical pdf and cdf and histogram data
  step1 <- median(diff(Grid.x[1,]))
  step2 <- median(diff(Grid.x[2,]))
  
  pdf1 <- step2*apply(P,1,sum)
  cdf1 <- approxfun(x=Grid.x[1,],y=step1*cumsum(pdf1),yleft=0,yright=1)
  hist1 <- hist(X1.time1000,seq(-10,10,by=0.05),plot=FALSE)
  xlim1.min <- hist1$breaks[min(which(hist1$counts!=0))]
  xlim1.max <- hist1$breaks[max(which(hist1$counts!=0))+1]
  pdf2 <- step1*apply(P,2,sum)
  cdf2 <- approxfun(x=Grid.x[2,],y=step2*cumsum(pdf2),yleft=0,yright=1)
  hist2 <- hist(X2.time1000,seq(-10,10,by=0.05),plot=FALSE)
  xlim2.min <- hist2$breaks[min(which(hist2$counts!=0))]
  xlim2.max <- hist2$breaks[max(which(hist2$counts!=0))+1]
  xlim.min <- min(xlim1.min,xlim2.min)
  xlim.max <- max(xlim1.max,xlim2.max)
  ylim.max <- max(c(hist1$density,hist2$density))
  
  # Plot histogram and pdf
  x1.text = expression(italic("X")[1]*"(1000)")
  x2.text = expression(italic("X")[2]*"(1000)")
  
  # X1
  par(mar=c(3,3,3,3)+c(0,.1,0,0),mgp=c(2,.5,0),cex=2/3,lwd=0.5)
  hist(X1.time1000,seq(-10,10,by=0.05),prob=TRUE,xlim=c(xlim.min,xlim.max),ylim=c(0,1.05*ylim.max),
       xaxs="i",xlab = x1.text, ylab = "Density", main="OU-WVAG",yaxs="i")
  box(lwd=1)
  lines(Grid.x[1,],pdf1,type="l",col="blue",lwd=1)
  # legend(x="topright",legend="Theoretical density",lty=1,col="blue")
  # X2
  par(mar=c(3,3,3,3)+c(0,.1,0,0),mgp=c(2,.5,0),cex=2/3,lwd=0.5)
  hist(X2.time1000,seq(-10,10,by=0.05),prob=TRUE,xlim=c(xlim.min,xlim.max),ylim=c(0,1.05*ylim.max),
       xaxs="i",xlab = x2.text, ylab = "Density", main="",yaxs="i")
  box(lwd=1)
  lines(Grid.x[2,],pdf2,type="l",col="blue",lwd=1)
  # legend(x="topright",legend="Theoretical density",lty=1,col="blue")
  
  # KS test table
  kstest.table <- rbind(c(ks.test(X1.time1000[1:100],cdf1)$p.val,
                          ks.test(X1.time1000[1:1000],cdf1)$p.val,
                          ks.test(X1.time1000[1:10000],cdf1)$p.val),
                        c(ks.test(X2.time1000[1:100],cdf2)$p.val,
                          ks.test(X2.time1000[1:1000],cdf2)$p.val,
                          ks.test(X2.time1000[1:10000],cdf2)$p.val))
  print(round(kstest.table,4))
}


################################################################################
# Check 2: Autocorrelation function

check2.bdlp <- function(a,alpha1,alpha2,mu1,mu2,S11,S22,S12,eta1,eta2,lambda,X.sample.all){
  X1 <- X.sample.all[[1]][,1]
  X2 <- X.sample.all[[1]][,2]
  xx <- seq(0,30,len=1000)
  x1.text = expression(paste("ACF of ",italic("X")[1]))
  x2.text = expression(paste("ACF of ",italic("X")[2]))
  
  par(mar=c(3,3,3,3)+c(0,.1,0,0),mgp=c(2,.5,0),cex=2/3)
  acf(X1,ylab=x1.text,xlab="Lag",main="")
  lines(xx,exp(-lambda*xx),col="red")
  title(main="OU-WVAG")
  # legend(x="topright",legend="Theoretical ACF",lty=1,col="red")
  par(mar=c(3,3,3,3)+c(0,.1,0,0),mgp=c(2,.5,0),cex=2/3)
  acf(X2,ylab=x2.text,xlab="Lag",main="")
  lines(xx,exp(-lambda*xx),col="red")
  # legend(x="topright",legend="Theoretical ACF",lty=1,col="red")
}

################################################################################
# Check 3: Cov(X_1(1000),X_2(1000))

check3.bdlp <- function(a,alpha1,alpha2,mu1,mu2,S11,S22,S12,eta1,eta2,lambda,X1.time1000,X2.time1000){
  # Theoretical cov
  E2 <- (exp(2*lambda*Delta.t)-1)/2
  b <- exp(-lambda*Delta.t)
  true.cov <- E2*a*(min(alpha1,alpha2)*S12+alpha1*alpha2*mu1*mu2)*b^2/(1-b^2)
  # Sample cov and bootstrap CI
  set.seed(572957)
  sample.cov <- cov(X1.time1000,X2.time1000)
  boot.cov <- boot(cbind(X1.time1000,X2.time1000), cov.boot, R=10000)
  boot.cov.ci<- boot.ci(boot.cov,type="perc")
  
  print(round(c(true.cov,sample.cov,boot.cov.ci$percent[4:5]),4))
}

###############################################################################

###### SD

load("SD Spec Sample.Rdata") # Loads the required variables and functions

set.seed(890836)

# Copy sample code, simulate more sample paths: 10000

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

h2 <- lambda* Delta.t

# Sample Z
size2 <- 10000
X.sample.all <- list()

# Simulate Y
return.process <- sim.vag(a,alpha1,alpha2,mu1,mu2,S11,S22,rho,t.max=size2)
Y1 <- return.process[,1]+eta1
Y2 <- return.process[,2]+eta2

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

# Obtain X(1000)
X1.time1000.sd <- unlist(lapply(X.sample.all, function(x){x[1000,1]}))
X2.time1000.sd <- unlist(lapply(X.sample.all, function(x){x[1000,2]}))
X.sample.all.sd <- X.sample.all

###### BDLP

load("BDLP Spec Sample.Rdata") # Loads the required variables and functions

set.seed(328536)

# Copy sample code, simulate more sample paths: 10000

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


# Sample Y
h <- lambda*Delta
size2 <- 10000
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

# Obtain X(1000)
X1.time1000.bdlp <- unlist(lapply(X.sample.all, function(x){x[1000,1]}))
X2.time1000.bdlp <- unlist(lapply(X.sample.all, function(x){x[1000,2]}))
X.sample.all.bdlp <- X.sample.all

################################################################################

# Run checks and get output

pdf("plot3.pdf",width = 6, height = 6)
layout.matrix <- matrix(c(1,2,3,4), nrow = 2, ncol = 2)
layout(mat = layout.matrix, heights = c(1,1), widths=c(1,1))
check1.sd(a,alpha1,alpha2,0,0,S11,S22,S12,eta1,eta2,lambda,X1.time1000.sd,X2.time1000.sd,X.sample.all.sd)
check1.bdlp(a,alpha1,alpha2,mu1,mu2,S11,S22,S12,eta1,eta2,lambda,X1.time1000.bdlp,X2.time1000.bdlp,X.sample.all.bdlp)
dev.off()

pdf("plot4.pdf",width = 6, height = 6)
layout.matrix <- matrix(c(1,2,3,4), nrow = 2, ncol = 2)
layout(mat = layout.matrix, heights = c(1,1), widths=c(1,1))
check2.sd(a,alpha1,alpha2,0,0,S11,S22,S12,eta1,eta2,lambda,X.sample.all.sd)
check2.bdlp(a,alpha1,alpha2,mu1,mu2,S11,S22,S12,eta1,eta2,lambda,X.sample.all.bdlp)
dev.off()

check3.sd(a,alpha1,alpha2,0,0,S11,S22,S12,eta1,eta2,lambda,X1.time1000.sd,X2.time1000.sd)
check3.bdlp(a,alpha1,alpha2,mu1,mu2,S11,S22,S12,eta1,eta2,lambda,X1.time1000.bdlp,X2.time1000.bdlp)




