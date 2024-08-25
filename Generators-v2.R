# # # # CHAPTER 3: ESTIMATION OF THE GENERATOR


# # Preps: loading libraries

library(copula) # for claytonCopula


# # Simulating values from Archimedean copulas

m = 1000

W = 1:m/m #seq(from=1/n,to=1,by=0.01) # grid for the x-axis

CopulaName <- "Clayton" # one of c("Indep", "Clayton", "AMH", "Frank", Gumbel)
  
if(CopulaName == "Indep"){
  
  # Independence copula
  
  cop <- indepCopula(dim=2)
  Phi_true = -log(W)
  
}

if(CopulaName == "Clayton"){
  
  # Clayton copula
  
  theta = 1
  cop <- claytonCopula(param=theta,dim=2)
  Phi_true = (W**(-theta)-1)/theta
  
}

if(CopulaName == "AMH"){
  
  # AMH copula
  
  theta = 0.5
  cop <- amhCopula(param=theta,dim=2)
  Phi_true = log((1-theta*(1-W))/W)/(1-theta)
  
}

if(CopulaName == "Frank"){
  
  # Frank
  
  theta = -1
  cop <- frankCopula(param=theta,dim=2)
  Phi_true = -((exp(theta)-1)/theta)*log((exp(-theta*W)-1)/(exp(-theta)-1))
  
}

if(CopulaName == "Gumbel"){
  
  # Gumbel (counter-example)
  
  theta = 10
  cop <- gumbelCopula(param=theta,dim=2)
  Phi_true = (-log(W))**theta
  
}


# # Sampling X and Y

n = 500 # size of sample

alpha = 2 # Weibull distribution shape parameter
beta = 2 # Weibull distribution scale parameter

MyCopula <- mvdc(copula=cop, # copula for (F(X), F(Y))
                 margins=c("weibull","weibull"), # Weibull distribution for margins X and Y
                 paramMargins=list(shape=alpha,scale=beta)) # alpha:shape, beta:scale

set.seed(1) # 1, 3
X1X2 <- rMvdc(n=n,MyCopula)
X1 <- X1X2[,1]
X2 <- X1X2[,2]

# # Pseudo-sample

V = rep(NA,times=n)

for(i in 1:n){
  V[i] = (sum((X1 <= X1[i])*(X2 <= X2[i])))/(n-1)
}

# # Empirical distribution

Kn = ecdf(V)
Vn = V[order(V)]

# # # Estimator 1: differential equation

Phi_est = rep(0,times=n)

for(i in 1:n){
  if(Vn[i] < 1-1/n){
    T = unique(c(Vn[(Vn[i] <= Vn) & (Vn <= 1-1/n)],1-1/n))
    ProdTerm = rep(NA,times=length(T)-1)
    ProdTerm = (Kn(T[-length(T)])-T[-length(T)])/(Kn(T[-length(T)])-T[-1])
    Phi_est[i] = (1/n)*abs(prod(ProdTerm))
  }
}

Phi_est_diff = approx(x=c(0,Vn,1), y=c(Phi_est[1],Phi_est,0), xout=W, method="constant", ties = "ordered")$y 

# plot(Phi_est_diff,type="l")

# # # Estimator 2: integral equation

h = rep(NA,times=n)
h[n] = 1 

for(i in (n-1):1){
  h[i] = (1/(n*Kn(i/n)-i))*sum(h[(i+1):n])
}

# Step interpolation

cond = ((h>=0) & is.finite(h)) # may have 0's, when Kn(i/n)-i/n <= 0

h.interp = approx(x=c(0,Vn[cond],1), y=c(h[cond][1],h[cond],1), xout=W, method="constant", ties = mean)$y 

Phi_est_int = (sum(h.interp) - cumsum(c(0,h.interp[-m])))/m


# # # Graph: comparison Phi_est_int, Phi_est_diff, Phi_true

Ylim = range(c(Phi_est_int,Phi_est_diff,Phi_true),na.rm = TRUE)

Xlabels = seq(from=0,to=1,by=0.2)
Ylabels = format(seq(from=0,to=Ylim[2],length.out=5),digits=2)#scientific=TRUE,

plot(W,Phi_true,type='l',col='blue',ylim=Ylim,lwd = 2,xlab="",ylab="",xaxt="none",yaxt="none")

lines(W,Phi_est_int,col='red',lwd = 2)
lines(W,Phi_est_diff,col='orange',lwd = 2)

axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
mtext(side=1, line=2, "t", font=2,cex=1.5)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2,hadj=1,padj=0)

legend("topright",legend=c("true density","integral estimator","differential estimator"),lwd = 4,col=c("blue","red", "orange"),lty=1,cex=1.25,bty="n")


# Write to csv for LaTex

write.csv(x=data.frame(cbind(W,Phi_true,Phi_est_int,Phi_est_diff)),
          file="C:/Users/mloudegu/Documents/Thesis/phi_estim_frank_500.csv",
          row.names = FALSE)



# # # # # # # # #  Limiting behavior of Estimator 1 (differential) # # # # # # # # # # # # # 


# # Simulating values from the copulas

SimuCopula <- function(CopulaName,W){
  
  # # Description
  # Function to simulate values of Archimedean copula and of the generator. 
  # Limited to the independence, Clayton, Frank and AMH copulas
  
  # # Inputs
  # CopulaName: name of the copula. One of 'Indep','Clayton','AMH' or 'Frank'
  # W: vector of values in (0,1)
  
  # # Outputs
  # A list with two elements: 
  # * cop: copula 
  # * gen_true: generator (phi)
  
  if(CopulaName == 'Indep'){
    # Independence copula
    
    cop <- indepCopula(dim=2)
    Phi_true = -log(W)
  }
  
  if(CopulaName == 'Clayton'){
    # Clayton copula
    
    theta = 1
    cop <- claytonCopula(param=theta,dim=2)
    Phi_true = (W**(-theta)-1)/theta
  }
  
  if(CopulaName == 'AMH'){
    # AMH copula
    
    theta = 0.5
    cop <- amhCopula(param=theta,dim=2)
    Phi_true = log((1-theta*(1-W))/W)/(1-theta)
  }
  
  if(CopulaName == 'Frank'){
    # Frank copula
    
    theta = -1
    cop <- frankCopula(param=theta,dim=2)
    Phi_true = -((exp(theta)-1)/theta)*log((exp(-theta*W)-1)/(exp(-theta)-1))
  }
  return(list(cop=cop,gen_true=Phi_true))
}

# Grid for the W values
W = seq(from=0.05,to=0.95,by=0.01) # seq(from=0.2,to=0.9,by=0.2) 

# CopulaNames <- c('Indep','Clayton','AMH','Frank')

CopulaName <- 'Clayton' # select the copula from the list c('Indep','Clayton','AMH','Frank')
cop = SimuCopula(CopulaName,W)$cop
Phi_true = SimuCopula(CopulaName,W)$gen_true

# # Sampling X and Y

N <- seq(from=100,to=3000,by=500) # the 0 will be ignored. N should have at least 2 elements.

ConvergenceMAtrix <- matrix(rep(0,times=length(W)*length(N)),ncol=length(W)) # colnames = w, rownames = n
colnames(ConvergenceMAtrix) <- paste('w',1:length(W),sep='')

# Building the copula object

alpha = 2 # Weibull distribution shape parameter
beta = 2 # Weibull distribution scale parameter

MyCopula <- mvdc(copula=cop, # copula for (F(X), F(Y))
                 margins=c("weibull","weibull"), # Weibull distribution for margins X and Y
                 paramMargins=list(shape=alpha,scale=beta)) # alpha:shape, beta:scale

for(n.index in 1:length(N)){
  
  #n.index <- 1
  n <- N[n.index]
  
  set.seed(5) # 5 for Clayton, AMH, indep
  X1X2 <- rMvdc(n=n,MyCopula)
  X1 <- X1X2[,1]
  X2 <- X1X2[,2]
  
  # # Pseudo-sample
  
  V = rep(NA,times=n)
  
  for(i in 1:n){
    V[i] = (sum((X1 <= X1[i])*(X2 <= X2[i])))/(n-1)
  }
  
  # # Empirical distribution
  
  Kn = ecdf(V)
  Vn = V[order(V)]
  
  # # # Estimator 1: differential equation
  
  Phi_est = rep(0,times=n)
  
  for(i in 1:n){
    if(Vn[i] < 1-1/n){
      T = unique(c(Vn[(Vn[i] <= Vn) & (Vn <= 1-1/n)],1-1/n))
      ProdTerm = rep(NA,times=length(T)-1)
      ProdTerm = (Kn(T[-length(T)])-T[-length(T)])/(Kn(T[-length(T)])-T[-1])
      Phi_est[i] = (1/n)*abs(prod(ProdTerm))
    }
  }
  
  Phi_est_diff = approx(x=c(0,Vn,1), y=c(Phi_est[1],Phi_est,0), xout=W, method="constant", ties = "ordered")$y 
  ConvergenceMAtrix[n.index,] <- (sqrt(n)/log(1/n, base = exp(1)))*(Phi_est_diff-Phi_true)/Phi_true
}

# ConvergenceMAtrix may contain infinite values. Replacing such values.
for(w.index in 1:length(W)){
  if(sum(is.infinite(ConvergenceMAtrix[,w.index])) != 0){
    ConvergenceMAtrix[,w.index][is.infinite(ConvergenceMAtrix[,w.index])] <- min(ConvergenceMAtrix[,w.index][is.finite(ConvergenceMAtrix[,w.index])])-0.25
  }
}


# colors to be used to graph
MyColors <- c('aquamarine4','bisque1','blue','blueviolet','brown1','burlywood3','darkgreen','dodgerblue4','deeppink2','gray1')

# axes ranges
#Xlim = c(0,500)
Ylim = range(ConvergenceMAtrix,na.rm = TRUE)

# axis labels
Xlabels = seq(from=0.05,to=0.95,by=0.1) # N[-1] #seq(from=0,to=500,by=50)
Ylabels = format(seq(from=Ylim[1],to=Ylim[2],length.out=5),digits=2)

#plot(N[-1],ConvergenceMAtrix[,1],type='l',col=MyColors[1],ylim=Ylim,lwd = 2,xlab="",ylab="",xaxt="none",yaxt="none")
plot(W,ConvergenceMAtrix[1,],type='l',col=MyColors[1],ylim=Ylim,lwd = 2,xlab="",ylab="",xaxt="none",yaxt="none")

for(n.index in 2:length(N)){
  #lines(N[-1],ConvergenceMAtrix[,w.index],col=MyColors[w.index],lwd = 2)
  lines(W,ConvergenceMAtrix[n.index,],col=MyColors[n.index],lwd = 2)
}

axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
mtext(side=1, line=3, "w", font=2,cex=1.25)#
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2,hadj=1,padj=0)

legend("bottom",legend=paste('n=',N,sep=''),lwd = 3,col=MyColors[1:length(W)],lty=1,cex=1.25,bty="n")

# Write to csv for LaTex

write.csv(x=data.frame(cbind(W,t(ConvergenceMAtrix))),
          file=paste0("/Users/magloireloudeguidjimdou/Documents/Articles/PhD - Chapter 3/phi_conv_w_",tolower(CopulaName),".csv"),
          row.names = FALSE,na="")



# # # # # # # # # # # # # # Confidence band for Estimator 1 # # # # # # # # # # # # # # 

alpha = 0.05

epsilon_n = 1/n
Phi_est_diff_At_epsilon_n = approx(x=c(0,Vn,1), y=c(Phi_est[1],Phi_est,0), xout=(1-epsilon_n), method="constant", ties = "ordered")$y 

#plot(W,Phi_est_diff,type="l")
#points(x=(1-epsilon_n),y=Phi_est_diff_At_epsilon_n,col="red")

# # Kn'(.) is approximated by 1/n

KnPrime_1 = function(w){return(1/n)}

VarphiN1Diff2 = function(KnPrime,z){
  VarphiN = (Phi_est_diff_At_epsilon_n * KnPrime(w=z))/(Kn(z)-(z))**2
  return(VarphiN)
  }

VarphiN1Diff2_1 = VarphiN1Diff2(KnPrime_1,z=1-epsilon_n)
# Phi_est_diff_At_epsilon_n * KnPrime_1(w=1-epsilon_n)/(Kn(1-epsilon_n)-(1-epsilon_n))**2

Sigma0Hat_1 = sqrt(VarphiN1Diff2_1/2)

LowerBound_1 = apply(X=cbind(rep(0,times=m),Phi_est_diff*(1-qnorm(p=1-alpha/2)*abs(log(epsilon_n)/sqrt(n))*Sigma0Hat_1)),MARGIN=1,FUN=max)
UpperBound_1 = Phi_est_diff*(1+qnorm(p=1-alpha/2)*abs(log(epsilon_n)/sqrt(n))*Sigma0Hat_1)

# # Kernel smoother for Kn'(.)

bn = 0.75
KnPrime_2 = function(w,axis=Vn,bandwidth=bn){
  fhat = (sum(dnorm(x=(w-axis)/bn))/bn)*(1/length(axis))
  return(fhat)
}

# KnPrime = apply(X=matrix(W,ncol=1),FUN=L,MARGIN=1)
# plot(W,KnPrime,type="l",col="blue")

VarphiN1Diff2_2 = VarphiN1Diff2(KnPrime_2,z=1-epsilon_n)
# Phi_est_diff_At_epsilon_n * KnPrime_2(w=(1-epsilon_n))/(Kn(1-epsilon_n)-(1-epsilon_n))**2

Sigma0Hat_2 = sqrt(VarphiN1Diff2_2/2)

LowerBound_2 = apply(X=cbind(rep(0,times=m),Phi_est_diff*(1-qnorm(p=1-alpha/2)*abs(log(epsilon_n)/sqrt(n))*Sigma0Hat_2)),MARGIN=1,FUN=max)
UpperBound_2 = Phi_est_diff*(1+qnorm(p=1-alpha/2)*abs(log(epsilon_n)/sqrt(n))*Sigma0Hat_2)

# # Graphs

p=10

Ylim = range(c(
               #LowerBound_1[-(1:p)],
               LowerBound_2[-(1:p)],
               #UpperBound_1[-(1:p)],
               UpperBound_2[-(1:p)],
               Phi_est_diff[-(1:p)],
               Phi_true[-(1:p)]))

plot(W[-(1:p)],Phi_true[-(1:p)],type="l",col="blue",ylim=Ylim, ylab="",xlab="w",lwd = 2) # 
lines(W[-(1:p)],Phi_est_diff[-(1:p)],col="green",lwd = 2)
#lines(W[-(1:p)],LowerBound_1[-(1:p)],col="red",lty="dashed",lwd = 2)
#lines(W[-(1:p)],UpperBound_1[-(1:p)],col="red",lty="dashed",lwd = 2)
lines(W[-(1:p)],LowerBound_2[-(1:p)],col="orange",lty="dashed")
lines(W[-(1:p)],UpperBound_2[-(1:p)],col="orange",lty="dashed")
legend("topright",legend=c(expression(varphi),expression(varphi[n1]),"confidence band"),lwd = c(2,2,2),col=c("blue", "green","red"),lty=c(1,1,2),cex=1.25,bty="n")


# # Save for Overleaf

DonneesOverleaf = data.frame(cbind(W,Phi_true,Phi_est_diff,LowerBound_1,UpperBound_1,LowerBound_2,UpperBound_2))
names(DonneesOverleaf) = c("W","Phi_true","Phi_est_diff","LowerBound_1","UpperBound_1","LowerBound_2","UpperBound_2")

write.csv(x=DonneesOverleaf,
          file=paste0("/Users/magloireloudeguidjimdou/Documents/Articles/PhD - Chapter 3/confidence_bands_",tolower(CopulaName),".csv"),
          row.names = FALSE,na="")
