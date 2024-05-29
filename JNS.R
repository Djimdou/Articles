# # # # CHAPTER 3: ESTIMATION OF THE GENERATOR

# # For submission to JNS

# # Preparations: loading libraries

# install.packages("xts")
library(xts)
# install.packages("sp")
library(sp)
# install.packages("CASdatasets", repos = "http://cas.uqam.ca/pub/", type="source")
library(CASdatasets)

# data(lossalae)
data(lossalaefull)
lossalae_uncensored = subset(lossalaefull,Censored==0) # restricting to the uncensored cases
attach(lossalae_uncensored)



# # Description

# Kendall tau
res = cor.test(Loss,ALAE, method="kendall")

# scatter plot
plot(log(Loss),log(ALAE))

# Write to csv for LaTex

lossalae_foroverleaf = data.frame(cbind(log(Loss),log(ALAE)))
names(lossalae_foroverleaf) = c("logLoss","logALAE")

write.csv(x=lossalae_foroverleaf,
          file="~/Documents/Articles/PhD - Chapter 3/JNS/lossalae_uncensored.csv",
          row.names = FALSE)

# # # # # Estimation of phi # # # # #

m = 1000

W = 1:m/m #seq(from=1/n,to=1,by=0.01) # grid for the x-axis

# # Pseudo-sample

n = dim(lossalae_uncensored)[1]

V = rep(NA,times=n)
# X1: Loss
# X2: ALAE
for(i in 1:n){
  V[i] = (sum((Loss <= Loss[i])*(ALAE <= ALAE[i])))/(n-1)
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

# plot(Phi_est_int,type="l")


# # # Graph: comparison Phi_est_int, Phi_est_diff, Phi_true

Ylim = range(c(Phi_est_diff,Phi_est_int),na.rm = TRUE)

Xlabels = seq(from=0,to=1,by=0.2)
Ylabels = format(seq(from=0,to=Ylim[2],length.out=5),digits=2)#scientific=TRUE,

plot(W,Phi_est_diff,type='l',col='green',ylim=Ylim,lwd = 2,xlab="",ylab="",xaxt="none",yaxt="none")

lines(W,Phi_est_int,col='red',lwd = 2)

axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
mtext(side=1, line=2, "t", font=2,cex=1.5)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2,hadj=1,padj=0)

legend("topright",legend=c("true density","integral estimator","differential estimator"),lwd = 4,col=c("blue","red", "orange"),lty=1,cex=1.25,bty="n")


# Write to csv for LaTex

phi_estim_ALAE = data.frame(cbind(W,Phi_est_int,Phi_est_diff))
names(phi_estim_ALAE) = c("W","Phi_est_int","Phi_est_diff")

write.csv(x=phi_estim_ALAE,
          file="~/Documents/Articles/PhD - Chapter 3/JNS/phi_estim_ALAE.csv",
          row.names = FALSE)



# # # # # Confidence band for Estimator 1 (differential) # # # # #

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

VarphiN1Diff2_2 = VarphiN1Diff2(KnPrime_2,z=1-epsilon_n)
# Phi_est_diff_At_epsilon_n * KnPrime_2(w=(1-epsilon_n))/(Kn(1-epsilon_n)-(1-epsilon_n))**2

Sigma0Hat_2 = sqrt(VarphiN1Diff2_2/2)

LowerBound_2 = apply(X=cbind(rep(0,times=m),Phi_est_diff*(1-qnorm(p=1-alpha/2)*abs(log(epsilon_n)/sqrt(n))*Sigma0Hat_2)),MARGIN=1,FUN=max)
UpperBound_2 = Phi_est_diff*(1+qnorm(p=1-alpha/2)*abs(log(epsilon_n)/sqrt(n))*Sigma0Hat_2)


# # Graphs

p=2

Ylim = range(c(
  LowerBound_1[-(1:p)],
  LowerBound_2[-(1:p)],
  UpperBound_1[-(1:p)],
  UpperBound_2[-(1:p)],
  Phi_est_diff[-(1:p)]
  ))

plot(W[-(1:p)],Phi_est_diff[-(1:p)],type="l",col="green",ylim=Ylim, ylab="",xlab="w",lwd = 2)
#lines(W[-(1:p)],Phi_est_diff[-(1:p)],col="green",lwd = 2)
lines(W[-(1:p)],LowerBound_1[-(1:p)],col="red",lty="dashed")
lines(W[-(1:p)],UpperBound_1[-(1:p)],col="red",lty="dashed")
lines(W[-(1:p)],LowerBound_2[-(1:p)],col="red",lty="dashed")
lines(W[-(1:p)],UpperBound_2[-(1:p)],col="red",lty="dashed")
legend("topright",legend=c(expression(varphi[n1]),"confidence band"),lwd = c(2,2),col=c("green","red"),lty=c(1,2),cex=1.25,bty="n")


# # Save for Overleaf

confidence_bands_foroverleaf = data.frame(cbind(W,Phi_est_diff,LowerBound_1,UpperBound_1,LowerBound_2,UpperBound_2))
names(confidence_bands_foroverleaf) = c("W","Phi_est_diff","LowerBound_1","UpperBound_1","LowerBound_2","UpperBound_2")

write.csv(x=confidence_bands_foroverleaf,
          file="~/Documents/Articles/PhD - Chapter 3/JNS/confidence_bands_ALAE.csv",
          row.names = FALSE,na="")
