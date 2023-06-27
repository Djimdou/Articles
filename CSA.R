# # # # CSA paper: Distribution of the joint survival function of an Archimedean copula

# # # Preliminary

#install.packages("GoFKernel")
library(GoFKernel) # for inverse
#install.packages("copula")
library(copula) # for claytonCopula

# # Clayton

# Simulating values

# m = 1000

# W = 1:m/m # grid for the x-axis

Theta0 = 1

Theta = seq(from=0.1,to=5,by=0.1) # c(0.1,1,3,5)  # 7 seems not good a value for the numerical approximations later

Z = seq(from=0.01, to=0.99,by=0.005) # 
A = Z[Z<=1/2] # seq(from=0.01,to=1/2,by=0.005)  # 
M = 10**15

cop <- claytonCopula(param=Theta0,dim=2)

# Sampling X and Y

# n = 1000 # size of sample

alpha = 2 # Weibull distribution shape parameter
beta = 2 # Weibull distribution scale parameter

MyCopula <- mvdc(copula=cop, # copula for (F(X), F(Y))
                  margins=c("weibull","weibull"), # Weibull distribution for margins X and Y
                  paramMargins=list(shape=alpha,scale=beta)) # alpha:shape, beta:scale

set.seed(1)

XY <- rMvdc(n=length(Z),MyCopula)
X <- XY[,1]
Y <- XY[,2]

U1 = pweibull(X, shape=alpha, scale = beta, lower.tail = TRUE, log.p = FALSE)
U2 = pweibull(Y, shape=alpha, scale = beta, lower.tail = TRUE, log.p = FALSE)

V = Cn_star = rep(NA,times=length(Z))

for(i in 1:length(Z)){
  V[i] = (sum((U1 <= U1[i])*(U2 <= U2[i])))/(length(Z)-1)
  #Cn_star[i] = (sum((U1 > U1[i])*(U2 > U2[i])))/(length(Z))
}

# Empirical distribution of Cn_star

# Kn_star from simulated values

density_Cn_star <- matrix(NA,ncol=length(Z),nrow=length(Theta))

# for(t in 1:length(Theta)){
#   S <- rbeta(n=length(Z),shape1=1,shape2=1/Theta[t])
#   T <- (S**(-Theta[t])-1)/Theta[t]
#   Cn_star <- 1-(Theta[t]*W*T+1)**(-1/Theta[t])-(Theta[t]*(1-W)*T+1)**(-1/Theta[t])+(Theta[t]*T+1)**(-1/Theta[t])
#   #density_C_star[t,] <- density(C_star,from=0,to=1,n=length(Z))$y
# }

W <- runif(n=length(Z), min = 0, max = 1)
S = rep(NA,times=length(Z))

# seems good if Theta[t] >= 0.5

for(t in 1:length(Theta)){
  
  Const = 1+Theta[t]
  i=1
  
  while(i <= length(Z)){
    
    y = rbeta(n=1,shape1=1,shape2=1/Theta[t])
    ratio <- (1-y**Theta[t])/(Const*(1-y)**(1/Theta[t]-1))
    
    while(runif(n=1, min = 0, max = 1) >= ratio)
    {
      y = rbeta(n=1,shape1=1,shape2=1/Theta[t])
      ratio <- (1-y**Theta[t])/(Const*(1-y)**(1/Theta[t]-1))
    }
    
    S[i]=y
    i=i+1
    
  }
  
  T <- (S**(-Theta[t])-1)/Theta[t]
  
  Cn_star <- 1-(Theta[t]*W*T+1)**(-1/Theta[t])-(Theta[t]*(1-W)*T+1)**(-1/Theta[t])+(Theta[t]*T+1)**(-1/Theta[t])
  density_Cn_star[t,] <- density(Cn_star,from=0,to=1,n=length(Z))$y
}

Kn_star = ecdf(Cn_star)

Kn = ecdf(V)

# # Distance

# Cn_star = Cn_star[order(Cn_star)]

# Kstar

inv_a_z = integrand = diff_a_z = matrix(NA,ncol=length(Z),nrow=length(A)) 

k_star = matrix(NA,nrow=length(Theta),ncol=length(Z)) # density of C_star

for(t in 1:length(Theta)){
  for(i in 1:length(A)){
    psi_inverse = inverse(function (z) {1-(Theta[t]*A[i]*z+1)**(-1/Theta[t])-(Theta[t]*(1-A[i])*z+1)**(-1/Theta[t])+(Theta[t]*z+1)**(-1/Theta[t])}, 0, M)
    for(j in 1:length(Z)){
      inv_a_z[i,j] = unlist(psi_inverse(Z[j]))
      diff_a_z[i,j] = 1/(A[i]*(Theta[t]*A[i]*inv_a_z[i,j]+1)**(-1/Theta[t]-1)+(1-A[i])*(Theta[t]*(1-A[i])*inv_a_z[i,j]+1)**(-1/Theta[t]-1)-(Theta[t]*inv_a_z[i,j]+1)**(-1/Theta[t]-1))
      integrand[i,j] = (Theta[t]+1)*inv_a_z[i,j]*(Theta[t]*inv_a_z[i,j]+1)**(-2-1/Theta[t]) * diff_a_z[i,j]
    }
  }
  
  for(j in 1:length(Z)){
    k_star[t,j] = 2*integrate(approxfun(x=A,y=integrand[,j]), lower = min(A), upper = max(A))[[1]]
  }
}

K_star = matrix(NA,nrow=length(Theta),ncol=length(Z))  # distribution of C_star
K_star[,1] = rep(0,times=length(Theta))

for(t in 1:length(Theta)){
  for(j in 2:length(Z)){
    K_star[t,j] = integrate(approxfun(x=Z,y=k_star[t,]), lower = min(Z), upper = Z[j])[[1]]
  }
}

# K

K = matrix(NA,nrow=length(Theta),ncol=length(Z))

for(t in 1:length(Theta)){
  K[t,] = Z*(1+(1/Theta[t])*(1-Z**Theta[t]))
}

dn = dn_star = rep(NA,times=length(Theta))

for(t in 1:length(Theta)){
  dn_star[t] = max(abs(Kn_star(Z) - K_star[t,]))
  dn[t] = max(abs(Kn(Z)-K[t,]))
}

Thetan = Theta[which.min(dn)]
Thetan_star = Theta[which.min(dn_star)]
# plot(Z,Kn_star(Z),type="l",col="blue");lines(Z,K_star[10,],col="green")

Ylim = range(c(dn_star,dn))

plot(x=Theta,y=dn_star,type="l",col="blue",ylim=Ylim,xlab=expression(theta),ylab="",lwd = 3)
lines(x=Theta,y=dn,col="green",lwd = 3)
abline(v=Theta0,col='red',lwd = 2, lty=2)
legend("topright",legend=c(expression(d[n]^"*"),expression(d[n]),bquote(theta[0] == .(Theta0))),lwd = c(3,3,2),col=c("blue", "green","red"),lty=c(1,1,2),cex=1.25,bty="n")

# Densities

#density_Cn_star <- density(Cn_star,from=0,to=1,n=length(Z))$y

# plot(Z,density_Cn_star,type="l",col='red')
# lines(Z,k_est[which.min(dn_star),],col="blue")

t=which(mapply(function(x, y) {isTRUE(all.equal(x, y))}, Theta, Theta0))
Ylim = range(c(k_star[t,],density_Cn_star[t,]),xlab="z",ylab="density")

plot(Z,k_star[t,],type="l",col="blue")
lines(Z,density_Cn_star[t,],col='red')


# Write to csv for LaTex plots

write.csv(x=data.frame(cbind(Theta,dn_star,dn,Theta0=rep(Theta0,times=length(Theta)),Thetan=rep(Thetan,times=length(Theta)),Thetan_star=rep(Thetan_star,times=length(Theta)))),
          file=paste0("/Users/magloireloudeguidjimdou/Documents/Articles/PhD - Chapter 2/KSStatistics-2.csv"),
          row.names = FALSE)

# set.seed(NULL)






Densities = read.csv(file="/Users/magloireloudeguidjimdou/Documents/Articles/PhD - Chapter 2/clayton_density_simu.csv",header = TRUE)

max(abs(Densities[,"V2"]-Densities[,"X1"]))
max(abs(Densities[,"V3"]-Densities[,"X2"]))

Ylim=range(Densities[,c("V2","X1")])
plot(Densities[,c("Z","V2")],type="l",ylim=Ylim)
lines(Densities[,c("Z","X1")])

Distances1 = read.csv(file="/Users/magloireloudeguidjimdou/Documents/Articles/PhD - Chapter 2/KSStatistics-1_25.csv",header = TRUE)
Distances2 = read.csv(file="/Users/magloireloudeguidjimdou/Documents/Articles/PhD - Chapter 2/KSStatistics-3.csv",header = TRUE)

Ylim=range(Distances1[,c("dn_star","dn")])
plot(Distances1[,c("Theta","dn_star")],type="l",ylim=Ylim)
lines(Distances1[,c("Theta","dn")])

Distances1[which.min(Distances1[,"dn_star"]),"Theta"]
Distances1[which.min(Distances1[,"dn"]),"Theta"]

Distances2[which.min(Distances2[,"dn_star"]),"Theta"]
Distances2[which.min(Distances2[,"dn"]),"Theta"]



