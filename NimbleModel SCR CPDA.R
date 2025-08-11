NimModel <- nimbleCode({
  lambda.N ~ dunif(0,1000) #expected abundance
  p0 ~ dunif(0,1) #baseline detection probability
  sigma ~ dunif(0,10) #detection spatial scale parameter
  N ~ dpois(lambda.N) #realized abundance
  for(i in 1:M){
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    pd[i,1:J] <- GetDetectionProb(s=s[i,1:2],X=X[1:J,1:2],J=J,sigma=sigma,p0=p0,z=z[i])
    y[i,1:J] ~ dBernoulliVector(pd=pd[i,1:J],K1D=K1D[1:J],z=z[i]) #vectorized obs mod
  }
  z.up <- 1 #used to record adaptation if using zSampler3
})
#custom Metropolis-Hastings update for N/z