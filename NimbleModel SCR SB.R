NimModel <- nimbleCode({
  lambda.N ~ dunif(0,1000) #expected abundance
  p0 ~ dunif(0,1) #baseline detection probability
  sigma ~ dunif(0,10) #detection spatial scale parameter
  N ~ dpois(lambda.N) #realized abundance
  # comb <- lfactorial(M)-lfactorial(N)+lfactorial((N-n)*step(N-n)) #log-combinatorial terms for z prior
  comb <- lfactorial(M)-lfactorial(N)+lfactorial(N-n) #don't actually need step(N-n) that S&B included
  zerouse ~ dpois(comb) #zeros trick to convert log-combinatorial terms to log(1/combinatorial terms)
  #constrain N between n and M
  dummy.data1 ~ dconstraint(N >= n)
  dummy.data2 ~ dconstraint(N <= M)
  for(i in 1:M){
    z[i] <- step(N-i)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    pd[i,1:J] <- GetDetectionProb(s=s[i,1:2],X=X[1:J,1:2],J=J,sigma=sigma,p0=p0,z=z[i])
    y[i,1:J] ~ dBernoulliVector(pd=pd[i,1:J],K1D=K1D[1:J],z=z[i]) #vectorized obs mod
  }
})
