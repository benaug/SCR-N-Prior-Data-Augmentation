e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.SCR <-
  function(N=NA,lam0=NA,p0=NA,sigma=NA,theta=NA,lambda=NA,K=NA,X=X,buff=NA){
    #simulate a population of activity centers
    xlim <- c(min(X[,1]) - buff,max(X[,1]) + buff)
    ylim <- c(min(X[,2]) - buff,max(X[,2]) + buff)
    s <- cbind(runif(N, xlim[1],xlim[2]), runif(N,ylim[1],ylim[2]))
    D <- e2dist(s,X)
    J <- nrow(X)
    #Capture individuals
    y <- array(0,dim=c(N,J,K))
    pd <- p0*exp(-D*D/(2*sigma*sigma))
    for(i in 1:N){
      for(j in 1:J){
        for(k in 1:K){
          y[i,j,k] <- rbinom(1,1,pd[i,j])
        }
      }
    }
    #discard uncaptured inds
    caught <- which(apply(y,c(1),sum)>0)
    y <- y[caught,,]
    if(K==1){
      y <- array(y,dim=c(dim(y),1))
    }
    n <- length(caught)
    out <- list(y=y,n=n,X=X,K=K,buff=buff,xlim=xlim,ylim=ylim)
    return(out)
  }