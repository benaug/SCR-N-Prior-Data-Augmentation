#This script allows you to compare 3 versions of the random-z approach

library(nimble)
library(coda)
source("sim.SCR.R")
source("NimbleModel SCR NPDA.R")
source("NimbleFunctions SCR NPDA.R")
source("sSampler.R")

#simulate some data
N <- 50
p0 <- 0.25
sigma <- 0.5
K <- 5
buff <- 2 #state space buffer. Should be at least 3 sigma.
X <- as.matrix(expand.grid(3:11,3:11))

#Simulate some data. setting seed so we can run 3 different zSamplers and confirm they produce same posterior
set.seed(123930963)
data <- sim.SCR(N=N,p0=p0,sigma=sigma,K=K,X=X,buff=buff)

#collapse 3D data to 2D for efficiency
y2D <- apply(data$y,c(1,2),sum)

#trap operation matrix
J <- nrow(X)
K1D <- rep(data$K,J) #assuming perfect operation here

#Augment data (collapsing to 2D for efficiency)
M <- 100
y2D <- matrix(0,M,J)
y2D[1:data$n,] <- apply(data$y,c(1,2),sum)

#inits for nimble
#must init N=sum(z)
z.init <- rep(0,M)
z.init[1:data$n] <- 1
N.init <- sum(z.init)
s.init <- cbind(runif(M, data$xlim[1],data$xlim[2]), runif(M,data$ylim[1],data$ylim[2])) #assign random locations
idx <- which(rowSums(y2D)>0) #switch for those actually caught
for(i in idx){
  trps <- matrix(data$X[y2D[i,]>0,1:2],ncol=2,byrow=FALSE)
  if(nrow(trps)>1){
    s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
  }else{
    s.init[i,] <- trps
  }
}

Niminits <- list(z=z.init,N=N.init, #must initialize N to be the sum of z init
                 lambda.N=N.init, #initializing lambda.N to be consistent with N.init
                 s=s.init,
                 p0=runif(1,0.1,0.9),sigma=runif(1,0.25,1))

#constants for Nimble
constants <- list(M=M,J=J,K1D=K1D,xlim=data$xlim,ylim=data$ylim)

#supply data to nimble
Nimdata <- list(y=y2D,X=data$X)

#set parameters to monitor
parameters <- c('lambda.N','p0','sigma','N','z.up') #z.up is to record the number of z's we propose to update at once in zSampler3 only
nt <- 1 #thinning rate

#Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
config.nodes <- c('lambda.N','p0','sigma') #we will add samplers for N/z and s below
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,nodes=config.nodes) 

zSampler.use <- 2 #which zSampler do we use?

if(zSampler.use==1){
  #N/z sampler Approach 1: Update 1 z at a time, but multiple updates per iteration.
  #add update, select one of M-N z=0 inds, 
  #subtract update, select from N z=1 inds, autoreject if select a captured individual.
  #This is what I have used historically, proposal probs cancel with z prior, 1/choose(M,N)
  z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal, setting to 25% of M here.
  #nodes used for update
  y.nodes <- Rmodel$expandNodeNames(paste("y[1:",M,",1:",J,"]"))
  pd.nodes <- Rmodel$expandNodeNames(paste("pd[1:",M,",1:",J,"]"))
  N.node <- Rmodel$expandNodeNames(paste("N"))
  z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
  calcNodes <- c(N.node,pd.nodes,y.nodes)
  ind.detected <- 1*(rowSums(y2D>0)>0) #use these to know when to reject setting z=0 if detected
  conf$addSampler(target = c("N"),
                  type = 'zSampler',control = list(z.ups=z.ups,M=M,J=J,ind.detected=ind.detected,
                                                   y.nodes=y.nodes,pd.nodes=pd.nodes,
                                                   N.node=N.node,z.nodes=z.nodes,
                                                   calcNodes=calcNodes),silent = TRUE)
}else if(zSampler.use==2){
  #N/z sampler Approach 2: Update 1 z at a time, but multiple updates per iteration.
  #add update, select one of M-N z=0 inds,
  #subtract update, select from N-n.det z=1 inds, so detected individuals are not selected
  #since they are always rejected. Asymmetric proposal probs need to be accounted for here.
  #Approach 2 produces a greater effective sample size per unit time than Approach 1
  z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal, setting to 25% of M here.
  #nodes used for update
  y.nodes <- Rmodel$expandNodeNames(paste("y[1:",M,",1:",J,"]"))
  pd.nodes <- Rmodel$expandNodeNames(paste("pd[1:",M,",1:",J,"]"))
  N.node <- Rmodel$expandNodeNames(paste("N"))
  z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
  calcNodes <- c(N.node,pd.nodes,y.nodes)
  ind.detected <- 1*(rowSums(y2D>0)>0) #use these to know when to reject setting z=0 if detected
  conf$addSampler(target = c("N"),
                  type = 'zSampler2',control = list(z.ups=z.ups,M=M,J=J,ind.detected=ind.detected,
                                                    y.nodes=y.nodes,pd.nodes=pd.nodes,
                                                    N.node=N.node,z.nodes=z.nodes,
                                                    calcNodes=calcNodes),silent = TRUE)
  
}else if(zSampler.use==3){
  #N/z sampler Approach 3: Same as Approach 2, but we update multiple z's at a time and this number is tuned.
  #We can also do multiple multi z updates per iteration
  #add update, select one of M-N z=0 inds, 
  #subtract update, select from N-n.det z=1 inds, so detected individuals are not selected
  #since they are always rejected. Asymmetric proposal probs need to be accounted for here.
  #z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal, setting to 25% of M here.
  #reversible jump approaches typically only update N once per iteration. This is not optimal,
  #particularly with large individual heterogeneity in detection probability as is typical in SCR.
  z.ups <- 1
  #nodes used for update
  y.nodes <- Rmodel$expandNodeNames(paste("y[1:",M,",1:",J,"]"))
  pd.nodes <- Rmodel$expandNodeNames(paste("pd[1:",M,",1:",J,"]"))
  N.node <- Rmodel$expandNodeNames(paste("N"))
  z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
  calcNodes <- c(N.node,pd.nodes,y.nodes)
  ind.detected <- 1*(rowSums(y2D>0)>0) #use these to know when to reject setting z=0 if detected
  conf$addSampler(target = c("N"),
                  type = 'zSampler3',control = list(z.ups=z.ups,M=M,J=J,ind.detected=ind.detected,
                                                    z.up.start=3, #How many z's do we propose to update at once? Starting value
                                                    z.up.max=10, #upper limit on the number updated at once
                                                    tune.interval=100, #tuning interval
                                                    tune.stop.iter=5000, #iteration where we stop tuning
                                                    y.nodes=y.nodes,pd.nodes=pd.nodes,
                                                    N.node=N.node,z.nodes=z.nodes,
                                                    calcNodes=calcNodes),silent = TRUE)
}else{print("zSampler.use must be 1, 2, or 3")}

for(i in 1:M){
  s.nodes <- Rmodel$expandNodeNames(paste("s[",i,", 1:2]", sep="")) #used when z=0 when proposing from prior
  calcNodes <- Rmodel$getDependencies(paste("s[",i,", 1:2]", sep="")) #used when z=1
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=data$xlim,ylim=data$ylim,
                                                 s.nodes=s.nodes,calcNodes=calcNodes,scale=1),silent = TRUE)
}

#can add block update for  correlated posteriors
conf$addSampler(target = c("p0","sigma"),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)


# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(5000,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[250:nrow(mvSamples),])) 
#z.up records number of z's updated at once. 
#Only relevant for zSampler3. Stops updating after tune.stop.iter iterations

#can run one chain with each N/z update and check Rhats to see they produce the same posterior distributions
#removing z.up field
mvSamples1 <- mvSamples[,-which(colnames(mvSamples)=="z.up")] #z1
mvSamples2 <- mvSamples[,-which(colnames(mvSamples)=="z.up")] #z2
mvSamples3 <- mvSamples[,-which(colnames(mvSamples)=="z.up")] #z3

# cor(mcmc(mvSamples[250:nrow(mvSamples),]))

tmp <- mcmc.list(mcmc(mvSamples1[250:nrow(mvSamples1),]),
                mcmc(mvSamples2[250:nrow(mvSamples2),]),
                mcmc(mvSamples3[250:nrow(mvSamples3),]))
gelman.diag(tmp)
plot(tmp)

#effective sample sizes for each sampler. Need to consider run time for ESS/unit time
effectiveSize(mcmc(mvSamples1[250:nrow(mvSamples1),]))
effectiveSize(mcmc(mvSamples2[250:nrow(mvSamples2),]))
effectiveSize(mcmc(mvSamples3[250:nrow(mvSamples3),]))

