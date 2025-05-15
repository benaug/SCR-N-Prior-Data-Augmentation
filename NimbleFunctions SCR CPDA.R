GetDetectionProb <- nimbleFunction(
  run = function(s=double(1), p0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- (s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2
      ans <- p0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

dBernoulliVector <- nimbleFunction(
  run = function(x=double(1), pd=double(1), K1D=double(1), z=double(0),
                 log=integer(0)){
    returnType(double(0))
    if(z==0){ #never propose to set z=0 for individuals with detections.
      return(0) #otherwise, need to return -Inf if these are proposed
    }else{
      logProb <- sum(dbinom(x, prob = pd, size = K1D, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rBernoulliVector <- nimbleFunction(
  run = function(n = integer(0),pd = double(1),K1D=double(1), z = double(0)){
    returnType(double(1))
    J <- nimDim(pd)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

#Required custom update for N/z
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    M <- control$M
    z.ups <- control$z.ups
    y.nodes <- control$y.nodes
    pd.nodes <- control$pd.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
    ind.detected <- control$ind.detected
  },
  run = function(){
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        # find all z's currently on
        z.on <- which(model$z==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        #prereject turning off detected individuals
        if(ind.detected[pick]==1){#is this an individual with captured?
          reject <- TRUE
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])
          # N.initial <- model$N[1]
          
          #propose new N/z
          model$N[1] <<-  model$N[1] - 1
          model$z[pick] <<- 0
          
          #turn pd off
          model$calculate(pd.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          # log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y) +
          #   log(N.initial) - log(M - N.initial + 1)
          accept <- decide(log_MH_ratio)
          
          if(accept){
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["pd",1][pick,] <<- model[["pd"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["pd"]][pick,] <<- mvSaved["pd",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0
          # N.initial <- model$N[1]
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1
          
          #turn pd on
          model$calculate(pd.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          # log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y) +
          #   log(M - N.initial) - log(N.initial + 1)
          accept <- decide(log_MH_ratio)
          if(accept){
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["pd",1][pick,] <<- model[["pd"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["pd"]][pick,] <<- mvSaved["pd",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function (){} )
)

zSampler2 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    M <- control$M
    z.ups <- control$z.ups
    y.nodes <- control$y.nodes
    pd.nodes <- control$pd.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
    ind.detected <- control$ind.detected
    n.det <- sum(control$ind.detected)
  },
  run = function(){
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){ #subtract
        # find all z's currently on
        # z.on <- which(model$z==1)
        z.on <- which(model$z == 1 & ind.detected==0)
        n.z.on <- length(z.on)
        if(n.z.on==0){
          reject <- TRUE
        }
        if(!reject){
          pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
          pick <- z.on[pick]
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])
          N.initial <- model$N[1]
          
          #propose new N/z
          model$N[1] <<-  model$N[1] - 1
          model$z[pick] <<- 0
          
          #turn pd off
          model$calculate(pd.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0
          
          #MH step
          # log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y) +
          #   log(N.initial - n.det) - log(M - N.initial)
          # log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y) +
          #   log(N.initial - n.det) - log(M - N.initial + 1)
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y) +
            log(N.initial - n.det) - log(N.initial)
          accept <- decide(log_MH_ratio)
          
          if(accept){
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["pd",1][pick,] <<- model[["pd"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["pd"]][pick,] <<- mvSaved["pd",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0
          N.initial <- model$N[1]
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1
          
          #turn pd on
          model$calculate(pd.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])
          
          #MH step
          # log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y) +
          #   log(M - N.initial) - log(N.initial + 1 - n.det)
          # log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y) +
          #   log(M - N.initial) - log(N.initial + 1)
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y) +
            log(N.initial + 1) - log(N.initial + 1 - n.det)
          accept <- decide(log_MH_ratio)
          if(accept){
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["pd",1][pick,] <<- model[["pd"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["pd"]][pick,] <<- mvSaved["pd",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function (){} )
)


zSampler3 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    M <- control$M
    z.ups <- control$z.ups
    y.nodes <- control$y.nodes
    pd.nodes <- control$pd.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
    ind.detected <- control$ind.detected
    n.det <- sum(control$ind.detected)
    z.up.start <- control$z.up.start  # Initial target number of z's to update
    z.up <- z.up.start
    z.up.max <- control$z.up.max  # Maximum z.up
    tune.interval <- control$tune.interval  # Tuning frequency
    tune.stop.iter <- control$tune.stop.iter  # Iteration to stop tuning
    # Variables for tuning
    n.proposals <- 0
    n.accept <- 0
    iter <- 0
  },
  run = function(){
    for(up in 1:z.ups){
      iter <<- iter + 1
      updown <- rbinom(1, 1, 0.5)
      reject <- FALSE
      if(updown == 0){ # subtract
        z.on <- which(model$z == 1 & ind.detected == 0)
        n.z.on <- length(z.on)
        if(n.z.on == 0){  # No non-detected individuals
          reject <- TRUE
        }
        if(!reject){
          # Set k to min(n.z.on, z.up)
          k <- min(n.z.on, z.up)
          # Select k indices without replacement
          picks <- numeric(k)
          indices <- z.on  # Copy of available indices
          n <- n.z.on
          for(j in 1:k){
            idx <- floor(runif(1, 1, n + 1))  # Random index from 1 to n
            picks[j] <- indices[idx]
            indices[idx] <- indices[n]
            n <- n - 1
          }
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- 0
          for(j in 1:k){
            lp.initial.y <- lp.initial.y + model$getLogProb(y.nodes[picks[j]])
          }
          N.initial <- model$N[1]
          # Propose new N/z
          model$N[1] <<- model$N[1] - k
          model$z[picks] <<- 0
          model$calculate(pd.nodes[picks])
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- 0
          for(j in 1:k){
            lp.proposed.y <- lp.proposed.y + model$calculate(y.nodes[picks[j]])
          }
          # MH ratio
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          #add z prior and asymmetric proposal prob terms that don't cancel
          for(j in 0:(k-1)){
            log_MH_ratio <- log_MH_ratio + log(N.initial - n.det - j) - log(N.initial - j)
          }
          n.proposals <<- n.proposals + 1
          accept <- decide(log_MH_ratio)
          if(accept){
            n.accept <<- n.accept + 1
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["pd",1][picks,] <<- model[["pd"]][picks,]
            mvSaved["z",1][picks] <<- model[["z"]][picks]
          } else {
            model[["N"]] <<- mvSaved["N",1][1]
            model[["pd"]][picks,] <<- mvSaved["pd",1][picks,]
            model[["z"]][picks] <<- mvSaved["z",1][picks]
            for(j in 1:k){
              model$calculate(y.nodes[picks[j]])
            }
            model$calculate(N.node)
          }
        }
      }else{ # add
        z.off <- which(model$z == 0)
        n.z.off <- length(z.off)
        if(n.z.off == 0){  # No available individuals or N at max
          reject <- TRUE
        }
        if(!reject){
          # Set k to min(n.z.off, z.up, M - N.initial)
          k <- min(n.z.off, z.up)
          if(k == 0){
            reject <- TRUE
          }
          if(!reject){
            # Select k indices without replacement
            picks <- numeric(k)
            indices <- z.off  # Copy of available indices
            n <- n.z.off
            for(j in 1:k){
              idx <- floor(runif(1, 1, n + 1))  # Random index from 1 to n
              picks[j] <- indices[idx]
              indices[idx] <- indices[n]
              n <- n - 1
            }
            lp.initial.N <- model$getLogProb(N.node)
            lp.initial.y <- 0
            for(j in 1:k){
              lp.initial.y <- lp.initial.y + model$getLogProb(y.nodes[picks[j]])
            }
            N.initial <- model$N[1]
            model$N[1] <<- model$N[1] + k
            model$z[picks] <<- 1
            model$calculate(pd.nodes[picks])
            lp.proposed.N <- model$calculate(N.node)
            lp.proposed.y <- 0
            for(j in 1:k){
              lp.proposed.y <- lp.proposed.y + model$calculate(y.nodes[picks[j]])
            }
            log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
            for(j in 1:k){
              log_MH_ratio <- log_MH_ratio + log(N.initial + j) - log(N.initial + j - n.det)
            }
            n.proposals <<- n.proposals + 1
            accept <- decide(log_MH_ratio)
            if(accept){
              n.accept <<- n.accept + 1
              mvSaved["N",1][1] <<- model[["N"]]
              mvSaved["pd",1][picks,] <<- model[["pd"]][picks,]
              mvSaved["z",1][picks] <<- model[["z"]][picks]
            }else{
              model[["N"]] <<- mvSaved["N",1][1]
              model[["pd"]][picks,] <<- mvSaved["pd",1][picks,]
              model[["z"]][picks] <<- mvSaved["z",1][picks]
              for(j in 1:k){
                model$calculate(y.nodes[picks[j]])
              }
              model$calculate(N.node)
            }
          }
        }
      }
      # Tune z.up every tune.interval iterations
      if((iter %% tune.interval == 0) & (iter <= tune.stop.iter)){
        if(n.proposals > 0){
          accept.rate <- n.accept / n.proposals
          if(accept.rate > 0.4 & z.up < z.up.max){
            z.up <<- z.up + 1
          } else if(accept.rate < 0.2 & z.up > 1){
            z.up <<- z.up - 1
          }
        }
        # Reset counters
        n.proposals <<- 0
        n.accept <<- 0
      }
    }
    model$z.up[1] <<- z.up
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function(){
    n.proposals <<- 0
    n.accept <<- 0
    iter <<- 0
    z.up <<- z.up.start  # Reset to initial z.up
  })
)