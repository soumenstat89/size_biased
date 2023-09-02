#######===================# DATA SIMULATION & NIMBLE MODEL FITTING SCRIPT ----------------########

{ ##--DO ALL
  rm(list=ls())
  cat("\014")
  gc()
  options(scipen = 999) #  no. of printed digits = 999+5 = 1004
  ts = format(Sys.time(), "%y%m%d_%H%M%S")
 ##------ IMPORT REQUIRED LIBRARIES 
  library(extraDistr)
  library(coda)               # Import the MCMC diagnostic tools
  library(nimble)             # Import the NIMBLE subroutines
  
  myVars <- list(WD = getwd(),
                 FunD = getwd(),
                 modelName = "ISROData",
                 shape.S = 0.5, 
                 rate.S = 0.01, 
                 M = 200,
                 J = 30,
                 epsVec =  c(25,50,75, 100, 150, 200),
                 ninput_extra_set  = c(25, 50, 75, 100, 150, 200, 250, 300),
                 nBatch = 500,
                 numIters = 10000, 
                 numBurnIn = 5000,
                 numChains = 3,
                 numThin = 1)
  
  source("dbin_geom_ISRO.R")
  
  source("simulateFromPosterior.R")

  getU <- nimbleFunction(
    run = function(x = double(1),
                   n = double(0)) {
    u <- nimNumeric(length = n, value = 0, init = TRUE)
    if(x[1] != -1){
      this.phase <- x[2]
      u[this.phase:n] <- 1
    }
    returnType(double(1))
    return(u)
  })


}##-- DO ALL

#######===================# GETTING INPUT DATA FOR MODEL SETUP ----------------########
{##-- DO ALL
  if(is.null(myVars$modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
  if(is.null(myVars$WD))stop("YOU SHOULD PROBABLY CHOOSE A WORKING DIRECTORY FOR THIS ANALYSIS/MODEL")
  if(!dir.exists(file.path(myVars$WD))){dir.create(file.path(myVars$WD))}
  if(!dir.exists(file.path(myVars$WD,myVars$modelName))){dir.create(file.path(myVars$WD, myVars$modelName))}
  
  outfolder <- file.path(myVars$WD,myVars$modelName,paste0("ISRO_OutFiles_Groupism_shape.S",
                                                           myVars$shape.S, "_", ts))
  if(!dir.exists(outfolder)){dir.create(outfolder,recursive=TRUE)}
  shape.S = myVars$shape.S 
  rate.S = myVars$rate.S 
  M <- myVars$M 
  
 
  bugd0_MT <- read.csv("Stages_A_B_03Oct2021_DETECTIONS_clean.csv",header = T, sep = ",")
  bugd0 <- read.csv("Stage_C_03Oct2021_DETECTIONS_clean.csv",header = T, sep = ",")
  inputs0 <- read.csv("Stage_C_03Oct2021_INPUTS_clean.csv", header = T, sep = ",")
  
  ##--- BUG ID ---##
  bugids <- c(paste0("MT_", bugd0_MT[, "BUG_ID"]), paste0("ST_", bugd0[, "BUG_ID"])) # n.detected x 1
  bugids
  
  ##---  NUMBER OF UNIQUE BUG DETECTED ---##
  n.detected <- length(bugids)
  n.bugs.detected <- sum(bugd0$COUNT) + sum(bugd0_MT$COUNT)
  
  
  ##---  NUMBER OF UNIQUE BUG DETECTED ---##
  n.detected.Software3.CI <- 15 # NUMBER OF UNIQUE BUG DETECTED DURING CODE INSPECTION OF SOFTWARE 3
  n.detected.Software4.CI <- 18 # NUMBER OF UNIQUE BUG DETECTED DURING CODE INSPECTION OF SOFTWARE 4
  n.detected.MT <- 27 # NUMBER OF UNIQUE BUG DETECTED DURING MODULE TESTING
  n.detected.ST <- 34 # NUMBER OF UNIQUE BUG DETECTED DURING SIMULATION TESTING
  n.detected.total <- n.detected.Software3.CI + n.detected.Software4.CI + n.detected.MT + n.detected.ST
  
  ##--- NUMBER OF INPUTS IN EACH PHASES ---##
  
  inputs0[nrow(inputs0),2] <- 1489 
  inputs_ST <- inputs0[, - c(1)]
  inputs_MT <- cbind(bugd0_MT[, "INPUTS"], matrix(1, nrow(bugd0_MT), ncol(inputs_ST)-1))
  colnames(inputs_MT) <- colnames(inputs_ST)
  inputs <- rbind(inputs_MT, inputs_ST[1:(nrow(inputs_ST)-1),])                
 
  inputs.aug <- matrix(as.numeric(inputs_ST[nrow(inputs_ST), ]),
                       nrow = M - n.detected, ncol = ncol(inputs),
                       byrow = T)
  colnames(inputs.aug) <- colnames(inputs)
  inputs.aug <- rbind(inputs, inputs.aug)

  
  ##--- NUMBER OF PHASES (including one phase corresp to module testing and 7 phases of simulations testing) ---## 
  nphases <- ncol(inputs)  # first col of inputs0 had missiona names
 
  
  ##--- CAPTURE-RECAPTURE DATA SET IN SPARSE FORMAT ---##
  bugd <- cbind(c(bugd0_MT$COUNT, bugd0$COUNT), c(bugd0_MT$PHASE +1, bugd0$PHASE+1))
  colnames(bugd) <- c("Detection", "Phase")
  sparseData <- rbind(bugd, matrix(-1, M - n.detected, ncol(bugd)))
  
  z <- c(rep(1, n.detected), rep(NA, M - n.detected))
  u <- do.call(rbind, lapply(1:M, function(i){getU(sparseData[i,], nphases)}))
  
  
  ### ====    3.8.SET THE INPUT FOR NIMBLE ====  
  
  
  ### ====      3.8.1. Define the constants to be used in the model ====
  
  nimConstants <- list(
    nphases = nphases,
    n.detected = n.detected,
    n.detected.Software3.CI = n.detected.Software3.CI,
    n.detected.Software4.CI = n.detected.Software4.CI,
    M = M, 
    shape.S = shape.S,
    rate.S = rate.S
  )
  ### ====      3.8.2. Define the data to be used in the model ====
  
  nimData <- list(
    inputs = inputs.aug, # inputs, #
    sparseData = sparseData,
    z = z
  )

  ### ====      3.8.3. Define the initial values to be used in the model  ====
  z.init <- ifelse(!is.na(z), NA, 1)
  z.init[!is.na(z.init)] <- rbinom(sum(!is.na(z.init)), size = 1, prob = 0.1)
  lambda.S.init <- rep(NA, M)
  Svec.init <- rep(NA, M)
  
  i <- 1
  for(i in 1:n.detected){
    lambda.S.init[i] <- 50 #sparseData[i,2] +5 # y[i,]
    Svec.init[i] <- rtpois(1, lambda.S.init[i], a = 10, b = 100)
   }
  lambda.S.init[is.na(Svec.init)] <- rgamma(M-n.detected, shape = 10, rate = 2)
  Svec.init[is.na(Svec.init)] <- rtpois(M-n.detected, lambda.S.init[is.na(Svec.init)], a = 5, b = 50)
  
  
  nimInits <- list(
    z = z.init,           
    psi = runif(1,0,1),
    logitr = logit(runif(1, 10^(-6), 10^(-5))),
   lambda.S = lambda.S.init,
    Svec = Svec.init
    
  )

  
  ### ====      3.8.4. Define the NIMBLE model code ====
  
  modelCode <- nimbleCode({
    
    ##---- LATENT PROCESS
    psi ~ dunif(0, 1)
    logitr ~ dunif( -20, 20)
    r <- ilogit(logitr)
    for (i in 1:M){
      
      z[i] ~ dbern(psi)
      lambda.S[i] ~ dgamma(shape = shape.S, rate = rate.S)
      Svec[i] ~ dpois(lambda = lambda.S[i])
      pvec[i] <-  1 - ((1-r)^(Svec[i]))
      sparseData[i, 1:2] ~ dbin_geom2(inputs = inputs[i, 1:nphases],
                                                       p = pvec[i],
                                                       indicator = z[i])
      
    }#i
    
    ##---- DERIVED QUANTITIES
    nGroup <- sum(z[1:M]) + n.detected.Software3.CI + n.detected.Software4.CI
    
  })#modelCode
  
  
  ### ====      3.8.5. Define the parameters to be monitored ====
  nimParams <- c("nGroup","psi", "z", "Svec", "lambda.S", "r", "logitr")
  
  inname <- file.path(outfolder, paste("InputForISROData_", ts, ".RData", sep = ""))
  
  save(myVars, sparseData, u, z, 
       bugd, bugd0, bugids, n.detected, inputs, M, n.detected, 
       nimConstants, nimData, nimInits, nimParams, modelCode,
       n.detected.Software3.CI,  n.detected.Software4.CI, n.detected.MT, 
       n.detected.ST, n.detected.total,
       file = inname)
}##--DO ALL

#######===================# MCMC # MCMC # MCMC # MCMC # MCMC  #===================########
{##--DO ALL
  
    model <- nimbleModel( code = modelCode,
                          constants = nimConstants,
                          data = nimData,
                          inits = nimInits,
                          check = F,
                          calculate = F)
    # CALCULATE LIKELIHOOD OF UNCOMPILED MODEL
    print(model$calculate())
    
    cmodel <- compileNimble(model)
    ## this should return a value. If NA: likely missing an inits value
    #                              if -inf something wrong with the inits valus and/or the model
    print(cmodel$calculate())
    
    # }##--DO ALL
    MCMCconf <- configureMCMC(model = model,
                              monitors = nimParams,
                              control = list(reflective = TRUE),
                              thin = 1)
    
    MCMC <- buildMCMC(MCMCconf)
    cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
    
    ### ====    6.3.RUN THE MCMC 
    {##-- DO ALL
     
      burnin = 0
      niter =  10000 
      nchains =  myVars$numChains 
      
      Runtime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC,
                                                        nburnin = burnin,
                                                        niter = niter,
                                                        nchains = nchains, 
                                                        samplesAsCodaMCMC = TRUE))
      
    }##-- DO ALL
  
  outname <- paste("OutputForISROData_nGroup.r_niter", niter, "_", ts, ".RData", sep = "")
  fname = file.path(outfolder, outname)
   
  save(myNimbleOutput, Runtime,
       niter, burnin, nchains,
       file=fname)
  
}##--DO ALL


#######===================# PREDICT N # PREDICT N # PREDICT N  #==================#################
   {##--DO ALL
    burnin <- myVars$numBurnIn 
    n.iterations <- nrow(myNimbleOutput[[1]])
    myNimbleOutput.subset <-  as.mcmc.list(lapply(myNimbleOutput,
                                                  function(x){mcmc(x[(burnin+1):n.iterations,])}))
    
    if(length(myNimbleOutput.subset) > 1){posteriorSamples <- do.call(rbind, myNimbleOutput.subset)} 
    if(length(myNimbleOutput.subset) == 1){posteriorSamples <- as.matrix(myNimbleOutput.subset)} 
    print(dim(posteriorSamples)[1])
    
    
    nBatch <- 100 
    seq.mat <- matrix(1:nrow(posteriorSamples), nrow = nBatch, byrow = T)
    
    ## SET-UP THE FUNCTION 
    simFromPost <- simulateFromPosterior(model = model, nodes = dimnames(posteriorSamples)[[2]]) # Just using the model and the tracked nodes
    
    ## COMPILE THE FUNCITON (building c++ equivalent version of the code)
    csimFromPost <- compileNimble(simFromPost, project = cmodel, resetFunctions = TRUE)
    batch = 1
    N.chains <- c()
    Nextra.chains <- c()
    for(batch in 1:nBatch){
      
      if (batch%%50 == 0)  cat('..... batch ', batch, '\n', sep = "")
      posteriorSamples.temp <- posteriorSamples[seq.mat[batch,],]
      
      csimFromPost$run(posteriorSamples = posteriorSamples.temp)
      
      sparseDataRep <- csimFromPost$mv[["sparseData"]]
      
      znames <- colnames(posteriorSamples)[grep("z", colnames(posteriorSamples))]
      
      iter <- 2
      N.chains.batch <- unlist(lapply(1:dim(posteriorSamples.temp)[1], function(iter){
        sparseDataRep.iter <- sparseDataRep[[iter]][,1]
        z <-   posteriorSamples.temp[iter,znames]
        
        undetectedGroups <- c(rep(F, n.detected), rep(T, length(z) - n.detected))
        bug.count <- n.bugs.detected + n.detected.Software3.CI +
          n.detected.Software4.CI + sum(sparseDataRep.iter[sparseDataRep.iter != -1 & undetectedGroups & z == 1])
        
        return(bug.count)
      }))
      N.chains <- c(N.chains, N.chains.batch)
      
      Nextra.chains.batch <- unlist(lapply(1:dim(posteriorSamples.temp)[1], function(iter){
        sparseDataRep.iter <- sparseDataRep[[iter]][,1]
        z <-   posteriorSamples.temp[iter,znames]
        undetectedGroups <- c(rep(F, n.detected), rep(T, length(z) - n.detected))
        bug.count <- sum(sparseDataRep.iter[sparseDataRep.iter != -1 & undetectedGroups & z == 1])
        return(bug.count)
      }))
      Nextra.chains <- c(Nextra.chains, Nextra.chains.batch)
      
    }#batch
    
    N.chains.mat <- matrix(N.chains, nrow = nrow(myNimbleOutput.subset[[1]]), ncol = nchains, byrow = F)
    N.chains.list <- as.mcmc.list(lapply(1:nchains, function(x){mcmc(N.chains.mat[,x])}))
   
    Nextra.chains.mat <- matrix(Nextra.chains, nrow = nrow(myNimbleOutput.subset[[1]]), ncol = nchains, byrow = F)
    Nextra.chains.list <- as.mcmc.list(lapply(1:nchains, function(x){mcmc(Nextra.chains.mat[,x])}))
    
    focal.parms <- c("nGroup","psi", "r", "logitr")
    
    myNimbleOutput.subset.new <- as.mcmc.list(lapply(1:nchains,function(x){
      out <- cbind(myNimbleOutput.subset[[x]][,focal.parms], N.chains.list[[x]], Nextra.chains.list[[x]])
      colnames(out) <- c(focal.parms, "predictedN", "predictedNextra")
      return(mcmc(out))
    }))
    
     
    ts <- strsplit(paste(strsplit(outname, "_")[[1]][4:5], collapse="_"), ".RData")[[1]][1]
    fname <- file.path(outfolder,
                       paste("OutputForISROData_Nbugs_niter", niter, "_", ts, ".RData", sep = ""))
    
    
    save(N.chains.list, Nextra.chains.list, myNimbleOutput.subset.new, 
         n.bugs.detected, n.detected.total,
         burnin, n.iterations, nchains,
         file = fname)
    
    
  }##--DO ALL

#######===================# RELIABILITY # RELIABILITY # RELIABILITY #===========############
{##--DO ALL
 
  ninput_extra_set <- myVars$ninput_extra_set
  epsVec <-  myVars$epsVec 
  rel.arr <- array(0, dim = c(length(epsVec), 2, length(ninput_extra_set)))
  dimnames(rel.arr) <- list(paste0("eps", epsVec),
                            paste0("phase", c(nimConstants$nphases, nimConstants$nphases + 1)),
                            paste0("extraTestCase", ninput_extra_set))
  cc <- 0
  ninput_extra <- ninput_extra_set[1]
  
  for(ninput_extra in ninput_extra_set){
    cat('------ extra input ', ninput_extra, '\n', sep = "")
        
        load( file.path(outfolder, paste("InputForISROData_", ts, ".RData", sep = "")))
        load( file.path(outfolder,list.files(outfolder, pattern = "OutputForISROData")[1]))
        burnin <- myVars$numBurnIn 
        n.iterations <- nrow(myNimbleOutput[[1]])
        myNimbleOutput.subset <-  as.mcmc.list(lapply(myNimbleOutput,
                                                      function(x){mcmc(x[(burnin+1):n.iterations,])}))
        
        if(length(myNimbleOutput.subset) > 1){posteriorSamples <- do.call(rbind, myNimbleOutput.subset)} 
        if(length(myNimbleOutput.subset) == 1){posteriorSamples <- as.matrix(myNimbleOutput.subset)} 
        
        M <- nimConstants$M
        
        J <- nimConstants$nphases + 1 
        phaseVec <- nimConstants$nphases:J
        
        nphases.initial <- nimConstants$nphases
        nimConstants$nphases <- J # The total number of possible phases (only first "nphases" were performed)
       
        inputs.initial <- nimData$inputs 
        inputs.aug <- matrix(ninput_extra, nrow = M, ncol = J-nphases.initial)
        nimData$inputs <- cbind(inputs.initial, inputs.aug) 
        nimData <- nimData[-which(names(nimData) %in% c("sparseData"))]
        model <- nimbleModel( code = modelCode,
                              constants = nimConstants,
                              data = nimData,
                              inits = nimInits,
                              check = FALSE,       
                              calculate = FALSE)  
        ## CALCULATE LIKELIHOOD OF UNCOMPILED MODEL
        # model$calculate() 
        cmodel <- compileNimble(model)
        ## this should return a value. If NA: likely missing an inits value
        #                              if -inf something wrong with the inits valus and/or the model
        # cmodel$calculate() 
        
        ## SET-UP THE FUNCTION 
        simFromPost <- simulateFromPosterior(model = model, nodes = dimnames(posteriorSamples)[[2]]) ##--Just using the model and the tracked nodes
        
        ## COMPILE THE FUNCITON (building c++ equivalent version of the code)
        csimFromPost <- compileNimble(simFromPost, project = cmodel, resetFunctions = TRUE)
        
        ## RUN THE COMPILED VERSION (much faster) (Note that, here we are using the already obtained MCMC samples of the tracked model parameters)
        
        nBatch <- myVars$nBatch 
        seq.mat <- matrix(1:dim(posteriorSamples)[1], nrow = nBatch, byrow = T)
        
        Bvec.chains <- NULL
        kvec.chains <- NULL
        
        batch <-  1
        
        for(batch in 1:nBatch){
          
          if (batch%%50 == 0)  cat('..... batch ', batch, '\n', sep = "")
          posteriorSamples.temp <- posteriorSamples[seq.mat[batch,],]
          
          csimFromPost$run(posteriorSamples = posteriorSamples.temp) # Using the MCMC samples for all the chains
          
          sparseDataRep <- csimFromPost$mv[["sparseData"]] 
          
          znames <- colnames(posteriorSamples)[grep("z", colnames(posteriorSamples))]
          SvecNames <- colnames(posteriorSamples)[grep("Svec", colnames(posteriorSamples))]
          
          iter <- 2
          moreSamples <- lapply(1:dim(posteriorSamples.temp)[1], function(iter){
            sparseDataRep.iter <- sparseDataRep[[iter]]
            
            z <-   posteriorSamples.temp[iter,znames]
            Svec <- posteriorSamples.temp[iter,SvecNames]
            
            # uij = 1 if i-th bug is detected on or before j-th phase
            u <- do.call(rbind, lapply(1:M, function(i){getU(sparseDataRep.iter[i,], J)}))
            
            
            undetectedGroups <- c(rep(F, n.detected), rep(T, length(z) - n.detected))
            indices <- sparseDataRep.iter[,1] != -1 & undetectedGroups & z == 1
            if(length(indices) > 0) undetectedBug.freq <- sparseDataRep.iter[indices, 1]
            if(length(indices) == 0) undetectedBug.freq <- rep(0, length(z))
            
            # B_i = remaining eventual bug size
            Bvec <- unlist(lapply(phaseVec, function(j) {sum(Svec[1:M]*z[1:M]*(1-u[1:M,j])*undetectedBug.freq)}))
            
            
            l <- 1
            kvec <- unlist(lapply(1:length(epsVec), function(l){
              k <- phaseVec[which(Bvec<= epsVec[l])[1]]
              if(is.na(k)) k <- J+1 
              return(k)
            }))
            
            out <- list(Bvec = Bvec,
                        kvec = kvec)
            return(out)
            
          }) 
          
          Bvec.chains <- rbind(Bvec.chains, do.call(rbind, lapply(moreSamples, function(x) x$Bvec)))
          kvec.chains <- rbind(kvec.chains, do.call(rbind, lapply(moreSamples, function(x) x$kvec)))
          
        }#batch
        
        niter.per.chain <- nrow(myNimbleOutput.subset[[1]])
        knames <-paste0("k_eps", epsVec)
        
       myNimbleOutput.new <- as.mcmc.list(lapply(c(1:nchains), function(chain){
          
          this.chain.index <- (chain -1)*niter.per.chain + 1:niter.per.chain
          a <-  cbind(kvec.chains[this.chain.index,], Bvec.chains[this.chain.index,])
          dimnames(a)[[2]] <- c(knames, paste0("B_", phaseVec)) 
          a <- mcmc(a)
          return(a)
        }))  
        
        Bnames <- paste0("B_", phaseVec)
        colnames(Bvec.chains) <- paste0("B_", phaseVec)
        
        l <- 1
        
        rel.mat <- do.call(rbind, lapply(1:length(epsVec), function(l){
          
          out <- unlist(lapply(Bnames, function(this.B){
            rel <- mean(Bvec.chains[, this.B]<= epsVec[l])
            if(is.na(rel)) rel <- 0 
            return(rel)
          }))
          return(out)
        }))
        cc <- cc + 1
        rel.arr[,,cc] <- rel.mat
        
        
        outnamekB <- paste("OutputForBugDatakB_J", J,"_ninputExtra", ninput_extra,"_", ts, ".RData", sep = "")
        fname = file.path(outfolder, outnamekB)
        
        save(myNimbleOutput.new, 
             knames, Bnames,
             myCodaOutput.new,
             Runtime,
             n.iterations, burnin, nchains,
             J, epsVec, ninput_extra_set,
             file=fname)
  }#ninput_extra
  

  dimnames(rel.arr) <- list(paste0("eps",epsVec), Bnames, ninput_extra_set) 
  
  outnamekB <- paste("ReliabilityEstimate_J", J,"_", ts, ".RData", sep = "")
  fname = file.path(outfolder, outnamekB)
  
  save(rel.arr,
       Bnames,epsVec,ninput_extra_set,
       file=fname)
  
}##--DO ALL

