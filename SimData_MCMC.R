##############################################################################################
#=========================# SIMULATION, MCMC AND PROCESSING SCRIPT   #=======================#
##############################################################################################
{##--DO ALL
  # }  ##--DO ALL
  rm(list=ls())
  cat("\014")
  gc()
  options(scipen = 999) #  no. of printed digits = 999+5 = 1004
  ts = format(Sys.time(), "%y%m%d_%H%M%S")
  
  ## ------ SOURCE THE REQUIRED LIBRARIES ------
  library(extraDistr)
  library(coda)                      # Import the MCMC diagnostic tools
  library(nimble)                    # Import the NIMBLE subroutines
  

  getU <- nimbleFunction(run = function(x = double(1), n = double(0)) {
    u <- nimNumeric(length = n, value = 0, init = TRUE)
    if(x[1] != -1){
      this.phase <- x[2]
      u[this.phase:n] <- 1
    }
    returnType(double(1))
    return(u)
  })
  

  ## ==============================================================================================
  ## ===== 1.SET ANALYSIS CHARACTERISTICS =====
  ## ==============================================================================================
  ### ====    1.1.GENERAL VARIABLES DECLARATION ====
  shape.S = 0.5
  myVars <- list( 
    WD = getwd(),
    FunD = getwd(),
    modelName = paste("size.biased_simulatedData_", ts, sep = ""),
    TRIALS = list(nphases = 5, 
                  J = 30), 
    BUGSIZE = list(shape.S = 0.5, 
                   rate.S = 0.01), 
    POPULATION = list( N = c(200),  #---TRUE NUMBER OF BUGS
                       M = 400), #---size of the augmented no. of bugs
    DETECTIONS = list( 
      epsVec =  c(25,50,75, 100, 150, 200),
                       ninput = c(1000,2000), 
                       ninput_extra = 2000, 
                       r = c( 0.75*10^(-5), 1.5*10^(-5)), 
      nrep = 50
                       ),
    MCMC = list( numIters = 10000, 
                 numBurnIn = 5000,
                 numChains = 3,
                 numThin = 1,
                 nBatch = 500),
    ## MISCELLANEOUS
    plot.check = TRUE)
  
  ### ==== 1.2. CREATE FOLDERS ====
  if(is.null(myVars$modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
  if(is.null(myVars$WD))stop("YOU SHOULD PROBABLY CHOOSE A WORKING DIRECTORY FOR THIS ANALYSIS/MODEL")
  if(!dir.exists(file.path(myVars$WD))){dir.create(file.path(myVars$WD))}
  if(!dir.exists(file.path(myVars$WD,myVars$modelName))){dir.create(file.path(myVars$WD, myVars$modelName))}
  
  infolder <- file.path(myVars$WD,myVars$modelName,"NimbleInFiles")
  if(!dir.exists(infolder)){dir.create(infolder,recursive=TRUE)}
  infolder_proc <- file.path(myVars$WD,myVars$modelName,"NimbleInFiles_Proc")
  if(!dir.exists(infolder_proc)){dir.create(infolder_proc,recursive=TRUE)}
  outfolder <- file.path(myVars$WD,myVars$modelName,"NimbleOutFiles")
  if(!dir.exists(outfolder)){dir.create(outfolder,recursive=TRUE)}
  figfolder <- file.path(myVars$WD,myVars$modelName,"FIGURES")
  if(!dir.exists(figfolder)){dir.create(figfolder,recursive=TRUE)}
  
  
  nphases <- myVars$TRIALS$nphases
  J <- myVars$TRIALS$J
  shape.S <- myVars$BUGSIZE$shape.S
  rate.S <- myVars$BUGSIZE$rate.S
  N <- myVars$POPULATION$N
  M <- myVars$POPULATION$M
  
  epsVec <- myVars$DETECTIONS$epsVec
  
  ## ==============================================================================================
  ## ------ 2. DEFINE THE PARAMETER SPACE -----
  ## ==============================================================================================
  
  ### --- 2.1. DEFINE VARYING PARAMETERS ACROSS SIMULATIONS ---
   
  param.space <- MakeParameterSpace( list.param = list(ninput= myVars$DETECTIONS$ninput,
                                                       r = myVars$DETECTIONS$r
  ),
  n.rep=myVars$DETECTIONS$nrep)
  
  dim(param.space)
  

  ### --- 2.2. ADD CONSTANT PARAMETERS ACROSS SIMULATIONS ---
  ## ADD THE CONSTANT PARAMATERS (from myVars) TO THE PARAMETER LIST
  temp <- unlist(myVars)
  sim.names <- names(myVars)
  temp.names <- names(temp)
  temp <- as.numeric(temp)
  
  names(temp) <- temp.names
  temp <- na.omit(temp)
  temp.names <- names(temp)
  for(sn in sim.names)  temp.names <- gsub(paste(sn, ".", sep=""),"",temp.names)
  names(temp) <- temp.names
  temp.names <- temp.names[!temp.names%in%names(param.space)]
  
  
  param.space <- data.frame(param.space,t(temp[temp.names]))
  summary(param.space)
  
  param.space$sim.id <- paste("Set", param.space$set_ID, "Rep", param.space$rep_ID,".RData", sep="")
  
  path <- file.path(myVars$WD, myVars$modelName, "param.space1.RData")
  save(param.space, file = path)
  dim(param.space)
}#-- DO ALL
## ==============================================================================================
## ===== 3. SIMULATE INPUT DATA  ====
## ==============================================================================================
{#-- DO ALL
  
  i <- 1
  
  for(i in 1:dim(param.space)[1]){#i
    
    print(param.space$sim.id[i])
    
    inputs <- rep(param.space$ninput[i], nphases)
  
   lambda.S <- rgamma(N, shape = shape.S, rate = rate.S)
  Svec <- rpois(N, lambda.S)
 
  r <- param.space$r[i]
  
  pvec <- 1 - (1-r)^Svec

  y <- matrix(0, N, nphases)
  u <- matrix(0, N, nphases)
  
  for (ii in 1:N){ 
    for (jj in 1:nphases){
      y[ii,jj] <- rbinom(1, inputs[jj],  pvec[ii]) 
      if(sum(y[ii,1:jj]) > 0){ 
        u[ii, jj:nphases] <- 1 
        break
      }#if 
    }#j
  }#i
  
  
  detected <- apply(y,1, function(x) sum(x)>0)
  n.detected <- sum(detected)
  
  Bvec <- apply(u, 2, function(x){sum(Svec*(1-x))})
  # Bvec
  
  
  kvec <- unlist(lapply(1:length(epsVec), function(l){
    k <- which(Bvec<= epsVec[l])[1]
    if(is.na(k)) k <- J+1 #0
    return(k)
  }))
  kvec
  
  
  ### ====    3.5.SUBSET TO INDIVIDUALS DETECTED AT LEAST ONE YEAR/DATA SET ==== 
  
  y.orig <- y
  
  y <- y[detected, ]
  
  
  ### ====    3.6.AUGMENT DATA ====
  
  #---augmentation factor that ensures that the total number of individuals 
  # (alive + available) is always the same, regardless of the simulation
  this.AF <- M/n.detected-1
  
  #---check it: 
  n.detected * (1+ this.AF) ==  M
  
  y <- MakeAugmentation(y = y, aug.factor = this.AF, replace.value = 0)
  z <- MakeAugmentation( y = rep(1, n.detected),
                         aug.factor = this.AF, replace.value = NA)
  
  detected <- apply(y,1, function(x) sum(x)>0)
  sum(detected)
  
  # # NUMBER OF DETECTIONS FOR EACH ID
  
  detNums <- apply(y, 1, function(x) length(which(x>0))) ## GET NO. OF CYCLES/PHASES WHERE BUG i IS DETECTED
  nMaxDetectors = max(detNums)
  detIndices <- matrix(-1, M, nMaxDetectors)
  ySparse <-  matrix(-1, M, nMaxDetectors)
  for(ii in 1:M){
    if(detNums[ii]>0){
      # GET WHERE (DETECTOR ID) DETECTIONS OCCUR
      detIndices[ii, 1:detNums[ii]] <- which(y[ii,]>0)
      #     # GET NUMBE OF DETECTIONS 
      ySparse[ii, 1:detNums[ii]] <- y[ii, detIndices[ii, 1:detNums[ii]]]
    }#if
  }#i
  sparseData <- cbind(ySparse, detIndices)
  
  SparseY <- list(ySparse = ySparse,
                  detIndices = detIndices,
                  detNums = detNums,
                  nMaxDetectors = nMaxDetectors)
  
  ### ====    3.7.TRANSFORM Y TO SPARSE MATRICES =====
  
  ##############################################################################################
  ##############################################################################################
  ### ====    3.8.SET THE INPUT FOR NIMBLE   ====
  ### ====      3.8.1.Define the nimble model ====
  
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
      sparseData[i, 1:2] ~ dbin_geom(inputs = inputs[1:nphases],
                                                       p = pvec[i],
                                                       indicator = z[i])
      
     
    }#i
    
    ##---- DERIVED QUANTITIES
    N <- sum(z[1:M])
    
  })#modelCode
  
  ##############################################################################################
  
  ### ====      3.8.2.Define the constants to be used in the model ====
  # Set the list of model constants
  nimConstants <- list(
    nphases = nphases,
    M = dim(y)[1],
    shape.S = shape.S,
    rate.S = rate.S#,
  )
  ### ====      3.8.3.Define the data to be used in the model ====
  # Set the list of data
  nimData <- list(
    inputs = inputs,
    sparseData = sparseData,
    z = z
  )
  
  
  ### ====      3.8.4.Define the initial values to be used in the model  ====
  z.init <- ifelse(!is.na(z), NA, 1)
  z.init[!is.na(z.init)] <- rbinom(sum(!is.na(z.init)), size = 1, prob = 0.1)
  
  lambda.S.init <- rep(NA, M)
  Svec.init <- rep(NA, M)
  
  ii <- 1
  for(ii in 1:n.detected){
    lambda.S.init[ii] <- sparseData[ii,2] +5 # y[i,]
    Svec.init[ii] <- rtpois(1, lambda.S.init[ii], a = 1, b = 15)
    
  }
  lambda.S.init[is.na(Svec.init)] <- rgamma(M-n.detected, shape = 200, rate = 2)
  Svec.init[is.na(Svec.init)] <- rtpois(M-n.detected, lambda.S.init[is.na(Svec.init)], a = 5, b = 1000)
  
  nimInits <- list(z = z.init,           
                   psi = runif(1,0,1),
                   logitr = logit(runif(1, 10^(-6), 10^(-5))),
                   lambda.S = lambda.S.init,
                   Svec = Svec.init
                   
  )
  
  
  ### ====      3.8.5.Define the parameters to be monitored ====
  nimParams <- c("N","psi", "z", "Svec", "lambda.S", "r", "logitr")
  
  paramspace <- param.space[i, ]
  
  ### ====    0.2.CREATE FOLDERS ====
 file.name <- paste0( "NimbleInFile", 
                      "_Set", param.space$set_ID[i],
                      "Rep", param.space$rep_ID[i],
                      ".RData")
  save( modelCode,
        nimConstants,
        nimData,
        nimParams,
        nimInits,
        sparseData,
        SparseY,
        myVars,
        y.orig, 
        Svec,
        Bvec,
        kvec,
        paramspace,
        detected, 
        file = file.path(infolder, file.name))
  
  
  # print(file.path(myVars$WD,myVars$modelName))
  }# row in parameter space
  
}##--DO ALL


## ==============================================================================================
## ===== 4. ADD BASIC SIM SUMMARIES TO PARAMETER TABLE =====
## ==============================================================================================

{##--DO ALL
  
  load(file.path(myVars$WD, myVars$modelName, "param.space1.RData"))
  
  filelist <- list.files(infolder)
  
  out2 <- lapply(1:dim(param.space)[1], function(i){
    this.pattern <- param.space$sim.id[i]
    this.path <- file.path(infolder, filelist[grep(this.pattern,filelist)][1]) ##--need only one chain, data are identical
    
    ## PLACE HOLDER FOR MISSING SIMS
    out <- data.frame(
      n.super = NA,
      n.detected = NA,
      n.detections = NA#,
    )
    
    ## EXTRACT NUMBER OF DETECTORS, SUPER POPULATION, INDIVIDUALS DETECTED & DETECTIONS
    try({
      load(this.path)
      out <- data.frame(
        n.super = nimConstants$M,
        n.detected = sum(detected),
        n.detections = sum(y.orig, na.rm=TRUE)
     )
    }, silent = TRUE)
    return(out)
  })
  
  ## CONVERT BACK TO A DATAFRAME
  out2 <- do.call(rbind, out2)
  
  ## ADD TO THE EXISTING PARAMETER SPACE DATAFRAME
  param.space <- cbind(param.space, out2)
  head(param.space)
  
  ## SAVE THE COMPLETE PARAMETER SPACE
  save(param.space, file = file.path(myVars$WD, myVars$modelName, "param.space.RData"))
}##--DO ALL

param.space$n.detections
range(param.space$n.detections)

param.space$n.detected
range(param.space$n.detected)

## ==============================================================================================
## =========================# MCMC # MCMC # MCMC # MCMC # MCMC # MCMC #========================#
## ==============================================================================================

{##--DO ALL
  

  load(file.path(myVars$WD,myVars$modelName,"param.space.RData"))
  file.list <- list.files(infolder)
  
   
    while(length(file.list)>0  ){
      this.infile <- sample(file.list, 1)
      
    print(this.infile)
    
    this.pattern <- strsplit(this.infile, split = "_")[[1]][2]
    this.pattern
    ## ----- LOAD THE FILE -----
    load(file.path(infolder, this.infile))
    file.rename(from = file.path(infolder, this.infile), to = file.path(infolder_proc, this.infile))
    
  ##RUN THE MODEL 
  model <- nimbleModel( code = modelCode,
                        constants = nimConstants,
                        data = nimData,
                        inits = nimInits,
                        check = F,       
                        calculate = F)  
  #CALCULATE LIKELIHOOD OF UNCOMPILED MODEL
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
  
  ### ====    6.3.RUN THE MCMC ====
  {##--DO ALL
    
    burnin = 0
    niter =  myVars$MCMC$numIters
    nchains = myVars$MCMC$numChains 
    
    
    Runtime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC, 
                                                      nburnin = burnin, 
                                                      niter = niter, 
                                                      nchains = nchains, 
                                                      samplesAsCodaMCMC = TRUE))
    
    
  }##--DO ALL
  

  burnin <- myVars$MCMC$numBurnIn 
  n.iterations <- nrow(myNimbleOutput[[1]])
  myNimbleOutput.subset <-  as.mcmc.list(lapply(myNimbleOutput,
                                                function(x){mcmc(x[(burnin+1):n.iterations,])}))
  
  if(length(myNimbleOutput.subset) > 1){posteriorSamples <- do.call(rbind, myNimbleOutput.subset)} 
  if(length(myNimbleOutput.subset) == 1){posteriorSamples <- as.matrix(myNimbleOutput.subset)} 
  print(dim(posteriorSamples)[1])
  
  M <- nimConstants$M
  epsVec <-  myVars$DETECTIONS$epsVec
  
  this.index <- which(param.space$sim.id == this.pattern)
  ninput_extra <- param.space$ninput[this.index] #nimData$inputs[1]
  J <- myVars$TRIALS$J 
  
  nphases.initial <- nimConstants$nphases
  nimConstants$nphases <- J # The total number of possible phases (only first "nphases" were performed)
  
  detected <- apply(y.orig,1, function(x) sum(x)>0)
  n.detected <- sum(detected)
  n.detected
  
  inputs.initial <- nimData$inputs 
  nimData$inputs <- c(inputs.initial, rep(myVars$DETECTIONS$ninput_extra, J-nphases.initial))
  sparseData.obs <- nimData$sparseData
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
  
  
  nBatch <- myVars$MCMC$nBatch 
  seq.mat <- matrix(1:dim(posteriorSamples)[1], nrow = nBatch, byrow = T)
  
  Bvec.chains <- NULL
  kvec.chains <- c()
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
      sparseData <- sparseDataRep[[iter]]
      z <-   posteriorSamples.temp[iter,znames]
      Svec <- posteriorSamples.temp[iter,SvecNames]
      u <- do.call(rbind, lapply(1:M, function(i){getU(sparseData[i,], J)}))
      Bvec <- unlist(lapply(1:J, function(j) {sum(Svec[1:M]*z[1:M]*(1-u[1:M,j]))}))
     
      kvec <- unlist(lapply(1:length(epsVec), function(l){
        k <- which(Bvec<= epsVec[l])[1]
        if(is.na(k)) k <- J+1 
        return(k)
      }))
      
      out <- list(Bvec = Bvec,
                  kvec = kvec)
      return(out)
      
    }) 
    # moreSamples[[1]]
    
    Bvec.chains <- rbind(Bvec.chains, do.call(rbind, lapply(moreSamples, function(x) x$Bvec)))
    kvec.chains <- rbind(kvec.chains, do.call(rbind, lapply(moreSamples, function(x) x$kvec)))
    
  }#batch
  
  niter.per.chain <- nrow(myNimbleOutput.subset[[1]])
  knames <-paste0("k_eps", epsVec)
  Bnames <- paste0("B_", 1:J)
  # colSums(sparseDataRep[[iter]]))
  myNimbleOutput.kvec.Bvec <- as.mcmc.list(lapply(c(1:nchains), function(chain){
    this.chain.index <- (chain -1)*niter.per.chain + 1:niter.per.chain
    a <-  cbind(kvec.chains[this.chain.index,], Bvec.chains[this.chain.index,])
    dimnames(a)[[2]] <- c(knames, Bnames) #c('k', paste0("B_", 1:J))
    a <- mcmc(a)
    return(a)
  }))  
  
  
  outname <- paste0("NimbleOutfile_", this.pattern)
  
  save(myNimbleOutput, Runtime, 
       niter, burnin, nchains,
       myNimbleOutput.kvec.Bvec, knames, Bnames,
       file=file.path(outfolder, outname))
  
  file.list <- list.files(infolder)#, pattern = "SARE")
  
  }#i
  
}##--DO ALL

