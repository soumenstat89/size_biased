############################################################################################
####### --------------- DATA SIMULATION & NIMBLE MODEL FITTING SCRIPT ----------------########
####### -------------- BAYESIAN SIZE BIASED MODEL FOR SOFTWARE TESTING  ---------------- ########
############################################################################################

{ ##--DO ALL
  rm(list=ls())
  cat("\014")
  gc()
  options(scipen = 999) #  no. of printed digits = 999+5 = 1004
  ts = format(Sys.time(), "%y%m%d_%H%M%S")
  ## ------ IMPORT REQUIRED LIBRARIES ------
  library(extraDistr)
  library(coda)               # Import the MCMC diagnostic tools
  library(nimble)             # Import the NIMBLE subroutines
  
 
  myVars <- list(WD = getwd(),
                 FunD = getwd(),
                 modelName = "SoftwareBugData",
                 shape.S = 0.5, 
                 rate.S = 0.01,
                 M = 500,
                 J = 50,
                 epsVec =  c(25,50,75, 100, 150, 200),
                 ninput_extra = 3000, 
                 ninput_extra_set  = c(1000, 2000, 3000, 4000, 5000),
                 nBatch = 500,
                 numIters = 10000, 
                 numBurnIn = 5000,
                 numChains = 3,
                 numThin = 1)
  
  source("dbin_geom.R")

  source("simulateFromPosterior.R")
  
  getU <- nimbleFunction(run = function(x = double(1), n = double(0)) {
    u <- nimNumeric(length = n, value = 0, init = TRUE)
    if(x[1] != -1){
      this.phase <- x[2]
      u[this.phase:n] <- 1
    }
    returnType(double(1))
    return(u)
  })

}##-- DO ALL

# ##-----------------
{##-- DO ALL
  if(is.null(myVars$modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
  if(is.null(myVars$WD))stop("YOU SHOULD PROBABLY CHOOSE A WORKING DIRECTORY FOR THIS ANALYSIS/MODEL")
  if(!dir.exists(file.path(myVars$WD))){dir.create(file.path(myVars$WD))}
  if(!dir.exists(file.path(myVars$WD,myVars$modelName))){dir.create(file.path(myVars$WD, myVars$modelName))}
  
  outfolder <- file.path(myVars$WD,myVars$modelName,paste0("Bug_OutFiles_shape.S", myVars$shape.S,"_ninputExtra", myVars$ninput_extra, "_", ts))
  if(!dir.exists(outfolder)){dir.create(outfolder,recursive=TRUE)}
  
  shape.S = myVars$shape.S 
  rate.S = myVars$rate.S
  M <- myVars$M 
  
  bugd0 <- read.csv("Software_Testing_Data_clean_SD.csv", header = T, sep = ",")
  bugd <- bugd0[, !colnames(bugd0) %in% c("Defect.Id")]
  
  
  ##--- NUMBER OF PHASES ---##
  nphases <- length(unique(bugd$Cycle))
 
  ##--- NUMBER OF INPUTS IN EACH PHASES ---##
  inputs <- rep(NA, nphases)
  for(j in 1:nphases){
    a1 <- bugd[,1]
    inputs[j] <- length(a1[bugd[,"Cycle"] == paste("Cycle ",j, sep= "")])
  }
  
  inputs
  
  n.inputs <- sum(inputs, na.rm = T)
  n.inputs
  
  
  ##--- DETECTION-NONDETECTION DATA MATRIX ---##
  detections <- bugd[, !colnames(bugd) %in% c("Test.Case.id", "Cycle")] # n.inputs x max no. of detection in a single input 
  detections[detections == ""] <- NA
 
  totaldet <- sum(!is.na(detections))
  totaldet
  
  ##--- BUG ID ---##
  bugids <- unique(detections[!is.na(detections)]) # n.detected x 1
  bugids
  
  ##---  NUMBER OF UNIQUE BUG DETECTED ---##
  n.detected <- length(bugids)
  n.detected 
  
  
  ##---  NO. OF TIMES EACH BUG WAS DETECTED ---##
  detfreq <- unlist(lapply(1:n.detected, function(ii){ # n.detected x 1
    sum(detections == bugids[ii], na.rm = T)
  }))
  detfreq 
  
  
  ##--- OBTAINING THE PHASES WHERE EACH BUG IS DETECTED ---##
  id <- 215
  detphases.list <- lapply(1:n.detected, function(id){
    detections2 <- detections
    detections2[is.na(detections2)] <- 10^5
    rows <- which(rowSums(detections2 == bugids[id]) >0)
    out <- sort(unique(unlist(lapply(rows, function(this.row){
      as.numeric(strsplit(bugd[this.row,"Cycle"], split = " ")[[1]][2])
    }))))
    return(out)
  })
  
  ##-- ONLY TAKING THE LAST PHASE WHERE A BUG IS DETECTED -- #
  detphases_clean <- unlist(lapply(detphases.list, function(x){x[length(x)]}))
  
  ##--- CAPTURE-RECAPTURE DATA SET IN SPARSE FORMAT ---##
  this.data <- cbind(detfreq, detphases_clean)
  
  indices.list <- lapply(1:nphases, function(i) c(1:length(detphases_clean))[which(detphases_clean==i)])
  n.detectedPerPhase <- unlist(lapply(1:nphases, function(i) length(detphases_clean[indices.list[[i]]])))
  n.detectionsPerPhase <- unlist(lapply(1:nphases, function(i) sum(detfreq[indices.list[[i]]])))

  sparseData <- rbind(this.data, matrix(-1, M - n.detected, 2))
  
  z <- c(rep(1, n.detected), rep(NA, M - n.detected))
  u <- do.call(rbind, lapply(1:M, function(i){getU(sparseData[i,], nphases)}))
  
 
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
    M = M, 
    shape.S = shape.S,
    rate.S = rate.S
  )
  ### ====      3.8.3.Define the data to be used in the model ====
  
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
  
  i <- 1
  for(i in 1:n.detected){
    lambda.S.init[i] <- sparseData[i,2] +5 
    Svec.init[i] <- rtpois(1, lambda.S.init[i], a = 1, b = 15)
    
  }
  lambda.S.init[is.na(Svec.init)] <- rgamma(M-n.detected, shape = 200, rate = 2)
  Svec.init[is.na(Svec.init)] <- rtpois(M-n.detected, lambda.S.init[is.na(Svec.init)], a = 5, b = 1000)
  
  nimInits <- list(
    z = z.init,           
    psi = runif(1,0,1),
    logitr = logit(runif(1, 10^(-6), 10^(-5))),
    lambda.S = lambda.S.init,
    Svec = Svec.init
    
  )
  
  
  ### ====      3.8.5.Define the parameters to be monitored ====
  nimParams <- c("N","psi", "z", "Svec", "lambda.S", "r", "logitr")
  
 
  inname <- file.path(outfolder, paste("InputForBugData_", ts, ".RData", sep = ""))
  
  save(myVars, sparseData, u, z, detections, detphases.list, detphases_clean, detfreq,
       bugd, bugd0, bugids, n.detected, totaldet, inputs,  M, n.detected, 
       nimConstants, nimData, nimInits, nimParams, modelCode, 
       file = inname)
  
  save(n.detectedPerPhase, n.detectionsPerPhase,
       file = file.path(outfolder,"DetectionsPerPhase.RData"))
  
}##--DO ALL

##############################################################################################
##############################################################################################
#=========================# MCMC # MCMC # MCMC # MCMC # MCMC # MCMC #========================#
##############################################################################################
##############################################################################################

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
  
  ### ====    6.3.RUN THE MCMC ====
  {##-- DO ALL
    
    burnin = 0
    niter =   myVars$numIters 
    nchains =  myVars$numChains 
    
    Runtime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC, 
                                                      nburnin = burnin, 
                                                      niter = niter, 
                                                      nchains = nchains, #6, 
                                                      samplesAsCodaMCMC = TRUE))
    
  }##-- DO ALL
  
  outname <- paste("OutputForBugData_Nr_niter", niter, "_", ts, ".RData", sep = "")
  fname = file.path(outfolder, outname)
  
  save(myNimbleOutput, Runtime, 
       niter, burnin, nchains,
       file=fname)
  
}##--DO ALL


##############################################################################################
##############################################################################################
#========================# RELIABILITY # RELIABILITY # RELIABILITY # #=======================#
##############################################################################################
##############################################################################################
{##--DO ALL
  
  ninput_extra_set <- myVars$ninput_extra_set
  epsVec <-  myVars$epsVec 

  rel.arr <- array(0, dim = c(length(epsVec),  myVars$J, length(ninput_extra_set)))
  dimnames(rel.arr) <- list(paste0("eps", epsVec),
                            paste0("phase", 1: myVars$J), 
                            paste0("extraTestCase", ninput_extra_set))

  cc <- 0
  ninput_extra <- ninput_extra_set[2]
  
  for(ninput_extra in ninput_extra_set){ 
    cat('------ extra input ', ninput_extra, '\n', sep = "")
    
    load( file.path(outfolder, paste("InputForBugData_", ts, ".RData", sep = "")))
    load( file.path(outfolder,list.files(outfolder, pattern = "OutputForBugData_Nr")[1]))
    burnin <- myVars$numBurnIn 
    n.iterations <- nrow(myNimbleOutput[[1]])
    nchains <- length(myNimbleOutput)
    
  myNimbleOutput.subset <-  as.mcmc.list(lapply(myNimbleOutput,
                                                function(x){mcmc(x[(burnin+1):n.iterations,])}))
  
  if(length(myNimbleOutput.subset) > 1){posteriorSamples <- do.call(rbind, myNimbleOutput.subset)} 
  if(length(myNimbleOutput.subset) == 1){posteriorSamples <- as.matrix(myNimbleOutput.subset)}
  print(dim(posteriorSamples)[1])
  
  
  M <- nimConstants$M
  J <- myVars$J 
  epsVec <-  myVars$epsVec 
  
  nphases.initial <- nimConstants$nphases
  nimConstants$nphases <- J # The total number of possible phases (only first "nphases" were performed)
  
  inputs.initial <- nimData$inputs 
  nimData$inputs <- c(inputs.initial, rep(ninput_extra, J-nphases.initial))
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
  
  ## RUN THE COMPILED VERSION (much faster) 
  
  nBatch <- myVars$nBatch 
  seq.mat <- matrix(1:dim(posteriorSamples)[1], nrow = nBatch, byrow = T)
  Bvec.chains <- NULL
  kvec.chains <- c()
  batch <-  1
  
  for(batch in 1:nBatch){
    
    if (batch%%50 == 0)  cat('..... batch ', batch, '\n', sep = "")
    posteriorSamples.temp <- posteriorSamples[seq.mat[batch,],]
    
    csimFromPost$run(posteriorSamples = posteriorSamples.temp) 
    
    sparseDataRep <- csimFromPost$mv[["sparseData"]] 
   
    colSums(sparseDataRep[[5]])
    
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
    
    Bvec.chains <- rbind(Bvec.chains, do.call(rbind, lapply(moreSamples, function(x) x$Bvec)))
    kvec.chains <- rbind(kvec.chains, do.call(rbind, lapply(moreSamples, function(x) x$kvec)))
    
  }#batch
  
  niter.per.chain <- nrow(myNimbleOutput.subset[[1]])
  knames <-paste0("k_eps", epsVec)
  
 
  myNimbleOutput.new <- as.mcmc.list(lapply(c(1:nchains), function(chain){
    
    this.chain.index <- (chain -1)*niter.per.chain + 1:niter.per.chain
    a <-  cbind(kvec.chains[this.chain.index,], Bvec.chains[this.chain.index,])
    dimnames(a)[[2]] <- c(knames, paste0("B_", 1:J)) 
    a <- mcmc(a)
    return(a)
  }))  
  
  Bnames.full <- paste0("B_", 1: dim(Bvec.chains)[2])
  colnames(Bvec.chains) <- Bnames.full
  
  rel.mat <- do.call(rbind, lapply(1:length(epsVec), function(l){
    
    out <- unlist(lapply(Bnames.full, function(this.B){
      rel <- mean(Bvec.chains[, this.B]<= epsVec[l])
      if(is.na(rel)) rel <- 0 
      return(rel)
    }))
    return(out)
  }))
  cc <- cc + 1
  rel.arr[,,cc] <- rel.mat
 
  
  outnamekB <- paste("OutputForBugDatakB_J", J, "_ninputExtra", ninput_extra,"_", ts, ".RData", sep = "")
  fname = file.path(outfolder, outnamekB)
  
  save(myNimbleOutput.new, 
       knames, Bnames,
       rtime, Runtime,
       n.iterations, burnin, nchains,
       J, epsVec, ninput_extra, ninput_extra_set,
       file=fname)
  
  try(rm(myNimbleOutput.new), silent = TRUE)
 
  try(rm(Bvec.chains), silent = TRUE)
  try(rm(kvec.chains), silent = TRUE)
  
  }#ninput_extra
  
  dimnames(rel.arr) <- list(paste0("eps",epsVec), Bnames.full, ninput_extra_set) 
  
  outname.rel <- paste("ReliabilityEstimate_J", J,"_", ts, ".RData", sep = "")
  
  save(rel.arr, mean.out, ninput_extra_set,
       Bnames, Bnames.full, epsVec, ninput_extra_set,
       file = file.path(outfolder, outname.rel))
  
}##--DO ALL

# ###############################################################################

