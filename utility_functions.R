##=============================================================================
## dbin_geom: Function to create a NIMBLE custom distribution for faster model runs
##=============================================================================
#' @title Function to create a NIMBLE custom distribution for faster model runs.
#'
#' @description
#' \code{dbin_geom} returns the likelihood of a given bug & binary detection history y[i,1:n.detectors] 
#' 
#' @param x \code{Vector} containing detection/non detections 
#' @param p \code{Numeric} variable denoting the detection probability of each bug.
#' @param inputs A \code{Vector}  of dimensions n.phases with the number of test cases for each phase.
#' @param indicator A \code{Logical}; if == 0 no detections possible; saves time. 
#' @param log A \code{integer} indicates whether the log-probability or the probabilty should be returned (default to normal)


#'
#' @examples
#' y[i,1:n.phases] ~ dbin_geom(inputs, p,  z[i]==1)

#### 1.Density function ####
dbin_geom <- nimbleFunction(run = function( x = double(1)
                                            , inputs = double(1)
                                            , p = double(0)
                                            , indicator = double(0, default = 1.0)
                                            , log = integer(0, default = 0)){
  # RETURN TYPE DECLARATION
  returnType(double(0))
  
  ySparse <- x[1]
  detIndices <- x[2]
  ## CHECK INPUT DIMENSIONS
  if(detIndices == -1) detNums <- 0
  else detNums <- 1
  
  ## GET NECESSARY INFO
  nphases <- length(inputs)
  
  ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
  if(indicator == 0){
    if(detNums == 0){
      if(log == 0) return(1.0)
      else return(0.0)
    }else{
      if(log == 0) return(0.0)
      else return(-Inf)
    }
  }
  
  # ## CALCULATE THE LIKELIHOOD OF THE DETECTION VECTOR
  if (detNums > 0) {
    if(sum(inputs) == nphases){
      logProb <- log(p) + (detIndices-1)*log(1-p)  
    }
    if(sum(inputs) > nphases){
      y <- nimNumeric(length = nphases, value = 0, init = TRUE)
      y[detIndices] <- ySparse
      if(detIndices ==1){
        logProb <- dbinom(y[detIndices], prob = p, size = inputs[detIndices], log = TRUE)
      }
      if(detIndices >1){
        logProb <- sum(inputs[1:(detIndices-1)])*log(1-p)+
          dbinom(y[detIndices], prob = p, size = inputs[detIndices], log = TRUE)
        
      }
    }
  }
  # ## CALCULATE THE LIKELIHOOD OF THE DETECTION VECTOR
  if (detNums == 0) {
    if(sum(inputs) == nphases){
      logProb <- nphases*log(1-p)  
    }
    
    if(sum(inputs) > nphases){
      logProb <- sum(inputs[1:nphases])*log(1-p)
    }
  }
  
  if(log)return(logProb)
  return(exp(logProb))
})

#### 2.Sampling function ####
rbin_geom <- nimbleFunction(run = function( n = integer(0)
                                            , inputs = double(1)
                                            , p = double(0)
                                            , indicator = double(0, default = 1.0)){
  # # Return type declaration
  returnType(double(1))
  
  ## GET NECESSARY INFO
  nphases <- length(inputs)
  
  ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
  if(indicator == 0){return(c(-1,-1))}
  
  ## INITIALIZE THE OUTPUT VECTOR OF DETECTIONS
  detectOut <- rep(0, nphases)
  
  ## SAMPLE THE DETECTION HISTORY (FOR RELEVANT DETECTORS ONLY)
  for(j in 1:nphases){
    # Draw the observation at detector j from a binomial distribution with probability p
    detectOut[j] <- rbinom(1, size = inputs[j], prob = p)
    if(detectOut[j] > 0) return(c(detectOut[j],j))
  }#j
  if(sum(detectOut) == 0){return(c(-1,-1))}
  
  
})

#### 3.Registration ####
registerDistributions(list(
  dbin_geom = list(
    BUGSdist = "dbin_geom(inputs, p, indicator)",
    types = c( "value = double(1)"
               , "inputs = double(1)"
               , "p = double(0)"
               , "indicator = double(0)"),
    pqAvail = FALSE,
    mixedSizes = TRUE)))

##====================================
## Simulation of Posterior replicates
##====================================
simulateFromPosterior <- nimbleFunction(
  setup = function(model, nodes) { # 'setup' function  is written is R language
    mv <- modelValues(model)
    allVars <- model$getVarNames()
    downNodes <- model$getDependencies(nodes = nodes, self = FALSE, includeData = TRUE)
  },
  run = function(posteriorSamples = double(2)) { # 'run' function is written is nimble language
    nIter <- dim(posteriorSamples)[1]
    nNodes <- dim(posteriorSamples)[2]
    resize(mv, nIter)
    for(ii in 1:nIter){
      if (ii%%500 == 0)  cat('..... drawing sample #', ii, '\n')
      values(model, nodes) <<- posteriorSamples[ii,1:nNodes]
      model$simulate(nodes = downNodes, includeData = TRUE) # Specified 'nodes' (by MCMC samples) should not get updated 
      model$calculate(nodes = downNodes) # Log-likelihood p(yrep | theta_iter) of the simulated data 
      copy(from = model, nodes = allVars, to = mv, rowTo = ii, logProb = TRUE)
    }#ii
  }) 

##====================================
## Calculation 'u' variable from data
##====================================
getU <- nimbleFunction(run = function(x = double(1), n = double(0)) {
  u <- nimNumeric(length = n, value = 0, init = TRUE)
  if(x[1] != -1){
    this.phase <- x[2]
    u[this.phase:n] <- 1
  }
  returnType(double(1))
  return(u)
})

##=============================================================================
## Function to Generate a dataframe of parameter space for simulation purposes
##=============================================================================

MakeParameterSpace <- function(list.param = list.param
                               ,n.rep = 10 
){
  param.df   <- expand.grid(list.param)
  n.comb <- nrow(param.df)
  param.df <- param.df[rep(seq_len(nrow(param.df)), n.rep), ]
  param.df$set_ID <- rep(seq_len(n.comb), n.rep)
  
  param.df <- param.df[sort(param.df$set_ID),]
  param.df$rep_ID <- rep(c(1:n.rep), n.comb)
  
  return(param.df)
}  

##=============================================================================
## Function to increase the dimensions of an object
##=============================================================================

MakeAugmentation <- function( y,
                              aug.factor= NULL,
                              aug.years = NULL,
                              replace.value = NA){
  ## Vector Data augmentation
  if(is.vector(y)){
    if(is.null(names(y))){ 
      names(y) <- 1:length(y) 
    }
    
    if(!is.null(aug.factor)){
      y.aug <- c(y, rep(replace.value, round(length(y) * aug.factor)))
      names(y.aug) <- c(names(y), rep("Augmented", round(length(y) * aug.factor)))
      y <- y.aug
    }
    
    if(!is.null(aug.years)){
      y.aug <- c(y, rep(replace.value, aug.years))
      maxYear <- as.numeric(names(y)[length(names(y))])
      augmentedYears <- (maxYear + 1) : (maxYear + aug.years)
      names(y.aug) <- c(names(y), augmentedYears)
      y <- y.aug
    }
  }
  
  ## Matrix augmentation
  if(is.matrix(y)){
    if(is.null(dimnames(y))){
      dimnames(y) <- list(1:dim(y)[1], 1:dim(y)[2])
    }
    
    if(!is.null(aug.factor)){
      n.tot <- round(dim(y)[1]*(1 + aug.factor))
      y.aug <- matrix(replace.value, n.tot, dim(y)[2])
      y.aug[1:dim(y)[1], ] <- y
      dimnames(y.aug) <- list(c( dimnames(y)[[1]], rep("Augmented", n.tot - dim(y)[1])),
                              dimnames(y)[[2]])
      y <- y.aug
    }
    
    if(!is.null(aug.years)){
      y.aug <- matrix(replace.value, dim(y)[1], dim(y)[2] + aug.years)
      y.aug[ ,1:dim(y)[2]] <- y
      maxYear <- as.numeric(dimnames(y)[[2]][length(dimnames(y)[[2]])])
      augmentedYears <- (maxYear + 1) : (maxYear + aug.years)
      dimnames(y.aug) <- list( dimnames(y)[[1]],
                               c(dimnames(y)[[2]], augmentedYears))
      y <- y.aug
    }
  }  
  
  
  return (y)
}


