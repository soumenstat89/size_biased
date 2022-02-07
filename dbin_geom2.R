#' @title Function to create a NIMBLE custom distribution for faster SCR model runs.
#'
#' @description
#' \code{dbin_geom2} returns the likelihood of a given individual & spatial binary detection history y[i,1:n.detectors] 
#' 
#' @param x \code{Vector} containing detection/non detections 
#' @param p \code{Numeric} variable denoting the detection probability of each bug.
#' @param inputs A \code{Vector}  of dimensions n.phases with the number of test cases for each phase.
#' @param indicator A \code{Logical}; if == 0 no detections possible; saves time. 
#' @param log A \code{integer} indicates whether the log-probability or the probabilty should be returned (default to normal)


#'
#' @examples
#' y[i,1:n.phases] ~ dbin_geom2(inputs, p,  z[i]==1)

#### 1.Density function ####
dbin_geom2 <- nimbleFunction(run = function( x = double(1)
                                                , inputs = double(1)
                                                , p = double(0)
                                                , indicator = double(0, default = 1.0)
                                                , log = integer(0, default = 0)){
  # RETURN TYPE DECLARATION
  returnType(double(0))
  
  ySparse <- x[1]
  this.phase <- x[2]
  ## CHECK INPUT DIMENSIONS
  if(this.phase == -1) detCount <- 0
  else detCount <- 1
  
  ## GET NECESSARY INFO
  nphases <- length(inputs)
  
  ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
  if(indicator == 0){
    if(detCount == 0){
      if(log == 0) return(1.0)
      else return(0.0)
    }else{
      if(log == 0) return(0.0)
      else return(-Inf)
    }
  }

  ## RECREATE Y

  # ## CALCULATE THE LIKELIHOOD OF THE DETECTION VECTOR
  if (detCount > 0) {
    if(sum(inputs) == nphases){
      logProb <- log(p) + (this.phase-1)*log(1-p)  
    }# if(sum(inputs) == nphases)
    
    if(sum(inputs) > nphases){
      y <- nimNumeric(length = nphases, value = 0, init = TRUE)
      y[this.phase] <- ySparse
      if(this.phase ==1){
        logProb <- dbinom(y[this.phase], prob = p, size = inputs[this.phase], log = TRUE)
      }# if(this.phase ==1)
      if(this.phase >1){
        if(sum(inputs[1:(this.phase-1)]) == 0){ 
             logProb <- dbinom(y[this.phase], prob = p, size = inputs[this.phase], log = TRUE)
          }
        if(sum(inputs[1:(this.phase-1)]) > 0){  
             logProb <- sum(inputs[1:(this.phase-1)])*log(1-p)+
                        dbinom(y[this.phase], prob = p, size = inputs[this.phase], log = TRUE)
          }
        
      }# if(this.phase >1)
    }#if(sum(inputs) > nphases)
  }#if (detCount > 0) 
    if (detCount == 0) {
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
rbin_geom2 <- nimbleFunction(run = function( n = integer(0)
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
  dbin_geom2 = list(
    BUGSdist = "dbin_geom2(inputs, p, indicator)",
    types = c( "value = double(1)"
               , "inputs = double(1)"
               , "p = double(0)"
               , "indicator = double(0)"),
    pqAvail = FALSE,
    mixedSizes = TRUE)))


