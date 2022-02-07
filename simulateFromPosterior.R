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

###==== Fast posterior predictive sampling using a nimbleFunction
ppSamplerNF <- nimbleFunction(
  setup = function(model, mcmc) {
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)
    cat("Stochastic parents of data are:", paste(parentNodes, collapse = ','), ".\n")
    simNodes <- model$getDependencies(parentNodes, self = FALSE)
    vars <- mcmc$mvSamples$getVarNames()  # need ordering of variables in mvSamples / samples matrix
    cat("Using posterior samples of:", paste(vars, collapse = ','), ".\n")
    n <- length(model$expandNodeNames(dataNodes, returnScalarComponents = TRUE))
  },
  run = function(samples = double(2)) {
    nSamp <- dim(samples)[1]
    ppSamples <- matrix(nrow = nSamp, ncol = n)
    for(i in 1:nSamp) {
      values(model, vars) <<- samples[i, ]
      model$simulate(simNodes, includeData = TRUE)
      ppSamples[i, ] <- values(model, dataNodes)
    }
    returnType(double(2))
    return(ppSamples)
  })

# ## Create the sampler for this model and this MCMC.
# ppSampler <- ppSamplerNF(model, mcmc)
# ## Stochastic parents of data are: mu,sigma .
# ## Using posterior samples of: mu,sigma .
# cppSampler <- compileNimble(ppSampler, project = model)
# 
# ## Check ordering of variables is same in 'vars' and in 'samples'.
# colnames(samples)
# ## [1] "mu"    "sigma"
# identical(colnames(samples), model$expandNodeNames(mcmc$mvSamples$getVarNames()))
# ## [1] TRUE
# set.seed(1)
# system.time(ppSamples_via_nf <- cppSampler$run(samples))


