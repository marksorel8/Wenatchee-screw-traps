TimeStratPetersenDiagError_2<-function (title, prefix, time, n1, m2, u2, jump.after = NULL, 
          logitP.cov, n.chains = 3, n.iter = 2e+05, n.burnin = 1e+05, 
          n.sims = 2000, tauU.alpha = 1, tauU.beta = 0.05, taueU.alpha = 1, 
          taueU.beta = 0.05, mu_xiP = logit(sum(m2, na.rm = TRUE)/sum(n1, 
                                                                      na.rm = TRUE)), tau_xiP = 1/var(logit((m2 + 0.5)/(n1 + 
                                                                                                                          1)), na.rm = TRUE), tauP.alpha = 0.001, tauP.beta = 0.001, 
          debug = FALSE, debug2 = FALSE, engine = c("jags", "openbugs")[1], 
          InitialSeed) 
{
  set.seed(InitialSeed)
  working.directory <- getwd()
  model.file <- file.path(working.directory, "model.txt")
  data.file <- file.path(working.directory, "data.txt")
  init.files <- file.path(working.directory, paste("inits", 
                                                   1:n.chains, ".txt", sep = ""))
  sink(model.file)
  cat("\nmodel{\n# Time Stratified Petersen with Diagonal recapture (no spillover in subsequent weeks or marked fish)\n#    and allowing for error in the smoothed U curve.\n\n# Refer to Bonner (2008) Ph.D. thesis from Simon Fraser University available at\n#     http://www.stat.sfu.ca/people/alumni/Theses/Bonner-2008.pdf\n# The model is in Appendix B. The discussion of the model is in Chapter 2.\n\n#  Data input:\n#      Nstrata - number of strata\n#      n1         - number of marked fish released\n#      m2         - number of marked fish recaptured\n#      u2         - number of unmarked fish captured (To be expanded to population).\n#      logitP.cov   - covariates for logitP\n#      NlogitP.cov  - number of logitP covariates\n#      SplineDesign- spline design matrix of size [Nstrata, maxelement of n.b.notflat]\n#                   This is set up prior to the call.\n#      b.flat   - vector of strata indices where the prior for the b's will be flat.\n#                 this is normally the first two of each spline segment\n#      n.b.flat - number of b coefficients that have a flat prior\n#      b.notflat- vector of strata indices where difference in coefficients is modelled\n#      n.b.notflat- number of b coefficients that do not have a flat prior\n#      tauU.alpha, tauU.beta - parameters for prior on tauU\n#      taueU.alpha, taueU.beta - parameters for prior on taueU\n#      mu_xiP, tau_xiP  - parameters for prior on mean logit(P)'s [The intercept term]\n#                         mu_xiP = mean; tau_xiP = 1/variance\n#                       - the other beta terms are given a prior of a N(mu=0, variance=1000)\n#      tauP.alpha, tauP.beta - parameter for prior on tauP (residual variance of logit(P)'s after adjusting for\n#                         covariates)\n#\n#  Parameters of the model are:\n#      p[i]\n#       logitP[i]  = logit(p[i]) = logitP.cov*beta.logitP\n#         The first beta.logitP has a prior from N(xiP, tauP)\n#            and xiP and tauP are currently set internally\n#         The remaining beta's are assigned a wider prior N(mu=0, var=1000).\n#      U[i]\n#       etaU[i]  = log(U[i])\n#         which comes from spline with parameters bU[1... Knots+q]\n#         + error term eU[i]\n\n   ##### Fit the spline and specify hierarchial model for the logit(P)'s ######\n   for(i in 1:Nstrata){\n        logUne[i] <- inprod(SplineDesign[i,1:n.bU],bU[1:n.bU])  # spline design matrix * spline coeff\n", 
      fill = TRUE)
  sink()
  if (tolower(engine) == "jags") {
    sink("model.txt", append = TRUE)
    cat("\n        etaU[i] ~ dnorm(logUne[i], taueU)T(,20)    # add random error\n   ", 
        fill = TRUE)
    sink()
  }
  if (tolower(engine) %in% c("openbugs")) {
    sink("model.txt", append = TRUE)
    cat("\n        etaU[i] ~ dnorm(logUne[i], taueU)C(,20)    # add random error\n   ", 
        fill = TRUE)
    sink()
  }
  sink("model.txt", append = TRUE)
  cat("\n        eU[i] <- etaU[i] - logUne[i]\n        mu.logitP[i] <- inprod(logitP.cov[i,1:NlogitP.cov], beta.logitP[1:NlogitP.cov])\n\n#        logitPu[i] ~ dnorm(mu.logitP[i],tauP)       # uncontrained logitP value\n#        logitP [i] <- max(-10,min(10, logitPu[i]))  # keep the logit from wandering too far off\n\n        ## Matt's fix to improve mixing. Use u2copy to break the cycle (this doesn't work??)\n        mu.epsilon[i] <- mu.logitP[i] - log(u2copy[i] + 1) + etaU[i]   \n        epsilon[i] ~ dnorm(mu.epsilon[i],tauP)                     \n        logitP[i] <- max(-10,min(10,log(u2copy[i] + 1) - etaU[i] + epsilon[i]))  \n   }\n\n   ##### Hyperpriors #####\n   ## Run size - flat priors\n   for(i in 1:n.b.flat){\n", 
      fill = TRUE)
  sink()
  if (tolower(engine) == "jags") {
    sink("model.txt", append = TRUE)
    cat("\n      bU[b.flat[i]] ~ dnorm(0.0,1.0E-6) \n   ", 
        fill = TRUE)
    sink()
  }
  if (tolower(engine) %in% c("openbugs")) {
    sink("model.txt", append = TRUE)
    cat("\n      bU[b.flat[i]] ~ dflat()\n   ", fill = TRUE)
    sink()
  }
  sink("model.txt", append = TRUE)
  cat("\n   }\n   ## Run size - priors on the difference\n   for(i in 1:n.b.notflat){\n      xiU[b.notflat[i]] <- 2*bU[b.notflat[i]-1] - bU[b.notflat[i]-2]\n      bU [b.notflat[i]] ~ dnorm(xiU[b.notflat[i]],tauU)\n   }\n   tauU ~ dgamma(tauU.alpha,tauU.beta)  # Notice reduction from .0005 (in thesis) to .05\n   sigmaU <- 1/sqrt(tauU)\n   taueU ~ dgamma(taueU.alpha,taueU.beta) # dgamma(100,.05) # Notice reduction from .0005 (in thesis) to .05\n   sigmaeU <- 1/sqrt(taueU)\n\n   ## Capture probabilities. The logit(p[i]) are n(logitP.cov*beta.logitP.cov, sigmaP**2)\n   beta.logitP[1] ~ dnorm(mu_xiP,tau_xiP) # first term is usually an intercept\n   for(i in 2:NlogitP.cov){\n      beta.logitP[i] ~ dnorm(0, .01)   # rest of beta terms are normal 0 and a large variance\n   }\n   beta.logitP[NlogitP.cov+1] ~ dnorm(0, .01) # dummy so that covariates of length 1 function properly\n   tauP ~ dgamma(tauP.alpha,tauP.beta)\n   sigmaP <- 1/sqrt(tauP)\n\n   ##### Likelihood contributions #####\n   for(i in 1:Nstrata){\n      logit(p[i]) <- logitP[i]       # convert from logit scale\n      U[i]   <- round(exp(etaU[i]))       # convert from log scale\n      m2[i] ~ dbin(p[i],n1[i])     # recovery of marked fish\n      u2[i] ~ dbin(p[i],U [i])      # capture of newly unmarked fish\n   }\n\n   ##### Derived Parameters #####\n   Utot <- sum( U[1:Nstrata])          # Total number of unmarked fish\n   Ntot <- sum(n1[1:Nstrata]) + Utot  # Total population size including those fish marked and released\n} # end of model\n\n", 
      fill = TRUE)
  sink()
  Nstrata <- length(n1)
  ext.jump <- c(0, jump.after, Nstrata)
  SplineDesign <- matrix(0, nrow = 0, ncol = 0)
  SplineDegree <- 3
  b.flat <- NULL
  b.notflat <- NULL
  all.knots <- NULL
  for (i in 1:(length(ext.jump) - 1)) {
    nstrata.in.set <- ext.jump[i + 1] - ext.jump[i]
    if (nstrata.in.set > 7) {
      knots <- seq(5, nstrata.in.set - 1, 4)/(nstrata.in.set + 
                                                1)
    }
    else {
      knots <- 0.5
    }
    all.knots <- c(all.knots, knots)
    z <- bs((1:nstrata.in.set)/(nstrata.in.set + 1), knots = knots, 
            degree = SplineDegree, intercept = TRUE, Boundary.knots = c(0, 
                                                                        1))
    b.flat <- c(b.flat, ncol(SplineDesign) + (1:2))
    b.notflat <- c(b.notflat, ncol(SplineDesign) + 3:(ncol(z)))
    SplineDesign <- cbind(SplineDesign, matrix(0, nrow = nrow(SplineDesign), 
                                               ncol = ncol(z)))
    SplineDesign <- rbind(SplineDesign, cbind(matrix(0, nrow = nrow(z), 
                                                     ncol = ncol(SplineDesign) - ncol(z)), z))
  }
  n.b.flat <- length(b.flat)
  n.b.notflat <- length(b.notflat)
  n.bU <- n.b.flat + n.b.notflat
  logitP.cov <- as.matrix(logitP.cov)
  NlogitP.cov <- ncol(as.matrix(logitP.cov))
  
  ################# Modifications from Thomas Buehrens ######################
 
  
  #new code from TB
   u2_approx<-log(u2+1) #create log-transformed u2+1
  if(is.na(u2_approx[1])){u2_approx[1]<-u2_approx[!is.na(u2_approx)][1]} #if first value is missing, replace with first observed count
  if(is.na(u2_approx[length(u2_approx)])){u2_approx[length(u2_approx)]<-u2_approx[!is.na(u2_approx)][length(u2_approx[!is.na(u2_approx)])]}#if las value is missing, replace with first observed count
  #endof new code by TB
 
  u2copy <- exp(spline(x = 1:Nstrata, y = (u2_approx), xout = 1:Nstrata)$y)-1 #modified by TB to fit a spline to the log-transformed data, and then back transform to untransformed space, then make values less than 1 equal to 1.
  #u2copy[u2copy<=1]<-1;
  ####################### end of TB modification ########################### 
  
  
  u2copy <- round(u2copy)
  
  
  
  
  datalist <- list("Nstrata", "n1", "m2", 
                   "u2", "u2copy", "logitP.cov", "NlogitP.cov", 
                   "SplineDesign", "b.flat", "n.b.flat", 
                   "b.notflat", "n.b.notflat", "n.bU", 
                   "tauU.alpha", "tauU.beta", "taueU.alpha", 
                   "taueU.beta", "mu_xiP", "tau_xiP", 
                   "tauP.alpha", "tauP.beta")
  avgP <- sum(m2, na.rm = TRUE)/sum(n1, na.rm = TRUE)
  Uguess <- pmax((u2 + 1) * (n1 + 2)/(m2 + 1), u2/avgP, 1, 
                 na.rm = TRUE)
  Uguess[which(is.na(Uguess))] <- mean(Uguess, na.rm = TRUE)
  init.bU <- lm(log(Uguess + 1) ~ SplineDesign - 1)$coefficients
  if (debug2) {
    cat("compute init.bU \n")
    browser()
  }
  logitPguess <- pmax(-10, pmin(10, logit((m2 + 1)/(n1 + 1))))
  init.beta.logitP <- as.vector(lm(logitPguess ~ logitP.cov - 
                                     1)$coefficients)
  if (debug2) {
    cat(" obtained initial values of beta.logitP\n")
    browser()
  }
  pdf(file = paste(prefix, "-initialU.pdf", sep = ""))
  plot(time, log(Uguess), main = paste(title, "\nInitial spline fit to estimated U[i]"), 
       ylab = "log(U[i])", xlab = "Stratum")
  lines(time, SplineDesign %*% init.bU)
  dev.off()
  parameters <- c("logitP", "beta.logitP", "tauP", 
                  "sigmaP", "bU", "tauU", "sigmaU", 
                  "eU", "taueU", "sigmaeU", "Ntot", 
                  "Utot", "logUne", "etaU", "U")
  if (any(is.na(m2))) {
    parameters <- c(parameters, "m2")
  }
  if (any(is.na(u2))) {
    parameters <- c(parameters, "u2")
  }
  init.vals <- genInitVals(model = "TSPDE", n1 = n1, 
                           m2 = m2, u2 = u2, logitP.cov = logitP.cov, SplineDesign = SplineDesign, 
                           n.chains = n.chains)
  data.list <- list()
  for (i in 1:length(datalist)) {
    data.list[[length(data.list) + 1]] <- get(datalist[[i]])
  }
  names(data.list) <- as.list(datalist)
  results <- run.MCMC(modelFile = model.file, dataFile = data.file, 
                      dataList = data.list, initFiles = init.files, initVals = init.vals, 
                      parameters = parameters, nChains = n.chains, nIter = n.iter, 
                      nBurnin = n.burnin, nSims = n.sims, overRelax = FALSE, 
                      initialSeed = InitialSeed, working.directory = working.directory, 
                      engine = engine, debug = debug)
  return(results)
}