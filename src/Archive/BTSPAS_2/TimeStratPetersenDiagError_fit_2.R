function (title = "TSDPE", prefix = "TSPDE-", time, 
          n1, m2, u2, sampfrac, jump.after = NULL, bad.n1 = c(), bad.m2 = c(), 
          bad.u2 = c(), logitP.cov = rep(1, length(n1)), n.chains = 3, 
          n.iter = 2e+05, n.burnin = 1e+05, n.sims = 2000, tauU.alpha = 1, 
          tauU.beta = 0.05, taueU.alpha = 1, taueU.beta = 0.05, mu_xiP = logit(sum(m2, 
                                                                                   na.rm = TRUE)/sum(n1, na.rm = TRUE)), tau_xiP = 1/var(logit((m2 + 
                                                                                                                                                  0.5)/(n1 + 1)), na.rm = TRUE), tauP.alpha = 0.001, tauP.beta = 0.001, 
          run.prob = seq(0, 1, 0.1), debug = FALSE, debug2 = FALSE, 
          engine = c("jags", "openbugs")[1], InitialSeed = ceiling(runif(1, 
                                                                         min = 0, max = if (engine == "jags") {
                                                                           1e+06
                                                                         } else {
                                                                           14
                                                                         }))) 
{
  version <- "2014-09-01"
  options(width = 200)
  time <- as.vector(time)
  n1 <- as.vector(n1)
  m2 <- as.vector(m2)
  u2 <- as.vector(u2)
  sampfrac <- as.vector(sampfrac)
  if (var(c(length(n1), length(m2), length(u2), length(sampfrac), 
            length(time))) > 0) {
    cat("***** ERROR ***** Lengths of n1, m2, u2, sampfrac, time must all be equal. They are:", 
        length(n1), length(m2), length(u2), length(sampfrac), 
        length(time), "\n")
    return()
  }
  if (length(logitP.cov)%%length(n1) != 0) {
    cat("***** ERROR ***** Dimension of covariate vector doesn't match length of n1 etc They are:", 
        length(n1), length(logitP.cov), dim(logitP.cov), 
        "\n")
    return()
  }
  if (any(m2 > n1, na.rm = TRUE)) {
    cat("***** ERROR ***** m2 must be <= n1. The arguments are \n n1:", 
        n1, "\n m2:", m2, "\n")
    return()
  }
  if (!all(bad.n1 %in% time, na.rm = TRUE)) {
    cat("***** ERROR ***** bad.n1 must be elements of strata identifiers. You entered \n bad.n1:", 
        bad.n1, "\n Strata identifiers are \n time:", 
        time, "\n")
    return()
  }
  if (!all(bad.m2 %in% time, na.rm = TRUE)) {
    cat("***** ERROR ***** bad.m2 must be elements of strata identifiers. You entered \n bad.m2:", 
        bad.m2, "\n Strata identifiers are \n time:", 
        time, "\n")
    return()
  }
  if (!all(bad.u2 %in% time, na.rm = TRUE)) {
    cat("***** ERROR ***** bad.u2 must be elements of strata identifiers. You entered \n bad.u2:", 
        bad.u2, "\n Strata identifiers are \n time:", 
        time, "\n")
    return()
  }
  if (!all(jump.after %in% time, na.rm = TRUE)) {
    cat("***** ERROR ***** jump.after must be elements of strata identifiers. You entered \n jump.after:", 
        jump.after, "\n Strata identifiers are \n time:", 
        time, "\n")
    return()
  }
  results.filename <- paste(prefix, "-results.txt", sep = "")
  sink(results.filename)
  cat(paste("Time Stratified Petersen with Diagonal recaptures and error in smoothed U - ", 
            date()))
  cat("\nVersion:", version, "\n\n")
  cat("\n\n", title, "Results \n\n")
  cat("*** Raw data *** \n")
  temp <- cbind(time, n1, m2, u2, round(sampfrac, digits = 2), 
                logitP.cov)
  colnames(temp) <- c("time", "n1", "m2", 
                      "u2", "SampFrac", paste("logitPcov[", 
                                              1:ncol(as.matrix(logitP.cov)), "]", sep = ""))
  print(temp)
  cat("\n\n")
  cat("Jump point are after strata: ", jump.after)
  if (length(jump.after) == 0) 
    cat("none - A single spline is fit")
  cat("\n\nValues of bad.n1 are : ", bad.n1, ". The value of n1 will be set to 1 and m2 to NA for these strata")
  if (length(bad.n1) == 0) 
    cat("none.")
  cat("\nValues of bad.m2 are : ", bad.m2, ". The value of m2 will be set to NA for these strata")
  if (length(bad.m2) == 0) 
    cat("none.")
  cat("\nValues of bad.u2 are : ", bad.u2, ". The value of u2 will be set to NA for these strata")
  if (length(bad.u2) == 0) 
    cat("none.")
  cat("\n\n*** Pooled Petersen Estimate based on pooling over ALL strata")
  cat("\nValues of u2 are adjusting for sampling fraction \n\n")
  cat("Total n1=", sum(n1, na.rm = TRUE), ";  m2=", 
      sum(m2, na.rm = TRUE), ";  u2=", sum(u2/sampfrac, 
                                           na.rm = TRUE), "\n\n")
  pp <- SimplePetersen(sum(n1, na.rm = TRUE), sum(m2, na.rm = TRUE), 
                       sum(u2/sampfrac, na.rm = TRUE))
  cat("Est U(total) ", format(round(pp$est), big.mark = ","), 
      "  (SE ", format(round(pp$se), big.mark = ","), 
      ")\n\n\n")
  temp.n1 <- n1
  temp.n1[match(bad.n1, time)] <- NA
  temp.m2 <- m2
  temp.m2[match(bad.m2, time)] <- NA
  temp.u2 <- u2
  temp.u2[match(bad.u2, time)] <- NA
  select <- (n1 > 0) & (!is.na(temp.n1)) & (!is.na(temp.m2)) & 
    (!is.na(temp.u2))
  cat("\n\n*** Pooled Petersen Estimate after EXCLUDING strata with missing value or flagged as bad.n1, bad.m2 or bad.m2. ")
  cat("\nValues of u2 are adjusted for sampling fraction\n\n")
  cat("The following strata are excluded because n1=0, NA values in m2 or u2, or flagged by bad.n1, bad.m2 or bad.u2:", 
      time[!select], "\n\n")
  temp.n1 <- n1[select]
  temp.m2 <- m2[select]
  temp.u2 <- u2[select]
  temp.sampfrac <- sampfrac[select]
  cat("Total n1=", sum(temp.n1), ";  m2=", sum(temp.m2), 
      ";  u2=", sum(temp.u2/temp.sampfrac), "\n\n")
  pp <- SimplePetersen(sum(temp.n1), sum(temp.m2), sum(temp.u2/temp.sampfrac))
  cat("Est U(total) ", format(round(pp$est), big.mark = ","), 
      "  (SE ", format(round(pp$se), big.mark = ","), 
      ")\n\n\n")
  cat("*** Stratified Petersen Estimator for each stratum PRIOR to removing strata with bad.n1, bad.m2, or bad.u2 values.")
  cat("\n    Values of u2 are adjusted for sampling fraction ***\n\n")
  temp.n1 <- n1
  temp.m2 <- m2
  temp.u2 <- u2/sampfrac
  sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
  temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$est), 
                round(sp$se))
  colnames(temp) <- c("time", "n1", "m2", 
                      "u2", "U[i]", "SE(U[i])")
  print(temp)
  cat("\n")
  cat("Est U(total) ", format(round(sum(sp$est, na.rm = TRUE)), 
                              big.mark = ","), "  (SE ", format(round(sqrt(sum(sp$se^2, 
                                                                               na.rm = TRUE))), big.mark = ","), ")\n\n\n")
  cat("*** Stratified Petersen Estimator for each stratum EXCLUDING strata with n1=0, NA values, or flagged by bad.n1, bad.m2, or bad.u2 values ***")
  cat("\n    Values of u2 are adjusted for sampling fraction ***\n\n")
  temp.n1 <- n1
  temp.n1[match(bad.n1, time)] <- NA
  temp.n1[match(bad.m2, time)] <- NA
  temp.n1[match(bad.u2, time)] <- NA
  temp.n1[temp.n1 == 0] <- NA
  temp.m2 <- m2
  temp.m2[match(bad.n1, time)] <- NA
  temp.m2[match(bad.m2, time)] <- NA
  temp.m2[match(bad.u2, time)] <- NA
  temp.u2 <- u2
  temp.u2[match(bad.n1, time)] <- NA
  temp.u2[match(bad.m2, time)] <- NA
  temp.u2[match(bad.u2, time)] <- NA
  temp.u2 <- temp.u2/sampfrac
  sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
  temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$est), 
                round(sp$se))
  colnames(temp) <- c("time", "n1", "m2", 
                      "u2", "U[i]", "SE(U[i])")
  print(temp)
  cat("\n")
  cat("Est U(total) ", format(round(sum(sp$est, na.rm = TRUE)), 
                              big.mark = ","), "  (SE ", format(round(sqrt(sum(sp$se^2, 
                                                                               na.rm = TRUE))), big.mark = ","), ")\n\n\n")
  cat("*** Test if pooled Petersen is allowable on strata without problems in n1 or m2. [Check if marked fractions are equal] ***\n\n")
  select <- (temp.n1 > 0) & (!is.na(temp.n1)) & (!is.na(temp.m2))
  temp.n1 <- n1[select]
  temp.m2 <- m2[select]
  test <- TestIfPool(temp.n1, temp.m2)
  cat("(Large Sample) Chi-square test statistic ", test$chi$statistic, 
      " has p-value", test$chi$p.value, "\n\n")
  temp <- cbind(time[select], test$chi$observed, round(test$chi$expected, 
                                                       1), round(test$chi$residuals^2, 1))
  colnames(temp) <- c("time", "n1-m2", "m2", 
                      "E[n1-m2]", "E[m2]", "X2[n1-m2]", "X2[m2]")
  print(temp)
  cat("\n Be cautious of using this test in cases of small expected values. \n\n")
  new.n1 <- rep(0, max(time) - min(time) + 1)
  new.m2 <- rep(0, max(time) - min(time) + 1)
  new.u2 <- rep(0, max(time) - min(time) + 1)
  new.sampfrac <- rep(0, max(time) - min(time) + 1)
  new.logitP.cov <- matrix(NA, nrow = max(time) - min(time) + 
                             1, ncol = ncol(as.matrix(logitP.cov)))
  new.time <- min(time):max(time)
  new.n1[time - min(new.time) + 1] <- n1
  new.m2[time - min(new.time) + 1] <- m2
  new.m2[match(bad.m2, new.time)] <- NA
  new.u2[time - min(new.time) + 1] <- u2
  new.u2[match(bad.u2, new.time)] <- NA
  new.sampfrac[time - min(new.time) + 1] <- sampfrac
  new.logitP.cov[time - min(new.time) + 1, ] <- as.matrix(logitP.cov)
  new.logitP.cov[is.na(new.logitP.cov[, 1]), 1] <- 1
  new.logitP.cov[is.na(new.logitP.cov)] <- 0
  new.m2[new.n1 == 0] <- NA
  new.n1[new.n1 == 0] <- 1
  new.n1[match(bad.n1, new.time)] <- 1
  new.m2[match(bad.n1, new.time)] <- NA
  new.u2 <- round(new.u2/new.sampfrac)
  new.u2[new.sampfrac < 0.001] <- NA
  jump.indicator <- rep("   ", max(time) - min(time) + 
                          1)
  jump.indicator[jump.after - min(time) + 1] <- "***"
  cat("\n\n*** Revised data *** \n")
  temp <- data.frame(time = new.time, n1 = new.n1, m2 = new.m2, 
                     u2 = new.u2, sampfrac = round(new.sampfrac, digits = 2), 
                     new.logitP.cov = new.logitP.cov, jump.indicator = jump.indicator)
  print(temp)
  cat("\n\n")
  cat("\n\n*** Information on priors *** \n")
  cat("   Parameters for prior on tauU (variance in spline coefficients: ", 
      tauU.alpha, tauU.beta, " which corresponds to a mean/std dev of 1/var of:", 
      round(tauU.alpha/tauU.beta, 2), round(sqrt(tauU.alpha/tauU.beta^2), 
                                            2), "\n")
  cat("   Parameters for prior on taueU (variance of log(U) about spline: ", 
      taueU.alpha, taueU.beta, " which corresponds to a mean/std dev of 1/var of:", 
      round(taueU.alpha/taueU.beta, 2), round(sqrt(taueU.alpha/taueU.beta^2), 
                                              2), "\n")
  cat("   Parameters for prior on beta.logitP[1] (intercept) (mean, 1/var):", 
      round(mu_xiP, 3), round(tau_xiP, 5), " which corresponds to a median P of ", 
      round(expit(mu_xiP), 3), "\n")
  cat("   Parameters for prior on tauP (residual variance of logit(P) after adjusting for covariates: ", 
      tauP.alpha, tauP.beta, " which corresponds to a mean/std dev of 1/var of:", 
      round(tauP.alpha/tauP.beta, 2), round(sqrt(tauP.alpha/tauP.beta^2), 
                                            2), "\n")
  cat("\n\n*** Initial seed for this run is: ", InitialSeed, 
      "\n")
  sink()
  if (debug2) {
    cat("\nprior to formal call to TimeStratPetersenDiagError\n")
    browser()
  }
  if (debug) {
    results <- TimeStratPetersenDiagError(title = title, 
                                          prefix = prefix, time = new.time, n1 = new.n1, m2 = new.m2, 
                                          u2 = new.u2, jump.after = jump.after - min(time) + 
                                            1, logitP.cov = new.logitP.cov, n.chains = 3, 
                                          n.iter = 10000, n.burnin = 5000, n.sims = 500, tauU.alpha = tauU.alpha, 
                                          tauU.beta = tauU.beta, taueU.alpha = taueU.alpha, 
                                          taueU.beta = taueU.beta, debug = debug, debug2 = debug2, 
                                          engine = engine, InitialSeed = InitialSeed)
  }
  else {
    results <- TimeStratPetersenDiagError(title = title, 
                                          prefix = prefix, time = new.time, n1 = new.n1, m2 = new.m2, 
                                          u2 = new.u2, jump.after = jump.after - min(time) + 
                                            1, logitP.cov = new.logitP.cov, n.chains = n.chains, 
                                          n.iter = n.iter, n.burnin = n.burnin, n.sims = n.sims, 
                                          tauU.alpha = tauU.alpha, tauU.beta = tauU.beta, taueU.alpha = taueU.alpha, 
                                          taueU.beta = taueU.beta, debug = debug, debug2 = debug2, 
                                          engine = engine, InitialSeed = InitialSeed)
  }
  plot_logU <- function(title, time, n1, m2, u2, results) {
    Nstrata <- length(n1)
    Uguess <- (u2 + 1) * (n1 + 2)/(m2 + 1)
    min_logU <- min(log(Uguess), na.rm = TRUE)
    max_logU <- max(log(Uguess), na.rm = TRUE)
    results.row.names <- rownames(results$summary)
    etaU.row.index <- grep("etaU", results.row.names)
    etaU <- results$summary[etaU.row.index, ]
    min_logU <- min(c(min_logU, etaU[, "mean"]), na.rm = TRUE)
    max_logU <- max(c(max_logU, etaU[, "mean"]), na.rm = TRUE)
    min_logU <- min(c(min_logU, etaU[, "2.5%"]), na.rm = TRUE)
    max_logU <- max(c(max_logU, etaU[, "2.5%"]), na.rm = TRUE)
    min_logU <- min(c(min_logU, etaU[, "97.5%"]), na.rm = TRUE)
    max_logU <- max(c(max_logU, etaU[, "97.5%"]), na.rm = TRUE)
    plot(time, log(Uguess), type = "p", main = paste(title, 
                                                     "\nFitted spline curve to raw U[i] with 95% credible intervals"), 
         sub = "Open/closed circles - initial and final estimates", 
         ylab = "log(U[i])", xlab = "Time Index", 
         ylim = c(min_logU, max_logU))
    points(time, etaU[, "mean"], type = "p", 
           pch = 19)
    lines(time, etaU[, "mean"])
    segments(time, etaU[, "2.5%"], time, etaU[, "97.5%"])
    logUne.row.index <- grep("logUne", results.row.names)
    logUne <- results$summary[logUne.row.index, "mean"]
    points(time, logUne, type = "p", pch = 20)
    lines(time, logUne, lty = 2)
  }
  pdf(file = paste(prefix, "-logU.pdf", sep = ""))
  plot_logU(title = title, time = new.time, n1 = new.n1, m2 = new.m2, 
            u2 = new.u2, results = results)
  dev.off()
  logitP.plot <- plot_logitP(title = title, time = new.time, 
                             n1 = new.n1, m2 = new.m2, u2 = new.u2, logitP.cov = new.logitP.cov, 
                             results = results)
  ggsave(plot = logitP.plot, filename = paste(prefix, "-logitP.pdf", 
                                              sep = ""), height = 6, width = 10, units = "in")
  results$plots$logitP.plot <- logitP.plot
  pdf(file = paste(prefix, "-Utot-acf.pdf", sep = ""))
  acf(results$sims.matrix[, "Utot"], main = paste(title, 
                                                  "\nAutocorrelation function for U total"))
  dev.off()
  pdf(file = paste(prefix, "-Ntot-posterior.pdf", sep = ""))
  plot(x = density(as.vector(results$sims.array[, , "Ntot"])), 
       main = paste(title, "\nPosterior density plot of N-total"), 
       sub = "Vertical lines mark 2.5th and 97.5th percentile")
  abline(v = results$summary["Ntot", c("2.5%", 
                                       "97.5%")])
  dev.off()
  pdf(file = paste(prefix, "-Utot-posterior.pdf", sep = ""))
  plot(x = density(as.vector(results$sims.array[, , "Utot"])), 
       main = paste(title, "\nPosterior density plot of U-total"), 
       sub = "Vertical lines mark 2.5th and 97.5th percentile")
  abline(v = results$summary["Utot", c("2.5%", 
                                       "97.5%")])
  dev.off()
  pdf(file = paste(prefix, "-GOF.pdf", sep = ""))
  discrep <- PredictivePosterior.TSPDE(new.n1, new.m2, new.u2, 
                                       expit(results$sims.list$logitP), round(results$sims.list$U))
  gof <- PredictivePosteriorPlot.TSPDE(discrep)
  dev.off()
  varnames <- names(results$sims.array[1, 1, ])
  pdf(file = paste(prefix, "-trace-logitP.pdf", sep = ""))
  parm.names <- varnames[grep("^logitP", varnames)]
  trace_plot(title = title, results = results, parms_to_plot = parm.names, 
             panels = c(3, 2))
  dev.off()
  pdf(file = paste(prefix, "-trace-logU.pdf", sep = ""))
  parm.names <- varnames[c(grep("Utot", varnames), grep("Ntot", 
                                                        varnames), grep("^etaU", varnames))]
  trace_plot(title = title, results = results, parms_to_plot = parm.names, 
             panels = c(3, 2))
  dev.off()
  sink(results.filename, append = TRUE)
  cat("\n\n*** Summary of MCMC results *** \n\n")
  print(results, digits.summary = 3)
  cat("\n\n*** Alternate DIC computation based on p_D = var(deviance)/2 \n")
  results.row.names <- rownames(results$summary)
  deviance.row.index <- grep("deviance", results.row.names)
  deviance <- results$summary[deviance.row.index, ]
  p.D <- deviance["sd"]^2/2
  dic <- deviance["mean"] + p.D
  cat("    D-bar: ", deviance["mean"], ";  var(dev): ", 
      deviance["sd"]^2, "; p.D: ", p.D, "; DIC: ", 
      dic)
  cat("\n\n\n\n*** Summary of Unmarked Population Size ***\n")
  temp <- results$summary[grep("Utot", rownames(results$summary)), 
                          ]
  old.Rhat <- temp["Rhat"]
  temp <- formatC(temp, big.mark = ",", format = "d")
  temp["Rhat"] <- formatC(old.Rhat, digits = 2, format = "f", 
                          flag = "#")
  print(temp, quote = FALSE)
  cat("\n\n*** Summary of Total Population Size *** \n")
  temp <- results$summary[grep("Ntot", rownames(results$summary)), 
                          ]
  old.Rhat <- temp["Rhat"]
  temp <- formatC(temp, big.mark = ",", format = "d")
  temp["Rhat"] <- formatC(old.Rhat, digits = 2, format = "f", 
                          flag = "#")
  print(temp, quote = FALSE)
  cat("\n\n\n\n*** Summary of Quantiles of Run Timing *** \n")
  cat("    This is based on the sample weeks provided and the U[i] values \n")
  q <- RunTime(time = time, U = results$sims.list$U, prob = run.prob)
  temp <- rbind(apply(q, 2, mean), apply(q, 2, sd))
  rownames(temp) <- c("Mean", "Sd")
  print(round(temp, 2))
  cat("\n\n")
  cat(paste("*** end of fit *** ", date()))
  sink()
  results$data <- list(time = time, n1 = n1, m2 = m2, u2 = u2, 
                       sampfrac = sampfrac, jump.after = jump.after, bad.n1 = bad.n1, 
                       bad.m2 = bad.m2, bad.u2 = bad.u2, logitP.cov = logitP.cov, 
                       version = version, date_run = date(), title = title)
  results$gof <- gof
  return(results)
}