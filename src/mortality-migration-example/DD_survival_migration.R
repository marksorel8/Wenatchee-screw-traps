library(deSolve)
options(device = windows)

#-------------------------------------------------------
# Logistic density-dependent survival
# Density-independent migration with periodic forcing
#-------------------------------------------------------

## @knitr survive_migrate
survive_migrate <- function(t, x, parms) 
{
  with(as.list(c(x, parms)), {
    r <- log(s0)
    K <- Rmax*(s0 - 1)/s0
    dRdt <- r*R*(1 - R/K) - dMdt(t, R, m, tau)
    dMdt <- dMdt(t, R, m, tau)
    
    return(list(c(dRdt,dMdt)))
  })
}

#-------------------------------------------------------
# Density-independent instantaneous migration rate
#-------------------------------------------------------

## @knitr dMdt
dMdt <- function(t, R, m, tau)
{
  R*m*t*exp(sin(2*pi*t/tau))
}
## @knitr

# Solution
## @knitr solve_ode
times <- seq(0, 1, length = 100)
parms <- c(s0 = 0.5, Rmax = 100, m = 2, tau = max(times)/3.75)
RMt <- ode(func = survive_migrate, parms = parms, y = c(R = 100, M = 0), times = times)
RMt <- data.frame(RMt)
## @knitr

#-------------------------------------------------------
# Depensatory instantaneous migration rate
#-------------------------------------------------------

## @knitr dMdt_dep
dMdt <- function(t, R, m, tau)
{
  R^2*m*t*exp(sin(2*pi*t/tau))
}
## @knitr

# Solution
## @knitr solve_ode_dep
times <- seq(0, 1, length = 100)
parms <- c(s0 = 0.5, Rmax = 100, m = 0.05, tau = max(times)/3.75)
RMt <- ode(func = survive_migrate, parms = parms, y = c(R = 100, M = 0), times = times)
RMt <- data.frame(RMt)
## @knitr

#-------------------------------------------------------
# Plot trajectories
#-------------------------------------------------------

dev.new()
## @knitr plot_solution
par(mfrow = c(2,1), mar = c(4.5, 4.5, 0.5, 1))

# total abundance
plot(R ~ time, data = RMt, type = "l", col = "salmon", lwd = 3, las = 1,
     cex.axis = 1.2, cex.lab = 1.5, ylim = range(RMt[,2:3]), xaxs = "i", yaxs = "i", 
     xlab = "", ylab = "Abundance")
lines(M ~ time, data = RMt, col = "darkblue", lwd = 3)
text(rep(max(RMt[,"time"]), 2), RMt[nrow(RMt),c("R","M")], c("Resident", "Migrant"), 
     cex = 1.5, adj = c(1.2,-0.5), col = c("salmon","darkblue"))

# migration rate
plot(c(NA, diff(M)) ~ time, data = RMt, type = "l", col = "darkblue", lwd = 3, 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, xaxs = "i", yaxs = "i", 
     xlab = "Time", ylab = "Migration rate")
abline(v = seq(0.75*parms["tau"], max(times) - parms["tau"], parms["tau"]), col = "lightgray")
## @knitr

#-------------------------------------------------------
# Now solve over range of R[t=0]
#-------------------------------------------------------

## @knitr solve_ode_R0
R0 <- seq(1, parms["Rmax"]*10, length = 100)
R0Mt <- do.call(rbind, lapply(R0, function(R) {
  out <- ode(func = survive_migrate, parms = parms, y = c(R = R, M = 0), 
             times = c(0, seq(0.75*parms["tau"], max(times), parms["tau"])))
  dplyr::mutate(data.frame(R0 = R, out), dM = c(NA, diff(M)), alrM = log(dM / tail(dM,1)))
}))
## @knitr

#-------------------------------------------------------
# Plot R[t=1] + M[t=1] vs R[t=0]
# If m = 0, this is the Beverton-Holt
#-------------------------------------------------------

dev.new()
## @knitr plot_DD_total
par(mar = c(4.5, 4.5, 0.5, 1))
plot(R + M ~ R0, data = R0Mt, subset = time == max(time), 
     type = "l", col = "salmon", lwd = 3, las = 1, cex.axis = 1.2, cex.lab = 1.5, 
     xlab = bquote(italic(R)[0]), ylab = bquote(italic(R)[italic(T)] + italic(M)[italic(T)]))
## @knitr

#-------------------------------------------------------
# Plot total migrants in each pulse vs R[t=0]
#-------------------------------------------------------

dev.new(width = 12, height = 3)
## @knitr plot_DD_Mt
tt <- seq(0.75*parms["tau"], max(times), parms["tau"])
par(mfrow = c(1,4), mar = c(1,3,2,1), oma = c(3,3,0,0))
for(i in 1:length(tt)) 
  plot(M ~ R0, data = R0Mt, subset = time == tt[i], type = "l", col = "darkblue", lwd = 3, 
       las = 1, cex.axis = 1.5, cex.main = 1.8, main = paste("stage", i))
mtext(bquote(italic(R)[0]), cex = 1.8*par("cex"), side = 1, line = 2, outer = TRUE)
mtext(bquote(Delta * italic(M)[italic(t)]), cex = 1.8*par("cex"), side = 2, line = 0.5, outer = TRUE)
## @knitr

#-------------------------------------------------------
# Plot log ratio of proportion of migrants in each pulse vs R[t=0]
# (reference class is final pulse)
#-------------------------------------------------------

dev.new(width = 9, height = 3)
## @knitr plot_DD_alrMt
tt <- seq(0.75*parms["tau"], max(times) - parms["tau"], parms["tau"])
par(mfrow = c(1,3), mar = c(1,3,2,1), oma = c(3,3,0,0))
for(i in 1:length(tt)) 
  plot(alrM ~ R0, data = R0Mt, subset = time == tt[i], type = "l", col = "darkblue", lwd = 3, 
       las = 1, cex.axis = 1.5, cex.main = 1.8, main = paste("stage", i))
mtext(bquote(italic(R)[0]), cex = 1.8*par("cex"), side = 1, line = 2, outer = TRUE)
mtext(bquote(alr(Delta * italic(M)[italic(t)])), cex = 1.8*par("cex"), side = 2, line = 0.5, outer = TRUE)
## @knitr


