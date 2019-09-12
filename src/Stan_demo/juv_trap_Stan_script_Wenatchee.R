if(.Platform$OS.type == "windows") options(device=windows)

library(dataRetrieval)
library(matrixStats)
library(lubridate)
library(rstan)
library(loo)
library(shinystan)
library(yarrr)
library(corrplot)
library(here)

if(file.exists(here("src","Stan_demo","juv_trap_fit.RData"))) 
  load(here("src","Stan_demo","juv_trap_fit.RData"))
if(file.exists(here("src","Stan_demo","juv_trap_my_fit.RData"))) 
  load(here("src","Stan_demo","juv_trap_my_fit.RData"))

# Convenience functions to simplify extracting posterior means
# or single parameters from stanfit objects
stan_mean <- function(object, pars)
{
  mm <- get_posterior_mean(object, pars)
  return(mm[,ncol(mm)])
}

extract1 <- function(object, par)
{
  extract(object, par)[[1]]
}

#===========================================================================
# DATA
#===========================================================================

source(here("src","Load Screw Trap Data.R"))

screw_trap_dat<-load_dat()
screw_trap_dat$chiw$chiw_catch<-droplevels(screw_trap_dat$chiw$chiw_catch[!is.na(screw_trap_dat$chiw$chiw_catch$allSubs),])



min_DOY<-min(as.numeric(screw_trap_dat$chiw$chiw_catch$DOY))
screw_trap_dat$chiw$chiw_catch$DOY<-as.numeric(screw_trap_dat$chiw$chiw_catch$DOY)-min_DOY+1

# add week 
screw_trap_dat$chiw$chiw_catch$week<-ceiling(as.numeric(screw_trap_dat$chiw$chiw_catch$DOY)/7)
#min_week<-min(screw_trap_dat$chiw$chiw_catch$week)
#screw_trap_dat$chiw$chiw_catch$week<-screw_trap_dat$chiw$chiw_catch$week
#add elapsed time
screw_trap_dat$chiw$chiw_catch$elapsed_time<-1

trap_catch_all <-screw_trap_dat$chiw$chiw_catch[,c("year","DOY","week","allSubs","elapsed_time")]
trap_catch_all$year<-as.integer(trap_catch_all$year)
trap_catch_all$DOY<-as.integer(trap_catch_all$DOY)


#Mark Recapture data
screw_trap_dat$chiw$chiw_effic$Date2<-as.Date(c(as.Date(screw_trap_dat$chiw$chiw_effic$Date[1:88],format="%d-%b-%y"),as.Date(screw_trap_dat$chiw$chiw_effic$Date[89:143],format="%m/%d/%y"),rep(NA, times=length(144:329))))
screw_trap_dat$chiw$chiw_effic$year<-as.integer(format(screw_trap_dat$chiw$chiw_effic$Date2,"%Y"))#add year
screw_trap_dat$chiw$chiw_effic$week<-ceiling((as.numeric(format(screw_trap_dat$chiw$chiw_effic$Date2,"%j"))-min_DOY+1)/7)#add week




rel<-reshape::melt( tapply(screw_trap_dat$chiw$chiw_effic$rel, screw_trap_dat$chiw$chiw_effic[,c("year","week")], sum))
rec<-reshape::melt( tapply(screw_trap_dat$chiw$chiw_effic$rec, screw_trap_dat$chiw$chiw_effic[,c("year","week")], sum))

MR_all <-cbind(rel[,1:3],rec[,3])[!is.na(rec[,3]),]
colnames(MR_all)[3:4]<-c("mark","recap")

# Download USGS daily mean discharge data
Chiwawa<-"12456500"
flow_all <- renameNWISColumns(readNWISdv(siteNumbers = Chiwawa, parameterCd = "00060", statCd = "00003",
                                         startDate = paste(min(trap_catch_all$year) , "01", "01", sep = "-"),
                                         endDate = paste(max(trap_catch_all$year), "12", "31", sep = "-")))
names(flow_all)[3:4] <- c("date","flow")
write.csv(flow_all, here("src","Stan_demo","flow_all.csv"), row.names = FALSE)


# Align flow data to trapping dates and merge with trap and MR data
dates <- data.frame(year = rep(unique(trap_catch_all$year), 
                                     tapply(trap_catch_all$DOY, trap_catch_all$year, max)),
                    day = unlist(sapply(tapply(trap_catch_all$DOY, trap_catch_all$year, max),
                                        function(x) 1:x)))
row.names(dates) <- NULL

dates$date <- as.Date(paste0(dates$year , "-02-22"), format = "%Y-%m-%d") + dates$day-1
dates$week <- ceiling( dates$day/7)

flow_daily <- merge(dates, flow_all[,c("date","flow")], all.x = TRUE)

flow_weekly <- aggregate(flow ~ week + year, data = flow_daily, mean)



# Stan data for single-year model
year <- 2011
trap_catch <- na.omit(trap_catch_all[trap_catch_all$year==year,])
trap_catch <- trap_catch[trap_catch$elapsed_time != 0,]
MR <- na.omit(MR_all[MR_all$year==year,])
X_p <- cbind(1, scale(flow_weekly$flow[flow_weekly$year==year]))
X_M <- matrix(1,max(trap_catch$DOY))

stan_dat <- list(N_MR = nrow(MR),
                 MR_week = MR$week,
                 mark = MR$mark,
                 recap = MR$recap,
                 NX_p = ncol(X_p),
                 X_p = X_p,
                 N_trap = nrow(trap_catch),
                 trap_day = trap_catch$DOY,
                 trap_week = trap_catch$week,
                 NX_M = ncol(X_M),
                 X_M = X_M,
                 C = trap_catch$allSubs,
                 elapsed_time = round(trap_catch$elapsed_time),
                 Use_NB=1)

# Stan data for multi-year model
trap_catch_my <- na.omit(trap_catch_all[is.element(trap_catch_all$year, MR_all$year),])
trap_catch_my <- trap_catch_my[trap_catch_my$elapsed_time != 0,]
trap_catch_my$year_f <- as.numeric(factor(trap_catch_my$year))
MR_my <- na.omit(MR_all[is.element(MR_all$year, trap_catch_all$year),])
MR_my$year_f <- as.numeric(factor(MR_my$year))
X_p_my <- matrix(1,max(trap_catch_my$year_f)*max(trap_catch_my$week),1)
X_M_my <- matrix(1,max(trap_catch_my$year_f)*max(trap_catch_my$DOY),1)

stan_dat_my <- list(N_MR = nrow(MR_my),
                    MR_year = MR_my$year_f,
                    MR_week = MR_my$week,
                    mark = MR_my$mark,
                    recap = MR_my$recap,
                    NX_p = ncol(X_p_my),
                    X_p = X_p_my,
                    N_trap = nrow(trap_catch_my),
                    trap_year = trap_catch_my$brood_year_f,
                    trap_week = trap_catch_my$week,
                    trap_day = trap_catch_my$day,
                    NX_M = ncol(X_M_my),
                    X_M = X_M_my,
                    C = trap_catch_my$catch,
                    elapsed_time = round(trap_catch_my$elapsed_time))


#===========================================================================
# CALL STAN TO FIT MODELS
#===========================================================================

#---------------------------
# Single year
#---------------------------

# Function to generate initial values
stan_init <- function(data, chains)
{
  with(data,
       {
         return(lapply(1:chains, function(i)
           list(beta_p = array(rnorm(NX_p, c(qlogis(0.2), rep(0, NX_p - 1)), 0.5), dim = NX_p),
                sigma_p = runif(1,0.1,2),
                beta_M = array(rnorm(NX_M, c(log(20), rep(0, NX_M - 1)), 0.5), dim = NX_M),
                phi_M = runif(1,0,0.9),
                sigma_M = runif(1,0.5,2),
                phi_obs=array(runif(stan_dat$Use_NB,.2,2000)))))
       })
}

# Call Stan to fit model
juv_trap_fit2 <- stan(file = here::here("src","Stan_demo","juv_trap.stan"),
                     data = stan_dat, 
                     init = stan_init(stan_dat,3), 
                     pars = c("beta_M","phi_M","sigma_M",
                              "beta_p","sigma_p","p",
                              "M_hat","M","M_tot","C_hat","phi_obs"),
                     chains = 3, iter = 1500, warmup = 500, thin = 1, cores = 3,
                     control = list(adapt_delta = 0.99, max_treedepth = 13))



# Print fitted model
print(juv_trap_fit2, pars = c("phi_M","sigma_M",
                             "beta_p","sigma_p"
                             ,"M_tot","phi_obs"), include =T, probs = c(0.05,0.5,0.95))

# Check it out in Shinystan
launch_shinystan(juv_trap_fit)

# Save stanfit
save(juv_trap_fit, file = here("src","Stan_demo","juv_trap_fit.RData"))

# Effective sample size per time
runtime <- mean(rowSums(get_elapsed_time(juv_trap_fit)))
n_eff <- summary(juv_trap_fit)$summary[,"n_eff"]
min(n_eff)/runtime
mean(n_eff)/runtime

#---------------------------
# Multi-year
#--------------------------

# Function to generate initial values
stan_init_my <- function(data, chains)
{
  with(data,
       {
         N_year <- max(trap_year)
         N_week <- max(trap_week)
         
         return(lapply(1:chains, function(i)
           list(beta_p = array(rnorm(NX_p, c(qlogis(0.2), rep(0, NX_p - 1)), 0.5), dim = NX_p),
                sigma_p_year = runif(1,0.1,2),
                sigma_p_week = runif(1,0.1,2),
                # beta_M = array(rnorm(NX_M, c(log(20), rep(0, NX_M - 1)), 0.5), dim = NX_M),
                mu_M = array(rnorm(N_year, log(20), 0.5), dim = N_year),
                phi_M = runif(N_year,0,0.9),
                sigma_M = runif(N_year,0.5,2))))
       })
}

# Call Stan to fit model
juv_trap_my_fit <- stan(file = here("src","Stan_demo","juv_trap_multiyear.stan"),
                        data = stan_dat_my, 
                        init = stan_init_my(stan_dat_my,3), 
                        pars = c("mu_M","phi_M","sigma_M","Q_M",
                                 "beta_p","sigma_p_year","sigma_p_week","p",
                                 "M_hat","M","M_tot","C_hat"),
                        chains = 3, iter = 1500, warmup = 500, thin = 1, cores = 3,
                        control = list(adapt_delta = 0.99, max_treedepth = 13))

# Print fitted model
print(juv_trap_my_fit, pars = c("M_hat","M","p","C_hat","Q_M"), include = F, probs = c(0.05,0.5,0.95))

# Check it out in Shinystan
launch_shinystan(juv_trap_my_fit)

# Save stanfit
save(juv_trap_my_fit, file = here("src","Stan_demo","juv_trap_my_fit.RData"))


#===========================================================================
# FIGURES
#===========================================================================

#---------------------------
# Single year
#---------------------------

# Time series of obs and predicted catch, capture probability, 
# and predicted true outmigrants
C_hat <- extract1(juv_trap_fit,"C_hat")
p <- extract1(juv_trap_fit,"p")
M <- extract1(juv_trap_fit,"M")
p_fall <- rowSums(M[,1:125])/rowSums(M)
M_tot <- extract1(juv_trap_fit,"M_tot")

dev.new(width = 7, height = 14)
par(mfrow = c(3,1), mar = c(4.5,4.5,1,1), oma = c(0,0,3,0))

with(stan_dat, {
  c1 <- transparent("blue", 0.7)
  plot(trap_day, colMeans(C_hat), type = "l", lwd = 2, col = "blue",
       ylim = c(0, max(C, colQuantiles(C_hat, probs = 0.975))),
       xlab = "Day", ylab = "Catch", las = 1, cex.lab = 1.5, cex.axis = 1.2)
  polygon(c(trap_day, rev(trap_day)),
          c(colQuantiles(C_hat, probs = 0.025), rev(colQuantiles(C_hat, probs = 0.975))),
          col = c1, border = NA)
  points(trap_day, C, pch = 1)
  mtext(side = 3, line = 1, paste("Chiwawa year", year))
  
  plot(1:max(trap_week), colMeans(p), type = "l", lwd = 2, col = "blue",
       ylim = c(0, 1), xlab = "Week", ylab = "Capture probability", 
       las = 1, cex.lab = 1.5, cex.axis = 1.2)
  polygon(c(1:max(trap_week), max(trap_week):1),
          c(colQuantiles(p, probs = 0.025), rev(colQuantiles(p, probs = 0.975))),
          col = c1, border = NA)
  points(MR_week, recap/mark, pch = 1)
  
  plot(1:max(trap_day), colMeans(M), type = "l", lwd = 2, col = "blue",
       ylim = c(0, max(colQuantiles(M, probs = 0.975))),
       xlab = "Day", ylab = "Outmigrants", las = 1, cex.lab = 1.5, cex.axis = 1.2)
  polygon(c(1:max(trap_day), max(trap_day):1),
          c(colQuantiles(M, probs = 0.025), rev(colQuantiles(M, probs = 0.975))),
          col = c1, border = NA)
  legend("topleft", 
         paste("total = ", round(mean(M_tot), 1), 
               " (", round(quantile(M_tot, 0.025), 1), ", ",
               round(quantile(M_tot, 0.975), 1), ")"#,
               #"proportion fall = ", round(quantile(p_fall, 0.975), 2), 
              # " (", round(quantile(p_fall, 0.025), 2), ", ",
               #round(quantile(p_fall, 0.975), 2), ")", sep = ""
              ), 
         bty = "n")
})


#---------------------------
# Multi-year
#---------------------------

# Time series of obs and predicted catch
C_hat <- extract1(juv_trap_my_fit_2,"C_hat")

dev.new(width = 10, height = 14)
par(mfcol = c(5,2), mar = c(3,4.5,1,0.5), oma = c(1,0,1,0))

with(stan_dat_my, {
  c1 <- transparent("blue", 0.3)
  for(i in 1:max(trap_year))
  {
    plot(trap_day[trap_year==i], C[trap_year==i], pch = 1,
         xlim = c(1, max(trap_day)), ylim = c(0, max(apply(C_hat, 2, quantile, 0.975))),
         xlab = "", ylab = "", las = 1, cex.lab = 1.2, cex.axis = 1, xaxt = "n")
    axis(1, at = axTicks(1), labels = format(axTicks(1) + as.Date("09-30", format = "%m-%d"), "%m/%d"), 
         cex.axis = 1)
    lines(trap_day[trap_year==i], colMeans(C_hat[,trap_year==i]), lwd = 1, col = "blue")
    polygon(c(trap_day[trap_year==i], rev(trap_day[trap_year==i])),
            c(apply(C_hat[,trap_year==i], 2, quantile, 0.025), 
              rev(apply(C_hat[,trap_year==i], 2, quantile, 0.975))),
            col = c1, border = NA)
    if(par("mfg")[1]==par("mfg")[3]) mtext(side = 1, line = 2.5, "Date") 
    if(par("mfg")[2]==1) mtext(side = 2, line = 3, "Catch")
    mtext(side = 3, line = 0.1, sort(unique(trap_catch_my$year))[i])
  }
})

# Annual time series of predicted true outmigrants
M <- extract1(juv_trap_my_fit_2,"M")
M_tot <- extract1(juv_trap_my_fit_2,"M_tot")

dev.new(width = 10, height = 14)
par(mfcol = c(5,2), mar = c(3,5,1,0.5), oma = c(1,0,1,0))

with(stan_dat_my, {
  c1 <- transparent("blue", 0.3)
  for(i in 1:max(trap_year))
  {
    plot(1:max(trap_day), colMeans(M[,i,]), type = "l", lwd = 2, col = "blue",
         ylim = c(0, max(apply(M, 2:3, quantile, 0.975))),
         xlab = "", ylab = "", las = 1, cex.lab = 1.2, cex.axis = 1, xaxt = "n")
    polygon(c(1:max(trap_day), max(trap_day):1),
            c(colQuantiles(M[,i,], probs = 0.025), rev(colQuantiles(M[,i,], probs = 0.975))),
            col = c1, border = NA)
    axis(1, at = axTicks(1), labels = format(axTicks(1) + as.Date("09-30", format = "%m-%d"), "%m/%d"), 
         cex.axis = 1)
    if(par("mfg")[1]==par("mfg")[3]) mtext(side = 1, line = 2.5, "Date") 
    if(par("mfg")[2]==1) mtext(side = 2, line = 3.5, "Outmigrants")
    mtext(side = 3, line = 0.1, sort(unique(trap_catch_my$year))[i])
    p_fall <- rowSums(M[,i,1:125])/rowSums(M[,i,])
    legend("topleft", 
           paste("total = ", round(mean(M_tot[,i]), 1), 
                 " (", round(quantile(M_tot[,i], 0.025), 1), ", ",
                 round(quantile(M_tot[,i], 0.975), 1), ") \n",
                 "proportion fall = ", round(quantile(p_fall, 0.975), 2), 
                 " (", round(quantile(p_fall, 0.025), 2), ", ",
                 round(quantile(p_fall, 0.975), 2), ")", sep = ""), 
           bty = "n")
  }
})


# Annual time series of standardized process errors in true outmigrants
log_M_hat_z <- matrix(stan_mean(juv_trap_my_fit_2,"log_M_hat_z"),
                      nrow = max(stan_dat_my$trap_year), ncol = max(stan_dat_my$trap_day),
                      byrow = TRUE)

dev.new(width = 10, height = 14)
pdf("Stan_chiwawa_subyearling_process_error.pdf" )
par(mfcol = c(7,2), mar = c(3,2,1,0.5), oma = c(1,1.5,1,0))

with(stan_dat_my, {
  c1 <- transparent("blue", 0.3)
  for(i in 1:max(trap_year))
  {
    plot(1:max(trap_day), log_M_hat_z[i,], type = "l", col = "blue",
         xlab = "", ylab = "", las = 1, cex.lab = 1.2, cex.axis = 1, xaxt = "n")
    axis(1, at = axTicks(1), labels = format(axTicks(1) + as.Date("09-30", format = "%m-%d"), "%m/%d"), 
         cex.axis = 1)
    if(par("mfg")[1]==par("mfg")[3]) mtext(side = 1, line = 2.5, "Date") 
    # if(par("mfg")[2]==1) mtext(side = 2, line = 3.5, "Process error")
    mtext(side = 3, line = 0.1, sort(unique(trap_catch_my$year))[i])
  }
  mtext("Process error", side = 2, outer = TRUE)
})

dev.off()
# Time series of total predicted true outmigrants and proportion fall by brood year
M <- extract1(juv_trap_my_fit_2,"M")/1000
M_tot <- extract1(juv_trap_my_fit_2,"M_tot")/1000
p_fall <- apply(M[,,1:125], 1:2, sum)/M_tot
y <- sort(unique(trap_catch_my$brood_year))
c1 <- transparent("blue", 0.7)

dev.new()
par(mfcol = c(2,1), mar = c(3,5,1,0.5), oma = c(2,0,0,0))

plot(y, colMeans(M_tot), type = "l", lwd = 2, col = "blue",
     ylim = c(0, max(colQuantiles(M_tot, probs = 0.975))),
     xlab = "", ylab = "Outmigrants (x 1000)", 
     las = 1, cex.lab = 1.2, cex.axis = 1)
rug(y[y %% 2 != 0], side = 1, ticksize = -0.04)
polygon(c(y, rev(y)),c(colQuantiles(M_tot, probs = 0.025), rev(colQuantiles(M_tot, probs = 0.975))),
        col = c1, border = NA)

plot(y, colMeans(p_fall), type = "l", lwd = 2, col = "blue",
     ylim = c(0,max(colQuantiles(p_fall, probs = 0.975))), 
     xlab = "", ylab = "Proportion fall", las = 1, cex.lab = 1.2, cex.axis = 1)
rug(y[y %% 2 != 0], side = 1, ticksize = -0.04)
polygon(c(y, rev(y)), c(colQuantiles(p_fall, probs = 0.025), rev(colQuantiles(p_fall, probs = 0.975))),
        col = c1, border = NA)
mtext(side = 1, line = 0, "Brood year", cex = 1.2, outer = TRUE)

# Interannual correlations among daily process error innovations
dev.new()
with(stan_dat_my, {
  Q_M <- matrix(stan_mean(juv_trap_my_fit,"Q_M"), max(trap_year), max(trap_year))
  dimnames(Q_M) <- list(sort(unique(trap_catch_my$brood_year)),
                        sort(unique(trap_catch_my$brood_year)))
  corrplot(Q_M, method = "ellipse", diag = FALSE)
})




