require(data.table)
require(ProbCast)
require(ggplot2)
require(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

zone <- 1


## load persistence model forecasts 
load(file=paste0("../persistence/zone", zone, "_forecasts.rda"))
horizons <- unique(persistence_fcs[[2]]$Horizon)

## load VAR forecasts
load(file=paste0("../VAR/zone",zone,"_qs_constvar_with24.rda"))
VAR_fcs_with24 <- copy(VAR_fcs)
rm(VAR_fcs)
load(file=paste0("../VAR/zone",zone,"_qs_constvar.rda"))

## load Markov Chain forecasts 
load(file=paste0("../MC/test_forecasts_zone",zone,".rda"))

## check that none of the forecasts are NA
# na_persistence <- which(!complete.cases(persistence_fcs[[1]]))
# na_var <- which(!complete.cases((VAR_fcs[[1]])))
# na_mc <- which(!complete.cases((MC_fcs[[1]])))

# get the skill scores
model_ss <- CJ(Horizon = horizons, Model=c("VAR", "VAR_24"))
persistence_rel <- data.table()
var_rel <- data.table()
var_24_rel <- data.table()
mc_rel <- data.table()
emd_rel <- data.table()


for (h in horizons){
  persistence_times <- persistence_fcs[[2]][Horizon==h,issue_time]# issue_times that there is a persistence forecast for
  VAR_times <- VAR_fcs[[2]][Horizon==h, issue_time]
  VAR_24_times <- VAR_fcs_with24[[2]][Horizon==h, issue_time]
  ## now can get the issue_times that have both persistence and VAR forecasts:
  common_times <- as.POSIXct(Reduce(intersect, list(persistence_times, VAR_times, VAR_24_times)), tz="UTC", origin="1970-01-01 00:00.00 UTC")
  persistence_indices <- persistence_fcs[[2]][(Horizon==h) & (issue_time %in% common_times), which=TRUE]
  VAR_indices <- VAR_fcs[[2]][(Horizon==h) & (issue_time %in% common_times), which=TRUE]
  VAR_24_indices <- VAR_fcs_with24[[2]][(Horizon==h) & (issue_time %in% common_times), which=TRUE]
  
  ## first, MAE and RMSE with q50
  p_errors <- persistence_fcs[[1]][persistence_indices, q50] - persistence_fcs[[2]][persistence_indices, ActualPower]
  p_MAE <- mean(abs(p_errors))
  p_RMSE <- sqrt(mean(p_errors^2))
  
  v_errors <- VAR_fcs[[1]][VAR_indices, q50] - VAR_fcs[[2]][VAR_indices, ActualPower]
  v_MAE <- mean(abs(v_errors))
  v_RMSE <- sqrt(mean(v_errors^2))
  model_ss[(Model=="VAR")&(Horizon==h), MAE_ss := (1-v_MAE/p_MAE)]
  model_ss[(Model=="VAR")&(Horizon==h), RMSE_ss := (1-v_RMSE/p_RMSE)]
  
  VAR_24_errors <- VAR_fcs_with24[[1]][VAR_24_indices, q50] - VAR_fcs_with24[[2]][VAR_24_indices, ActualPower]
  v24_MAE <- mean(abs(VAR_24_errors))
  v24_RMSE <- sqrt(mean(VAR_24_errors^2))
  model_ss[(Model=="VAR_24")&(Horizon==h), MAE_ss := (1-v24_MAE/p_MAE)]
  model_ss[(Model=="VAR_24")&(Horizon==h), RMSE_ss := (1-v24_RMSE/p_RMSE)]
  
  
  
  ## and also pinball skill score.
  p_pinball <- pinball(persistence_fcs[[1]][persistence_indices], realisations = persistence_fcs[[2]][persistence_indices, ActualPower])
  v_pinball <- pinball(VAR_fcs[[1]][VAR_indices], realisations = VAR_fcs[[2]][VAR_indices, ActualPower])
  v24_pinball <- pinball(VAR_fcs_with24[[1]][VAR_24_indices], realisations = VAR_fcs_with24[[2]][VAR_24_indices, ActualPower])
  # pinball_ss <- 1-v_pinball$Loss/p_pinball$Loss # if you uwant pinball skill score per quantile.
  # plot(v_pinball$Quantile, pinball_ss)
  
  model_ss[(Model=="VAR") & (Horizon==h), Pinball_ss := 1 - mean(v_pinball$Loss)/mean(p_pinball$Loss)]
  model_ss[(Model=="VAR_24") & (Horizon==h), Pinball_ss := 1 - mean(v24_pinball$Loss)/mean(p_pinball$Loss)]
  
  
  p_h_rel <- reliability(persistence_fcs[[1]][persistence_indices], realisations = persistence_fcs[[2]][persistence_indices, ActualPower])
  p_h_rel$Horizon <- h
  persistence_rel <- rbind(persistence_rel, p_h_rel)
  
  var_h_rel <- reliability(VAR_fcs[[1]][VAR_indices], realisations = VAR_fcs[[2]][VAR_indices, ActualPower])
  var_h_rel$Horizon <- h
  var_rel <- rbind(var_rel, var_h_rel)
  
  v24_h_rel <- reliability(VAR_fcs_with24[[1]][VAR_24_indices], realisations = VAR_fcs_with24[[2]][VAR_24_indices, ActualPower])
  v24_h_rel$Horizon <- h
  var_24_rel <- rbind(var_24_rel, v24_h_rel)
  
}

ggplot(data=model_ss, aes(x=Horizon, y=MAE_ss, colour=Model)) +
  geom_line()

ggplot(data=model_ss, aes(x=Horizon, y=RMSE_ss, colour=Model)) +
  geom_line()

ggplot(data=model_ss, aes(x=Horizon, y=Pinball_ss, colour=Model)) +
  geom_line()


ggplot(data=var_rel, aes(x=Nominal, y=Empirical, colour=Horizon)) +
  geom_line()

h <- 6
plot(var_rel[Horizon==h, Nominal], var_rel[Horizon==h,Empirical], type='l', col='red')
lines(var_24_rel[Horizon==h, Nominal], var_24_rel[Horizon==h, Empirical], type='l', col='blue')
