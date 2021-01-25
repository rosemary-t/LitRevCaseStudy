require(data.table)
require(ProbCast)
require(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

zone <- 4
horizon <- 2
models <- c("Persistence", "VAR", "MC", "EMD")
Nboots <- 500 # number of times to bootstrap resample


## load all the different model forecasts for this zone
load(file=paste0("../persistence/zone", zone, "_forecasts.rda"))
load(file=paste0("../VAR/final_quantiles_zone",zone,".rda"))

## how about just the constant variance model?? 
# load(file=paste0("../VAR/zone",zone,"_qs_constvar.rda"))
# load(file=paste0("../VAR/zone",zone,"_qs_expsmoothvar.rda"))

load(file=paste0("../MC/test_forecasts_zone",zone,".rda"))
load(file=paste0("../EMD/final_quantiles_zone",zone,".rda"))

times1 <- merge(persistence_fcs[[2]][Horizon==horizon, .(issue_time, target_time)], VAR_fcs[[2]][Horizon==horizon, .(issue_time, target_time)], by=c('issue_time', 'target_time'))
times2 <- merge(MC_fcs[[2]][Horizon==horizon, .(issue_time, target_time)], EMD_fcs[[2]][Horizon==horizon, .(issue_time, target_time)])
mutual_times <- merge(times1, times2)

get_reliability <- function(forecasts, horizon, timestamps, nb, model){
  ## timestamps are the target times we want to evaluate for
  rowinds <- forecasts[[2]][(Horizon==horizon) & (target_time %in% timestamps), which=TRUE]
  rel <- reliability(qrdata=forecasts[[1]][rowinds,], realisations = forecasts[[2]][rowinds, ActualPower], bootstrap=nb, plot.it=T)
  rel <- as.data.table(rel)
  rel[, flat_empirical := Empirical-Nominal]
  rel[, flat_upper := upper-Nominal]
  rel[, flat_lower := lower-Nominal]
  rel[, Model := model]
  
  return(rel)
}

p_rel <- get_reliability(persistence_fcs, horizon, mutual_times$target_time, Nboots, "Persistence")
var_rel <- get_reliability(VAR_fcs, horizon, mutual_times$target_time, Nboots, "VAR")
mc_rel <- get_reliability(MC_fcs, horizon, mutual_times$target_time, Nboots, "MC")
emd_rel <- get_reliability(EMD_fcs, horizon, mutual_times$target_time, Nboots, "EMD")
reliability_dt <- rbind(p_rel, var_rel, mc_rel, emd_rel)
fwrite(reliability_dt, file=paste0("./zone",zone, "_h",horizon,"_reliability.csv"))

