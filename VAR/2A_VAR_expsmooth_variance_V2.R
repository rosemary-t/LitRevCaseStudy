require(data.table)
require(forecast)
require(rstudioapi)
require(doParallel)
require(ProbCast)

setwd(dirname(getActiveDocumentContext()$path))

zones <- c(1:10)
horizons <- c(1:6) # forecast 1 to 6 hours ahead
eps <- 0.01
LL <- log(eps/(1-eps)) # lower and upper limits in transformed space.
UL <- log((1-eps)/eps)
qs <- seq(0.05, 0.95, 0.05) # quantiles we will produce forecasts for.

################ function to get exponential smoothing forecasts of variance ################
horizon_var_fcs <- function(zonedata, H, splitdate){
  require(data.table)
  require(forecast)
  
  zone_h_ts <- zonedata[Horizon==H]
  zone_h_ts <- zone_h_ts[order(issue_time)]
  zone_h_ts[, errors := approx(.I, errors, .I)$y] ## interpolate missing values in errors
  horizon_fcs <- copy(zone_h_ts[issue_time>splitdate, .(issue_time,target_time)])
  test_times <- horizon_fcs$target_time
  
  for (it in test_times){
    # fit exponential smoothing to all data where we have forecast and actual power before it (ie before target_time<=it)
    data <-zone_h_ts[target_time <= it,] 
    modelfit <- ses(ts(zone_h_ts[target_time <= it, errors]), h=1, inital='optimal')
    horizon_fcs[(issue_time==it), "varfc" := modelfit$mean]
  }
  return(horizon_fcs)
}

#############################################################################################

load("./mean_forecasts_trainset.rda")
load("./mean_forecasts_testset.rda")
split_date <- max(training_mean_fcs$issue_time)
all_mean_fcs <- rbind(training_mean_fcs, testing_mean_fcs)
load("./all_sites_power.rda") # this is power in transformed space already
setnames(alldata, 'timestamp', 'target_time') # so this matches with training_mean_fcs colnames.

## want to get residuals for each zone.
power_fcts <- alldata[target_time %in% all_mean_fcs$target_time] # keep only the timestamps we also have forecasts for



for (z in zones){
  print (z)
  ## find the variance of the residuals from the training set (in transformed space)
  zone_data <- merge(power_fcts[,.(target_time, power=get(paste0("zone",z,"power")))], 
                           all_mean_fcs[, .(issue_time, target_time, Horizon, meanfcs=get(paste0("zone",z,"power_meanfcs")))])
  zone_data[,errors :=(power-meanfcs)]
  
  cores <- detectCores()
  cl <- makeCluster(cores-2)
  registerDoParallel(cl)
  start_time <- Sys.time()
  variance_fcs <- foreach(H = horizons) %dopar% {
    horizon_var_fcs(zone_data, H, split_date)
  }
  print (Sys.time() - start_time)
  stopCluster(cl)
  
  variance_fcs <- rbindlist(variance_fcs)
  meanvardt <- merge(zone_data, variance_fcs, by=c("issue_time", "target_time")) # means and variances in one data.table.
  meanvardt[, varfc := abs(varfc)] # standard deviations must be positive
  meanvardt[, ActualPower := 1/(1+exp(-power))]
  meanvardt <- na.omit(meanvardt) ## delete rows with NAs
  
  ## now we have a mean and a variance for every time point, make the quantile forecasts.
  z_quantiles <- data.table()
  for (q in qs){
    z_quantiles[, paste0('q',q*100) := qnorm(p=q, mean=meanvardt$meanfcs, sd=meanvardt$varfc)]
  }
  
  class(z_quantiles) <- c("MultiQR", class(z_quantiles)) # make a multiQR object so we can sort the quantiles.
  z_quantiles <- SortQuantiles(z_quantiles, Limits=list(U=UL, L=LL)) # make sure quantiles don't cross and crop at eps (in transformed space)
  z_quantiles <- 1/(1+exp(-z_quantiles)) ## now transform back to 'power' space. However, this operation returns a data.frame object so:
  class(z_quantiles) <- c("MultiQR", "data.table", class(z_quantiles)) # make a multiQR object again (SortQuantiles just returns a data.table).
  other_info <- meanvardt[, .(issue_time, target_time, Horizon, ActualPower)] # all the other information needed to evaluate the quantiles.
  
  ## make a list with the MultiQR object and the time info data.table
  VAR_fcs <-list(z_quantiles,other_info) 
  save(VAR_fcs, file=paste0("./zone",z,"_qs_expsmoothvar.rda"))
}
