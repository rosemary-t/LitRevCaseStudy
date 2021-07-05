require(data.table)
require(ProbCast)
require(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

zones <- c(1:10)
horizons <- c(1:6) # forecast 1 to 6 hours ahead
eps <- 0.01
LL <- log(eps/(1-eps)) # lower and upper limits in transformed space.
UL <- log((1-eps)/eps)
qs <- seq(0.05, 0.95, 0.05) # quantiles we will produce forecasts for.


load("./mean_forecasts_trainset.rda") # find variance from the mean forecasts made over the training set.
load("./mean_forecasts_testset.rda") # these will be used to generate the final quantile forecasts for the test set.
# load("./VAR_fold_settings.rda")
load("./all_sites_power.rda") # this is power in transformed space already.....
setnames(alldata, 'timestamp', 'target_time') # so this matches with training_mean_fcs colnames.
power_fcts <- alldata[target_time %in% training_mean_fcs$target_time] # the power realisations for the same time points as the training data mean forecasts



for (z in zones){
  print (z)
  ## find the variance of the residuals from the training set
  zone_train_data <- merge(power_fcts[,.(target_time, power=get(paste0("zone",z,"power")))], 
                          training_mean_fcs[, .(issue_time, target_time, Horizon, meanfcs=get(paste0("zone",z,"power_meanfcs")))])
  zone_train_data[,errors :=(power-meanfcs)]
  
  ## if want to look at ACF/PACF of residuals
  # hist(zone_train_data[Horizon==6, errors],breaks=50)
  # plot(zone_train_data[Horizon==6, target_time], zone_train_data[Horizon==6, errors], type='l')
  # require(forecast)
  # acf(zone_train_data[Horizon==6, errors])

  
  zone_sd <- sd(zone_train_data$errors) # standard deviation of the residuals
  
  ## mean forecasts in 'power' (untransformed) space, needed to generate quantile forecasts
  zone_test_data <- testing_mean_fcs[, .(issue_time, target_time, Horizon, mean_T_fcs=get(paste0("zone",z,"power_meanfcs")))]
  
  # load the original ActualPower (unaffected by clipping, transforming)
  load(paste0('../Data/data_zone', z,'.rda'))
  assign("zdata", get(paste0("data_zone", z)))
  rm(list=paste0("data_zone", z))
  zdata <- zdata[,.(TARGETVAR, timestamp)]
  setnames(zdata, c('TARGETVAR', 'timestamp'), c('ActualPower', 'target_time'))
  zone_test_data <- merge(zone_test_data, zdata, by='target_time')
  
  ## now make quantile forecasts (in transformed space first) for this zone
  z_quantiles <- data.table()
  for (q in qs){
    z_quantiles[, paste0('q',q*100) := qnorm(p=q, mean=zone_test_data$mean_T_fcs, sd=zone_sd)]
  }
  
  class(z_quantiles) <- c("MultiQR", class(z_quantiles)) # make a multiQR object so we can sort the quantiles.
  z_quantiles <- SortQuantiles(z_quantiles, Limits=list(U=UL, L=LL)) # make sure quantiles don't cross and crop at eps (in transformed space)
  z_quantiles <- 1/(1+exp(-z_quantiles)) ## now transform back to 'power' space. However, this operation returns a data.frame object so:
  class(z_quantiles) <- c("MultiQR", "data.table", class(z_quantiles)) # make a multiQR object again (SortQuantiles just returns a data.table).
  
  other_info <- zone_test_data[, .(issue_time, target_time, Horizon, ActualPower)] # all the other information needed to evaluate the quantiles.

  ## make a list with the MultiQR object and the time info data.table
  VAR_fcs <-list(z_quantiles,other_info) 
  save(VAR_fcs, file=paste0("./zone",z,"_qs_constvar.rda"))
  
}




