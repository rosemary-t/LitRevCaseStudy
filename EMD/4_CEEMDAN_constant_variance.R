require(data.table)
require(ProbCast)
require(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

zones <- c(1:10)
horizons <- c(1:6) # forecast 1 to 6 hours ahead
split_date <- as.POSIXct("2012-12-31 12:00") # split the training/testing data here
eps <- 0.01
LL <- log(eps/(1-eps)) # lower and upper limits in transformed space.
UL <- log((1-eps)/eps)
qs <- seq(0.05, 0.95, 0.05) # quantiles we will produce forecasts for.

load("../VAR/all_sites_power.rda") # this is power in transformed space already
setnames(alldata, 'timestamp', 'target_time') # so this matches with training_mean_fcs colnames.


for (z in c(1:10)){
  print (z)
  
  ## load the site mean forecasts
  load(paste0("./mean_forecasts_zone",z,".rda"))
  
  ## want to get residuals for each zone.
  power_fcts <- alldata[target_time %in% all_mean_fcs$target_time] # keep only the timestamps we also have forecasts for
  
  
  ## find the residuals 
  zone_data <- merge(power_fcts[,.(target_time, power=get(paste0("zone",z,"power")))], 
                     all_mean_fcs[, .(issue_time, target_time, Horizon, meanfcs=mean_fc)])
  zone_data[,errors :=(power-meanfcs)]
  
  zone_train_data <- zone_data[issue_time <= split_date]
  sdevs <- zone_train_data[, .(sd = sd(errors)), by=Horizon]
  
  zone_test_data <- zone_data[issue_time > split_date,]
  zone_test_data <- merge(zone_test_data, sdevs, by="Horizon")
  zone_test_data[, ActualPower := 1/(1+exp(-power))]
  
  ## now make quantile forecasts (in transformed space first) for this zone
  z_quantiles <- data.table()
  for (q in qs){
    z_quantiles[, paste0('q',q*100) := qnorm(p=q, mean=zone_test_data$meanfcs, sd=zone_test_data$sd)]
  }
  
  class(z_quantiles) <- c("MultiQR", class(z_quantiles)) # make a multiQR object so we can sort the quantiles.
  z_quantiles <- SortQuantiles(z_quantiles, Limits=list(U=UL, L=LL)) # make sure quantiles don't cross and crop at eps (in transformed space)
  z_quantiles <- 1/(1+exp(-z_quantiles)) ## now transform back to 'power' space. However, this operation returns a data.frame object so:
  class(z_quantiles) <- c("MultiQR", "data.table", class(z_quantiles)) # make a multiQR object again (SortQuantiles just returns a data.table).
  other_info <- zone_test_data[, .(issue_time, target_time, Horizon, ActualPower)] # all the other information needed to evaluate the quantiles.
  
  ## make a list with the MultiQR object and the time info data.table
  EMD_fcs <-list(z_quantiles,other_info) 
  save(EMD_fcs, file=paste0("./zone",z,"_qs_constvar.rda"))
}




