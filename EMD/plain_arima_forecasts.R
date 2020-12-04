require(data.table)
require(forecast)
require(doParallel)
require(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

z <- 1 # which area we are forecasting for
horizons <- c(1:6) # number of steps ahead we are forecasting for
split_date <- as.POSIXct("2012-12-31 12:00") # split the training/testing data here
eps <- 0.01 # cutoff point at boundaries

opt_params <- fread(file=paste0("./opt_parameters.csv"))
opt_windowlen <- opt_params[zone==z, opt_wl]

# now load the data and prepare it
load(paste0("../../Data/GEFcom2014data/data_zone",z,".rda"))
assign('zdata', get(paste0('data_zone', z)))
rm(list=paste0("data_zone",z))

zdata[, TARGETVAR := approx(.I, TARGETVAR, .I)$y] ## interpolate missing values in  TARGETVAR
zdata <- zdata[!is.na(zdata$TARGETVAR), ] ## if the last row in the column is NA, needs to be removed as can't be interpolated

## lognormal transform: 
zdata[, T_targetvar := lapply(.SD, function(x){ifelse(x <= eps,  eps, ifelse(x >= (1-eps), (1-eps), x))}), .SDcols=c("TARGETVAR")] # clip to [eps, (1-eps)] range
zdata[, T_targetvar := log(T_targetvar/(1-T_targetvar))]

## need to produce mean forecasts now for sliding windows using opt_windowlen, for ALL possible time points 
## then the ones in the training set can be used to find the variance of the residuals
## and the test set is used to evaluate the forecasts out of sample.
# start_ts <- zdata[opt_windowlen, timestamp] # first possible issue time 
# end_ts <- zdata[(dim(zdata)[1] - max(horizons)), timestamp]
# allposs_window_times <- zdata[(timestamp >= start_ts) & (timestamp <= end_ts), timestamp] # all forecast issue times

load(paste0("./train_mean_forecasts_zone",z,".rda"))
train_times <- unique(Train_forecasts$issue_time)


## function to return forecasts
arima_forecasts <- function(ts_datatable, WT, WL, H, ar_d){
  require(data.table)
  require(forecast)
  
  horizon_forecasts <- data.table(issue_time=WT, target_time=WT+c(1:H)*60*60)
    
  window_start <- WT - (WL-1)*60*60
  windowdata <- ts_datatable[(timestamp >= window_start) & (timestamp <= WT),]
  if (ar_d==0){
    fit <- auto.arima(ts_datatable$T_targetvar, max.p=6, max.q=4, max.order=8, stepwise=FALSE)
  }else{
    fit <- auto.arima(ts_datatable$T_targetvar, d=ar_d, max.p=6, max.q=4, max.order=8, stepwise=FALSE)
  }
  
  horizon_forecasts[, mean_ar_fc := as.numeric(forecast(fit, h=H)[["mean"]])]
  return (horizon_forecasts)
}


## forecasts for arima with d=0
cl <- makeCluster(6)
registerDoParallel(cl)
start_time <- Sys.time()
horizon_forecasts <- foreach(wt=train_times) %dopar% {
  arima_forecasts(zdata, wt, opt_windowlen, max(horizons), 0)
}
print (Sys.time() - start_time)
stopCluster(cl)

arima0_forecasts <- rbindlist(horizon_forecasts)
save(arima0_forecasts, file=paste0("./arima0_mean_forecasts_zone",z,".rda"))


## and for d=1
cl <- makeCluster(6)
registerDoParallel(cl)
start_time <- Sys.time()
horizon_forecasts <- foreach(wt=train_times) %dopar% {
  arima_forecasts(zdata, wt, opt_windowlen, max(horizons), 1)
}
print (Sys.time() - start_time)
stopCluster(cl)

arima1_forecasts <- rbindlist(horizon_forecasts)
save(arima1_forecasts, file=paste0("./arima1_mean_forecasts_zone",z,".rda"))

  
