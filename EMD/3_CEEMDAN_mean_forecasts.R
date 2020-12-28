require(Rlibeemd)
require(data.table)
require(forecast)
require(doParallel)
require(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

z <- 2 # which area we are forecasting for
horizons <- c(1:6) # number of steps ahead we are forecasting for
split_date <- as.POSIXct("2012-12-31 12:00") # split the training/testing data here
eps <- 0.01 # cutoff point at boundaries

opt_params <- fread(file="./opt_parameters.csv")
opt_nimfs <- opt_params[Zone==z, Opt_nimf]
opt_windowlen <- opt_params[Zone==z, Opt_wl]

# now load the data and prepare it
load(paste0("../Data/data_zone", z, ".rda"))
assign('zdata', get(paste0('data_zone', z)))
rm(list=paste0("data_zone",z))

zdata[, TARGETVAR := approx(.I, TARGETVAR, .I)$y] ## interpolate missing values in  TARGETVAR
zdata <- zdata[!is.na(zdata$TARGETVAR), ] ## if the last row in the column is NA, needs to be removed as can't be interpolated

## lognormal transform: 
zdata[, T_targetvar := lapply(.SD, function(x){ifelse(x <= eps,  eps, ifelse(x >= (1-eps), (1-eps), x))}), .SDcols=c("TARGETVAR")] # clip to [eps, (1-eps)] range
zdata[, T_targetvar := log(T_targetvar/(1-T_targetvar))]

## also load the table specifying what model for each IMF:
opt_model <- fread(file=paste0("./whichmodels_errorinfo_zone",z,".csv"))
imfnames <- c(paste("IMF", c(1:(opt_nimfs-1))), paste("Residual", opt_nimfs))

## need to produce mean forecasts now using opt_nimfs and opt_windowlen, for ALL possible time points 
## then the ones in the training set can be used to find the variance of the residuals
## and the test set is used to evaluate the final probabilistic forecasts out of sample.
start_ts <- zdata[opt_windowlen, timestamp] # first possible issue time 
end_ts <- zdata[(dim(zdata)[1] - max(horizons)), timestamp]
allposs_window_times <- zdata[(timestamp >= start_ts) & (timestamp <= end_ts), timestamp] # all forecast issue times


## first put everything in a function that makes the forecasts for all horizons 1:H, for a single window time. 
produce_forecasts <- function(ts_datatable, WT, WL, nimfs, H, opt_model){
  require(Rlibeemd)
  require(data.table)
  require(forecast)
  
  ## set up data.table to put forecasts into
  horizon_forecasts <- data.table(issue_time=WT, Horizon=c(1:H), target_time=WT+c(1:H)*60*60)
  imfnames <- c(paste("IMF", c(1:(nimfs-1))), paste("Residual", nimfs))
  
  window_start <- WT - (WL-1)*60*60
  windowdata <- ts_datatable[(timestamp >= window_start) & (timestamp <= WT), T_targetvar]
 
  imfs <- ceemdan(windowdata,num_imfs=nimfs) # create the imfs for this rolling window.
  colnames(imfs) <- imfnames
  imflist <- matrix(nrow=nimfs, ncol=H) # empty list to put each IMF forecast in, for each horizon.
  
  for (i in c(1:opt_nimfs)){
    imfname <- imfnames[i]
    opt_arma <- opt_model[imf_name == imfname, `best_model (value for d)`]
    
    if (opt_arma == 0){
      # fit arma with d==0 (in reality seems to be faster to not specify d)
      fit <- auto.arima(imfs[,imfname], max.p=6, max.q=4, max.order=8, stepwise = FALSE)
    }else if(opt_arma==1){
      # fit arma with d=1.
      fit <- auto.arima(imfs[,imfname], d=1, max.p=6, max.q=4, max.order=8, stepwise = FALSE)
    }else{cat("wt=",wt,"," ,imfname,  'not selected a model to fit!')}
    
    imflist[i,] <-  as.numeric(forecast(fit, h=H)[["mean"]]) # save the forecast
  }
  
  horizon_forecasts[, mean_fcs := colMeans(imflist)]
  return (horizon_forecasts)
}

# test <- produce_forecasts(zdata, allposs_window_times[1], opt_windowlen, opt_nimfs, max(horizons), opt_model)

## first do for 500 windows in the training set, to get standard deviation of residuals from
## find the window times used already....
load(paste0("./zone_", z, "_issuetimesdt.rda"))
poss_train_times <- allposs_window_times[(!allposs_window_times %in% issue_timesdt$timestamps) & (allposs_window_times <= split_date)] 
training_times <- sample(poss_train_times, size=500, replace=FALSE)

cores <- detectCores()
cl <- makeCluster(cores-2)
registerDoParallel(cl)
start_time <- Sys.time()
train_horizon_forecasts <- foreach(wt =  training_times) %dopar% {
  produce_forecasts(zdata, wt, opt_windowlen, opt_nimfs, max(horizons), opt_model)
}
print (Sys.time() - start_time)
stopCluster(cl)

Train_forecasts <- rbindlist(train_horizon_forecasts)
save(Train_forecasts, file=paste0("./train_mean_forecasts_zone",z,".rda"))

## now get mean forecasts for all windows in the test set
test_times <- allposs_window_times[allposs_window_times > split_date]
second_test_times <- test_times[!test_times %in% unique(Test_forecasts$issue_time)]

cl <- makeCluster(cores-2)
registerDoParallel(cl)
start_time <- Sys.time()
test_horizon_forecasts <- foreach(wt =  second_test_times) %dopar% {
  produce_forecasts(zdata, wt, opt_windowlen, opt_nimfs, max(horizons), opt_model)
}
print (Sys.time() - start_time)
stopCluster(cl)

Test1_forecasts <- rbindlist(test_horizon_forecasts)
save(Test_forecasts, file=paste0("./test_mean_forecasts_zone",z,".rda"))
load(paste0("./test_mean_forecasts_zone",z,".rda"))
#Test_forecasts <- rbind(Test_forecasts, Test1_forecasts)
