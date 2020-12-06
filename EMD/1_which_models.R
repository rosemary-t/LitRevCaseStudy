require(data.table)
require(doParallel)
require(rstudioapi)
require(forecast)
require(frbs)
require(Rlibeemd)
require(tseries)
## also loads forecast Rlibeemd and frbs packages within make_forecasts function.

setwd(dirname(getActiveDocumentContext()$path))



### Code to run several different models on each IMF (and the residuals from creating different numbers of IMFs)
### so we can evaluate error and see which is the best model for each decomposed time series.

########################### function to fit and predict from ANFIS model ######################

anfis_predict <- function(window_imf, nlags, hor, eps){
  ## window_imf is a single time series object.
  ## hor= horizon in number of steps ahead
  
  lags <- c(1:nlags) # specify which lags to use as ANFIS model inputs.
  
  ## get the lags needed for the actual forecast
  testinds <- length(window_imf) - rev(lags) + 1 # reverse lags so these indices are in the same order as the columns in the training data.
  testdata <- window_imf[testinds] # these are the inputs to generate a forecast for this rolling window.
  
  ## now make the training data
  shifts <- lags + hor - 1
  traindata <- as.data.table(window_imf) # data.table we can add the lags to.
  
  ## now create the lags
  for (s in shifts){
    traindata[, paste0('lag',s) := shift(x,n=s, fill=NA, type="lag")]
  }
  
  ## need to reorder the columns so x (the forecast output) is the last column.
  colnames <- paste0('lag',shifts) # a list with the names of the lag columns (i.e. the forecast inputs)
  colnames <- append(colnames, 'x') # add the forecast output column to the END of this list.
  setcolorder(traindata,colnames)
  
  traindata <- traindata[complete.cases(traindata),] # remove rows with any NA
  
  ## finally, train ANFIS model
  anfis <- frbs.learn(traindata, method.type="ANFIS")
  
  ans <- predict(anfis, testdata)
  return (ans)
}

############################################################################################

for(zone in 5:10){ # which area we are forecasting for
  split_date <- as.POSIXct("2012-12-31 12:00", tz="UTC") # split the training/testing data here
  horizon <- 1 # number of steps ahead we are forecasting for
  eps <- 0.01 # cutoff point at boundaries
  wl <- 1000 # window length to use for this test
  max_nimfs <- 10
  min_nimfs <- 3 #(this is 2 IMF series plus a residual)
  
  # now load the data and prepare it
  load(paste0('../Data/data_zone', zone,'.rda'))
  assign('zdata', get(paste0('data_zone', zone)))
  rm(list=paste0("data_zone",zone))
  
  
  zdata[, TARGETVAR := approx(.I, TARGETVAR, .I)$y] ## interpolate missing values in  TARGETVAR
  zdata <- zdata[!is.na(zdata$TARGETVAR), ] ## if the last row in the column is NA, needs to be removed as can't be interpolated
  
  ## lognormal transform: 
  zdata[, T_targetvar := lapply(.SD, function(x){ifelse(x <= eps,  eps, ifelse(x >= (1-eps), (1-eps), x))}), .SDcols=c("TARGETVAR")] # clip to [eps, (1-eps)] range
  zdata[, T_targetvar := log(T_targetvar/(1-T_targetvar))]
  
  
  ## split into training and testing data (at same point as for other models)
  z_train <- copy(zdata[timestamp <= split_date,])
  # z_test <- copy(zdata[timestamp > split_date,])
  
  start_ts <- z_train[wl, timestamp]
  end_ts <- z_train[(dim(z_train)[1] - horizon), timestamp]
  poss_window_times <- z_train[(timestamp >= start_ts) & (timestamp <= end_ts), timestamp] # the timestamps for the end of all possible sliding windows in the training data.
  issue_times <- sample(poss_window_times, 100, replace=FALSE) # forecast issue times for 100 random windows.
  save(issue_times, file=paste0("./zone_",zone,"_issuetimes.rda")) # keep a record of which training windows have been used.
  
  ## define a function to return the forecast for every IMF, for every model.
  make_forecasts <- function(IT, zonedata, h, min_nimfs, max_nimfs, WL){
    require(Rlibeemd)
    require(data.table)
    require(forecast)
    require(frbs)
    ## IT = issue time of this window, h=horizon, WL=window length.
    models <- c("arma", "anfis", "arma1")
    window_forecasts <- data.table(issue_time=IT, model=models) #  data.table to put forecasts in
    
    window_start <- IT - (WL-1)*60*60
    windowdata <- zonedata[(timestamp >= window_start) & (timestamp<= IT),]
    imfs <- ceemdan(windowdata$T_targetvar,num_imfs=max_nimfs) # create the imfs for this rolling window.
    
    for (imfi in c(1:max_nimfs)){
      fit <- auto.arima(imfs[,imfi], max.p=6, max.q=4, max.order=8, stepwise = FALSE) # this is fitting ARIMA for this imf
      window_forecasts[model=="arma", colnames(imfs)[imfi] := as.numeric(forecast(fit, h=horizon)[["mean"]])] # save arma forecast
      
      ## and for the arma(d=1)
      fit1 <- auto.arima(imfs[,imfi],d=1, max.p=6, max.q=4, max.order=8, stepwise = FALSE)
      window_forecasts[model=="arma1", colnames(imfs)[imfi] := as.numeric(forecast(fit1, h=horizon)[["mean"]])]
      
      ## and finally the ANFIS model.
      ans <- anfis_predict(imfs[,imfi], nlags=2, hor=horizon, eps)
      window_forecasts[model=="anfis", colnames(imfs)[imfi] := ans]
    }
    
    setnames(window_forecasts, 'Residual', paste('Residual', max_nimfs))
    
    range_nimfs <- c(min_nimfs: (max_nimfs-1)) # the different numbers of IMFs we need the residual for (we're already got the residual for nimfs=10)
    
    for (nimfs in range_nimfs){
      imfs <- ceemdan(windowdata$T_targetvar,num_imfs=nimfs)
      fit <- auto.arima(imfs[,"Residual"], max.p=6, max.q=4, max.order=8, stepwise = FALSE)
      window_forecasts[model=="arma", paste("Residual", nimfs) := as.numeric(forecast(fit, h=horizon)[["mean"]])]
      
      ## and for the arma(d=1) and arma(d=2) too
      fit1 <- auto.arima(imfs[,"Residual"],d=1, max.p=6, max.q=4, max.order=8, stepwise = FALSE)
      window_forecasts[model=="arma1", paste("Residual", nimfs) := as.numeric(forecast(fit1, h=horizon)[["mean"]])]
      
      ## and finally the ANFIS model.
      ans <- anfis_predict(imfs[,"Residual"], nlags=2, hor=horizon, eps)
      window_forecasts[model=="anfis", paste("Residual", nimfs) := ans]
      
    }
    return (window_forecasts)
  }
  
  ## now run the forecasts for each window, in parallel
  cores <- detectCores()
  cl <- makeCluster(cores-2)
  registerDoParallel(cl)
  start_time <- Sys.time()
  model_forecasts <- foreach(it = issue_times) %dopar% {
    make_forecasts(it, z_train, horizon, min_nimfs, max_nimfs, wl)
  }
  print (Sys.time() - start_time)
  stopCluster(cl)
  
  model_forecasts <- rbindlist(model_forecasts)
  
  
  ## now convert forecasts into errors, for model evaluation
  ## first need the 'actual values' to forecast:
  rw_actualvals <- data.table('issue_time'=issue_times)
  range_nimfs <- c(min_nimfs: (max_nimfs-1)) 
  
  for (rwt in issue_times){
    window_start <- rwt - (wl-1)*60*60
    window_end <- rwt + horizon*60*60 # include the data from forecast issue time -> forecast target time as well
    windowdata <- z_train[(timestamp >= window_start) & (timestamp<= window_end),]
    imfs <- ceemdan(windowdata$T_targetvar,num_imfs=max_nimfs) # create the imfs for this rolling window.
    actual_imfs <- tail(imfs, 1) # the last value corresponds to the forecast target time
    rw_actualvals[issue_time == rwt, colnames(actual_imfs) := as.data.table(actual_imfs)]
    
    ## also need the residual series for different numbers of imfs
    for (nimfs in range_nimfs){
      imfs <- ceemdan(windowdata$T_targetvar,num_imfs=nimfs) # create the imfs for this rolling window and this nimfs
      actual_imfs <- tail(imfs, 1)[,'Residual']
      rw_actualvals[issue_time == rwt, paste("Residual", nimfs) := as.numeric(actual_imfs)]
    }
  }
  setnames(rw_actualvals,'Residual', 'Residual 10')
  series_names <- names(rw_actualvals)[2:18] # names of all the IMF/Residual columns we want to get errors for
  
  ## now make forecast errors..
  allmodel_errors <- data.table(issue_time = model_forecasts$issue_time, model=model_forecasts$model)
  for (sn in series_names){
    series <- merge(model_forecasts[, .(issue_time, model, modelimf = get(sn))], rw_actualvals[, .(issue_time, actualimf = get(sn))], by='issue_time')
    series[, paste0(sn) := abs(modelimf-actualimf)]
    seriescols <- c('issue_time', 'model', paste0(sn))
    allmodel_errors <- merge(allmodel_errors, series[,..seriescols], by=c('issue_time', 'model'))
  }
  
  models_mae <- allmodel_errors[, lapply(.SD, mean), .SDcols=series_names, by=model]
  models_mae[,type := "MAE"]
  
  
  ## now do bootstrap resampling to get a standard error for the error values too, since it's only on a small sample.
  require(tseries)
  model_names <- models_mae$model
  
  bs_ses <- data.table(model=model_names) # data.table to put bootstrapped standard errors into. 
  
  for (mn in model_names){
    for (sn in series_names){
      bs_error <- tsbootstrap(x=allmodel_errors[model==mn,get(sn)], nb=200, mean)
      bs_ses[model==mn, paste0(sn) := bs_error$se]
    }
  }
  
  bs_ses[, type:= 'SE']
  all_error_info <- rbind(models_mae, bs_ses)
  fwrite(all_error_info, file=paste0("./whichmodels_errorinfo_zone",zone,".csv"))
  
  
  ## now plot errors for the different imf:
  require(ggplot2)
  
  long_plot_data <- melt(all_error_info, id.vars=c("model", "type"), measure.vars = series_names)
  mae_data <- long_plot_data[type=="MAE", .(model, series=variable, mae=value)]
  se_data <- long_plot_data[type=="SE", .(model, series=variable, se=value)]
  plot_dt <- merge(mae_data,se_data, by=c("model", "series"))
  
  pdf(file=paste0("./which_models_plot_zone", zone, ".pdf"), width=20, height=14, family="Times")
  ggplot(data=plot_dt, aes(x=series, y=mae, colour=model))+
    geom_point() +
    geom_errorbar(aes(ymin=mae-se, ymax=mae+se), width=.2)
  dev.off()
  
}
## conclusions:
## anfis are always worse than the best performing arma model. (zone 1 at least)

