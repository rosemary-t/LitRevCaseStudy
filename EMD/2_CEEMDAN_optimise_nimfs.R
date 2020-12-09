require(Rlibeemd)
require(data.table)
require(forecast)
require(rstudioapi)
require(doParallel)
require(tseries)


setwd(dirname(getActiveDocumentContext()$path))

### Code to run the best model on each IMF (and the residuals from creating different numbers of IMFs)
### so we can evaluate error and find the optimal number of IMFs to use.
### and then the same for the optimal window length too.

for(zone in 5:10){ # which area we are forecasting for
  horizon <- 1 # number of steps ahead we are forecasting for
  split_date <- as.POSIXct("2012-12-31 12:00") # split the training/testing data here
  eps <- 0.01 # cutoff point at boundaries
  i_wl <- 1000 # initial window length 
  max_nimfs <- 10
  min_nimfs <- 3 #(this is 2 IMF series plus a residual)
  try_windowlens <- seq(1000,5000,500) # the set of window lengths to test over
  
  # now load the data and prepare it
  load(paste0("../Data/data_zone", zone, ".rda"))
  assign('zdata', get(paste0('data_zone', zone)))
  rm(list=paste0("data_zone",zone))
  
  
  zdata[, TARGETVAR := approx(.I, TARGETVAR, .I)$y] ## interpolate missing values in  TARGETVAR
  zdata <- zdata[!is.na(zdata$TARGETVAR), ] ## if the last row in the column is NA, needs to be removed as can't be interpolated
  
  ## lognormal transform: 
  zdata[, T_targetvar := lapply(.SD, function(x){ifelse(x <= eps,  eps, ifelse(x >= (1-eps), (1-eps), x))}), .SDcols=c("TARGETVAR")] # clip to [eps, (1-eps)] range
  zdata[, T_targetvar := log(T_targetvar/(1-T_targetvar))]
  
  
  ## split into validation and testing data
  z_train <- copy(zdata[timestamp<= split_date,])
  
  ## select more random windows, over which to find the optimal number of IMFs
  load(paste0("./zone_",zone, "_issuetimesdt.rda")) # get the issue times of the random windows already used to find the optimal model
  
  start_ts <- z_train[i_wl, timestamp]
  end_ts <- z_train[(dim(z_train)[1] - horizon), timestamp]
  allposs_window_times <- z_train[(timestamp >= start_ts) & (timestamp <= end_ts), timestamp] # all forecast issue times
  poss_window_times <- allposs_window_times[!allposs_window_times %in% issue_timesdt$timestamps] # minus the windows already used
  issue_times1 <- sample(poss_window_times, size=200, replace=FALSE) # pick 200 random windows to optimise nimfs over
  issue_times1dt <- data.table(timestamps=issue_times1, stage=2) # add issue_times1 to the list of used window times.
  
  
  ## load the table specifying what model for each IMF:
  opt_model <- fread(file=paste0("./whichmodels_errorinfo_zone",zone,".csv"))
  
  
  ## define function to generate forecasts, using the best model, for every IMF.
  series_forecasts <- function(IT, zonedata, hor, min_nimfs, max_nimfs, WL, opt_models){
    require(Rlibeemd)
    require(data.table)
    require(forecast)
    
    best_forecasts <- data.table(issue_time=IT, target_time=IT+hor*60*60) # data.table to store IMF forecasts
    range_nimfs <- c(min_nimfs: (max_nimfs-1))
    
    window_start <- IT - (WL-1)*60*60
    windowdata <- zonedata[(timestamp >= window_start) & (timestamp <= IT),]
    
    ## first forecast IMFs 1-9 and Residual 10
    imfs <- ceemdan(windowdata$T_targetvar,num_imfs=max_nimfs) # create the imfs for this rolling window.
    imfnames <- colnames(imfs)
    imfnames[imfnames=="Residual"] <- "Residual 10" # rename the residual as "Residual 10" to avoid confusion
    colnames(imfs) <- imfnames
    
    for (imfname in colnames(imfs)){
      opt_arma <- opt_model[model == "opt_model", get(imfname)]
      
      if (opt_arma == "arma"){
        # fit arma with d==0 (in reality seems to be faster to not specify d)
        fit <- auto.arima(imfs[,imfname], max.p=6, max.q=4, max.order=8, stepwise = FALSE)
      }else if(opt_arma=="arma1"){
        # fit arma with d=1.
        fit <- auto.arima(imfs[,imfname], d=1, max.p=6, max.q=4, max.order=8, stepwise = FALSE)
      }else{cat("IT=",IT,"," ,imfname,  'not selected a model to fit!')}
      
      best_forecasts[issue_time== IT, paste0(imfname) := as.numeric(forecast(fit, h=horizon)[["mean"]])] # save the forecast
    }
    
    for (nimfs in range_nimfs){
      imfs <- ceemdan(windowdata$T_targetvar,num_imfs=nimfs)
      imfname <- paste("Residual", nimfs)
      opt_arma <- opt_model[model == "opt_model", get(imfname)]
      
      if (opt_arma == "arma"){
        # fit arma with d==0 (in reality seems to be faster to not specify d)
        fit <- auto.arima(imfs[,"Residual"], max.p=6, max.q=4, max.order=8, stepwise = FALSE)
      }else if(opt_arma=="arma1"){
        # fit arma with d=1.
        fit <- auto.arima(imfs[,"Residual"], d=1, max.p=6, max.q=4, max.order=8, stepwise = FALSE)
      }else{cat("wt=",wt,"," ,imfname,  'not selected a model to fit!')}
      
      best_forecasts[issue_time== IT, paste0(imfname) := as.numeric(forecast(fit, h=horizon)[["mean"]])] # save the forecast
    }
    
    return(best_forecasts)
  }
  
  # test <- series_forecasts(issue_times1[2], z_train, horizon, min_nimfs, max_nimfs, i_wl, opt_model)
  
  ## now run the forecasts for each window
  cores <- detectCores()
  cl <- makeCluster(cores-2)
  registerDoParallel(cl)
  start_time <- Sys.time()
  model_forecasts <- foreach(it = issue_times1) %dopar% {
    series_forecasts(it, z_train, horizon, min_nimfs, max_nimfs, i_wl, opt_model)
  }
  print (Sys.time() - start_time)
  stopCluster(cl)
  
  best_forecasts <- rbindlist(model_forecasts)
  
  save(best_forecasts, file=paste0("./best_forecasts_zone",zone,".rda"))
  # load(paste0("./best_forecasts_zone",zone,".rda"))
  # issue_times1 <- best_forecasts$issue_time
  
  ## now need to make the 'final' summed forecasts for different numbers of IMFs.
  for (nimfs in c(min_nimfs:max_nimfs)){
    # first get a list of best_forecasts columns that must be summed for an overall forecast with nimfs (e.g. nimfs=3 means IMF 1+IMF 2+Residual 3)
    imf_names <- c(paste("IMF", c(1:(nimfs-1))), paste("Residual", nimfs))
    best_forecasts[, paste0("nimfs=",nimfs) := rowSums(.SD), .SDcols=imf_names]
  }
  
  final_fcs_cols <- paste0("nimfs=", c(min_nimfs:max_nimfs)) # names of the columns in best_forecasts containing the summed forecasts
  best_forecasts[, (final_fcs_cols) := lapply(.SD, function(x){1/(1+exp(-x))}), .SDcols = final_fcs_cols] # convert summed 'total' forecasts back to power space 
  select_cols <- c("issue_time", "target_time", final_fcs_cols)
  IMF_fcs <- best_forecasts[, ..select_cols]
  
  
  # get actual values for error evaluation 
  IMF_fcs <- merge(IMF_fcs, z_train[,.(TARGETVAR, target_time=timestamp)])
  IMF_errors <- IMF_fcs
  IMF_errors[, (final_fcs_cols) := lapply(.SD, function(x){abs(x-TARGETVAR)}), .SDcols=final_fcs_cols]
  IMF_MAEs <- data.table(nimfs =paste0(c(min_nimfs:max_nimfs)), MAE=as.numeric(IMF_errors[,lapply(.SD, mean), .SDcols=final_fcs_cols]))
  
  
  ## Calculate standard error for the error values too via bootstrapping, since it's only a small sample.
  for (ni in c(min_nimfs:max_nimfs)){
    bs_error <- tsbootstrap(x=IMF_errors[,get(paste0("nimfs=",ni))], nb=200, mean)
    IMF_MAEs[nimfs==ni, 'S.E.' := bs_error$se]
  }
  
  
  plot(IMF_MAEs$nimfs, IMF_MAEs$MAE)
  arrows(x0=as.numeric(IMF_MAEs$nimfs), y0=(IMF_MAEs$MAE-IMF_MAEs$S.E.), y1=(IMF_MAEs$MAE+IMF_MAEs$S.E.), angle=90, length=0.05, code=3)
  
  best_nimfs <- as.numeric(IMF_MAEs[ , .SD[which.min(MAE)]]$nimfs) # the optimal number of IMFs to use.
  ## if there's no significant difference in the plot, can manually choose the lowest number to keep the model simple.
  # best_nimfs <- 3
  opt_params <- fread(file="./opt_parameters.csv")
  opt_params[Zone==zone, Opt_nimf := best_nimfs] # input best_nimfs into opt_parameters csv
  
  ## now find the optimal window length
  ## only select from issue_times where there is enough previous data for the longest rolling window
  start_ts <- z_train[max(try_windowlens), timestamp]
  end_ts <- z_train[(dim(z_train)[1] - horizon), timestamp]
  allposs_window_times <- z_train[(timestamp >= start_ts) & (timestamp <= end_ts), timestamp] # all possible forecast issue times
  poss_window_times <- allposs_window_times[!allposs_window_times %in% c(issue_timesdt$timestamps, issue_times1)] # minus the windows already used
  issue_times2 <- sample(poss_window_times, size=200, replace=FALSE)
  issue_times2dt <- data.table(timestamps=issue_times2, stage=3)
  issue_timesdt <- rbind(issue_timesdt, issue_times1dt,issue_times2dt)
  save(issue_timesdt, file=paste0("./zone_",zone, "_issuetimesdt.rda"))
  
  ## define function to fit models to windows with different lengths
  WL_forecasts <- function(IT, zonedata, hor, opt_nimfs, WLs, opt_models){
    require(Rlibeemd)
    require(data.table)
    require(forecast)
    
    wl_forecasts <- data.table(issue_time=IT, target_time=IT+hor*60*60) # data.table to store IMF forecasts
    imfnames <- c(paste("IMF", c(1:(opt_nimfs-1))), paste("Residual", opt_nimfs))
    
    for (WL in WLs){
      window_start <- IT - (WL-1)*60*60
      windowdata <- zonedata[(timestamp >= window_start) & (timestamp <= IT),]
      
      ## first forecast IMFs 1-9 and Residual 10
      imfs <- ceemdan(windowdata$T_targetvar,num_imfs=opt_nimfs) # create the imfs for this rolling window.
      colnames(imfs) <- imfnames
      
      imf_fcs <- numeric(opt_nimfs) # empty vector to put individual IMF forecasts in.
      
      for (i in c(1:dim(imfs)[2])){
        imfname <- colnames(imfs)[i]
        opt_arma <- opt_model[model == "opt_model", get(imfname)]
        
        if (opt_arma == "arma"){
          # fit arma with d==0 (in reality seems to be faster to not specify d)
          fit <- auto.arima(imfs[,imfname], max.p=6, max.q=4, max.order=8, stepwise = FALSE)
        }else if(opt_arma=="arma1"){
          # fit arma with d=1.
          fit <- auto.arima(imfs[,imfname], d=1, max.p=6, max.q=4, max.order=8, stepwise = FALSE)
        }else{cat("IT=",IT,"," ,imfname,  'not selected a model to fit!')}
        
        imf_fcs[i] <- as.numeric(forecast(fit, h=horizon)[["mean"]])
      }
      
      wl_forecasts[issue_time==IT, paste0("wl=", WL) := sum(imf_fcs)]
    }
    return (wl_forecasts)
  }
  
  # test <- WL_forecasts(issue_times2[2], z_train, horizon, best_nimfs, try_windowlens, opt_model)
  
  cl <- makeCluster(cores-2)
  registerDoParallel(cl)
  start_time <- Sys.time()
  wl_forecasts <- foreach(it = issue_times2) %dopar% {
    WL_forecasts(it, z_train, horizon, best_nimfs, try_windowlens, opt_model)
  }
  print (Sys.time() - start_time)
  stopCluster(cl)
  
  wl_forecasts <- rbindlist(wl_forecasts)
  
  save(wl_forecasts, file=paste0("./wl_forecasts_zone",zone,".rda"))
  # load(file=paste0("./wl_forecasts_zone",zone,".rda"))
  # issue_times2 <- wl_forecasts$issue_time
  
  
  # get actual values for error evaluation 
  wl_fc_cols <- paste0("wl=", try_windowlens)
  wl_forecasts[, (wl_fc_cols) := lapply(.SD, function(x){1/(1+exp(-x))}), .SDcols = wl_fc_cols] # convert forecasts back to power space 
  wl_forecasts <- merge(wl_forecasts, z_train[,.(TARGETVAR, target_time=timestamp)]) # merge on target_time
  wl_errors <- wl_forecasts
  
  wl_errors[, (wl_fc_cols) := lapply(.SD, function(x){abs(x-TARGETVAR)}), .SDcols=wl_fc_cols]
  wl_MAEs <- data.table(windowlen =paste0(try_windowlens), MAE=as.numeric(wl_errors[,lapply(.SD, mean), .SDcols=wl_fc_cols]))
  
  ## bootstrap errors
  for (wl in try_windowlens){
    bs_error <- tsbootstrap(x=wl_errors[,get(paste0("wl=",wl))], nb=200, mean)
    wl_MAEs[windowlen==wl, 'S.E.' := bs_error$se]
  }
  
  plot(wl_MAEs$windowlen, wl_MAEs$MAE)
  arrows(x0=as.numeric(wl_MAEs$windowlen), y0=(wl_MAEs$MAE-wl_MAEs$S.E.), y1=(wl_MAEs$MAE+wl_MAEs$S.E.), angle=90, length=0.05, code=3)
  
  best_wl <- as.numeric(wl_MAEs[ , .SD[which.min(MAE)]]$windowlen) # the optimal window length.
  opt_params[Zone==zone, Opt_wl := best_wl]
  fwrite(opt_params, file="./opt_parameters.csv")
}
