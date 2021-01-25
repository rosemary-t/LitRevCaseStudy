require(Rlibeemd)
require(data.table)
require(doParallel)
require(rstudioapi)
require(ProbCast)

setwd(dirname(getActiveDocumentContext()$path))

## train models on portion of training set, and produce mean forecasts for rest of training set as well as all of test set

############# function to make input/output sets from sliding window time series of all IMFs ##############

make_xydfs <- function(iwt, zdata, wl, nimfs, nlags, hors){
  ## create data.table with all lags, and lags of the differenced series too to fit AR(d=1)
  ## including lags 12 and 24.
  
  window_XYs <- data.table() #  data.table to put X,Y inputs/outputs in
  hors_cols <- paste0('y_',hors)
  
  window_start <- iwt - (wl-1)*60*60
  windowdata <- zdata[(timestamp >= window_start) & (timestamp<= iwt),]
  
  ## make imfs from window time series
  imfts <- as.data.table(ceemdan(windowdata$T_targetvar,num_imfs=nimfs)) 
  setnames(imfts, 'Residual', paste('Residual', nimfs))
  
  
  ## now set up matrix of inputs (lags) and outputs, to train AR with, for each IMF separately.
  Ps <- c(c(1:nlags),12,24)
  lagcols <- paste0('lag',Ps)
  for (imfname in names(imfts)){
    dtx <- copy(imfts[, .(y=get(imfname), series=imfname, target_time=windowdata$timestamp)]) # initialise data.table
    dtx[, (hors_cols) := shift(y, (max(hors)-hors), type='lag', fill=NA)]
    dtx[, (lagcols) := shift(y, (max(hors)-1+Ps), type='lag', fill=NA)]
    dtx <- na.omit(dtx)
    window_XYs <- rbind(window_XYs, dtx)
  }
  return (window_XYs)
}

############# function to make input lags for a given issue time, for forecasts ##############

make_Xs <- function(iwt, zdata, wl, nimfs, nlags){
  require(data.table)
  require(Rlibeemd)
  
  window_x <- data.table() #  data.table to put X,Y inputs/outputs in

  window_start <- iwt - (wl-1)*60*60
  windowdata <- zdata[(timestamp >= window_start) & (timestamp<= iwt),]
  
  ## make imfs from window time series
  imfts <- as.data.table(ceemdan(windowdata$T_targetvar,num_imfs=nimfs)) 
  setnames(imfts, 'Residual', paste('Residual', nimfs))
  
  ## now set up matrix of inputs (lags) and outputs, to train AR with, for each IMF separately.
  Ps <- c(24,12, c(nlags:1))
  lagcols <- paste0('lag',Ps)
  for (imfname in names(imfts)){
    dtx <- data.table(issue_time=iwt, series=imfname)
    dtx[, (lagcols) := as.list(imfts[(windowlen-Ps+1), get(imfname)])]
    window_x <- rbind(window_x, dtx)
  }
  return(window_x)
}

####################### function to return forecasts for each horizon ########################

horizon_fcs <- function(all_xydt, all_x, hor, bestmodel){
  ## fit model on rows in all_xydt
  ## forecast for each row in all_x
  
  imfnames <- bestmodel$series
  all_imf_fcs <- data.table(issue_time=unique(all_x$issue_time))
  
  # fit and forecast each IMF
  for (imfname in imfnames){
    imf_model <- bestmodel[series==imfname, model]
    imf_p <- bestmodel[series==imfname, P]
    
    XY <- all_xydt[series==imfname]
    X <- all_x[series==imfname]
    
    if(imf_model=="AR"){
      lagcols <- paste0('lag', c(1:imf_p))
    }else if(imf_model=="AR_24"){
      lagcols <- paste0('lag', c(c(1:imf_p), 12,24))
    }else{print("different model is optimum")}
    
    form <- as.formula(paste0("y_",hor," ~ ", paste0(lagcols, collapse = "+")))
    fit <- lm(form, XY)
    forecasts <- data.table(issue_time=X$issue_time, fc=predict(fit, X))
    setnames(forecasts, "fc", paste0(imfname))
    
    all_imf_fcs <- merge(all_imf_fcs, forecasts, by="issue_time")
  }
  
  all_imf_fcs[, mean_fc := rowSums(.SD), .SDcols = imfnames]
  all_imf_fcs[, target_time := issue_time+hor*60*60]
  
  return (all_imf_fcs[, .(issue_time, target_time, mean_fc)])
}

##############################################################################################

z <- 10 # which area we are forecasting for
horizons <- c(1:6) # number of steps ahead we are forecasting for
split_date <- as.POSIXct("2012-12-31 12:00") # split the training/testing data here
eps <- 0.01 # cutoff point at boundaries for lognormal transformation
windowlen <- 500 # window length to use for EMD

## list of the model (and lags) to fit for each IMF
load(paste0("./model_errors_zone",z,".rda"))

## and get the number of IMFs to use.
opt_params <- fread(file="./opt_parameters.csv")
opt_nimfs <- opt_params[Zone==z, Opt_nimf]
imfnames <- c(paste("IMF", c(1:(opt_nimfs-1))),paste("Residual", opt_nimfs))
best_model <- model_error_dt[series %in% imfnames , .SD[which.min(MAE)], by = c("series")]
Nlags <- max(best_model$P)


# now load the data and prepare it
load(paste0('../Data/data_zone', z,'.rda'))
assign('zdata', get(paste0('data_zone', z)))
rm(list=paste0("data_zone",z))

zdata[, TARGETVAR := approx(.I, TARGETVAR, .I)$y] ## interpolate missing values in  TARGETVAR
zdata <- zdata[!is.na(zdata$TARGETVAR), ] ## if the last row in the column is NA, needs to be removed as can't be interpolated

## lognormal transform: 
zdata[, T_targetvar := lapply(.SD, function(x){ifelse(x <= eps,  eps, ifelse(x >= (1-eps), (1-eps), x))}), .SDcols=c("TARGETVAR")] # clip to [eps, (1-eps)] range
zdata[, T_targetvar := log(T_targetvar/(1-T_targetvar))]

## split into validation and testing data
z_train <- copy(zdata[timestamp<= split_date,])

## use 1/4 of training data for finding best model per imf, 1/4 for finding opt number of imfs, last 1/2 for fittingmodels and generating mean forecasts to fit variance of residuals to.
train_chunk_length <- as.integer(dim(z_train)[1]/4)
train_chunk_three <- z_train[(2*train_chunk_length):(4*train_chunk_length)]

## need to produce mean forecasts now using opt_nimfs and best_model
## mean forecasts for times in the training set are used to find the variance of the residuals
## and the test set is used to evaluate the forecasts out of sample.

## first, make XY 
start_ts <- train_chunk_three[windowlen, timestamp]
end_ts <- train_chunk_three[(dim(train_chunk_three)[1] - max(horizons)), timestamp]
issue_times <- seq(from=start_ts, to=end_ts, by=windowlen*60*60)

all_imf_xys <- data.table()
for (it in issue_times){
  it_dt <- make_xydfs(it, train_chunk_three, windowlen, opt_nimfs, Nlags, horizons)
  all_imf_xys <- rbind(all_imf_xys, it_dt)
}


end_test_ts <- zdata[(dim(zdata)[1] - max(horizons)), timestamp]
forecast_window_times <- zdata[(timestamp >= start_ts) & (timestamp <= end_test_ts), timestamp] # the timestamps for the end of all possible sliding windows in the training data.
# test <- make_Xs(forecast_window_times[10], zdata, windowlen, opt_nimfs, Nlags)

cores <- detectCores()
cl <- makeCluster(cores-2)
registerDoParallel(cl)
start_time <- Sys.time()
X_fc_rows <- foreach(wt =  forecast_window_times) %dopar% {
  make_Xs(wt, zdata, windowlen, opt_nimfs, Nlags)
}
print (Sys.time() - start_time)
stopCluster(cl)

all_X_rows <- rbindlist(X_fc_rows)

all_mean_fcs <- data.table()
for (h in horizons){
  hfcs <- horizon_fcs(all_imf_xys, all_X_rows, h, best_model)
  all_mean_fcs <- rbind(all_mean_fcs, hfcs)
}

all_mean_fcs[, Horizon := (target_time - issue_time)]

save(all_mean_fcs, file=paste0("./mean_forecasts_zone", z,".rda"))
