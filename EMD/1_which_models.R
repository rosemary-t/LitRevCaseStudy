require(data.table)
require(rstudioapi)
require(frbs)
require(Rlibeemd)
require(tseries)
require(caret)
require(doParallel)

setwd(dirname(getActiveDocumentContext()$path))



### Code to run several different models on each IMF (and the residuals from creating different numbers of IMFs)
### so we can evaluate error and see which is the best model for each decomposed time series.

############# function to make input/output sets from sliding window time series of all IMFs ##############

make_xydfs <- function(iwt, zdata, wl, min_nimfs, max_nimfs, nlags, hor){
  ## create data.table with all lags, and lags of the differenced series too to fit AR(d=1)
  ## including lags 12 and 24.
  
  window_XYs <- data.table() #  data.table to put X,Y inputs/outputs in
  
  window_start <- iwt - (wl-1)*60*60
  windowdata <- zdata[(timestamp >= window_start) & (timestamp<= iwt),]
  
  ## make imfs from window time series
  imfts <- as.data.table(ceemdan(windowdata$T_targetvar,num_imfs=max_nimfs)) 
  setnames(imfts, 'Residual', paste('Residual', max_nimfs))
  
  ## also need the residual series for all other numbers of imfs:
  range_nimfs <- c(min_nimfs: (max_nimfs-1))
  for (nimfs in range_nimfs){
    imfs <- ceemdan(windowdata$T_targetvar,num_imfs=nimfs)
    imfts[, paste("Residual", nimfs) := as.numeric(imfs[,"Residual"])]
  }
    
  ## now set up matrix of inputs (lags) and outputs, to train AR with, for each IMF separately.
  Ps <- c(c(1:nlags),12,24)
  lagcols <- paste0('lag',Ps)
  lagcols_diff <- paste0('lag', Ps, '_diff')
  for (imfname in names(imfts)){
    dtx <- copy(imfts[, .(y=get(imfname), series=imfname)]) # initialise data.table
    dtx[, (lagcols) := shift(y, (hor-1+Ps), type='lag', fill=NA)]
    dtx[, y_diff := y-shift(y,1)] # create differenced y column (to train arma_d1 on)
    dtx[, (lagcols_diff) := shift(y_diff, (hor-1+Ps), type='lag', fill=NA)] # and the differenced lags too
    dtx <- na.omit(dtx)
    window_XYs <- rbind(window_XYs, dtx)
  }
  return (window_XYs)
}

########################### function to fit and predict from ANFIS model ######################

anfis_cv <- function(XYdt, imfname, cv){
  require(data.table)
  require(frbs)
  require(tseries)
  
  ## always fits ANFIS model with 2 lags as inputs
  ## XYdt is data.table with model inputs(lags) and outputs (output must be last column when supplied to frbs.learn)
  ## imfname = name of IMF to fit to
  ## cv = number of folds to do cross validation with.
  
  error_dt <- data.table(series=imfname, model="ANFIS", P=2)
  
  lagcols <- paste0('lag', c(1:2))
  anfis_cols <- c(lagcols, "y") # y must be the last column in the list (anfis model uses last column as model output)
  anfis_xy <- XYdt[series==imfname, ..anfis_cols]

  anfis_xy <- anfis_xy[complete.cases(anfis_xy),] # remove rows with any NA
  fold_len <- as.integer(dim(anfis_xy)[1]/cv)
  rownums <- c(1:dim(anfis_xy)[1])
  anfis_resids <- numeric(length(rownums))
  
  for (i in c(1:cv)){
    testrows <- c(((i-1)*fold_len):(i*fold_len))
    trainrows <- setdiff(rownums, testrows)
    
    Xtest <- anfis_xy[testrows, !"y"]
    Ytest <- anfis_xy[testrows, y]
    XYtrain <- anfis_xy[trainrows]
    anfis <- frbs.learn(XYtrain, method.type="ANFIS")
    ans <- predict(anfis, Xtest)
    anfis_resids[testrows] <- abs(ans-Ytest)
  }
  
  bs_err_an <- tsbootstrap(x=anfis_resids, nb=200, mean)
  error_dt[,MAE := mean(abs(anfis_resids))]
  error_dt[,SE := bs_err_an$se]
  
  return (error_dt)
}

############################################################################################

for (zone in c(7:10)){
  print (zone)
  split_date <- as.POSIXct("2012-12-31 12:00", tz="UTC") # split the training/testing data here
  horizon <- 1 # number of steps ahead we are forecasting for
  eps <- 0.01 # cutoff point at boundaries
  windowlen <- 500 # window length to use for EMD
  Max_nimfs <- 10
  Min_nimfs <- 3 #(this is 2 IMF series plus a residual)
  Nlags <- 8 # max number of lags we're testing in AR models. Will also test adding in 12,24hr lags as well.
  opt_anfis_lags <- 2 # 8 lags gives a lower error for the case tested in 0_anfis_lags but takes too long to train 
  
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
  
  ## use 1/4 of training data for finding best model per imf, 1/4 for finding opt number of imfs, last 1/2 for fittingmodels and generating mean forecasts to fit variance of residuals to.
  train_chunk_length <- as.integer(dim(z_train)[1]/4)
  train_chunk_one <- z_train[c(1:train_chunk_length)] # just use the first quarter in this code.
  
  start_ts <- train_chunk_one[windowlen, timestamp]
  end_ts <- train_chunk_one[(dim(train_chunk_one)[1] - horizon), timestamp]
  poss_window_times <- train_chunk_one[(timestamp >= start_ts) & (timestamp <= end_ts), timestamp] # the timestamps for the end of all possible sliding windows in the training data.
  issue_times <- sample(poss_window_times, 100, replace=FALSE) # forecast issue times for 100 random windows.

  ## decompose all the sliding windows into their IMFs:
  all_imf_xys <- data.table()
  start_time <- Sys.time()
  for (it in issue_times){
    it_dt <- make_xydfs(it, train_chunk_one, windowlen, Min_nimfs, Max_nimfs, Nlags, horizon)
    all_imf_xys <- rbind(all_imf_xys, it_dt)
  }
  print (Sys.time() - start_time)
  
  imfseries <- unique(all_imf_xys$series)
  
  ## first, fit AR model options.
  models <- models <- c("AR", "AR_d1", "AR_24", "AR_d124")
  model_error_dt <- CJ(series=imfseries, model=models, P=c(2:Nlags))
  
  # fit models for every IMF, for every lag up to Nlags
  registerDoParallel(1)
  for (imfname in imfseries){
    print (imfname)
    for (p in unique(model_error_dt$P)){
      lagcols <- paste0('lag', c(1:p))
      diff_lagcols <- paste0('lag', c(1:p), '_diff')
      ARform <- as.formula(paste0("y ~ ", paste0(lagcols, collapse = "+")))
      ARd1form <- as.formula(paste0("y_diff ~", paste0(diff_lagcols, collapse="+")))
      
      ## fit AR(p) model with d=0
      imf_xy <- all_imf_xys[series==imfname]
      train.control <- trainControl(method = "cv", number = 5)
      model0 <- train(ARform, data = imf_xy, method = "lm", trControl = train.control)
      resids0 <- model0$finalModel$residuals
      bs_err <- tsbootstrap(x=resids0, nb=200, mean)
      model_error_dt[(series==imfname) & (model=="AR") & (P==p), MAE := mean(abs(resids0))]
      model_error_dt[(series==imfname) & (model=="AR") & (P==p), SE := bs_err$se]
      
      ## and AR(p) with d=1
      imf_diff_xy <- all_imf_xys[series==imfname]
      model1 <- train(ARd1form, data = imf_xy, method = "lm", trControl = train.control)
      pred_diffs <- model1$finalModel$fitted.values
      resids1 <- pred_diffs + imf_xy$lag1 - imf_xy$y_diff # forecast y= pred_diff+lag1, then residual by subtracting actual value.
      bs_err1 <- tsbootstrap(x=resids1, nb=200, mean)
      model_error_dt[(series==imfname) & (model=="AR_d1") & (P==p), MAE := mean(abs(resids1))]
      model_error_dt[(series==imfname) & (model=="AR_d1") & (P==p), SE := bs_err1$se]
      
      
      lagcols_24 <- paste0('lag', c(c(1:p), 12,24))
      diff_lagcols_24 <- paste0('lag', c(c(1:p),12,24), '_diff')
      ARform24 <- as.formula(paste0("y ~ ", paste0(lagcols_24, collapse = "+")))
      ARd1form24 <- as.formula(paste0("y_diff ~", paste0(diff_lagcols_24, collapse="+")))
      
      ## fit AR(p) model with d=0 and 12 and 24 hr lags
      model024 <- train(ARform24, data = imf_xy, method = "lm", trControl = train.control)
      resids024 <- model024$finalModel$residuals
      bs_err24 <- tsbootstrap(x=resids024, nb=200, mean)
      model_error_dt[(series==imfname) & (model=="AR_24") & (P==p), MAE := mean(abs(resids024))]
      model_error_dt[(series==imfname) & (model=="AR_24") & (P==p), SE := bs_err24$se]
      
      ## and AR(p) with d=1 and 12 and 24 hr lags
      model124 <- train(ARd1form24, data = imf_xy, method = "lm", trControl = train.control)
      pred_diffs <- model124$finalModel$fitted.values
      resids124 <- pred_diffs + imf_xy$lag1 - imf_xy$y_diff # forecast y= pred_diff+lag1, then residual by subtracting actual value.
      bs_err124 <- tsbootstrap(x=resids124, nb=200, mean)
      model_error_dt[(series==imfname) & (model=="AR_d124") & (P==p), MAE := mean(abs(resids124))]
      model_error_dt[(series==imfname) & (model=="AR_d124") & (P==p), SE := bs_err124$se]
    }
  }
  
  ## and now fit ANFIS for each IMF, in parallel
  # test <- anfis_cv(all_imf_xys, imfname, 2)

  cores <- detectCores()
  cl <- makeCluster(cores-2)
  registerDoParallel(cl)
  start_time <- Sys.time()
  anfis_errors <- foreach (imfname = imfseries) %dopar% {
    anfis_cv(all_imf_xys, imfname, 2)
  }
  print (Sys.time() - start_time)
  stopCluster(cl)
  
  model_error_dt <- rbind(model_error_dt, rbindlist(anfis_errors))
  save(model_error_dt, file=paste0("./model_errors_zone",zone,".rda"))
}
