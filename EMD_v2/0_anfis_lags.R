require(data.table)
require(ggplot2)
require(frbs)
require(Rlibeemd)
require(tseries)
require(doParallel)
require(rstudioapi)


setwd(dirname(getActiveDocumentContext()$path))

## to find the optimal inputs to ANFIS model
horizon <- 1 # just one step ahead forecast for now.


make_xydfs_short <- function(iwt, zdata, wl, max_nimfs, nlags, hor){
  window_XYs <- data.table() #  data.table to put X,Y inputs/outputs in
  
  window_start <- iwt - (wl-1)*60*60
  windowdata <- zdata[(timestamp >= window_start) & (timestamp<= iwt),]

  imfts <- as.data.table(ceemdan(windowdata$T_targetvar,num_imfs=max_nimfs))
  
  ## now set up matrix of inputs (lags) and outputs, to train AR with, for each IMF separately.
  Ps <- c(1:nlags)
  lagcols <- paste0('lag',Ps)
  for (imfname in names(imfts)){
    dtx <- copy(imfts[, .(y=get(imfname), series=imfname)]) # initialise data.table 
    dtx[, (lagcols) := shift(y, (hor-1+Ps), type='lag', fill=NA)]
    dtx <- na.omit(dtx)
    window_XYs <- rbind(window_XYs, dtx)
  }
  
  return (window_XYs)
}

anfis_cv <- function(XYdt, cv){
  ## XYdt is data.table with only model inputs(lags) and outputs (output must be last column)
  ## cv = number of folds to do cross validation with.
  
  XYdt <- XYdt[complete.cases(XYdt),] # remove rows with any NA
  keep_rows <- sample.int(dim(XYdt)[1], as.integer(dim(XYdt)[1]/2)) # only use half the rows - otherwise takes ages to train
  XYdt <- XYdt[keep_rows]
  fold_len <- as.integer(dim(XYdt)[1]/cv)
  rownums <- c(1:dim(XYdt)[1])
  residuals <- numeric(length(rownums))
  
  for (i in c(1:cv)){
    testrows <- c(((i-1)*fold_len):(i*fold_len))
    trainrows <- setdiff(rownums, testrows)
    
    Xtest <- XYdt[testrows, !"y"]
    Ytest <- XYdt[testrows, y]
    XYtrain <- XYdt[trainrows]
    start_time <- Sys.time()
    anfis <- frbs.learn(XYtrain, method.type="ANFIS")
    print (Sys.time()-start_time)
    ans <- predict(anfis, Xtest)
    residuals[testrows] <- abs(ans-Ytest)
  }
  
  return (residuals)
}

####

zone <- 1
split_date <- as.POSIXct("2012-12-31 12:00", tz="UTC") # split the training/testing data here
horizon <- 1 # number of steps ahead we are forecasting for
eps <- 0.01 # cutoff point at boundaries
windowlen <- 500 # window length to use for EMD
nimfs <- 10
Nlags <- 8 # max number of lags

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

## use 1/4 of training data for finding best model per imf, 1/4 for finding opt number of imfs, 1/4 for opt window length, 1/4 for generating mean forecasts to fit variance of residuals to.
train_chunk_length <- as.integer(dim(z_train)[1]/4)
train_chunk_one <- z_train[c(1:train_chunk_length)] # just use the first quarter in this code.

load(paste0("./zone_",zone,"_issuetimesdt.rda")) # keep a record of which training windows have been used.


## decompose all the sliding windows into their IMFs:
all_imf_xys <- data.table()
start_time <- Sys.time()
for (it in issue_timesdt$timestamps){
  it_dt <- make_xydfs_short(it, train_chunk_one, windowlen, nimfs, Nlags, horizon)
  all_imf_xys <- rbind(all_imf_xys, it_dt)
}
print (Sys.time() - start_time)



## function to find errors per lag, for IMF 4 only
imf4_errors <- function(p, all_imf_xys){
  require(data.table)
  require(frbs)
  require(tseries)

  error_dt <- data.table(series="IMF 4", lag=p)
  lagcols <- paste0('lags', c(1:p))
  anfis_cols <- c(lagcols, "y") # y must be the last column in the list (anfis model uses last column as model output)
  anfis_xy <- all_imf_xys[series=="IMF 4", ..anfis_cols]
  anfis_resids <- anfis_cv(anfis_xy, 2)
  bs_err_an <- tsbootstrap(x=anfis_resids, nb=200, mean)

  error_dt[lag==p, MAE := mean(anfis_resids)]
  error_dt[lag==p, SE := bs_err_an$se]

  return (error_dt)
}

## this takes more then 12 hours...
cores <- detectCores()
cl <- makeCluster(cores-2)
registerDoParallel(cl)
start_time <- Sys.time()
errorsdf <- foreach (P = c(2:Nlags)) %dopar% {
  imf4_errors(P, all_imf_xys)
}
print (Sys.time() - start_time)
stopCluster(cl)

error_dt <- rbindlist(errorsdf)
save(error_dt, file=('./anfis_lag_errors.rda'))
