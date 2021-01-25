require(Rlibeemd)
require(data.table)
require(rstudioapi)
require(tseries)
require(caret)

setwd(dirname(getActiveDocumentContext()$path))

### Code to run the best model on each IMF (and the residuals from creating different numbers of IMFs)
### so we can evaluate error and find the optimal number of IMFs to use.

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
  # lagcols_diff <- paste0('lag', Ps, '_diff')
  for (imfname in names(imfts)){
    dtx <- copy(imfts[, .(y=get(imfname), series=imfname, target_time=windowdata$timestamp)]) # initialise data.table
    dtx[, (lagcols) := shift(y, (hor-1+Ps), type='lag', fill=NA)]
    # dtx[, y_diff := y-shift(y,1)] # create differenced y column (to train arma_d1 on)
    # dtx[, (lagcols_diff) := shift(y_diff, (hor-1+Ps), type='lag', fill=NA)] # and the differenced lags too
    dtx <- na.omit(dtx)
    window_XYs <- rbind(window_XYs, dtx)
  }
  return (window_XYs)
}

############################################################################################

zone <- 1 # which area we are forecasting for
horizon <- 1 # number of steps ahead we are forecasting for
split_date <- as.POSIXct("2012-12-31 12:00") # split the training/testing data here
eps <- 0.01 # cutoff point at boundaries
windowlen <- 500 # window length to use for EMD
Max_nimfs <- 10
Min_nimfs <- 3 #(this is 2 IMF series plus a residual)

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

## use 1/4 of training data for finding best model per imf, 1/4 for finding opt number of imfs, last 1/2 for fittingmodels and generating mean forecasts to fit variance of residuals to.
train_chunk_length <- as.integer(dim(z_train)[1]/4)
train_chunk_two <- z_train[(1*train_chunk_length):(2*train_chunk_length)] # just use the first quarter in this code.

start_ts <- train_chunk_two[windowlen, timestamp]
end_ts <- train_chunk_two[(dim(train_chunk_two)[1] - horizon), timestamp]
issue_times <- seq(from=start_ts, to=end_ts, by=windowlen*60*60)

## load the table specifying what model for each IMF:
opt_model <- load(file=paste0("./model_errors_zone",zone,".rda"))
best_model <- model_error_dt[ , .SD[which.min(MAE)], by = c("series")]
Nlags <- max(best_model$P)

## make data table with forecast inputs and outputs
## decompose all the sliding windows into their IMFs:
all_imf_xys <- data.table()
start_time <- Sys.time()
for (it in issue_times){
  it_dt <- make_xydfs(it, train_chunk_two, windowlen, Min_nimfs, Max_nimfs, Nlags, horizon)
  all_imf_xys <- rbind(all_imf_xys, it_dt)
}
print (Sys.time() - start_time)

imfseries <- unique(all_imf_xys$series)
all_imf_fcs <- data.table(target_time=unique(all_imf_xys$target_time)) # data.table to put forecasts in, for all IMF series.

# fit and predict (by cross validation) 'best' model for each IMF
for (imfname in imfseries){
  imf_model <- best_model[series==imfname, model]
  imf_p <- best_model[series==imfname, P]
  
  if(imf_model=="AR"){
    lagcols <- paste0('lag', c(1:imf_p))
  }else if(imf_model=="AR_24"){
    lagcols <- paste0('lag', c(c(1:imf_p), 12,24))
  }
  
  form <- as.formula(paste0("y ~ ", paste0(lagcols, collapse = "+")))
  imf_xy <- all_imf_xys[series==imfname]
  train.control <- trainControl(method = "cv", number = 5)
  model <- train(form, data = imf_xy, method = "lm", trControl = train.control)
  imf_fcs <- data.table(target_time=imf_xy$target_time)
  imf_fcs[, paste0(imfname):=model$finalModel$fitted.values]
  
  all_imf_fcs <- merge(all_imf_fcs, imf_fcs, by="target_time")
}

actual_values <- train_chunk_two[, .(target_time=timestamp, TARGETVAR)]
all_imf_fcs <- merge(actual_values, all_imf_fcs, by="target_time")


## now need to make the 'final' summed forecasts for different numbers of IMFs.
for (nimfs in c(Min_nimfs:Max_nimfs)){
  # first get a list of best_forecasts columns that must be summed for an overall forecast with nimfs (e.g. nimfs=3 means IMF 1+IMF 2+Residual 3)
  imf_names <- c(paste("IMF", c(1:(nimfs-1))), paste("Residual", nimfs))
  all_imf_fcs[, paste0("nimfs=",nimfs) := rowSums(.SD), .SDcols=imf_names]
}

final_fcs_cols <- paste0("nimfs=", c(Min_nimfs:Max_nimfs)) # names of the columns in best_forecasts containing the summed forecasts
all_imf_fcs[, (final_fcs_cols) := lapply(.SD, function(x){1/(1+exp(-x))}), .SDcols = final_fcs_cols] # convert summed 'total' forecasts back to power space 
IMF_errors <- all_imf_fcs
IMF_errors[, (final_fcs_cols) := lapply(.SD, function(x){abs(x-TARGETVAR)}), .SDcols=final_fcs_cols]
IMF_MAEs <- data.table(nimfs =paste0(c(Min_nimfs:Max_nimfs)), MAE=as.numeric(IMF_errors[,lapply(.SD, mean), .SDcols=final_fcs_cols]))

## Calculate standard error for the error values too via bootstrapping, since it's only a small sample.
for (ni in c(Min_nimfs:Max_nimfs)){
  bs_error <- tsbootstrap(x=IMF_errors[,get(paste0("nimfs=",ni))], nb=200, mean)
  IMF_MAEs[nimfs==ni, 'S.E.' := bs_error$se]
}


plot(IMF_MAEs$nimfs, IMF_MAEs$MAE)
arrows(x0=as.numeric(IMF_MAEs$nimfs), y0=(IMF_MAEs$MAE-IMF_MAEs$S.E.), y1=(IMF_MAEs$MAE+IMF_MAEs$S.E.), angle=90, length=0.05, code=3)

best_nimfs <- as.numeric(IMF_MAEs[ , .SD[which.min(MAE)]]$nimfs) # the optimal number of IMFs to use.
## if there's no significant difference in the plot, can manually choose the lowest number near the minimum to keep the model simple.
# best_nimfs <- 5
# opt_params <- fread(file="./opt_parameters.csv")
# opt_params[Zone==zone, Opt_nimf := best_nimfs] # input best_nimfs into opt_parameters csv
# fwrite(opt_params, file="./opt_parameters.csv" )


