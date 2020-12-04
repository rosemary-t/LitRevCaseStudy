require(Rlibeemd)
require(data.table)
require(forecast)
require(rstudioapi)

## also loads forecast and Rlibeemd packages within make_forecasts function.

setwd(dirname(getActiveDocumentContext()$path))

### Code to run the best model on each IMF (and the residuals from creating different numbers of IMFs)
### so we can evaluate error and find the optimal number of IMFs to use.
### and then the same for the optimal window length too.

zone <- 3 # which area we are forecasting for
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
load(paste0("./zone_",zone, "_issuetimes.rda")) # get the issue times of the random windows already used to find the optimal model
used_times <- issue_times 
rm(issue_times)


start_ts <- z_train[i_wl, timestamp]
end_ts <- z_train[(dim(z_train)[1] - horizon), timestamp]
allposs_window_times <- z_train[(timestamp >= start_ts) & (timestamp <= end_ts), timestamp] # all forecast issue times
poss_window_times <- allposs_window_times[!allposs_window_times %in% used_times] # minus the windows already used
issue_times1 <- sample(poss_window_times, size=200, replace=FALSE) # pick 200 random windows to optimise nimfs over
issue_times <- c(used_times, issue_times1) # add issue_times1 to the list of used window times.
save(issue_times, file=paste0("./zone_",zone, "_issuetimes.rda"))


## load the table specifying what model for each IMF:
opt_model <- fread(file=paste0("./whichmodels_errorinfo_zone",zone,".csv"))


## define function to generate forecasts, using the best model, for every IMF.
series_forecasts <- function(IT, zonedata, h, min_nimfs, max_nimfs, WL, opt_models){
  require(Rlibeemd)
  require(data.table)
  require(forecast)
  
  best_forecasts <- data.table(issue_time=IT) # data.table to store IMF forecasts
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
    }else{cat("wt=",wt,"," ,imfname,  'not selected a model to fit!')}
    
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

test <- series_forecasts(issue_times1[2], z_train, horizon, min_nimfs, max_nimfs, i_wl, opt_model)

## now run the forecasts for each window
cl <- makeCluster(6)
registerDoParallel(cl)
start_time <- Sys.time()
model_forecasts <- foreach(it = issue_times1) %dopar% {
  series_forecasts(it, z_train, horizon, min_nimfs, max_nimfs, i_wl, opt_model)
}
print (Sys.time() - start_time)
stopCluster(cl)

# save(best_forecasts, file=paste0(path,"\\best_forecasts_zone",zone,".rda"))
# load(paste0(path,"\\best_forecasts_zone",zone,".rda"))
# issue_times1 <- best_forecasts$windownum

## now need to make the 'final' summed forecasts for different numbers of IMFs.
for (nimfs in c(min_nimfs:max_nimfs)){
  # first get a list of best_forecasts columns that must be summed for an overall forecast with nimfs (e.g. nimfs=3 means IMF 1+IMF 2+Residual 3)
  imf_names <- c(paste("IMF", c(1:(nimfs-1))), paste("Residual", nimfs))
  best_forecasts[, paste0("nimfs=",nimfs) := rowSums(.SD), .SDcols=imf_names]
}

final_fcs_cols <- paste0("nimfs=", c(min_nimfs:max_nimfs)) # names of the columns in best_forecasts containing the summed forecasts
best_forecasts[, (final_fcs_cols) := lapply(.SD, function(x){1/(1+exp(-x))}), .SDcols = final_fcs_cols] # convert summed 'total' forecasts back to power space 
select_cols <- c("issue_time", final_fcs_cols)
IMF_fcs <- best_forecasts[, ..select_cols]


# get actual values for error evaluation 
IMF_fcs[, target_time := issue_time + horizon*60*60]
IMF_fcs <- merge(IMF_fcs, z_train[,.(TARGETVAR, target_time=timestamp)])
IMF_errors <- IMF_fcs
IMF_errors[, (final_fcs_cols) := lapply(.SD, function(x){abs(x-TARGETVAR)}), .SDcols=final_fcs_cols]
MAEs <- data.table(nimfs =paste0(c(min_nimfs:max_nimfs)), MAE=as.numeric(IMF_errors[,lapply(.SD, mean), .SDcols=final_fcs_cols]))


## Calculate standard error for the error values too via bootstrapping, since it's only a small sample.
require(tseries)
for (ni in c(min_nimfs:max_nimfs)){
  bs_error <- tsbootstrap(x=IMF_errors[,get(paste0("nimfs=",ni))], nb=200, mean)
  MAEs[nimfs==ni, 'S.E.' := bs_error$se]
}


plot(MAEs$nimfs, MAEs$MAE, ylim=c(0.05,0.07))
arrows(x0=as.numeric(MAEs$nimfs), y0=(MAEs$MAE-MAEs$S.E.), y1=(MAEs$MAE+MAEs$S.E.), angle=90, length=0.05, code=3)

best_nimfs <- as.numeric(MAEs[ , .SD[which.min(MAE)]]$nimfs) # the optimal number of IMFs to use.
## if there's no significant difference in the plot, can manually choose the lowest number to keep the model simple.
# best_nimfs <- 3


## now find the optimal window length
## only select from issue_times where there is enough previous data for the longest rolling window
start_ts <- z_train[max(try_windowlens), timestamp]
end_ts <- z_train[(dim(z_train)[1] - horizon), timestamp]
allposs_window_times <- z_train[(timestamp >= start_ts) & (timestamp <= end_ts), timestamp] # all possible forecast issue times
poss_window_times <- allposs_window_times[!allposs_window_times %in% c(used_times, issue_times1)] # minus the windows already used
issue_times2 <- sample(poss_window_times, size=200, replace=FALSE)
wl_forecasts <- data.table(issue_time=issue_times2)
imfnames <- c(paste("IMF", c(1:(best_nimfs-1))), paste("Residual", best_nimfs))


i <- 1
for (wt in issue_times2){#for every sliding window
  print (i)
  i <- i+1
  for (wl in try_windowlens){#for every window length
    
    window_start <- wt - (wl-1)*60*60
    windowdata <- z_train[(timestamp >= window_start) & (timestamp <= wt),] 

    imfs <- ceemdan(windowdata$T_targetvar,num_imfs=best_nimfs) # create the imfs for this rolling window.
    imflist <- numeric(best_nimfs) # empty list to put each IMF forecast in.
    colnames(imfs) <- imfnames # so the residual column now matches its name in opt_model.
    
    for (i in c(1:best_nimfs)){
      imfname <- imfnames[i]
      opt_arma <- opt_model[imf_name == imfname, `best_model (value for d)`]
      
      if (opt_arma == 0){
        # fit arma with d==0 (in reality seems to be faster to not specify d)
        fit <- auto.arima(imfs[,imfname], max.p=6, max.q=4, max.order=8, stepwise = FALSE)
      }else if(opt_arma==1){
        # fit arma with d=1.
        fit <- auto.arima(imfs[,imfname], d=1, max.p=6, max.q=4, max.order=8, stepwise = FALSE)
      }else{cat("wt=",wt,"," ,imfname,  'not selected a model to fit!')}
      
      imflist[i] <- as.numeric(forecast(fit, h=horizon)[["mean"]]) # collect the forecast for each IMF.
    }
    
    wl_forecasts[issue_time== wt, paste0("wl=",wl) := sum(imflist)] # save the overall forecast
  }
}

# save(wl_forecasts, file=paste0(path,"\\wl_forecasts_zone",zone,".rda"))
# load(file=paste0(path,"\\wl_forecasts_zone",zone,".rda"))


# get actual values for error evaluation 
wl_fc_cols <- paste0("wl=", try_windowlens)
wl_forecasts[, target_time := issue_time + horizon*60*60]
wl_forecasts[, (wl_fc_cols) := lapply(.SD, function(x){1/(1+exp(-x))}), .SDcols = wl_fc_cols] # convert forecasts back to power space 
wl_forecasts <- merge(wl_forecasts, z_train[,.(TARGETVAR, target_time=timestamp)]) # merge on target_time
wl_errors <- wl_forecasts

wl_errors[, (wl_fc_cols) := lapply(.SD, function(x){abs(x-TARGETVAR)}), .SDcols=wl_fc_cols]
wl_MAEs <- data.table(windowlen =paste0(try_windowlens), MAE=as.numeric(wl_errors[,lapply(.SD, mean), .SDcols=wl_fc_cols]))
require(tseries)
for (wl in try_windowlens){
  bs_error <- tsbootstrap(x=wl_errors[,get(paste0("wl=",wl))], nb=200, mean)
  wl_MAEs[windowlen==wl, 'S.E.' := bs_error$se]
}

plot(wl_MAEs$windowlen, wl_MAEs$MAE, ylim=c(0.08, 0.13))
arrows(x0=as.numeric(wl_MAEs$windowlen), y0=(wl_MAEs$MAE-wl_MAEs$S.E.), y1=(wl_MAEs$MAE+wl_MAEs$S.E.), angle=90, length=0.05, code=3)

best_wl <- as.numeric(wl_MAEs[ , .SD[which.min(MAE)]]$windowlen) # the optimal window length.
