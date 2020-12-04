require(data.table)
require(glmnet)
require(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))


## define variables...
split_date <- as.POSIXct("2012-12-31 12:00") # every target time up to and including split_date is the training data.
horizons <- c(1:6) # forecast 1 to 6 hours ahead
eps <- 0.01
nfolds <- 5
Ps <- c(1:7) # the range of lags we want to test over
N_lambda <- 25 # number of lambda (regularisation parameters) to test over to get optimum
zones <- c(1:10)


make_lags <- function(dt, p, varcols, horizon){
  ## dt = the data.table with all sites' time series plus timestamp column.
  ## p = the maximum lag to include. Plus 12 and 24 hr lags are also always added as well.
  dtx <- copy(dt)
  #Ps <- c(c(1:p),12,24)
  Ps <- c(1:p)
  for (i in Ps){
    newcols <- paste0(varcols, '_lag', i)
    dtx[, (newcols) := lapply(.SD, function(x){shift(x, (horizon-1+i), type='lag', fill=NA)}), .SDcols=varcols]
  }
  dtx <- na.omit(dtx)
  setnames(dtx, 'timestamp', 'target_time') # timestamp is now the target time of the forecast.
  dtx[, 'issue_time' := target_time - horizon*60*60] # issue time is the time of the most recent (called lag 1) forecast inputs
  return (dtx)
}


if (file.exists("./all_sites_power.rda")){
  load("./all_sites_power.rda")
  varcols <- names(alldata)
  varcols <- varcols[varcols!= 'timestamp'] # names of all the columns with site power data.
}else{
  ## load all zones' data and put in one data table.
  for (zone in zones){
    load(paste0('../Data/data_zone', zone,'.rda'))
    assign("zdata", get(paste0("data_zone", zone)))
    rm(list=paste0("data_zone", zone))
    zdata <- zdata[,.(TARGETVAR, timestamp)] # just keep the power (not using the NWP forecasts)
    setnames(zdata, old="TARGETVAR",new=paste0("zone",zone, "power"))

    if (zone == 1){
      alldata <- copy(zdata)
    }else{
      alldata <- merge(alldata, zdata)
    }

  }
  rm(zdata)


  ## Lognormal transform.... first clip upper and lower values to threshold given by eps
  varcols <- names(alldata)
  varcols <- varcols[varcols!= 'timestamp'] # exclude the timestamp column from transformations!
  alldata[, (varcols) := lapply(.SD, function(x){ifelse(x <= eps,  eps, ifelse(x >= (1-eps), (1-eps), x))}), .SDcols=varcols] # clip to [eps, (1-eps)] range
  alldata[, (varcols) := lapply(.SD, function(x){log(x/(1-x))}), .SDcols=varcols] # lognormal transformation.
  save(alldata, file="./all_sites_power.rda")
}


## split into training and testing data
traindata <- copy(alldata[timestamp <= split_date,])
# testdata <- copy(alldata[timestamp > split_date,])

## create data.table to put cv results for different p (lags) in: glmnet optimises for lambda by cross validation.
cv_results <- CJ(pval=Ps, Horizons=horizons) # row for each value of p, for each horizon.
training_mean_fcs <- data.table() # need out-of-sample forecasts on the training set, to get the variance of the residuals
testing_mean_fcs <- data.table() # mean forecasts for the test set will be used to generate qunatile forecasts.

for (h in horizons){
  print (h)
  for (p in Ps){
    lag_dt = make_lags(traindata, p, varcols, horizon=h)
    FXp <- lag_dt[,.SD, .SDcols = !c('issue_time', 'target_time', varcols)] # our forecast inputs are all the lags, but not varcols (those are forecast outputs) and timestamps.
    FYp <- lag_dt[,..varcols] # the forecast output columns.
    
    #start_time = Sys.time()
    cvfit = cv.glmnet(as.matrix(FXp), as.matrix(FYp), type.measure = "mae", family="mgaussian", nfolds = 5, alpha=1, nlambda=N_lambda)
    #print (Sys.time() - start_time)
    opt_lambda <- cvfit$lambda.1se
    opt_lambda_index <- which(cvfit$lambda == opt_lambda)
    
    cv_results[(pval==p & Horizons==h), lambda := opt_lambda]
    cv_results[(pval==p & Horizons==h), MAE := cvfit$cvm[opt_lambda_index]]
    cv_results[(pval==p & Horizons==h), MAE_sd := cvfit$cvsd[opt_lambda_index]] # if you also want the standard deviation of the error
  }
   
  
  # plot(cv_results$pval, cv_results$MAE, ylim=c(9:10))
  # arrows(x0=as.numeric(cv_results$pval), y0=(cv_results$MAE-cv_results$MAE_sd), y1=(cv_results$MAE+cv_results$MAE_sd), angle=90, length=0.05, code=3)

  cv_opt <- cv_results[Horizons==h,.SD[which.min(MAE)]] # get the minimum values.
  p_opt <- cv_opt$pval
  
  ## now we have the optimal number of lags and regularisation on the training set, fit model
  bigX <- make_lags(alldata, p_opt, varcols, horizon=h) # make forecast inputs and outputs for all time points
  
  ## first, get out-of-sample forecasts for the training dates by cross validation (to use in fitting model for variance)
  traindata_matrix <- bigX[target_time <= split_date,] # forecast inputs and outputs, for this horizon and the optimal lag value
  print (dim(traindata_matrix))
  ## add a kfold label to traindata_matrix
  foldlen <- as.integer(dim(traindata_matrix)[1]/nfolds)
  traindata_matrix$kfold <- "fold 1"
  for (i in c(1:(nfolds-1))){
    traindata_matrix[c((i*foldlen+1):dim(traindata_matrix)[1]), kfold := paste("fold", i+1)]
  }
  
  for (k in unique(traindata_matrix$kfold)){
    xtrain <- as.matrix(traindata_matrix[kfold != k, .SD, .SDcols = !c('issue_time', 'target_time', 'kfold', varcols)])
    ytrain <- as.matrix(traindata_matrix[kfold != k, ..varcols])
    xtest <- as.matrix(traindata_matrix[kfold == k, .SD, .SDcols = !c('issue_time', 'target_time', 'kfold', varcols)])
    test_times <- traindata_matrix[kfold == k, c('issue_time','target_time')]
    
    model <- glmnet(xtrain, ytrain, family='mgaussian', lambda=cv_opt$lambda)
    fcs = as.data.table(drop(predict(model, newx=xtest)))
    names(fcs) <- paste0(names(fcs),"_meanfcs")
    fcs <- cbind(test_times, fcs)
    training_mean_fcs <- rbind(training_mean_fcs, fcs)
    
  }
  
  
  
  ## now fit model on all training data and use to forecast for the test set
  Xtrain <- as.matrix(bigX[target_time <= split_date, .SD, .SDcols = !c('issue_time', 'target_time', varcols)])
  Ytrain <- as.matrix(bigX[target_time <= split_date, ..varcols])
  Xtest <- as.matrix(bigX[target_time > split_date, .SD, .SDcols = !c('issue_time', 'target_time', varcols)])
  # Ytest <- as.matrix(bigX[target_time > split_date, ..varcols])
  test_times <- bigX[target_time > split_date, c('issue_time','target_time')]

  model <- glmnet(Xtrain, Ytrain, family='mgaussian', lambda=cv_opt$lambda)
  fcs = as.data.table(drop(predict(model, newx=Xtest)))
  names(fcs) <- paste0(names(fcs),"_meanfcs")
  fcs <- cbind(test_times, fcs)
  testing_mean_fcs <- rbind(testing_mean_fcs, fcs) # put this horizon's forecasts into overall data.table
}

training_mean_fcs[,'Horizon':= (target_time - issue_time)]
testing_mean_fcs[,'Horizon':= (target_time - issue_time)]


save(training_mean_fcs, file="./mean_forecasts_trainset.rda")
save(testing_mean_fcs, file="./mean_forecasts_testset.rda")
save(cv_results, file="./VAR_fold_settings.rda")



     