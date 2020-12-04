require(data.table)
require(ggplot2)
require(frbs)
require(stats)
require(tseries)


## to find the optimal learning rate (gradient step) and/or maximum number of iterations.



horizon <- 1 # just one step ahead forecast for now.


fileloc <- 'C:\\Users\\rosemaryt\\ShareFile\\Personal Folders\\reproducing short term forecasting papers\\zone1_w100_rollingwindows'
IMFs <- fread(paste0(fileloc,'\\window0.csv'))
#IMFs[,timestamp:=as.POSIXct(timestamp,format="%d/%m/%Y %H:%M",tz="UTC")]
IMFs[,timestamp:=as.POSIXct(V1,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
IMFs[,V1:= NULL]
nimfs <- dim(IMFs)[2] - 2 # minus two as the timestamp column doesn't count.

splitpt <- dim(IMFs)[1]/2 # halfway through the data points; use half for training and half for testing.
train_IMFs <- IMFs[c(1:splitpt)]
test_IMFs <- IMFs[c((splitpt+1):dim(IMFs)[1])]

maxs <- lapply(IMFs, max)
mins <- lapply(IMFs, min)

i <- 1 ## just use IMF 1.
imf_no <- paste0("imf", i)
## get the training/testing data and also normalise it.
train_imf <- (train_IMFs[,..imf_no] - mins[[imf_no]])/(maxs[[imf_no]]-mins[[imf_no]])
test_imf <- (test_IMFs[,..imf_no] - mins[[imf_no]])/(maxs[[imf_no]]-mins[[imf_no]])

## create anfis input arrays
lags <- c(1:3) # specify which lags to use as ANFIS model inputs.

## create the ANFIS inputs (from already normalised data)
createinputs <- function(input_ts, horizons){
  for (j in c(1:length(horizons))){
    ## for each  horizon, train a model, predict and save the prediction in 'data'
    h <- horizons[j]
    shifts <- lags + h - 1
    normdata <- input_ts #initialise the df where all the lags will be added.
    
    #shift the time series to create lags as model inputs
    for (s in shifts){
      lag <- shift(input_ts,n=s, fill=NA, type="lag", give.names=TRUE)
      lag_df <- as.data.table(lag)
      normdata <- cbind(lag_df, normdata) # put normdata second so the imfs remain in the last column (this is taken as the target (output) variable in the anfis model)
    }
  }
  return(normdata)
  
}


anfis_train <- createinputs(train_imf, horizon)
anfis_train <- anfis_train[complete.cases(anfis_train),] # remove rows with any NAs.
anfis_test <- createinputs(test_imf[,get(imf_no)], horizon)
anfis_test <- anfis_test[complete.cases(anfis_test),] # remove rows with any NAs.
actuals <- anfis_test[, normdata]
anfis_test[, normdata :=NULL]
anfis_test[, timestamp :=NULL]




max_iters <- seq(2,20,by=2)
## initialise empty lists to put results in.
bs_nmae <-numeric(length(max_iters)) ## bootstrapped nmae of anfis model
bs_sterr <- numeric(length(max_iters)) ## bootstrapped standard error of anfis nmae. 



for (mi in max_iters){

  ## now run anfis:
  start_time <- Sys.time()
  print ("ANFIS train")
  model <- ANFIS(anfis_train, max.iter=2, step.size=0.01) # this line won't run and frbs.predict won't accept a max.iter argument.
  anfis <- frbs.learn(anfis_train, method.type=ANFIS(anfis_train, max.iter=2))
  print (Sys.time()-start_time)
  
  ans <- predict(anfis, anfis_test)
  anfis_residuals <- actuals - ans
  anfis_residuals <- as.data.table(anfis_residuals)
  anfis_residuals <- anfis_residuals[,.(anfis_resid=V1)] # rename V1.
  anfis_residuals[,timestamp := anfis_test_ts]
  
  result <- merge.data.table(test_imf, anfis_residuals, by="timestamp")
  result[,diff := abs(result$arima_resid) - abs(result$anfis_resid)]
  
  sarima_nmae[i] <- mean(abs(result$arima_resid))
  anfis_nmae[i] <- mean(abs(result$anfis_resid))
  
  bs <- tsbootstrap(result$diff, nb=1000, statistic = mean, m=1)
  bs_nmae[i] <- mean(bs[["statistic"]])
  bs_sterr[i] <- as.numeric(bs[["se"]])
}
  
  
  
  
  
#write.csv(results, file="C:\\Users\\rosemaryt\\ShareFile\\Personal Folders\\reproducing short term forecasting papers\\zone1_w100_rollingwindows\\results.csv")
