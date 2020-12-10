require(data.table)
require(expm)
require(stats)
require(ProbCast)
require(doParallel)
require(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

## define all the model parameters
nlevels <- seq(80,120,10)   # possible number of power levels to test
cs <- seq(400, 1200, 200)   # weighting given to data (over prior)
split_date <- as.POSIXct("2012-12-31 12:00") # every target time up to and including split_date is the training data.
horizons <- c(1:6) # forecast 1 to 6 hours ahead (1:6 steps ahead since it's hourly resolution data)
training_len_options <- seq(5000,7000,500) # the set of window lengths to test over
max_wl = max(training_len_options) + max(horizons)  # the maximum window length we will test over.
Quantiles <- seq(0.05, 0.95, 0.05)
Qnames <- paste0("q",Quantiles*100)


## load the data
zone <- 10
load(paste0("../Data/data_zone", zone, ".rda"))
assign("zdata", get(paste0("data_zone", zone)))
rm(list=paste0("data_zone", zone))
zdata <- zdata[,.(TARGETVAR, timestamp)] # just keep the power (not using the NWP forecasts)


## define functions:
## one to do the prediction (return a probability for each power level)
prediction <- function(timeseries, nlevels, horizon, inputstate, c){
  ## timeseries must be a data.table with a single column named 'ts'.
  ## c is the strength of the counts information relative to the prior - i.e. c=2 means counts gets twice as much weight as prior.
  ## function returns a data.table with all possible power states and their forecast probabilities, given the input state.

  timeseries[,ts_shift := shift(ts, n=horizon, type="lead")]
  timeseries <- na.omit(timeseries) # delete rows with NAs
  
  ## we only care about the row in the transition matrix corresponding to the known forecast input for this particular rolling window:
  tm_row <- timeseries[ts==inputstate]
  
  counts <- tm_row[, .N, by=.(ts, ts_shift)]

  
  ## now all input states are present, check if there are any output states missing too.
  missing_outputstates <- setdiff(c(1:nlevels), counts[,ts_shift])
  if (length(missing_outputstates) != 0){
    extra_rows <- data.table(ts=inputstate, ts_shift = missing_outputstates, N=0) # N=0 means no transitions to these states.
    counts <- rbind(counts, extra_rows)
  }
  
  ## check the counts matrix contains a row for each output state now.
  if (dim(counts)[1] != nlevels){
    print ("counts data table doesn't have the right number of rows")
    return (NA)
  }else{
    counts[,alpha := nlevels-abs(counts$ts-counts$ts_shift)] # add in a column for the alpha parameters of dirichlet prior
    counts[,unnorm_probs := c*N+alpha-1] # the state probabilities are given by k*N+alpha-1 (before normalisation of the vector)
    counts[,probs := unnorm_probs/sum(unnorm_probs)] # normalise the probabilities.
    
    result <- counts[,.(ts_shift, probs)]
    result <- result[order(ts_shift)]
    
    return (result)
  }
}

## and one to return the pinball for a given set of parameters, from the i th row of MCsettings
avg_pinball <- function(MCsettings, i, powerdata, windowtimes, quantiles){
  require(data.table)
  require(expm)
  require(stats)
  require(ProbCast)
  
  qnames <- paste0("q",quantiles*100)
  nl <- MCsettings[i, Nlevels] # nlevels
  wlen <- MCsettings[i, windowlen] # window length
  C <- MCsettings[i, Cval] # relative weight for counts data.
  
  h <- 1 # just evaluate for a one step ahead horizon.
  levels <- seq(from=0, to=1, length.out=nl) # between 0 and 1 because TARGETVAR is already normalised
  powerdata[, ts := findInterval(powerdata[,TARGETVAR], levels)] # column 'ts' is the discretised time series.
  
  ## empty data.table to record quantile forecasts for each sliding window:
  window_q_forecasts <- data.table('issue_time' = windowtimes, 'target_time' = windowtimes + h*60*60) 
  
  ## iterate over each sliding window
  for (wt in windowtimes){
    window_start <- wt - (wlen-1)*60*60
    windowdata <- powerdata[(timestamp >= window_start) & (timestamp<= wt),]
    Inputstate <- powerdata[timestamp == wt, ts] # power level of the input to the markov chain.
    actual_value <- powerdata[timestamp == (wt + h*60*60), ts]
    
    if ((!is.na(Inputstate)) & (!is.na(actual_value))){
      # both the forecast input and the actual value are NOT NA, so we can generate a forecast.
      forecast_probs <- prediction(windowdata, nl, h, Inputstate, C) # calculate transition probabilities
      
      
      if (all(is.na(forecast_probs))){
        print ('All the quantile forecasts are NA')
      }else{
        forecast_probs[,powers := levels]
        forecast_probs[,CDF := cumsum(probs)]
        
        ## now just have to turn these into quantiles......
        quantile_fc <- approx(x=forecast_probs[,CDF], y=forecast_probs[,powers], xout=quantiles)
        
        ## fill in the NAN quantiles if there are any (high or low powers so quantiles are out of range of cdf points.)
        if (min(forecast_probs[,CDF]) > min(quantiles)){
          ## there will be NAs at low powers in the interpolation
          NonNAindex <- which(!is.na(quantile_fc$y))
          firstNonNA <- min(NonNAindex)
          quantile_fc$y[c(1:(firstNonNA-1))] <- 0 # replace the lower power NAs with 0
        }
        
        
        if (max(forecast_probs[,CDF]) < max(quantiles)){
          ## there will be NAs at high powers in the interpolation
          NonNAindex <- which(!is.na(quantile_fc$y))
          lastNonNA <- max(NonNAindex)
          quantile_fc$y[c((lastNonNA+1): length(quantile_fc$y))] <- 1
          
        }
      }
      
      
      window_q_forecasts[issue_time == wt, (qnames) := as.list(quantile_fc$y)]
      window_q_forecasts[issue_time == wt, 'actual' := powerdata[timestamp == wt+h*60*60, TARGETVAR]]
    } 
  }
  
  ## now get quantile forecasts into a MultiQR object and evaluate Pinball.
  
  if (!any(is.na(window_q_forecasts))){
    ## only make the pinball score if we have complete data (ie nothing has gone wrong)
    quantilevals <- window_q_forecasts[,..qnames]
    class(quantilevals) <- c("MultiQR", class(quantilevals))
    pbscore <- pinball(quantilevals, realisations=window_q_forecasts$actual, plot.it=F)
    return(mean(pbscore$Loss))
  }else{
    return (NA)
  }
}

#############################################################################################

## rolling windows are indexed by the forecast issue time (the timestamp for the last point in the training window)
train_data <- zdata[timestamp <= split_date]
first_window <- train_data[max(training_len_options), timestamp]
last_window <- train_data[(dim(train_data)[1] - max(horizons)), timestamp]
window_times <- train_data[(timestamp >= first_window) & (timestamp <= last_window), timestamp] # the timestamps for the end of all possible sliding windows in the training data.


## optimise (nlevels, c, window_len) simultaneously
h <- horizons[1] # optimise parameters for one step ahead forecast.

#test <- avg_pinball(MC_settings_results, 2, train_data, window_times, Quantiles)

## data.table to store pinball loss for each (nlevels, k) combination
MC_settings_results <- CJ(Cval=cs, Nlevels=nlevels, windowlen=training_len_options)

## test each parameter combination in parallel
cores <- detectCores()
cl <- makeCluster(cores-2)
registerDoParallel(cl)
start_time <- Sys.time()
pbscores <- foreach(i = c(1:dim(MC_settings_results)[1])) %dopar% {
  avg_pinball(MC_settings_results, i, train_data, window_times, Quantiles)
}
print (Sys.time() - start_time)
stopCluster(cl)

pbscores <- unlist(pbscores)
MC_settings_results[, Pinball := pbscores]

save(MC_settings_results, file=paste0("./MCsettings_zone",zone,".rda"))

opt_params <- MC_settings_results[ , .SD[which.min(Pinball)]]

plot(MC_settings_results[(Nlevels==opt_params$Nlevels) & (windowlen==opt_params$windowlen), Cval], MC_settings_results[(Nlevels==opt_params$Nlevels) & (windowlen==opt_params$windowlen),Pinball], main="Cval")
plot(MC_settings_results[(Cval==opt_params$Cval) & (windowlen==opt_params$windowlen), Nlevels], MC_settings_results[(Cval==opt_params$Cval) & (windowlen==opt_params$windowlen),Pinball], main="Nlevels")
plot(MC_settings_results[(Cval==opt_params$Cval) & (Nlevels==opt_params$Nlevels), windowlen], MC_settings_results[(Cval==opt_params$Cval) & (Nlevels==opt_params$Nlevels),Pinball], main="windowlen")
