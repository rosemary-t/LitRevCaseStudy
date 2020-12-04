require(data.table)
require(expm)
require(stats)
require(ProbCast)

## define all the model parameters
#nlevels <- seq(60,90,10)   # possible number of power levels to test
nlevels <- 70
cs <- seq(500, 1000, 100)
split_date <- as.POSIXct("2012-12-31 12:00") # every target time up to and including split_date is the training data.
horizons <- c(1:6) # forecast 1 to 6 hours ahead (1:6 steps ahead since it's hourly resolution data)
training_len_options <- 7000
#training_len_options <- seq(5500,7000,500) # the set of window lengths to test over
max_wl = max(training_len_options) + max(horizons)  # the maximum window length we will test over.
quantiles <- seq(0.05, 0.95, 0.05)
qnames <- paste0("q",quantiles*100)


## load the data
zone <- 1
load(paste0("C:\\Users\\rosemaryt\\ShareFile\\Personal Folders\\Data\\GEFcom2014data\\data_zone", zone, ".rda"))
assign("zdata", get(paste0("data_zone", zone)))
rm(list=paste0("data_zone", zone))
zdata <- zdata[,.(TARGETVAR, timestamp)] # just keep the power (not using the NWP forecasts)


## define functions:
## one to do the prediction
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

## and one to return the pinball for a given set of parameters
## the parameters are stored in a vector c(nlevels, window_len, C)
avg_pinball <- function(params, powerdata, windowtimes){
  nl <- params[1] # nlevels
  wlen <- params[2] # window length
  C <- params[3] # relative weight for counts data.
  
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
h <- 1 # for now.. we need to do this whole process for all horizons eventually.... or just use optimal settings from h=1 for all other horizons too


## data.table to store pinball loss for each (nlevels, k) combination
MC_settings_results <- CJ(Cval=cs, Nlevels=nlevels, windowlen=training_len_options)


for (nl in nlevels){
  ### discretise data
  print (nl)

  for (C in cs){
    for (WL in training_len_options){
      pbloss <- avg_pinball(c(nl, WL, C), train_data, window_times)
      MC_settings_results[(Cval==C)&(Nlevels==nl)&(windowlen==WL), Pinball := pbloss]
    }
  }
}



#save(MC_settings_results, file=paste0("C:\\Users\\rosemaryt\\ShareFile\\Personal Folders\\reproducing short term forecasting papers\\Markov Chain\\MCsettings_zone",zone,"param_range_2.rda"))
#MC_settings_results_1 <- copy(MC_settings_results)
#zone <- 1
load(file=paste0("C:\\Users\\rosemaryt\\ShareFile\\Personal Folders\\reproducing short term forecasting papers\\Markov Chain\\MCsettings_zone",zone,"param_range_2.rda"))
MC_settings_results <- rbind(MC_settings_results, MC_settings_results_1)

## param_range_1 is nlevels <- seq(20,50,10) ; cs <- seq(80, 160, 10) ; training_len_options <- seq(1000,5000,500)
## param_range_2 is nlevels <- seq(60,90,10) ; cs <- seq(170, 250, 10) ; training_len_options <- seq(5500,9500,500)
# opt_params <- MC_settings_results[ , .SD[which.min(Pinball)]]
# 
plot(MC_settings_results[(Nlevels==70) & (windowlen==7000), Cval], MC_settings_results[(Nlevels==70) & (windowlen==7000),Pinball])
# plot(MC_settings_results[(Cval==250) & (windowlen==7000), Nlevels], MC_settings_results[(Cval==250) & (windowlen==7000),Pinball])
# plot(MC_settings_results[(Cval==250) & (Nlevels==70), windowlen], MC_settings_results[(Cval==250) & (Nlevels==70),Pinball])
