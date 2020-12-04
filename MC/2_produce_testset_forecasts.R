require(data.table)
require(expm)
require(stats)
require(ProbCast)

## define all the model parameters
split_date <- as.POSIXct("2012-12-31 12:00") # every target time up to and including split_date is the training data.
horizons <- c(1:6) # forecast 1 to 6 hours ahead (1:6 steps ahead since it's hourly resolution data)
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
  ## timeseries must be a data.table with a column named 'ts'.
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

## load the table of pinball scores for different model parameters
load(file=paste0("C:\\Users\\rosemaryt\\ShareFile\\Personal Folders\\reproducing short term forecasting papers\\Markov Chain\\MCsettings_zone",zone,".rda"))
opt_params <- MC_settings_results[ , .SD[which.min(Pinball)]]  # the best (c, nlevels) combination
opt_nl <- opt_params$Nlevels
opt_c <- opt_params$Cval
opt_wl <- opt_params$windowlen

# make a list of windows whose issue_time (end point) falls in the test data set
first_window <- zdata[timestamp==(split_date+1*60*60), timestamp] # one after 'split date' as this was included in the training data.
last_window <- zdata[(dim(zdata)[1] - max(horizons)), timestamp]
window_times <- zdata[(timestamp >= first_window) & (timestamp <= last_window), timestamp] # the timestamps for the end of all possible sliding windows in the training data.

# discretise data into nlevels
levels <- seq(from=0, to=1, length.out=opt_nl) # between 0 and 1 because TARGETVAR is already normalised
zdata[, ts := findInterval(zdata[,TARGETVAR], levels)] # column 'ts' is the discretised time series.


all_test_forecasts <- data.table() # to store all quantile forecasts for all horizons.
other_info <- data.table() # to also store issue time, target time, actual power values.

for (h in horizons){
  horizon_forecasts <- data.table(issue_time=window_times, target_time=window_times+h*60*60, ActualPower=zdata[timestamp %in% (window_times + h*60*60), TARGETVAR])
  for (wt in window_times){
    window_start <- wt - (opt_wl-1)*60*60
    windowdata <- zdata[(timestamp >= window_start) & (timestamp<= wt)]
    Inputstate <- zdata[timestamp == wt, ts] # power level of the input to the markov chain.
    actual_state <- zdata[timestamp == (wt + h*60*60), ts]
    
    if ((!is.na(Inputstate)) & (!is.na(actual_state))){
      # both the forecast input and the actual value are NOT NA, so we can generate a forecast.
      forecast_probs <- prediction(windowdata, opt_nl, h, Inputstate, opt_c) # calculate transition probabilities
      
      
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
      
      
      horizon_forecasts[issue_time == wt, (qnames) := as.list(quantile_fc$y)]
    }
  }
  
  ## now horizon_forecasts is filled in, save quantile forecasts into all_test_forecasts and other data into other_info
  all_test_forecasts <- rbind(all_test_forecasts, horizon_forecasts[,..qnames])
  other_info <- rbind(other_info, horizon_forecasts[,.(issue_time, target_time, ActualPower)])
}

na_quantiles <- which(!complete.cases(all_test_forecasts))
na_powers <- which(is.na(other_info[, ActualPower]))
na_rows <- unique(c(na_quantiles, na_powers)) # indices of rows with NA either in quantile forecasts or in actual power.
complete_rows <- setdiff(c(1:dim(all_test_forecasts)[1]), na_rows)

all_test_forecasts <- all_test_forecasts[complete_rows]
other_info <- other_info[complete_rows]

other_info[, Horizon := (target_time-issue_time)]
class(all_test_forecasts) <- c("MultiQR", class(all_test_forecasts))
MC_fcs <- list(all_test_forecasts, other_info)
save(MC_fcs, file=paste0("C:\\Users\\rosemaryt\\ShareFile\\Personal Folders\\reproducing short term forecasting papers\\Markov Chain\\test_forecasts_zone",zone,".rda"))


