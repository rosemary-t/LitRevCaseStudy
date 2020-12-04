require(data.table)
require(ProbCast)
require(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

split_date <- as.POSIXct("2012-12-31 12:00") # every target time up to and including split_date is the training data.
horizons <- c(1:6) # forecast 1 to 6 hours ahead (1:6 steps ahead since it's hourly resolution data)
quantiles <- seq(0.05, 0.95, 0.05)
qnames <- paste0("q",quantiles*100)

## load the data
zone <- 1
load(paste0("../Data/data_zone", zone, ".rda"))
assign("zdata", get(paste0("data_zone", zone)))
rm(list=paste0("data_zone", zone))
zdata <- zdata[,.(TARGETVAR, timestamp)] # just keep the power (not using the NWP forecasts)
mean_forecasts <- data.table() # store the point persistence forecast for the test times, for each horizon


for (h in horizons){
  horizonforecasts <- zdata[,.(issue_time=timestamp, forecast=TARGETVAR)]
  horizonforecasts[, target_time := issue_time + h*60*60 ]
  horizonforecasts[, ActualPower := shift(forecast, n=h, type='lead')]
  horizonforecasts_train <- horizonforecasts[issue_time <= split_date,]
  horizonforecasts_train[, residual := forecast - ActualPower]
  horizonforecasts_test <- horizonforecasts[issue_time > split_date]
  horizonforecasts_test[, SD := sd(horizonforecasts_train$residual)]
  mean_forecasts <- rbind(mean_forecasts, horizonforecasts_test)

}

mean_forecasts[, Horizon := (target_time - issue_time)]

quantilefcs <- as.data.table(sapply(quantiles, function(x){qnorm(x, mean=mean_forecasts$forecast, sd=mean_forecasts$SD)}))
names(quantilefcs) <- qnames
quantilefcs[, (qnames) := lapply(.SD, function(x){ifelse(x <= 0,  0, ifelse(x >= 1, 1, x))}), .SDcols=qnames] # clip to [0,1] range

## drop rows with NA in either the quantiles or ActualPower
NAquantilerows <- which(!complete.cases(quantilefcs)) # get a list of rows where the quantiles are NA
NApowerrows <- which(is.na(mean_forecasts$ActualPower))
NArows <- unique(sort(c(NAquantilerows, NApowerrows)))
keeprows <- setdiff(c(1:dim(quantilefcs)[1]), NArows) # indices of rows we do want to keep.
quantilefcs <- quantilefcs[keeprows]
mean_forecasts <- mean_forecasts[keeprows]

class(quantilefcs) <- c("MultiQR", class(quantilefcs))

other_info <- mean_forecasts[, .(issue_time, target_time, Horizon, ActualPower)]
persistence_fcs <-list(quantilefcs,other_info) 
save(persistence_fcs, file=paste0("zone",zone,"_forecasts.rda"))
