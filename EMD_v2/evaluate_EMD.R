require(data.table)
require(ProbCast)
require(ggplot2)
require(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

## choose between constant variance and exponential smoothing model for the variance, 
## for all sites and all horizons

zone <- 1 # what zone do we want to evaluate forecasts for?
Nboots <- 1000 # number of times to bootstrap resample.
quantile_bounds <- c(0.025,0.975) # quantiles of bootstraped metric to plot. If (0,1) will use min and max observed pinball scores.


load(paste0("./zone",zone,"_qs_constvar.rda")) # load quantile forecasts from fitting a normal distribution
assign('const_var', EMD_fcs)
rm(EMD_fcs)

load(paste0("./zone",zone,"_qs_expsmoothvar.rda")) # load quantile forecasts from fitting a normal distribution
assign('expsm_var', EMD_fcs)
rm(EMD_fcs)

horizons <- unique(const_var[[2]]$Horizon)
## plot pinball loss per horizon for each model
# get the skill scores
model_scores <- CJ(Horizon = horizons, Model=c("constvar", "expsmoothvar"))
bootsamples <- list()

for (h in horizons){
  print (h)
  # vectors to store bootstrap sampled pinball scores in
  constvar_pbsamples <- numeric(Nboots)
  expvar_pbsamples <- numeric(Nboots)

  const_times <- sort(const_var[[2]][Horizon==h, target_time])
  exp_times <- sort(expsm_var[[2]][Horizon==h, target_time])

  mutual_times <- const_times[const_times %in% exp_times]

  for (nb in c(1:Nboots)){
    # print (nb)
    sample_times <- sample(mutual_times, length(mutual_times), replace=TRUE)
    const_indices <- const_var[[2]][(Horizon==h) & (target_time %in% sample_times), which=TRUE]
    exp_indices <- expsm_var[[2]][(Horizon==h) & (target_time %in% sample_times), which=TRUE]

    const_pb <- pinball(const_var[[1]][const_indices], realisations = const_var[[2]][const_indices, ActualPower], plot.it=FALSE)
    constvar_pbsamples[nb] <- mean(const_pb$Loss)

    exp_pb <- pinball(expsm_var[[1]][exp_indices], realisations = expsm_var[[2]][exp_indices, ActualPower], plot.it=FALSE)
    expvar_pbsamples[nb] <- mean(exp_pb$Loss)
  }

  ## now record the mean and quantiles of bootstrap samples
  model_scores[(Horizon==h) & (Model=="constvar"), avPinball := mean(constvar_pbsamples)]
  model_scores[(Horizon==h) & (Model=="expsmoothvar"), avPinball := mean(expvar_pbsamples)]

  if (quantile_bounds[1] == 0){
    ## just uses the min value for the error bar
    model_scores[(Horizon==h) & (Model=="constvar"), lower := min(constvar_pbsamples)]
    model_scores[(Horizon==h) & (Model=="expsmoothvar"), lower := min(expvar_pbsamples)]
  }else{
    model_scores[(Horizon==h) & (Model=="constvar"), lower := quantile(constvar_pbsamples, quantile_bounds[1])]
    model_scores[(Horizon==h) & (Model=="expsmoothvar"), lower := quantile(expvar_pbsamples, quantile_bounds[1])]
  }

  if (quantile_bounds[2] ==1){
    ## just uses the max observed value
    model_scores[(Horizon==h) & (Model=="constvar"), upper := max(constvar_pbsamples)]
    model_scores[(Horizon==h) & (Model=="expsmoothvar"), upper := max(expvar_pbsamples)]
  }else{
    model_scores[(Horizon==h) & (Model=="constvar"), upper := quantile(constvar_pbsamples, quantile_bounds[2])]
    model_scores[(Horizon==h) & (Model=="expsmoothvar"), upper := quantile(expvar_pbsamples, quantile_bounds[2])]
  }
}

ggplot(data=model_scores, aes(x=Horizon, y=avPinball, colour=Model)) +
  geom_line() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2)


## now we need to combine the quantile forecasts from the correct variance model for each horizon, to get a 'final' set:
# which_variance <- fread('../evaluate_models/const_vs_expsmooth.csv')
# const_horizons <- which_variance[(Zone==zone)&(EMD_opt=="const"), Horizon]
# exp_horizons <- which_variance[(Zone==zone)&(EMD_opt=="exp"), Horizon]
# 
# const_quantile_indices <- const_var[[2]][(Horizon %in% const_horizons), which=TRUE]
# exp_quantile_indices <- expsm_var[[2]][(Horizon %in% exp_horizons), which=TRUE]
# 
# new_quantiles <- rbind(const_var[[1]][const_quantile_indices], expsm_var[[1]][exp_quantile_indices])
# class(new_quantiles) <- c("MultiQR", class(new_quantiles))
# new_otherinfo <- rbind(const_var[[2]][const_quantile_indices], expsm_var[[2]][exp_quantile_indices])
# 
# EMD_fcs <- list(new_quantiles, new_otherinfo)
# save(EMD_fcs, file=paste0('./final_quantiles_zone',zone,'.rda'))
