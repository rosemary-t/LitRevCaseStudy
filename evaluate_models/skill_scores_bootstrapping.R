require(data.table)
require(ProbCast)
require(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

zones <- c(1:10)
model <- "EMD"
Nboots <- 1000 # number of times to bootstrap resample
bounds <- c(0.025,0.975) # quantiles of bootstraped metric to plot. If (0,1) will use min and max observed pinball scores.


## load persistence model forecasts
persistence_quantiles_list <- vector(mode = "list", length = length(zones))
persistence_details_list <- vector(mode = "list", length = length(zones))
for (i in zones){
  load(file=paste0("../persistence/zone", i, "_forecasts.rda"))
  persistence_fcs[[2]][, id := c(1:dim(persistence_fcs[[2]])[1])] # add an row index count
  persistence_quantiles_list[[i]] <- persistence_fcs[[1]]
  persistence_details_list[[i]] <- persistence_fcs[[2]]
  
  if (i == 1){
    ## initialise a record of the (issue_time, target_time) pairs that exist for ALL zones
    persistence_times_mutual <- persistence_fcs[[2]][,.(issue_time, target_time, Horizon)]
  }else{
    zone_rows <- persistence_fcs[[2]][,.(issue_time, target_time, Horizon)]
    persistence_times_mutual <- merge(persistence_times_mutual, zone_rows) # keep the rows from mutual_rows, that also appear in zone_rows
  }
  # print (dim(persistence_times_mutual))
}
horizons <- unique(persistence_details_list[[1]]$Horizon)


## load model forecasts
model_quantiles_list <- vector(mode = "list", length = length(zones))
model_details_list <- vector(mode = "list", length = length(zones))
if (model=="VAR"){
  for (i in zones){
    load(file=paste0("../VAR/final_quantiles_zone",i,".rda"))
    VAR_fcs[[2]][, id := c(1:dim(VAR_fcs[[2]])[1])] # add an row index count
    model_quantiles_list[[i]] <- VAR_fcs[[1]]
    model_details_list[[i]] <- VAR_fcs[[2]]
    
    if (i == 1){
      ## initialise a record of the (issue_time, target_time) pairs that exist for ALL zones
      model_times_mutual <- VAR_fcs[[2]][,.(issue_time, target_time, Horizon)]
    }else{
      zone_rows <- VAR_fcs[[2]][,.(issue_time, target_time, Horizon)]
      model_times_mutual <- merge(model_times_mutual, zone_rows) # keep the rows from mutual_rows, that also appear in zone_rows
    }
  }
  rm(VAR_fcs)
}else if(model=="EMD"){
  for (i in zones){
    load(file=paste0("../EMD_v2/final_quantiles_zone",i,".rda"))
    EMD_fcs[[2]][, id := c(1:dim(EMD_fcs[[2]])[1])] # add an row index count
    model_quantiles_list[[i]] <- EMD_fcs[[1]]
    model_details_list[[i]] <- EMD_fcs[[2]]
    
    if (i == 1){
      ## initialise a record of the (issue_time, target_time) pairs that exist for ALL zones
      model_times_mutual <- EMD_fcs[[2]][,.(issue_time, target_time, Horizon)]
    }else{
      zone_rows <- EMD_fcs[[2]][,.(issue_time, target_time, Horizon)]
      model_times_mutual <- merge(model_times_mutual, zone_rows) # keep the rows from mutual_rows, that also appear in zone_rows
    }
  }
  rm(EMD_fcs)
}else if(model=="MC"){
  for (i in zones){
    load(file=paste0("../MC/test_forecasts_zone",i,".rda"))
    MC_fcs[[2]][, id := c(1:dim(MC_fcs[[2]])[1])] # add an row index count
    model_quantiles_list[[i]] <- MC_fcs[[1]]
    model_details_list[[i]] <- MC_fcs[[2]]
    
    if (i == 1){
      ## initialise a record of the (issue_time, target_time) pairs that exist for ALL zones
      model_times_mutual <- MC_fcs[[2]][,.(issue_time, target_time, Horizon)]
    }else{
      zone_rows <- MC_fcs[[2]][,.(issue_time, target_time, Horizon)]
      model_times_mutual <- merge(model_times_mutual, zone_rows) # keep the rows from mutual_rows, that also appear in zone_rows
    }
  }
  rm(MC_fcs)
}else{print("no model forecasts loaded")}

# now get the (issue_time, target_time) pairs we have both persistence and model forecasts for, for all zones
mutual_times <- merge(persistence_times_mutual, model_times_mutual)

model_skillscores <- CJ(Horizon=horizons, type=c("mean", "lower", "upper"))
## we are ready to bootstrap.
for (h in horizons){
  print (h)
  horizon_times <- mutual_times[Horizon==h, target_time]
  h_mae_ss <- matrix(nrow=Nboots, ncol=length(zones))
  h_rmse_ss <- matrix(nrow=Nboots, ncol=length(zones))
  h_pb_ss <- matrix(nrow=Nboots, ncol=length(zones))
  
  for (nb in c(1:Nboots)){
    sample_times <- sample(horizon_times, length(horizon_times), replace=TRUE)
    
    ## now calculate skill score for each zone and save in relevant matrices
    for (z in zones){
      persistence_h_details <- persistence_details_list[[z]][Horizon==h]
      model_h_details <- model_details_list[[z]][Horizon==h]
      
      persistence_short_indices <- match(sample_times, persistence_h_details[, target_time])
      persistence_indices <- persistence_h_details[persistence_short_indices, id]
      
      model_short_indices <- match(sample_times, model_h_details[, target_time])
      model_indices <- model_h_details[model_short_indices, id]
      
      ## calculate MAE and RMSE with q50
      ## for persistence
      p_errors <- persistence_quantiles_list[[z]][persistence_indices, q50] - persistence_details_list[[z]][persistence_indices, ActualPower]
      p_MAE <- mean(abs(p_errors))
      p_RMSE <- sqrt(mean(p_errors^2))
      
      ## and the model
      model_errors <- model_quantiles_list[[z]][model_indices, q50] - model_details_list[[z]][model_indices, ActualPower]
      model_MAE <- mean(abs(model_errors), na.rm=TRUE)
      model_RMSE <- sqrt(mean(model_errors^2, na.rm=TRUE))
      
      ## enter skill score into matrices
      h_mae_ss[nb,z] <- 1-model_MAE/p_MAE
      h_rmse_ss[nb,z] <- 1-model_RMSE/p_RMSE
      
      ## and then pinball score too.
      p_pinball <- pinball(persistence_quantiles_list[[z]][persistence_indices], realisations = persistence_details_list[[z]][persistence_indices, ActualPower], plot.it=FALSE)
      model_pinball <- pinball(model_quantiles_list[[z]][model_indices], realisations = model_details_list[[z]][model_indices, ActualPower], plot.it=FALSE)
      h_pb_ss[nb,z] <- 1 - mean(model_pinball$Loss)/mean(p_pinball$Loss)
    }
  }
  
  ## now fill in entries in skill_scores data.table
  model_skillscores[(Horizon==h)&(type=="mean"), MAE_ss := mean(h_mae_ss)]
  model_skillscores[(Horizon==h)&(type=="mean"), RMSE_ss := mean(h_rmse_ss)]
  model_skillscores[(Horizon==h)&(type=="mean"), pb_ss := mean(h_pb_ss)]
  
  if (bounds[1] == 0){
    model_skillscores[(Horizon==h)&(type=="lower"), MAE_ss := min(h_mae_ss)]
    model_skillscores[(Horizon==h)&(type=="lower"), RMSE_ss := min(h_rmse_ss)]
    model_skillscores[(Horizon==h)&(type=="lower"), pb_ss := min(h_pb_ss)]
  }else{
    model_skillscores[(Horizon==h)&(type=="lower"), MAE_ss := quantile(h_mae_ss, bounds[1])]
    model_skillscores[(Horizon==h)&(type=="lower"), RMSE_ss := quantile(h_rmse_ss, bounds[1])]
    model_skillscores[(Horizon==h)&(type=="lower"), pb_ss := quantile(h_pb_ss, bounds[1])]
  }
  
  if (bounds[2] == 1){
    model_skillscores[(Horizon==h)&(type=="upper"), MAE_ss := min(h_mae_ss)]
    model_skillscores[(Horizon==h)&(type=="upper"), RMSE_ss := min(h_rmse_ss)]
    model_skillscores[(Horizon==h)&(type=="upper"), pb_ss := min(h_pb_ss)]
  }else{
    model_skillscores[(Horizon==h)&(type=="upper"), MAE_ss := quantile(h_mae_ss, bounds[2])]
    model_skillscores[(Horizon==h)&(type=="upper"), RMSE_ss := quantile(h_rmse_ss, bounds[2])]
    model_skillscores[(Horizon==h)&(type=="upper"), pb_ss := quantile(h_pb_ss, bounds[2])]
  }
}

model_skillscores[, Model := paste0(model)]  
fwrite(model_skillscores, file=paste0("./",model,"_skillscores_horizons.csv"))
# model_skillscores <- fread(file=paste0("./",model,"_skillscores_horizons.csv"))

require(ggplot2)
ggplot(data=model_skillscores[type=='mean'], aes(x=Horizon, y=pb_ss)) +
  geom_line()

# ggplot(data=model_scores, aes(x=Horizon, y=avPinball, colour=Model)) +
#   geom_line() +
#   geom_errorbar(aes(ymin=lower, ymax=upper), width=.2)
