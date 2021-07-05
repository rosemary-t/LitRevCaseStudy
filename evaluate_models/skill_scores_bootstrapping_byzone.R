require(data.table)
require(ProbCast)
require(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

zones <- c(1:10)
model <- "MC"
Nboots <- 1000 # number of times to bootstrap resample
bounds <- c(0.025,0.975) # quantiles of bootstraped metric to plot. If (0,1) will use min and max observed pinball scores.


for (i in zones){
  print (i)
  
  load(file=paste0("../persistence/zone", i, "_forecasts.rda"))
  persistence_fcs[[2]][, id := c(1:dim(persistence_fcs[[2]])[1])] # add an row index count
  persistence_quantiles <- persistence_fcs[[1]]
  persistence_details <- persistence_fcs[[2]]
  rm(persistence_fcs)
  
  mutual_times <- persistence_details[,.(issue_time, target_time, Horizon)]
  
  if (model=="VAR"){
    load(file=paste0("../VAR/final_quantiles_zone",i,".rda"))
    VAR_fcs[[2]][, id := c(1:dim(VAR_fcs[[2]])[1])] # add an row index count
    model_quantiles <- VAR_fcs[[1]]
    model_details <- VAR_fcs[[2]]
    rm(VAR_fcs)
  }else if (model=="EMD"){
    load(file=paste0("../EMD/final_quantiles_zone",i,".rda"))
    EMD_fcs[[2]][, id := c(1:dim(EMD_fcs[[2]])[1])] # add an row index count
    model_quantiles <- EMD_fcs[[1]]
    model_details <- EMD_fcs[[2]]
    rm(EMD_fcs)
  }else if (model=="MC"){
    load(file=paste0("../MC/test_forecasts_zone",i,".rda"))
    MC_fcs[[2]][, id := c(1:dim(MC_fcs[[2]])[1])] # add an row index count
    model_quantiles <- MC_fcs[[1]]
    model_details <- MC_fcs[[2]]
    rm(MC_fcs)
  }else{print("no model forecasts loaded")}
  
  horizons <- unique(persistence_details$Horizon)
  mutual_times <- merge(mutual_times, model_details[,.(issue_time, target_time,Horizon)])
  
  model_skillscores <- CJ(Horizon=horizons, type=c("mean", "lower", "upper"))
  
  for (h in horizons){
    print (h)
    horizon_times <- mutual_times[Horizon==h, target_time]
    h_mae_ss <- rep(0,Nboots)
    h_rmse_ss <- rep(0,Nboots)
    h_pb_ss <- rep(0,Nboots)
    
    for (nb in c(1:Nboots)){
      sample_times <- sample(horizon_times, length(horizon_times), replace=TRUE)
      
      ## now calculate skill score for each zone and save in relevant matrices
      persistence_h_details <- persistence_details[Horizon==h]
      model_h_details <- model_details[Horizon==h]
      
      persistence_short_indices <- match(sample_times, persistence_h_details[, target_time])
      persistence_indices <- persistence_h_details[persistence_short_indices, id]
      
      model_short_indices <- match(sample_times, model_h_details[, target_time])
      model_indices <- model_h_details[model_short_indices, id]
      
      ## calculate MAE and RMSE with q50
      ## for persistence
      p_errors <- persistence_quantiles[persistence_indices, q50] - persistence_details[persistence_indices, ActualPower]
      p_MAE <- mean(abs(p_errors))
      p_RMSE <- sqrt(mean(p_errors^2))
      
      ## and the model
      model_errors <- model_quantiles[model_indices, q50] - model_details[model_indices, ActualPower]
      model_MAE <- mean(abs(model_errors), na.rm=TRUE)
      model_RMSE <- sqrt(mean(model_errors^2, na.rm=TRUE))
      
      ## enter skill score into matrices
      h_mae_ss[nb] <- 1-model_MAE/p_MAE
      h_rmse_ss[nb] <- 1-model_RMSE/p_RMSE
      
      ## and then pinball score too.
      p_pinball <- pinball(persistence_quantiles[persistence_indices], realisations = persistence_details[persistence_indices, ActualPower], plot.it=FALSE)
      model_pinball <- pinball(model_quantiles[model_indices], realisations = model_details[model_indices, ActualPower], plot.it=FALSE)
      h_pb_ss[nb] <- 1 - mean(model_pinball$Loss)/mean(p_pinball$Loss)
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
  model_skillscores[, Zone := paste0(i)]
  fwrite(model_skillscores, file=paste0("./ss_by_zone/",model,"_zone",i,"skillscores_horizons.csv"))


}

# require(ggplot2)
# ggplot(data=model_skillscores[type=='mean'], aes(x=Horizon, y=pb_ss)) +
#   geom_line()

# ggplot(data=model_scores, aes(x=Horizon, y=avPinball, colour=Model)) +
#   geom_line() +
#   geom_errorbar(aes(ymin=lower, ymax=upper), width=.2)
