require(data.table)
require(ProbCast)
require(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

zones <- c(1:10)
model <- "VAR"
horizons <- c(1:6)

skill_scores <- CJ(Zone=zones, Horizon=c(1:6))

for (i in zones){
  load(file=paste0("../persistence/zone", i, "_forecasts.rda"))
  if (model=="VAR"){
    load(file=paste0("../VAR/final_quantiles_zone",i,".rda"))
    assign('model_fcs', VAR_fcs)
    rm(VAR_fcs)
  }else if(model=="MC"){
    load(file=paste0("../MC/test_forecasts_zone",i,".rda"))
    assign('model_fcs', MC_fcs)
    rm(MC_fcs)
  }else if(model=="EMD"){
    load(file=paste0("../EMD_v2/final_quantiles_zone",i,".rda"))
    assign('model_fcs', EMD_fcs)
    rm(EMD_fcs)
  }
  
  persistence_times <- persistence_fcs[[2]][,.(issue_time, target_time, Horizon)]
  model_times <- model_fcs[[2]][,.(issue_time, target_time, Horizon)]
  
  mutual_times <- merge(persistence_times, model_times)
  
  for (h in horizons){
    persistence_indices <- match(mutual_times[Horizon==h, target_time], persistence_times[Horizon==h, target_time])
  }
  
  
}


model_skillscores <- CJ(Horizon=horizons, type=c("mean", "min", "max"))
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
      persistence_indices <- match(sample_times, persistence_details_list[[z]][Horizon==h, target_time])
      model_indices <- match(sample_times, model_details_list[[z]][Horizon==h, target_time])
      
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
  model_skillscores[(Horizon==h)&(type=="min"), MAE_ss := min(h_mae_ss)]
  model_skillscores[(Horizon==h)&(type=="max"), MAE_ss := max(h_mae_ss)]
  
  model_skillscores[(Horizon==h)&(type=="mean"), RMSE_ss := mean(h_rmse_ss)]
  model_skillscores[(Horizon==h)&(type=="min"), RMSE_ss := min(h_rmse_ss)]
  model_skillscores[(Horizon==h)&(type=="max"), RMSE_ss := max(h_rmse_ss)]
  
  model_skillscores[(Horizon==h)&(type=="mean"), pb_ss := mean(h_pb_ss)]
  model_skillscores[(Horizon==h)&(type=="min"), pb_ss := min(h_pb_ss)]
  model_skillscores[(Horizon==h)&(type=="max"), pb_ss := max(h_pb_ss)]
}

model_skillscores[, Model := paste0(model)]  
fwrite(model_skillscores, file=paste0("./",model,"_skillscores_horizons.csv"))
# model_skillscores <- fread(file=paste0("./",model,"_skillscores_horizons.csv"))

require(ggplot2)
ggplot(data=model_skillscores[type=='mean'], aes(x=Horizon, y=pb_ss)) +
  geom_line()

ggplot(data=model_scores, aes(x=Horizon, y=avPinball, colour=Model)) +
  geom_line() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2)
