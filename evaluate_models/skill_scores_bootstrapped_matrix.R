require(data.table)
require(ProbCast)
require(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

zones <- c(1:10)
models <- c("Persistence", "VAR", "MC", "EMD")
Nboots <- 1000 # number of times to bootstrap resample


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


## load forecasts from each other model too
VAR_quantiles_list <- vector(mode = "list", length = length(zones))
VAR_details_list <- vector(mode = "list", length = length(zones))
MC_quantiles_list <- vector(mode = "list", length = length(zones))
MC_details_list <- vector(mode = "list", length = length(zones))
EMD_quantiles_list <- vector(mode = "list", length = length(zones))
EMD_details_list <- vector(mode = "list", length = length(zones))


for (i in zones){
  load(file=paste0("../VAR/final_quantiles_zone",i,".rda"))
  VAR_fcs[[2]][, id := c(1:dim(VAR_fcs[[2]])[1])] # add a row index count
  VAR_quantiles_list[[i]] <- VAR_fcs[[1]]
  VAR_details_list[[i]] <- VAR_fcs[[2]]
  
  load(file=paste0("../MC/test_forecasts_zone",i,".rda"))
  MC_fcs[[2]][, id := c(1:dim(MC_fcs[[2]])[1])] # add an row index count
  MC_quantiles_list[[i]] <- MC_fcs[[1]]
  MC_details_list[[i]] <- MC_fcs[[2]]
  
  load(file=paste0("../EMD_v2/final_quantiles_zone",i,".rda"))
  EMD_fcs[[2]][, id := c(1:dim(EMD_fcs[[2]])[1])] # add an row index count
  EMD_quantiles_list[[i]] <- EMD_fcs[[1]]
  EMD_details_list[[i]] <- EMD_fcs[[2]]
  
  if (i == 1){
    ## initialise a record of the (issue_time, target_time) pairs that exist for ALL zones
    VAR_times_mutual <- VAR_fcs[[2]][,.(issue_time, target_time, Horizon)]
    MC_times_mutual <- MC_fcs[[2]][,.(issue_time, target_time, Horizon)]
    EMD_times_mutual <- EMD_fcs[[2]][,.(issue_time, target_time, Horizon)]
  }else{
    varzone_rows <- VAR_fcs[[2]][,.(issue_time, target_time, Horizon)]
    VAR_times_mutual <- merge(VAR_times_mutual, varzone_rows) # keep the rows from mutual_rows, that also appear in zone_rows
    mczone_rows <- MC_fcs[[2]][,.(issue_time, target_time, Horizon)]
    MC_times_mutual <- merge(MC_times_mutual, mczone_rows)
    emdzone_rows <- VAR_fcs[[2]][,.(issue_time, target_time, Horizon)]
    EMD_times_mutual <- merge(EMD_times_mutual, emdzone_rows)
  }
}
rm(VAR_fcs)
rm(MC_fcs)
rm(EMD_fcs)


# now get the (issue_time, target_time) pairs we have forecasts for all models at, for all zones
mutual_times <- merge(merge(persistence_times_mutual, VAR_times_mutual), merge(MC_times_mutual, EMD_times_mutual))


model_errors <- CJ(Horizon=horizons, model=models)
## we are ready to bootstrap.
for (h in horizons){
  print (h)
  horizon_times <- mutual_times[Horizon==h, target_time]
  p_h_errs <- matrix(nrow=Nboots, ncol=3)
  var_h_errs <- matrix(nrow=Nboots, ncol=3)
  mc_h_errs <- matrix(nrow=Nboots, ncol=3)
  emd_h_errs <- matrix(nrow=Nboots, ncol=3)
  
  for (nb in c(1:Nboots)){
    sample_times <- sample(horizon_times, length(horizon_times), replace=TRUE)
    
    p_errs <- matrix(nrow=length(zones), ncol=3) # 1st col=MAE, 2nd col=RMSE, 3rd col=pinball.
    var_errs <- matrix(nrow=length(zones), ncol=3)
    mc_errs <- matrix(nrow=length(zones), ncol=3)
    emd_errs <- matrix(nrow=length(zones), ncol=3)
    
    ## now calculate errors for each zone and save in relevant matrices
    for (z in zones){
      
      ## first, persistence
      persistence_h_details <- persistence_details_list[[z]][Horizon==h]
      persistence_short_indices <- match(sample_times, persistence_h_details[, target_time])
      persistence_indices <- persistence_h_details[persistence_short_indices, id]
      p_errors <- persistence_quantiles_list[[z]][persistence_indices, q50] - persistence_details_list[[z]][persistence_indices, ActualPower]
      p_errs[z, 1] <- mean(abs(p_errors))
      p_errs[z, 2] <- sqrt(mean(p_errors^2))
      p_pinball <- pinball(persistence_quantiles_list[[z]][persistence_indices], realisations = persistence_details_list[[z]][persistence_indices, ActualPower], plot.it=FALSE)
      p_errs[z, 3] <- mean(p_pinball$Loss)
      
      
      ## now, VAR
      var_h_details <- VAR_details_list[[z]][Horizon==h]
      var_short_indices <- match(sample_times, var_h_details[, target_time])
      var_indices <- var_h_details[var_short_indices, id]
      var_errors <- VAR_quantiles_list[[z]][var_indices, q50] - VAR_details_list[[z]][var_indices, ActualPower]
      var_errs[z, 1] <- mean(abs(var_errors), na.rm=TRUE)
      var_errs[z, 2] <- sqrt(mean(var_errors^2, na.rm=TRUE))
      var_pinball <- pinball(VAR_quantiles_list[[z]][var_indices], realisations = VAR_details_list[[z]][var_indices, ActualPower], plot.it=FALSE)
      var_errs[z, 3] <- mean(var_pinball$Loss)
      
      ## MC
      mc_h_details <- MC_details_list[[z]][Horizon==h]
      mc_short_indices <- match(sample_times, mc_h_details[, target_time])
      mc_indices <- mc_h_details[mc_short_indices, id]
      mc_errors <- MC_quantiles_list[[z]][mc_indices, q50] - MC_details_list[[z]][mc_indices, ActualPower]
      mc_errs[z, 1] <- mean(abs(mc_errors), na.rm=TRUE)
      mc_errs[z, 2] <- sqrt(mean(mc_errors^2, na.rm=TRUE))
      mc_pinball <- pinball(MC_quantiles_list[[z]][mc_indices], realisations = MC_details_list[[z]][mc_indices, ActualPower], plot.it=FALSE)
      mc_errs[z, 3] <- mean(mc_pinball$Loss)
      
      ## and EMD.
      emd_h_details <- EMD_details_list[[z]][Horizon==h]
      emd_short_indices <- match(sample_times, emd_h_details[, target_time])
      emd_indices <- emd_h_details[emd_short_indices, id]
      emd_errors <- EMD_quantiles_list[[z]][emd_indices, q50] - EMD_details_list[[z]][emd_indices, ActualPower]
      emd_errs[z, 1] <- mean(abs(emd_errors), na.rm=TRUE)
      emd_errs[z, 2] <- sqrt(mean(emd_errors^2, na.rm=TRUE))
      emd_pinball <- pinball(EMD_quantiles_list[[z]][emd_indices], realisations = EMD_details_list[[z]][emd_indices, ActualPower], plot.it=FALSE)
      emd_errs[z, 3] <- mean(emd_pinball$Loss)
    }
    
    p_h_errs[nb,] <- colMeans(p_errs)
    var_h_errs[nb,] <- colMeans(var_errs)
    mc_h_errs[nb,] <- colMeans(mc_errs)
    emd_h_errs[nb,] <- colMeans(emd_errs)
    
  }
  
  cols <- c("MAE", "RMSE", "Pinball")
  model_errors[(Horizon==h) & (model=="Persistence"), (cols) := as.list(colMeans(p_h_errs))]
  model_errors[(Horizon==h) & (model=="VAR"), (cols) := as.list(colMeans(var_h_errs))]
  model_errors[(Horizon==h) & (model=="MC"), (cols) := as.list(colMeans(mc_h_errs))]
  model_errors[(Horizon==h) & (model=="EMD"), (cols) := as.list(colMeans(emd_h_errs))]
  
  
}

 
# fwrite(model_errors, file=paste0("./all_model_errors.csv"))
# model_errors <- fread(file=paste0("./all_model_errors.csv"))

