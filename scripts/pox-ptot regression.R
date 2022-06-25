###linear regression (predict Pox:Ptot as a function of Alox+Feox and PSI)

#load libraries
require(caret); require(data.table)

# set LOCATION OF DATA FILES
loc <- paste0("D:/OneDrive - SPRINGG/NMI-PROJ/JustPmaps/02 models/pox-ptot ratio/")

#constants
ptop2o5 <- 2.29
Pmolarmass <- 30.973762

#read NMI dataset
dt <- fread(paste0(loc,"data.csv"))

#regression function, returning lm.model, prediction plot and data

  #calc PSI
  dt[, psi := `p_ox (mmol/kg)` / (`al_ox (mmol/kg)` + `fe_ox (mmol/kg)`)]
  
  #only select when psi != NA
  dt <- dt[!is.na(psi)]
  
  #convert prt to mmol/kg
  dt[, `prt (mg p2o5/kg)` := `p_rt (mg p2o5/100g)` * 10]
  dt[, `prt (mg p/kg)` := `prt (mg p2o5/kg)` / ptop2o5]
  dt[, `prt (mmol/kg)` := `prt (mg p/kg)` / Pmolarmass]
  
  #only keep pox values <= prt
  dt <- dt[`p_ox (mmol/kg)` <= `prt (mmol/kg)`]
  
  #calc pox_ptot ratio
  dt[, pox_ptot := `p_ox (mmol/kg)` / `prt (mmol/kg)`]
  
  #calc alox + feox
  dt[, `alox_feox (mmol/kg)` := `al_ox (mmol/kg)` + `fe_ox (mmol/kg)`]
  
  # change names by Gerard to simplify use in models
  dt[,alfeox := `alox_feox (mmol/kg)`]
  dt[,alfeox_log10 := log10(alfeox)]
  dt[,pox_ptot_exp := exp(pox_ptot)]
  
  # create regression model for pox-ptot ratio
  model.lm <- lm(pox_ptot ~ alfeox_log10 + psi, data = dt) 
  model.lm <- lm(pox_ptot ~ alfeox_log10, data = dt) 
  
  #predict on training set
  dt$pox_ptot_predicted <- predict(model.lm, dt)
  
  #extract metrics
  metrics <- data.table(R2 = R2(dt$pox_ptot, dt$pox_ptot_predicted),
                        RMSE = RMSE(dt$pox_ptot, dt$pox_ptot_predicted),
                        formula = paste0(round(model.lm$coefficients[2], 3), "log(feox_alox) + ", round(model.lm$coefficients[3], 3), "psi ", round(model.lm$coefficients[1], 3)))
  metrics[, label := paste0("R2 = ", round(R2, 2), "\nRMSE = ", round(RMSE, 2), "\n", formula)]
  
  #plot
  plot1 <- ggplot(dt, aes(x = pox_ptot, y = pox_ptot_predicted)) +
    theme_bw() + geom_point() +
    scale_x_continuous(limits = c(0, 1.2), breaks = seq(0, 1.25, 0.25)) +
    scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1.25, 0.25)) +
    theme(text = element_text(size = 15)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "blue", size = 1) +
    annotate(geom = "text", x = 0.25, y = 1.1, label = metrics$label) +
    labs(x = "Measured Pox:Ptotal (ratio)", y = "Predicted Pox:Ptotal (ratio)")
  plot1
  #bind to list and return
  result <- list(model = model.lm,
                 modelsummary = summary(model.lm),
                 plot = plot1, 
                 data = dt)
  

saveRDS(result, "data/poxptot_lm.RDS")
