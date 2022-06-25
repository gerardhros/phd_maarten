###linear regression (predict Pox:Ptot as a function of Alox+Feox and PSI)

#load libraries
library(caret)
library(data.table)

#set wd
setwd("C:/Users/Maarten van Doorn/SPRINGG/Gerard Ros - NMI-PROJ/JustPmaps/02 models/pox-ptot ratio")

#constants
ptop2o5 <- 2.29
Pmolarmass <- 30.973762

#regression function, returning lm.model, prediction plot and data
regression_poxptot <- function(){
  #read NMI dataset
  dt <- fread("data.csv")
  
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
  
  #create regression model
  model.lm <- lm(pox_ptot ~ log10(`alox_feox (mmol/kg)`) + psi, data = dt) #pox_ptot
  
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
  
  #bind to list and return
  result <- list(model = model.lm,
                 modelsummary = summary(model.lm),
                 plot = plot1, 
                 data = dt)
  
  #return
  return(result)
}

#run regression and save as .RDS
result <- regression_poxptot()
saveRDS(result, "regression results.RDS")
