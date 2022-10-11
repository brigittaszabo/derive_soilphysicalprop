#########################################################
# Script for soil properties calculation using the format of SWAT soil table
# Author/developer: Brigitta Szabó, János Mészáros
# Last modification: 18/09/2022
#########################################################

# # installing mandatory libraries for euptfv2
# install.packages("devtools")
# install.packages("glue")
# install.packages("Rdpack")
# devtools::install_github("tkdweber/euptf2")


# install and import other necessary packages
packages = c("openxlsx", "stringr")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

library(euptf2)

options(warn = -1)

# defining working directory automatically to the path where R file placed
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# reading soil data
load("soil_prop_in.rdata") # adapt
summary(soil_prop)

# basic structure of the previous file to see variables, data structure etc.
str(soil_prop)

# creating a second soil data frame for later modifications and keep the original one intact
soil_prop_extended <- soil_prop
summary(soil_prop_extended)

# # replace 0 with NA if 0 means no data
# soil_prop_extended[soil_prop_extended == 0] <- NA
# summary(soil_prop_extended)


# compute bulk density with PTF of Alexander et al. (1980)
# Alexander E B. 1980. Bulk densities of California soils in relation to other soil properties. Soil Sci Soc Am J. 44: 689–692.
CBNx = names(soil_prop_extended)[(startsWith(names(soil_prop_extended), "SOL_CBN"))]
cbnn = length(CBNx)

for (j in 1:cbnn) {
  for (i in 1:nrow(soil_prop_extended)){
    soil_prop_extended[i, (paste0("BD", j))] <- 1.72-0.294*((soil_prop_extended[i,(paste0("SOL_CBN", j))])^0.5)
    }
}

summary(soil_prop_extended)

BDx = names(soil_prop_extended)[(startsWith(names(soil_prop_extended), "BD"))]
bdn = length(BDx)

for (j in 1:bdn) {
  for (i in 1:nrow(soil_prop)){
    if (soil_prop[i,(paste0("SOL_CBN", j))] > 1){
      soil_prop_extended[i, (paste0("SOL_BD", j))] <- soil_prop_extended[i,(paste0("BD", j))] + 0.009 * soil_prop_extended[i,(paste0("CLAY", j))]
    } else {
      soil_prop_extended[i,(paste0("SOL_BD", j))] <- soil_prop_extended[i,(paste0("BD", j))] + 0.005 * soil_prop_extended[i,(paste0("CLAY", j))] + 0.001 * soil_prop_extended[i,(paste0("SILT", j))]
    }
  }
}


summary(soil_prop_extended)

# computed 0 values mean NA -> rewrite those to NA
soil_prop_extended[soil_prop_extended == 0] <- NA
summary(soil_prop_extended)

# new DF validation 
soil_prop_extended


# compute soil hydraulic properties
# predict van Genuchten parameters
# use PTF07 of Szabó et al. (2021) if SOL_Z, CLAY, SILT, SAND, BD and SOL_CBN1 are available, else use the PTF recommended in Table 11 of https://doi.org/10.5194/gmd-14-151-2021 

SOL_Zx = names(soil_prop_extended)[(startsWith(names(soil_prop_extended), "SOL_Z"))]
SOL_Zx = SOL_Zx[SOL_Zx != "SOL_ZMX"]
zn = length(SOL_Zx)
soil_prop_extended[soil_prop_extended == 0] <- NA
soil_prop_extended$DEPTH_M1 <- soil_prop_extended$SOL_Z1/2/10

soil_prop2 <- cbind(soil_prop_extended, (((soil_prop_extended[, SOL_Zx[2:zn]])-(soil_prop_extended[,SOL_Zx[1:zn-1]]))/2 + soil_prop_extended[,SOL_Zx[1:zn-1]])/10)
names(soil_prop2)[(ncol(soil_prop2)-zn+2):ncol(soil_prop2)] <- paste("DEPTH_M", 2:zn, sep="")
soil_prop2$rownum <- 1:nrow(soil_prop2)

for (i in 1:zn){
  input <- soil_prop2[, c("rownum", paste0("DEPTH_M",i), paste0("SOL_BD",i), paste0("SOL_CBN",i), paste0("CLAY",i), paste0("SILT",i),paste0("SAND",i))]
  
  # input <- Filter(function(x)!all(is.na(x)), input)
  names(input)[1:7] <- c("rownum","DEPTH_M","BD", "OC", "USCLAY", "USSILT", "USSAND")
  
  pred_VG1 <- euptfFun(ptf = "PTF07", predictor = input, target = "VG", query = "predictions")
  names(pred_VG1)[2:6] <- c("THS","THR", "ALP", "N", "M")
  pred_VG <- merge(pred_VG1[, c(1:6,8,10:14)], input[,c(1,2)], by="rownum", all.y=TRUE)
  
  FC <- pred_VG$THR+(pred_VG$THS-pred_VG$THR)*((1+(((pred_VG$N-1)/pred_VG$N)^(1-2*pred_VG$N)))^((1-pred_VG$N)/pred_VG$N))
  nam <- paste("FC", i, sep = "")
  assign(nam, FC, envir=.GlobalEnv)
  
  WP <- pred_VG$THR+((pred_VG$THS-pred_VG$THR)/((1+pred_VG$THR*(15000^pred_VG$N))^(1-(1/pred_VG$N))))
  nam <- paste("WP", i, sep = "")
  assign(nam, WP, envir=.GlobalEnv)
  
  SOL_AWC <- FC-WP
  nam <- paste("SOL_AWC", i, sep = "")
  assign(nam, SOL_AWC, envir=.GlobalEnv)
  
  SOL_K <- (4.65 * (10^4) * pred_VG$THS * (pred_VG$ALP^2))*0.41666001
  nam <- paste("SOL_K", i, sep = "")
  assign(nam, SOL_K, envir=.GlobalEnv)
  
  # compute albedo
  # method of Gascoin et al. (2009) from Table 6 of Abbaspour, K.C., AshrafVaghefi, S., Yang, H. & Srinivasan, R. 2019. Global soil, landuse, evapotranspiration, historical and future weather databases for SWAT Applications. Scientific Data, 6:263.
  SOL_ALB <- 0.15+0.31*exp(-12.7*FC)
  nam <- paste("SOL_ALB", i, sep = "")
  assign(nam, SOL_ALB, envir=.GlobalEnv)
  
  # compute USLE erodibility K factor
  # method published in Sharpley and Williams (1990) based on Table 5 of Abbaspour, K.C., AshrafVaghefi, S., Yang, H. & Srinivasan, R. 2019. Global soil, landuse, evapotranspiration, historical and future weather databases for SWAT Applications. Scientific Data, 6:263.
  ES <- 0.2+0.3*exp(-0.0256*input$USSAND*(1-(input$USSILT/100)))
  ECT <- (input$USSILT/(input$USCLAY+input$USSILT))^0.3
  EOC <- 1-(0.25*input$OC/(input$OC+exp(3.72-2.95*input$OC)))
  EHS <- 1-(0.7*(1-input$USSAND/100)/((1-input$USSAND/100)+exp(-5.51+22.9*(1-input$USSAND/100))))
  USLE_K <- ES*ECT*EOC*EHS
  nam <- paste("USLE_K", i, sep = "")
  assign(nam, USLE_K, envir=.GlobalEnv)
}


num <- str_count(ls(), "SOL_AWC")
sum(num)
awcn <- sum(num)-1

pred_SOL_par <- data.frame(mget(c((paste0("SOL_AWC", 1:awcn)), (paste0("SOL_K", 1:awcn)), (paste0("SOL_ALB", 1:awcn)), (paste0("USLE_K", 1:awcn)))))

sol_bd_col <- names(soil_prop_extended)[(startsWith(names(soil_prop_extended), "SOL_BD"))]
SOL_BD <- soil_prop_extended[, c(sol_bd_col)]

patterns <- c("SOL_AWC","SOL_K", "SOL_ALB", "USLE_K", "BD", "DEPTH")
soil_prop_final1 <- soil_prop_extended[, -grep(paste(patterns, collapse="|"), colnames(soil_prop_extended))]

soil_prop_final <- cbind(soil_prop_final1, SOL_BD, pred_SOL_par)
soil_prop_final
write.xlsx(soil_prop_final, file = "user_soil_table.xlsx")


pred_FC_WP <- data.frame(mget(c((paste0("FC", 1:awcn)), (paste0("WP", 1:awcn)))))
pred_FC_WP
write.xlsx(pred_FC_WP, file = "pred_FC_WP.xlsx")
