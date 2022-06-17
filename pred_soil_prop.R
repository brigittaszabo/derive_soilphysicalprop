#########################################################
# Script for soil properties calculation using the format of SWAT soil table
# Author/developer: Brigitta Szabó, János Mészáros
# Last modification: 17/06/2022
#########################################################

# installing mandatory libraries for euptfv2
install.packages("devtools")
install.packages("glue")
install.packages("Rdpack")
devtools::install_github("tkdweber/euptf2")


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

# defining working directory automatically to the path where R file placed
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# reading soil data from CSV
soil_prop <- read.csv2("soil_prop_in.csv", dec = ".", sep = ";") # adapt
summary(soil_prop)

# basic structure of the previous file to see variables, data structure etc.
str(soil_prop)

# creating a second soil data frame for later modifications and keep the original one intact
soil_prop_extended <- soil_prop
summary(soil_prop_extended)

# numeric variables to numeric
soil_prop_extended[9:ncol(soil_prop_extended)] = lapply(soil_prop_extended[9:ncol(soil_prop_extended)], FUN = function(y){as.numeric(y)})
summary(soil_prop_extended)

# replace 0 with NA if 0 means no data
soil_prop_extended[soil_prop_extended == 0] <- NA
summary(soil_prop_extended)

# compute effective bulk density
# reference: Wessolek et al. (2009)
for (j in 1:10) {
  for (i in 1:nrow(soil_prop)){
    if (soil_prop[i,((j*12)+5)] > 1){
      soil_prop_extended[i,(ncol(soil_prop)+j)] <- soil_prop[i,((j*12)+2)] + 0.009 * soil_prop[i,((j*12)+6)]
    } else {
      soil_prop_extended[i,(ncol(soil_prop)+j)] <- soil_prop[i,((j*12)+2)] + 0.005 * soil_prop[i,((j*12)+6)] + 0.001 * soil_prop[i,((j*12)+7)]
    }
  }
}

summary(soil_prop_extended)

# computed 0 values mean NA -> rewrite those to NA
soil_prop_extended[soil_prop_extended == 0] <- NA
summary(soil_prop_extended)

# renaming effective bulk density columns
dfnames <- names(soil_prop)
for (i in 1:10) {
  dfnames[length(dfnames)+1] <- paste0("SOL_BD", i, "_eff")
}

# rewriting extended DF names
names(soil_prop_extended) <- dfnames


# new DF validation 
soil_prop_extended


# compute soil hydraulic properties
# predict van Genuchten parameters
# use PTF07 of Szabó et al. (2021) if SOL_Z, CLAY, SILT, SAND, BD and SOL_CBN1 are available, else use the PTF recommended in Table 11 of https://doi.org/10.5194/gmd-14-151-2021 

soil_prop$DEPTH_M1 <- soil_prop$SOL_Z1/2
soil_prop2 <- cbind(soil_prop, ((soil_prop[,c(25,37,49,61,73,85,97,109,121)])-(soil_prop[,c(13,25,37,49,61,73,85,97,109)]))/2 + soil_prop[,c(13,25,37,49,61,73,85,97,109)])
names(soil_prop2)[154:length(soil_prop2)] <- paste("DEPTH_M", 2:10, sep="")


for (i in 1:10){
  input <- soil_prop2[, c(paste0("DEPTH_M",i), paste0("SOL_BD",i), paste0("SOL_CBN",i), paste0("CLAY",i), paste0("SILT",i),paste0("SAND",i))]
  
  # if 0 means NA, define it
  input[input == 0] <- NA
  # input <- Filter(function(x)!all(is.na(x)), input)
  names(input)[1:6] <- c("DEPTH_M","BD", "OC", "USCLAY", "USSILT", "USSAND")
  
  pred_VG <- euptfFun(ptf = "PTF07", predictor = input, target = "VG", query = "predictions")
  names(pred_VG)[2:6] <- c("THS","THR", "ALP", "N", "M")
  
  FC <- pred_VG$THR+(pred_VG$THS-pred_VG$THR)*((1+(((pred_VG$N-1)/pred_VG$N)^(1-2*pred_VG$N)))^((1-pred_VG$N)/pred_VG$N))
  nam <- paste("FC", i, sep = "")
  assign(nam, FC, envir=.GlobalEnv)
  
  WP <- pred_VG$THR+((pred_VG$THS-pred_VG$THR)/((1+pred_VG$THR*(15000^pred_VG$N))^(1-(1/pred_VG$N))))
  nam <- paste("WP", i, sep = "")
  assign(nam, WP, envir=.GlobalEnv)
  
  SOL_AWC <- FC-WP
  nam <- paste("SOL_AWC", i, sep = "")
  assign(nam, SOL_AWC, envir=.GlobalEnv)
  
  SOL_K <- 4.65 * (10^4) * pred_VG$THS * (pred_VG$ALP^2)
  nam <- paste("SOL_K", i, sep = "")
  assign(nam, SOL_K, envir=.GlobalEnv)
  
  # compute albedo
  # Table 6 of Abbaspour, K.C., AshrafVaghefi, S., Yang, H. & Srinivasan, R. 2019. Global soil, landuse, evapotranspiration, historical and future weather databases for SWAT Applications. Scientific Data, 6:263.
  SOL_ALB <- 0.26+0.1068*exp(-4.9*FC)
  nam <- paste("SOL_ALB", i, sep = "")
  assign(nam, SOL_ALB, envir=.GlobalEnv)
  
  # compute USLE erodibility K factor
  # Table 5 of Abbaspour, K.C., AshrafVaghefi, S., Yang, H. & Srinivasan, R. 2019. Global soil, landuse, evapotranspiration, historical and future weather databases for SWAT Applications. Scientific Data, 6:263.
  ES <- 0.2+0.3*exp(-0.256*input$USSAND*(1-(input$USSILT/100)))
  ECT <- (input$USSILT/(input$USCLAY+input$USSILT))^0.3
  EOC <- 1-(0.25*input$OC/(input$OC+exp(0.72-2.95*input$OC)))
  EHS <- 1-(0.7*(1-input$USSAND/100)/((1-input$USSAND/100)+exp(-5.51+22.9*(1-input$USSAND/100))))
  USLE_K <- ES*ECT*EOC*EHS
  nam <- paste("USLE_K", i, sep = "")
  assign(nam, USLE_K, envir=.GlobalEnv)
}


num <- str_count(ls(), "SOL_AWC")
sum(num)

pred_SOL_par <- data.frame(mget(c((paste0("SOL_AWC", 1:(sum(num)-1))), (paste0("SOL_K", 1:(sum(num)-1))), (paste0("SOL_ALB", 1:(sum(num)-1))), (paste0("USLE_K", 1:(sum(num)-1))))))

SOL_BD <- soil_prop_extended[, c(153:162)]
names(SOL_BD) <- gsub("_eff", "", names(SOL_BD))

patterns <- c("SOL_AWC","SOL_K", "SOL_ALB", "USLE_K")
soil_prop_final1 <- soil_prop_extended[, -grep(paste(patterns, collapse="|"), colnames(soil_prop_extended))]

soil_prop_final <- cbind(soil_prop_final1[, 1:112], SOL_BD, pred_SOL_par)
soil_prop_final

write.xlsx(soil_prop_final, file = "soil_prop_out.xlsx")


pred_FC_WP <- data.frame(mget(c((paste0("FC", 1:(sum(num)-1))), (paste0("WP", 1:(sum(num)-1))))))
write.xlsx(pred_FC_WP, file = "pred_FC_WP.xlsx")
