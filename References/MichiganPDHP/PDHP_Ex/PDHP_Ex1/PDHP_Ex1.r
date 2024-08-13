#---------------------------------------------------
#      Application: Meyer, Viscusi and Durbin (1995, AER)
#---------------------------------------------------
#-----------------------------------------------------------------------------
# Startup - clear memory, load packages, set working directory, and import the dataset
# Clear memory
rm(list=ls())
library(wooldridge)
library(here)
library(fixest)

#-----------------------------------------------------------------------------
# import the data
injury <- wooldridge::injury
# Get unit's id (repeated cross section)
injury$id <- 1:dim(injury)[1]

#-----------------------------------------------------------------------------
# create two state subsets
# Kentucky subset
injury_ky <- subset(injury, injury$ky==1)
# Michigan subset
injury_mi <- subset(injury, injury$mi==1)
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#     A. Kentucky Sample
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Unconditional DiD analysis
#-----------------------------------------------------------------------------
# Duration of benefits
twfe_ky_dur <- fixest::feols(durat ~ highearn + afchnge + I(highearn * afchnge), 
                      data = injury_ky,
                      cluster = ~id)

summary(twfe_ky_dur)

# Log of duration of benefits
twfe_ky_ldur <- fixest::feols(ldurat ~ highearn + afchnge + I(highearn * afchnge), 
                             data = injury_ky,
                             cluster = ~id)


summary(twfe_ky_ldur)
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#     B. Michigan Sample
#-----------------------------------------------------------------------------
 #Your turnn!