
#This script extracts growth rates out of FATES output
from_new_data <- T

#load packages
library(tidyverse)
library(ncdf4)
library(lubridate)

#load supporting functions
source('utils/supporting_funcs_esm_tools.R')
source('utils/figure_formatting_esm_tools.R')
source('generalFunctions.R')

#set path to the VDM output files
data_path <- '~/cloud/gdrive/FATES/FATES_data'

RA_values <- c(0.009, 0.079, 0.1, 0.149, 0.219, 0.289, 0.359, 0.43)

cases <- c()
for(c in RA_values){
  cases <- append(cases,paste0('bci-1pft-fates-main-from-inventory.Ccd4bc7aba-Fa0568083.2021-09-29_RAis',c,'AllVars'))
}

case_aliases <- c()
for (c in RA_values){
  case_aliases <- append(case_aliases,paste0('RA=',c))
}

ModelDataComparisonVars <- c('PFTnindivs','DDBH_CANOPY_SCPF','DDBH_UNDERSTORY_SCPF')

#put all model data into a data frame
if(from_new_data == T){
  model_data <- tibble()
  j <- 0
  for(c in cases){
    j <- j + 1
    print(paste("processing case",case_aliases[j]))
    files <- getFilePaths(data_path,c)
    for(f in files){
      
      ba <- getMultiDimSiteVar(f,varName = 'BA_SCPF')
      growth_canopy <- getMultiDimSiteVar(f,varName = 'DDBH_CANOPY_SCPF')
      growth_und <- getMultiDimSiteVar(f,varName = 'DDBH_UNDERSTORY_SCPF')
      size_class <- getMultiDimSiteVar(f,varName = 'fates_levscls')
      n_per_size_class <- getMultiDimSiteVar(f,varName = 'NPLANT_SCPF')
      date.x <- getDateFromNC(f)
      mort_und <- getMultiDimSiteVar(f,varName = 'MORTALITY_UNDERSTORY_SCPF')
      mort_can <- getMultiDimSiteVar(f,varName = 'MORTALITY_CANOPY_SCPF')
      agb <- getMultiDimSiteVar(f,varName = 'AGB_SCPF')
      
      temp <- tibble(RA = RA_values[j], 
                     date = date.x,
                     sc = size_class, 
                     N = n_per_size_class,
                     gC = growth_canopy, 
                     gU = growth_und, 
                     ba = ba,
                     mort_und = mort_und,
                     mort_can = mort_can,
                     agb = agb)
      
      model_data <- rbind(model_data,temp)
    }
  }
  write_csv(model_data,"data/model_predictions_growth_varying_RA.csv")
}

if(from_new_data == F){
  model_data <- read_csv("data/model_predictions_growth_varying_RA.csv")
}


model_data2 <- model_data %>%
  
  mutate(growth_per_ind_und = gU / N) %>%
  mutate(growth_per_ind_can = gC / N) %>%
  group_by(sc,RA) %>%
  #mutate_at(.vars = c('DDBH_CANOPY_SCPF','DDBH_UNDERSTORY_SCPF'),.funs = function(x){x / 1e4}) %>% #convert from cm ha-1 yr-1 to cm m-2 yr-1
  # mutate(growth_per_ind_und = DDBH_UNDERSTORY_SCPF / PFTnindivs,
  #        growth_per_ind_can = DDBH_CANOPY_SCPF / PFTnindivs) %>%
  #select(-date) %>%
  #group_by(case) %>%
  summarise(gU = mean(growth_per_ind_und),
            gC = mean(growth_per_ind_can)) %>%
  ungroup()

RAvsGrowth <- model_data2 %>%
  ggplot(aes(RA,gC)) +
  geom_line() +
  facet_wrap(~sc,scales = "free") +
  adams_theme

makePNG(fig = RAvsGrowth,path_to_output.x = "figures/",file_name = "RAvsGrowthCanopy.png",width = 14,height = 14)
