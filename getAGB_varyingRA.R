#load packages
library(tidyverse)
library(ncdf4)
library(lubridate)

#load supporting functions
source('utils/supporting_funcs_esm_tools.R')
source('utils/figure_formatting_esm_tools.R')
source('generalFunctions.R')


from_new_data <- T




###################################
########model predictions##########
###################################

#set path to the VDM output files
data_path <- '~/cloud/gdrive/FATES/FATES_data'

RA_values <- c(0.009, 0.079, 0.149, 0.219, 0.289, 0.359, 0.43)

cases <- c()
for(c in RA_values){
  cases <- append(cases,paste0('bci-1pft-fates-main-from-inventory.Ccd4bc7aba-Fa0568083.2021-09-29_RAis',c,'HgtMinIs2.55'))
}

case_aliases <- c()
for (c in RA_values){
  case_aliases <- append(case_aliases,paste0('RA=',c))
}


ModelDataComparisonVars <- c('AGB','PFTbiomass')

#put all model data into a data frame
if(from_new_data == T){
  model_data <- tibble()
  j <- 0
  for(c in cases){
    j <- j + 1
    files <- getFilePaths(data_path,c)
    print(paste("processing case",case_aliases[j]))
    c_data <- standardize_units(getManySiteVarsOverTime(files,ModelDataComparisonVars),ModelDataComparisonVars,varUnitsMapFates)
    c_data <- c_data %>% add_column(case = case_aliases[j])
    model_data <- rbind(model_data,c_data)
  }
  write_csv(model_data,"data/AGB_predictions_varying_RA.csv")
}

AGB_accumulation <- model_data %>%
  ggplot(aes(date,AGB)) +
  geom_line() +
  facet_wrap(~case,scales = "free") +
  adams_theme

makePNG(fig = RAvsGrowth,path_to_output.x = "figures/",file_name = "RAvsGrowthCanopy.png",width = 14,height = 14)
