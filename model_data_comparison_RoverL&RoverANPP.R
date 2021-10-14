
#This script takes data from multiple model simulations and compares it with observations
from_new_data <- F

#load packages
library(tidyverse)
library(ncdf4)
library(lubridate)

#load supporting functions
source('utils/supporting_funcs_esm_tools.R')
source('utils/figure_formatting_esm_tools.R')
source('generalFunctions.R')

#load obs data
obs_data_RoL <- read_csv("data/RoverL_BCI_obs.csv") %>%
  mutate(var_new = case_when(
    var == "RoL" ~ "R/L",
    var == "RoANPP" ~ "R/ANPP",
    TRUE ~ var
  )) %>%
  mutate(var = var_new)


#placeholder data until Rachel gets this data
obs_data_RoANPP <- read_csv("data/RoverL_BCI_obs.csv") %>% 
  mutate_at(.vars = "value",.funs = function(x){x*0.5}) %>%
  mutate(var = 'RoANPP')
obs_data_recruitment <-  read_csv("data/Recruitment_BCI_obs.csv") %>%
  filter(simYr %in% c(6,7)) %>%
  group_by(case) %>% summarise(Mvalue = mean(value), sd = sd(value)) %>%
  rename(value = Mvalue)

#set path to the VDM output files
data_path <- '~/cloud/gdrive/FATES/FATES_data'

RA_values <- c(0.009, 0.079, 0.149, 0.219, 0.289, 0.359, 0.43)

cases <- c()
for(c in RA_values){
  cases <- append(cases,paste0('bci-1pft-fates-main-from-inventory.Ccd4bc7aba-Fa0568083.2021-09-29_RAis',c))
}

case_aliases <- c()
for (c in RA_values){
  case_aliases <- append(case_aliases,paste0('RA=',c))
}


ModelDataComparisonVars <- c('SEEDS_IN_LOCAL_ELEM','LEAF_LITTER_IN','NPP','RECRUITMENT','NPP_FNRT_SCPF', 'NPP_BGSW_SCPF','NPP_BGDW_SCPF')

#put all model data into a data frame
if(from_new_data == T){
  model_data <- tibble()
  j <- 0
  for(c in cases){
    j <- j + 1
    files <- getFilePaths(data_path,c)
    print(paste("processing case",case_aliases[j]))
    c_data <- standardize_units(getManySiteVarsOverTime(files,ModelDataComparisonVars),ModelDataComparisonVars,varUnitsMapFates)
    c_data <- c_data %>% mutate(ANPP = NPP - (NPP_FNRT_SCPF + NPP_BGSW_SCPF + NPP_BGDW_SCPF)) %>% add_column(case = case_aliases[j])
    model_data <- rbind(model_data,c_data)
  }
  write_csv(model_data,"data/model_predictions_varying_RA")
}

model_data <- read_csv("data/model_predictions_varying_RA")

#calculate annual means of R/L and R/ANPP by case
model_data_for_fig <- model_data %>% filter(date > as.Date("1907-01-01")) %>%
  mutate(RoL = SEEDS_IN_LOCAL_ELEM / LEAF_LITTER_IN,
         RoANPP = SEEDS_IN_LOCAL_ELEM / ANPP) %>%
  mutate_at(.vars = "RECRUITMENT",.funs = function(x){x * days_per_yr * m2_per_ha}) %>% #converting recruitment to ind. per ha per year
  select(simYr,date,RoL,RoANPP,RECRUITMENT,case) %>%
  toLongForm(vars = c('RoL','RoANPP',"RECRUITMENT")) %>%
  group_by(case,simYr,var) %>%
  summarise(value = mean(value)) %>%
  ungroup()
  


#create R/L figure comparison
litter_flux_fig <- model_data_for_fig %>%
  mutate(RA = as.numeric(substr(case,4,8))) %>%
  mutate(var_new = case_when(
    var == "RoL" ~ "R/L",
    var == "RoANPP" ~ "R/ANPP",
    TRUE ~ var
  )) %>%
  mutate(var = var_new) %>%
  filter(var %in% c('R/L','R/ANPP')) %>% 
  ggplot(aes(case,value, color = RA)) +
  scale_color_continuous() +
  facet_wrap(~var) +
  #scale_color_manual(values = c("blue","yellow")) +
  geom_boxplot() +
  geom_boxplot(data = obs_data_RoANPP %>%
                 mutate(var_new = case_when(
                   var == "RoL" ~ "R/L",
                   var == "RoANPP" ~ "R/ANPP",
                   TRUE ~ var
                 )) %>% 
                 mutate(var = var_new) %>% rbind(obs_data_RoL), mapping = aes(case,value), color ="black") +
  ylab("") +
  theme_minimal() +
  adams_theme +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())



makePNG(litter_flux_fig,"figures/","litter_flux_fig")







#create figure for recruitment comparison at BCI
recruitment_comp_df <- model_data_for_fig %>%
  filter(var == 'RECRUITMENT') %>%
  group_by(case) %>%
  summarise(Mvalue = mean(value), 
            sd = sd(value,na.rm = T)) %>%
  mutate(RA = as.numeric(substr(case,4,8))) %>%
  rename(value = Mvalue)



recruitment_comp_fig <- recruitment_comp_df %>%
  ggplot(aes(case,value,color = RA)) +
  geom_point(size = 7) +
  geom_errorbar(aes(ymin = value-sd,ymax = value+sd, width = 0)) +
  geom_point(data = obs_data_recruitment, mapping = aes(case,value), color ="black", shape = 2, size = 7) +
  geom_errorbar(data = obs_data_recruitment, mapping = aes(ymin = value-sd,ymax = value+sd, width = 0), color = "black") +
  scale_color_continuous() +
  rec.y.axis +
  scale_y_log10(breaks = c(100,1000,10000,32000),
                labels = c(100,1000,10000,32000),
                limits = c(50,32000)) +
  theme_minimal() +
  adams_theme +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())

makePNG(recruitment_comp_fig,"figures/","model-data-comp-recruitment",height = 4, width = 5)



