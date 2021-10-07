
#This script takes data from multiple model simulations and compares it with observations

#load packages
library(tidyverse)
library(ncdf4)
library(lubridate)

#load supporting functions
source('utils/supporting_funcs_esm_tools.R')
source('utils/figure_formatting_esm_tools.R')
source('generalFunctions.R')

#load obs data
obs_data_RoL <- read_csv("data/RoverL_BCI_obs.csv")
#placeholder data until Rachel gets this data
obs_data_RoANPP <- read_csv("data/RoverL_BCI_obs.csv") %>% 
  mutate_at(.vars = "value",.funs = function(x){x*0.5}) %>%
  mutate(var = 'RoANPP')
obs_data_recruitment <-  read_csv("data/Recruitment_BCI_obs.csv") 

#set path to the VDM output files
data_path <- '~/cloud/gdrive/FATES/FATES_data'

cases <- c('bci-1pft-fates-main-from-inventory.Ccd4bc7aba-Fa0568083.2021-09-29_RAis0.009','bci-1pft-fates-main-from-inventory.Ccd4bc7aba-Fa0568083.2021-09-29_RAis0.079')
case_aliases <- c('RA=0.009','RA=0.079')

ModelDataComparisonVars <- c('SEEDS_IN_LOCAL_ELEM','LEAF_LITTER_IN','NPP','RECRUITMENT','NPP_FNRT_SCPF', 'NPP_BGSW_SCPF','NPP_BGDW_SCPF')

#put all model data into a data frame
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
#could add recruitment easily into this figure
model_data_for_fig %>%
  mutate(RaVal = as.numeric(substr(case,4,8))) %>%
  filter(var %in% c('RoL','RoANPP')) %>% 
  ggplot(aes(case,value, color = RaVal)) +
  facet_wrap(~var,scales = "free") +
  scale_color_continuous() +
  #scale_color_manual(values = c("blue","yellow")) +
  geom_boxplot() +
  geom_boxplot(data = obs_data_RoANPP %>% 
                 rbind(obs_data_RoL), mapping = aes(case,value), color ="black") +
  ylab("") +
  theme_minimal() +
  adams_theme +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())

#create figure for recruitment comparison at BCI
rbind(model_data_for_fig) %>%
  mutate(RaVal = as.numeric(substr(case,4,8))) %>%
  filter(var == 'RECRUITMENT') %>% 
  ggplot(aes(case,value,color = RaVal)) +
  facet_wrap(~var,scales = "free") +
  geom_boxplot() +
  geom_boxplot(data = obs_data_recruitment, 
               mapping = aes(case,value), color ="black") +
  scale_color_continuous() +
  scale_y_log10() +
  theme_minimal() +
  adams_theme +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())





