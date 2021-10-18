
#This script takes data from multiple model simulations and compares it with observations
from_new_data <- F
ANPP_df <- read_csv('data/ANPP_ests_BCI_Luquillo.csv')

bci_ANPP <- ANPP_df %>% filter(Site == 'BCI') %>% pull(ANPP) %>% `/`(2) #dividing by 2 to convert biomass to carbon 

#load packages
library(tidyverse)
library(ncdf4)
library(lubridate)

#load supporting functions
source('utils/supporting_funcs_esm_tools.R')
source('utils/figure_formatting_esm_tools.R')
source('generalFunctions.R')

#load R/L obs
RoL_BCI_obs <- read_csv("data/RoverL_BCI_obs.csv") %>%
  mutate(var_new = case_when(
    var == "RoL" ~ "R/L",
    var == "RoANPP" ~ "R/ANPP",
    TRUE ~ var
  )) %>%
  mutate(var = var_new) %>%
  select(-var_new) %>%
  filter(var == "R/L")






RoANPP_BCI_obs <- read_csv("data/RoverL_BCI_obs.csv") %>%
  filter(var == 'Rgm2yr') %>%
  add_column(ANPP = bci_ANPP) %>%
  mutate(`R/ANPP` = value / ANPP) %>%
  mutate(var = "R/ANPP", value = `R/ANPP`) %>%
  select(case,simYr,var,value)

#load recruitment obs
recruitmentBCIobs <-  read_csv("data/Recruitment_BCI_obs.csv") %>%
  filter(simYr %in% c(6,7)) %>%
  group_by(case) %>% summarise(Mvalue = mean(value), sd = sd(value)) %>%
  rename(value = Mvalue) %>%
  mutate(var = "RECRUITMENT") %>%
  select(case,var,value,sd)

#join all obs.
#load R / ANPP obs
bci_all_obs <- read_csv('data/all_BCI_obs.csv')

BCI_obs <- bci_all_obs %>%
  select(-units) %>%
  drop_na(year) %>%
  spread(var,value) %>%
  mutate(`R/ANPP` = case_when(
    site == "BCI" ~ R / 1800,
    site == "Luquillo" ~ 51 / 1050
  )) %>%
  toLongForm(c('R','L','R/L','R/ANPP')) %>%
  select(-L) %>%
  rename(yr = year) %>%
  filter(var != "R") %>%
  mutate(site = "BCI obs.") %>% rename(case = site) %>%
  select(case,var,value) %>%
  group_by(case,var) %>% summarise(Mvalue = mean(value), sd = sd(value)) %>%
  rename(value = Mvalue) %>%
  ungroup() %>%
  rbind(recruitmentBCIobs)

#write_csv(BCI_obs,"data/bciObsClean.csv")

###################################
########model predictions##########
###################################

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
modelPredictionsBCI <- model_data %>% filter(date > as.Date("1907-01-01")) %>%
  mutate(RoL = SEEDS_IN_LOCAL_ELEM / LEAF_LITTER_IN,
         RoANPP = SEEDS_IN_LOCAL_ELEM / ANPP) %>%
  mutate_at(.vars = "RECRUITMENT",.funs = function(x){x * days_per_yr * m2_per_ha}) %>% #converting recruitment to ind. per ha per year
  select(simYr,date,RoL,RoANPP,RECRUITMENT,case) %>%
  toLongForm(vars = c('RoL','RoANPP',"RECRUITMENT")) %>%
  group_by(case,simYr,var) %>%
  summarise(value = mean(value)) %>%
  ungroup() %>%
  mutate(RA = as.numeric(substr(case,4,8))) %>%
  mutate(var_new = case_when(
    var == "RoL" ~ "R/L",
    var == "RoANPP" ~ "R/ANPP",
    TRUE ~ var
  )) %>%
  mutate(var = var_new) %>%
  filter(var %in% c('R/L','R/ANPP','RECRUITMENT')) %>%
  select(-var_new) %>%
  group_by(case,var) %>%
  summarise(Mvalue = mean(value), 
            sd = sd(value,na.rm = T)) %>%
  mutate(RA = as.numeric(substr(case,4,8))) %>%
  rename(value = Mvalue)



modelDataComp_RoL <- modelPredictionsBCI %>%
  filter(var == "R/L") %>%
  ggplot(aes(case,value,color = RA)) +
  geom_point(size = 7) +
  geom_errorbar(aes(ymin = value-sd,ymax = value+sd, width = 0)) +
  geom_point(data = BCI_obs %>% filter(var == "R/L"), 
             mapping = aes(case,value), color ="black", shape = 2, size = 5, stroke = 2) +
  geom_errorbar(data = BCI_obs %>% filter(var == "R/L"), 
                mapping = aes(ymin = value-sd,ymax = value+sd, width = 0), color = "black") +
  scale_color_continuous() +
  scale_y_continuous(limits = c(0,1.5), breaks = c(0,0.5,1.0,1.5)) +
  labs(title = "R/L") +
  ylab("") +
  adams_theme +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        legend.position = c(0.25,0.75)) +
  theme(legend.key.size = unit(0.8, 'cm'), #change legend key size
        legend.key.height = unit(0.6, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=16))



modelDataComp_RoANPP <- modelPredictionsBCI %>%
  filter(var == "R/ANPP") %>%
  ggplot(aes(case,value,color = RA)) +
  geom_point(size = 7) +
  geom_errorbar(aes(ymin = value-sd,ymax = value+sd, width = 0)) +
  geom_point(data = BCI_obs %>% filter(var == "R/ANPP"), 
             mapping = aes(case,value), color ="black", shape = 2, size = 5, stroke = 2) +
  geom_errorbar(data = BCI_obs %>% filter(var == "R/ANPP"), 
                mapping = aes(ymin = value-sd,ymax = value+sd, width = 0), color = "black") +
  scale_color_continuous() +
  scale_y_continuous(limits = c(0,0.4)) +
  ylab("") +
  labs(title = "R/ANPP") +
  adams_theme +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        legend.position = "none")


modelDataComp_Recruitment <- modelPredictionsBCI %>%
  filter(var == "RECRUITMENT") %>%
  ggplot(aes(case,value*0.19878,color = RA)) +
  geom_point(size = 7) +
  geom_errorbar(aes(ymin = value*0.19878-sd*0.19878,ymax = value*0.19878+sd*0.19878, width = 0)) +
  geom_point(data = BCI_obs %>% filter(var == "RECRUITMENT"), 
             mapping = aes(case,value), color ="black", shape = 2, size = 5, stroke = 2) +
  geom_errorbar(data = BCI_obs %>% filter(var == "RECRUITMENT"), 
                mapping = aes(ymin = value-sd,ymax = value+sd, width = 0), color = "black") +
  scale_color_continuous() +
  rec.y.axis +
  scale_y_log10(breaks = c(100,1000,10000),
                labels = c(100,1000,10000),
                limits = c(50,10000)) +
  labs(title = "Recruitment") +
  adams_theme +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        legend.position = "none")


BCI_model_data_comp_fig <- plot_grid(modelDataComp_RoL,modelDataComp_RoANPP,
                                     modelDataComp_Recruitment,
          nrow = 1, rel_widths = c(1,1,1.2), labels = c("(a)","(b)","(c)"),
          label_fontface = "bold", label_x = 0, label_size = 25)


makePNG(BCI_model_data_comp_fig,"figures/","BCI_model_data_comp_fig",width = 11,height = 4)




