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

RA_values <- c(0.079)

cases <- c('bci-1pft-fates-main-from-inventory.Ccd4bc7aba-Fa0568083.2021-09-29_RAis0.079','bci-1pft-fates-main-from-inventory.Ccd4bc7aba-Fa0568083.2021-09-29_RAis0.079HgtMinIs2.55')

case_aliases <- c("HgtMin=1.3","HgtMin=2.55")

ModelDataComparisonVars <- c('SEEDS_IN_LOCAL_ELEM','LEAF_LITTER_IN','NPP','RECRUITMENT','NPP_FNRT_SCPF', 'NPP_BGSW_SCPF','NPP_BGDW_SCPF')

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

HgtMinComp <- model_data %>% filter(date > as.Date("1907-01-01")) %>%
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
  rename(value = Mvalue) %>% filter(var == "RECRUITMENT") %>%
  pull(value)

HgtMinRatio <- HgtMinComp[2]/ HgtMinComp[1]



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


modelDataComp_Recruitment <- modelPredictionsBCI %>%
  filter(var == "RECRUITMENT") %>%
  ggplot(aes(case,value,color = RA)) +
  geom_point(size = 7) +
  geom_errorbar(aes(ymin = value-sd,ymax = value+sd, width = 0)) +
  geom_point(data = BCI_obs %>% filter(var == "RECRUITMENT"), 
             mapping = aes(case,value), color ="black", shape = 2, size = 5, stroke = 2) +
  geom_errorbar(data = BCI_obs %>% filter(var == "RECRUITMENT"), 
                mapping = aes(ymin = value-sd,ymax = value+sd, width = 0), color = "black") +
  scale_color_continuous() +
  rec.y.axis +
  scale_y_log10(breaks = c(100,1000,10000,32000),
                labels = c(100,1000,10000,32000),
                limits = c(50,32000)) +
  labs(title = "Recruitment") +
  adams_theme +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        legend.position = "none")



