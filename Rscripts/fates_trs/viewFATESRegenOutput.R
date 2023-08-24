

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
case_name <- 'bci-12pfts-fates-API20-TRS-12-17-2021.C2c05f5dc6-F6c38a6b9.2021-12-17_fates-trs'

#case_alias <- "Ris0.079"
case_alias <- "fates-trs"
files <- getFilePaths(data_path,case_name)


vars <- c('FATES_NPP_PF','FATES_RECRUITMENT_PF','FATES_NPLANT_PF','FATES_PATCHAREA_AP','FATES_MORTALITY_PF','FATES_NPATCHES','FATES_NCOHORTS','FATES_VEGC_PF','FATES_MORTALITY_CANOPY_SZAP','FATES_MORTALITY_USTORY_SZAP','FATES_SEED_BANK_TRS','FATES_SEEDLING_POOL_TRS','FATES_PARPROF_DIR_CLLL','FATES_PARPROF_DIF_CLLL','FATES_LAISUN_Z_CLLL','FATES_LAISHA_Z_CLLL','FATES_DDBH_CANOPY_SZAP','FATES_DDBH_USTORY_SZAP','FATES_SEED_PROD_USTORY_SZ','FATES_SEED_PROD_CANOPY_SZ','FATES_SEEDS_IN','FATES_SEEDS_IN_LOCAL')


# 'FATES_RECRUITMENT_PF','FATES_NPLANT_PF','FATES_PATCHAREA_AP','FATES_MORTALITY_PF','FATES_NPATCHES','FATES_NCOHORTS','FATES_VEGC_PF','FATES_MORTALITY_CANOPY_SZAP','FATES_MORTALITY_USTORY_SZAP','FATES_PARPROF_DIR_CLLL','FATES_PARPROF_DIF_CLLL','FATES_LAISUN_Z_CLLL','FATES_LAISHA_Z_CLLL','FATES_DDBH_CANOPY_SZAP','FATES_DDBH_USTORY_SZAP','FATES_SEED_PROD_USTORY_SZ','FATES_SEED_PROD_CANOPY_SZ','FATES_SEEDS_IN'


site_level_vars_index <- c()
j <- 0
for (var in vars){
  j <- j + 1
  site_level_vars_index[j] <- length(ncvar_get(nc = nc_open(files[1]),varid = var))
}
  
#process site-level variables
siteLevelVars <- vars[site_level_vars_index == 1]

#put 1D site vars into a df where all units are in g m-2 d-1##
#VDM_data <- getManySiteVarsOverTime(files,siteLevelVars)
#write_csv(VDM_data, "fates-trs-output-data.csv")
VDM_data <- read_csv("fates-trs-output-data.csv")

VDM_data_v2 <- VDM_data %>% select(date,FATES_NPP_PF,FATES_SEEDS_IN,
                                   FATES_SEED_BANK_TRS,FATES_SEEDLING_POOL_TRS,
                                   FATES_RECRUITMENT_PF) %>%
  mutate(model = "fates-trs")

TRS_data <- read_csv('~/cloud/gdrive/rec_submodel/regeneration_submodel/temp/forComparisonToFATES.csv') %>%
  rename(FATES_NPP_PF = NPP,FATES_SEEDS_IN = c_repro,
         FATES_SEED_BANK_TRS = seedbank,FATES_SEEDLING_POOL_TRS = seedpool,
         FATES_RECRUITMENT_PF = R) %>%
  mutate(model = "offline-trs") %>%
  mutate_at(.vars = "date",.funs = function(x){x - (365*100) + 6})

VDM_data_v3 <- rbind(VDM_data_v2,TRS_data) %>%
  filter(date < as.Date("1911-01-01"))


##############################
###FATES-TRS vs. offline TRS##
##############################

#Plot NPP
NPP <- VDM_data_v3 %>%
  filter(date > as.Date("1901-01-01")) %>%
  ggplot(aes(date,FATES_NPP_PF*sec_per_day*days_per_yr*m2_per_ha,color = model)) +
  geom_line(size = 1) +
  scale_x_date(date_labels = "%Y",date_breaks = "1 year") +
  ylab("NPP [kg C ha-1 yr-1]") +
  theme_minimal()

#Plot carbon allocated to reproduction
C_for_repro <- VDM_data_v3 %>%
  filter(date > as.Date("1901-01-01")) %>%
  ggplot(aes(date,FATES_SEEDS_IN*sec_per_day*days_per_yr*m2_per_ha,color = model)) +
  geom_line(size = 1) +
  scale_x_date(date_labels = "%Y",date_breaks = "1 year") +
  ylab("C for repro. [kg C ha-1 yr-1]") +
  theme_minimal()

seedbank <- VDM_data_v3 %>%
  filter(date > as.Date("1901-01-01")) %>%
  ggplot(aes(date,FATES_SEED_BANK_TRS*m2_per_ha,color = model)) +
  geom_line(size = 1) +
  scale_x_date(date_labels = "%Y",date_breaks = "1 year") +
  ylab("NPP [kg C ha-1]") +
  theme_minimal()
  
seedlingPool <- VDM_data_v3 %>%
  filter(date > as.Date("1901-01-01")) %>%
  ggplot(aes(date,FATES_SEEDLING_POOL_TRS*m2_per_ha,color = model)) +
  geom_line(size = 1) +
  scale_x_date(date_labels = "%Y",date_breaks = "1 year") +
  ylab("NPP [kg C ha-1]") +
  theme_minimal()

recruitment <- VDM_data_v3 %>%
  filter(date > as.Date("1901-01-01")) %>%
  ggplot(aes(date,FATES_RECRUITMENT_PF*m2_per_ha,color = model)) +
  geom_line(size = 1) +
  scale_x_date(date_labels = "%y",date_breaks = "1 year") +
  ylab("N recruits [ha-1 yr-1]") +
  xlab("simulation year") +
  geom_hline(aes(yintercept = 70, linetype = "BCI obs."),color = "black") +
  scale_color_manual(values = c("olivedrab","dodgerblue")) +
  adams_theme +
  scale_linetype_manual(name = "", values = c(2), 
                        guide = guide_legend(override.aes = list(color = c("black"))))
  
makePNG(fig = recruitment,path_to_output.x = 'figures/',file_name = "trs-fates-vs-offline.png")



##############################
###FATES-TRS vs. FATES-default
##############################

case_name <- 'bci-12pfts-fates-main-API20-12-16-2021.C2c05f5dc6-Fb2714927.2021-12-16_fates-default'
case_alias <- "fates-default"
files <- getFilePaths(data_path,case_name)

recruitment_default_FATES <- getManySiteVarsOverTime(files,varNames = c('FATES_VEGC_PF','FATES_RECRUITMENT_PF','FATES_SEEDS_IN'))
recruitment_default_FATES <- recruitment_default_FATES %>% add_column(model = "fates-default") %>%
  select(date,simYr,model,FATES_VEGC_PF,FATES_RECRUITMENT_PF,FATES_SEEDS_IN) 
  

FATES_TRS <- VDM_data %>%
  mutate(model = "fates-trs") %>%
  select(date,simYr,model,FATES_VEGC_PF,FATES_RECRUITMENT_PF,FATES_SEEDS_IN) 

FATES_comparison <- rbind(FATES_TRS,recruitment_default_FATES)


FATES_recruitment_comparison <- FATES_comparison %>%
  filter(date > as.Date("1901-01-01"), simYr < 11) %>%
  ggplot(aes(date,FATES_RECRUITMENT_PF*m2_per_ha,color = model)) +
  geom_line(size = 1) +
  scale_x_date(date_labels = "%y",date_breaks = "1 year") +
  scale_y_log10() +
  ylab("N recruits [ha-1 yr-1]") +
  xlab("simulation year") +
  geom_hline(aes(yintercept = 70, linetype = "BCI obs."),color = "black") +
  scale_color_manual(values = c("black","olivedrab")) +
  adams_theme +
  scale_linetype_manual(name = "", values = c(2), 
                        guide = guide_legend(override.aes = list(color = c("black"))))

makePNG(fig = FATES_recruitment_comparison,path_to_output.x = 'figures/',file_name = "trs-fates-vs-fates-default.png")

FATES_C_alloc_repro_comparison <- FATES_comparison %>%
  filter(date > as.Date("1901-01-01"), simYr < 11) %>%
  ggplot(aes(date,FATES_SEEDS_IN*sec_per_day*days_per_yr*m2_per_ha,color = model)) +
  geom_line(size = 1) +
  scale_x_date(date_labels = "%y",date_breaks = "1 year") +
  scale_y_log10() +
  ylab("C for repro. [kg C ha-1 yr-1]") +
  xlab("simulation year") +
  #geom_hline(aes(yintercept = 70, linetype = "BCI obs."),color = "black") +
  scale_color_manual(values = c("black","olivedrab")) +
  adams_theme +
  scale_linetype_manual(name = "", values = c(2), 
                        guide = guide_legend(override.aes = list(color = c("black"))))



getMultiDimSiteVar <- function(file,varName){ #input are netcdf files
  nc <- nc_open(file)
  r <- ncvar_get(nc,varName)
  #if(length(r) > 1){r <- sum(tail(r,3))} #temporary change
  nc_close(nc)
  return(r) #returns variable value of interest
}

getMultiDimSiteVar(files[100],"FATES_DDBH_USTORY_SZAP")

nc_open(files[1])

FATES_DDBH_USTORY_SZAP
#scratch
growth_rates_und <- getManySiteVarsOverTime(files,'FATES_DDBH_USTORY_SZAP')
ncvar_get(nc = nc_open(files[1]),varid = 'FATES_DDBH_USTORY_SZAP')
