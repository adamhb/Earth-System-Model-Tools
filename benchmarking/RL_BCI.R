
library(tidyverse)
library(lubridate)
source('utils/figure_formatting_esm_tools.R')
source('benchmarking/RL_functions.R')

##############################Ecosystem level R/L , all species###########################

data_path <- "data/"

BCI_mo_RL <- read_csv(paste(data_path,"BCI_mo_RL.csv", sep = ""))

BCI_mo_RL %>% 
  group_by(Year) %>% 
  summarise(Rm2yr = sum(Rm2), 
            Lm2yr = sum(Lm2),
            Sm2yr = sum(Sm2)) %>% 
  mutate(RL = Rm2yr/Lm2yr) %>% 
  ggplot() +
  geom_boxplot(aes(x = "", 
                   y = RL))+
  xlab("BCI R/L (all species)") +
  adams_theme

############################## Canopy sp  R/L##############################################


##load data
canopy_sp <- read_csv(paste(data_path, "canopy_sp_bci.csv", sep = "")) %>% 
  rename(sp = `unique(bcifull_clean$sp)`)

bci_traits <- read_csv("~/Desktop/NGEE-Tropics/BCI_data/BCITRAITS_20101220_updated3.21.19.csv") 


trap_data <- read_csv("~/Desktop/NGEE-Tropics/BCI_data/Hanbury_Brown_50ha_rawdata_20190404.csv") %>% 
  mutate(part = as_factor(part)) %>% 
  mutate(sp = tolower(sp6)) %>% 
  mutate(Date = as.Date(fecha, format = "%Y%m%d")) %>% 
  mutate(Year = lubridate::year(Date)) %>% 
  filter(sp %in% canopy_sp$sp) %>% # filter for only canopy sp
  filter(Year>2013 & Year < 2019) ## full years of data

fruit_mass_df <- bci_traits %>% 
  mutate(sp = tolower(`SP$`)) %>% 
  dplyr::select(sp, FRUIT_DRY, DSPR_DRY) %>%
  mutate(fruit_mass  = case_when(
    !is.na(FRUIT_DRY) ~ FRUIT_DRY, 
    is.na(FRUIT_DRY)  ~ DSPR_DRY
  ))

##designate part codes for reproductive and leaf material 
repro_parts <- c(0:10) #part 10 is not in metadtaa seed rain file but from pers. communication with Joe Wright, represents damaged fruit
leaf_parts <- 11
#################

##helper function to use with purr library functions 
purr_fun <- function(yrdata){
  yrdata$mass_m2_total
} 

## RL for canopy species with no herbivory correction 
R_yr <- 
  trap_data %>% 
  group_by(Year) %>%
  nest() %>% 
  mutate(R_yr = map(data, get_part_mass_total, parts = repro_parts)) %>% 
  transmute(Year, Rmass_m2 = map_dbl(R_yr, purr_fun))


L_yr <- 
  trap_data %>% 
  group_by(Year) %>% 
  nest() %>% 
  mutate(L_yr = map(data, get_part_mass_total, parts = leaf_parts)) %>% 
  transmute(Year, Lmass_m2 = map_dbl(L_yr, purr_fun))


## RL for canopy species with herbivory correction 
R_yr_corr_temp<- 
  trap_data %>% 
  group_by(Year) %>%
  nest() %>% 
  mutate(R_yr_corr = map(data, get_repro_total_m2)) %>% 
  transmute(Year, Rmass_m2_corr = map_dbl(R_yr_corr, purr_fun))

L_yr_corr <- L_yr %>% 
  mutate(Lmass_m2_corr = Lmass_m2*1.11)

##leaf herbivory correction source: https://daac.ornl.gov/NPP/guides/NPP_BRR.html
## see brr_npp_r1.txt, estsimated ~ 50g/m2yr lost to insect herbivory, ~30 g/m2yr to vert. herbiv.
## this is proprtional to ~ 11% of the observed leaf litterfall flux


## Join 
R_yr_corr <- R_yr %>%
  left_join(R_yr_corr_temp)

RL_df <- R_yr_corr %>% 
  left_join(L_yr_corr) %>% 
  mutate(RL = Rmass_m2 / Lmass_m2 ) %>% 
  mutate(RL_corr = Rmass_m2_corr / Lmass_m2_corr) %>% 
  select(Year, RL, RL_corr) %>% 
  pivot_longer(cols = c(RL, RL_corr), names_to = "Type", values_to = "RL") 


## Plot
RL_df %>% 
  ggplot() +
  geom_boxplot(aes(x = Type,
                   y = RL)) +
  adams_theme 




