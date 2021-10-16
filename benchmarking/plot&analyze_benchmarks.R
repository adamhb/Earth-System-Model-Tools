#This script is used to plot/analyze the recruitment data
library(ggplot2)
library(tidyverse)
source('utils/figure_formatting_esm_tools.R')
source('utils/supporting_funcs_esm_tools.R')
source('generalFunctions.R')

vitalRates_allSites <- read_csv("data/vital_rates_all_sites.csv")

str(vitalRates_allSites)

m2_per_ha <- 1e4
########################################################
#######Plotting recruitment data########################
########################################################

#################################
#plot 1. Mean recruitment  ######
#################################
obs_rec_rates_all_sites_df <- vitalRates_allSites %>% #turn this into a function?
  group_by(site, cen_interval) %>%
  summarise(R_ha_yr = sum(R_t1,na.rm = T) * m2_per_ha, #the sum of all species area-based recruitment rates
            M_ha_yr = sum(M_t1,na.rm = T) * m2_per_ha) #the sum of all species area-based mortality rates


obs_rec_rates_all_sites <- obs_rec_rates_all_sites_df %>%
  ggplot(aes(site,R_ha_yr)) +
  geom_boxplot() +
  #stat_summary(fun = median, fun.ymax = length,
  #             geom = "text", aes(label = ..ymax..), vjust = -1) +
  #stat_summary(fun.data = give.n, geom = "text", fun = median) +
  ylab(expression(paste("N recruits"," [ha"^"-1"," yr"^"-1","]"))) +
  scale_x_discrete(labels = c('bci\nn = 7','luq\nn = 5','scbi\nn = 2','serc\nn = 1')) +
  xlab(bquote('Site')) +
  scale_y_continuous(limits = c(0,150)) +
  theme_minimal() +
  adams_theme


makePNG(fig = obs_rec_rates_all_sites,path_to_figures, file_name = "obs_rec_rates_all_sites")

#interesting to see luquillo having a spike in recruitment at beginning
#was this right after a hurricane?
rec_benchmarks_plot_over_time


################################################################
#######Plot 2. Plotting litter flux data########################
################################################################
ANPP_df <- read_csv('data/ANPP_ests_BCI_Luquillo.csv')
bci_ANPP <- ANPP_df %>% filter(Site == 'BCI') %>% pull(ANPP) %>% `/`(2) #dividing by 2 to convert biomass to carbon 
luq_ANPP <- ANPP_df %>% filter(Site == 'Luquillo') %>% pull(ANPP) %>% mean()


#read in and clean data for each site
RoANPP_BCI_obs <- read_csv("data/RoverL_BCI_obs.csv") %>%
  filter(var == 'Rgm2yr') %>%
  add_column(ANPP = bci_ANPP) %>%
  mutate(`R/ANPP` = value / ANPP) %>%
  mutate(var = "R/ANPP", value = `R/ANPP`) %>%
  select(case,simYr,var,value) %>% add_column(site = "bci") %>%
  rename(yr = simYr) %>% select(site,yr,var,value)




RoL_BCI_obs <- read_csv("data/RoverL_BCI_obs.csv") %>%
  filter(var == 'RoL') %>% mutate(var = "R/L") %>% add_column(site = "bci") %>%
  rename(yr = simYr) %>% select(site,yr,var,value)

RoL_LUQ_obs <- read_csv("data/RoverL_Luquillo_obs.csv") %>% add_column(site = "luq") %>%
  filter(var == "RoL") %>% mutate(var = "R/L") %>%
  rename(yr = simYr) %>% select(site,yr,var,value)



#update this once I hear from Rachel. Use the new table here.
RoANPP_LUQ_obs <- read_csv("data/RoverL_Luquillo_obs.csv") %>% add_column(site = "luq") %>%
  filter(var == 'Rgm2yr') %>%
  add_column(ANPP = luq_ANPP) %>%
  mutate(`R/ANPP` = value / ANPP) %>%
  mutate(var = "R/ANPP", value = `R/ANPP`) %>%
  select(case,simYr,var,value)
  
RoL_SCBI <- read_csv("data/SCBI_RL.csv") %>% add_column(site = "scbi")
RoL_SERC <- read_csv("data/SERC_RL.csv") %>% add_column(site = "serc")

RoLSCBIandSERC <- rbind(RoL_SCBI,RoL_SERC) %>%
  toLongForm(vars = c("Ryr","Lyr","RL")) %>%
  rename(yr = Year) %>%
  select(site,yr,var,value) %>%
  filter(var == "RL") %>%
  mutate(var = "R/L")


#join and plot data
RoLandRoANPP_obs <- rbind(RoANPP_BCI_obs, RoL_BCI_obs, RoL_LUQ_obs,RoLSCBIandSERC) %>%
  filter(var == "R/L")

#plot R over L
RoL_obs_fig <- RoLandRoANPP_obs %>%
  ggplot(aes(site,value)) +
  geom_boxplot() +
  ylab(expression(paste('R [g m'^'-2','yr'^'-1',"] / ",'L [g m'^'-2','yr'^'-1',"]"))) +
  #stat_summary(fun.data = give.n, geom = "text", fun = median) +
  xlab(bquote('Site')) +
  scale_x_discrete(labels = c('bci\nn = 5','luq\nn = 17','scbi\nn = 5','serc\nn = 5')) +
  scale_y_continuous(limits = c(0,1.5)) +
  theme_minimal() +
  adams_theme

makePNG(RoL_obs_fig,'figures/',file_name = 'RoLandRoANPP_obs_fig',height = 4.5, width = 5)

#plot R over ANPP
rbind(RoANPP_BCI_obs, RoL_BCI_obs, RoL_LUQ_obs,RoLSCBIandSERC) %>%
  #add new Luq data for R here
  filter(var == "R/ANPP")






###############################
####Make table for SI Table S5#
###############################

serc_scbi <- rbind(RoL_SCBI,RoL_SERC) %>%
  rename(R = Ryr, L = Lyr, `R/L` = RL, year = Year) %>%
  select(-site) %>%
  toLongForm(vars = c("R","L","R/L")) %>%
  rename(site = Site)

bci <- read_csv("data/RoverL_BCI_obs.csv") %>% add_column(site = "bci") %>%
  rename(year = simYr) %>% select(-case) %>%
  select(site,year,var,value) %>%
  mutate(var_new = case_when(
    var == "RoL" ~ "R/L",
    var == "Rgm2yr" ~ "R",
    TRUE ~ var
  )) %>%
  mutate(var = var_new) %>%
  select(-var_new) 


luq <- read_csv("data/RoverL_Luquillo_obs.csv") %>% add_column(site = "luq") %>%
  rename(year = simYr) %>% select(-case) %>%
  select(site,year,var,value) %>%
  mutate(var_new = case_when(
    var == "RoL" ~ "R/L",
    var == "Rgm2yr" ~ "R",
    TRUE ~ var
  )) %>%
  mutate(var = var_new) %>%
  select(-var_new) 

SI_Table_S5 <- rbind(serc_scbi, bci, luq) %>%
  mutate(units = case_when(
    var == "R/L" ~ "-",
    var %in% c("R","L") ~ "g C m-2 yr-1",
    TRUE ~ var
  )) 
  
write_csv(x = SI_Table_S5,path = "data/Table_S5.csv")









#################################
#plot 3. Recruitment over time###
#################################
#NOTE: This should have date on the x axis. Need to calcuate mean cen int date for that.

#preparing recruitment data for plotting
df_plot_over_time <- vitalRates_allSites %>% #turn this into a function?
  group_by(site, cen_interval) %>%
  summarise(R_ha_yr = sum(R_t1,na.rm = T) * m2_per_ha, #the sum of all species area-based recruitment rates
            M_ha_yr = sum(M_t1,na.rm = T) * m2_per_ha) #the sum of all species area-based mortality rates

#plotting the recruitment observations
rec_benchmarks_plot_over_time <- df_plot_over_time %>% 
  ggplot(mapping = aes(x = cen_interval, y = M_ha_yr, color = site)) + #get the census interval year to line up on this plot instead of the census number
  geom_point(size = 5) +
  geom_line() +
  ylab(label = "recruitment rate (# ind. per ha yr)") +
  xlab(label = "census interval") +
  theme_minimal()

