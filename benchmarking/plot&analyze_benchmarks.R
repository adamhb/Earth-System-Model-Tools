#This script is used to plot/analyze the recruitment data
library(ggplot2)
library(tidyverse)
source('utils/figure_formatting_esm_tools.R')
source('utils/supporting_funcs_esm_tools.R')

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
  stat_summary(fun = median, fun.ymax = length,
               geom = "text", aes(label = ..ymax..), vjust = -1) +
  #stat_summary(fun.data = give.n, geom = "text", fun = median) +
  ylab(expression(paste("N recruits"," [ha"^"-1"," yr"^"-1","]")))+
  xlab(bquote('Site'))+
  theme_minimal() +
  adams_theme


makePNG(fig = obs_rec_rates_all_sites,path_to_figures, file_name = "obs_rec_rates_all_sites")


#################################
#plot 2. Recruitment over time###
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

#interesting to see luquillo having a spike in recruitment at beginning
#was this right after a hurricane?
rec_benchmarks_plot_over_time


########################################################
#######Plotting litter flux data########################
########################################################

#Do R/L and R/ANPP here.

#I already have r/l for BCI, need to get these metrics for other sites
read_csv('data/BCI_mo_RL.csv')
obs_data_RoL <- read_csv("data/RoverL_BCI_obs.csv") %>% add_column(site = "bci")



  

