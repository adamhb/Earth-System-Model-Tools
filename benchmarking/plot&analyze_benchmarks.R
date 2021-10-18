#This script is used to plot/analyze the recruitment data
library(ggplot2)
library(tidyverse)
source('utils/figure_formatting_esm_tools.R')
source('utils/supporting_funcs_esm_tools.R')
source('generalFunctions.R')

vitalRates_allSites <- read_csv("data/vital_rates_all_sites.csv")
spLatinMap <- read_csv('data/spLatinMap.csv')

vitalRates_allSites %>% left_join(spLatinMap,by = "sp") 

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
            M_ha_yr = sum(M_t1,na.rm = T) * m2_per_ha) %>% #the sum of all species area-based mortality rates
  group_by(site) %>%
  summarise(Rm = mean(R_ha_yr), sd = sd(R_ha_yr), n = length(R_ha_yr)) %>%
  rename(R_ha_yr = Rm) %>%
  ungroup() %>%
  mutate(se = sd / sqrt(n))





obs_rec_rates_all_sites <- obs_rec_rates_all_sites_df %>%
  ggplot(aes(site,R_ha_yr)) +
  geom_point(shape = 2, size = 5, stroke = 2) +
  geom_errorbar(aes(ymin = R_ha_yr-se,ymax = R_ha_yr+se, width = 0)) +
  #stat_summary(fun = median, fun.ymax = length,
  #             geom = "text", aes(label = ..ymax..), vjust = -1) +
  #stat_summary(fun.data = give.n, geom = "text", fun = median) +
  ylab(expression(paste("N recruits"," [ha"^"-1"," yr"^"-1","]"))) +
  scale_x_discrete(labels = c('bci\nn = 7','luq\nn = 5','scbi\nn = 2','serc\nn = 1')) +
  xlab(bquote('Site')) +
  scale_y_continuous(limits = c(0,100)) +
  adams_theme


makePNG(fig = obs_rec_rates_all_sites,path_to_figures, file_name = "obs_rec_rates_all_sites", width = 4.5, height = 3.5)

#interesting to see luquillo having a spike in recruitment at beginning
#was this right after a hurricane?
rec_benchmarks_plot_over_time


################################################################
#######Plot 2. Plotting litter flux data########################
################################################################
RoL_SCBI <- read_csv("data/SCBI_RL.csv") %>% add_column(site = "scbi")
RoL_SERC <- read_csv("data/SERC_RL.csv") %>% add_column(site = "serc")

RoLSCBIandSERC <- rbind(RoL_SCBI,RoL_SERC) %>%
  toLongForm(vars = c("Ryr","Lyr","RL")) %>%
  rename(yr = Year) %>%
  select(site,yr,var,value) %>%
  filter(var == "RL") %>%
  mutate(var = "R/L") 


bci_luq_all_obs <- read_csv('data/all_BCI_obs.csv') %>%
  rbind(read_csv('data/all_Luquillo_obs.csv'))

mean_sd_bci_luq <- bci_luq_all_obs %>%
  select(-units) %>%
  drop_na(year) %>%
  spread(var,value) %>%
  mutate(`R/ANPP` = case_when(
    site == "BCI" ~ R / 1800,
    site == "Luquillo" ~ 51 / 1050
  )) %>%
  toLongForm(c('R','L','R/L','R/ANPP')) %>%
  select(-L) %>%
  rename(yr = year)


RoLandRoANPP_obs <- RoLSCBIandSERC %>%
  rbind(mean_sd_bci_luq) %>%
  group_by(site,var) %>%
  summarise(M = mean(value), sd = sd(value), n = length(value)) %>%
  rename(value = M) %>%
  ungroup() %>%
  mutate(site2 = case_when(
    site == 'BCI' ~ "bci",
    site == "Luquillo" ~ "luq",
    TRUE ~ site
  )) %>% select(-site) %>% rename(site = site2) %>%
  mutate(site = factor(site, levels = c("bci",'luq','scbi','serc')),
         se = sd / sqrt(n)) 

RoL_obs_fig <- RoLandRoANPP_obs %>%
  filter(var == "R/L") %>%
  ggplot(aes(site,value)) +
  geom_point(shape = 2, size = 5, stroke = 2) +
  geom_errorbar(aes(ymin = value-se,ymax = value+se, width = 0)) +
  ylab(expression(paste('R / L', ' [g m'^'-2','yr'^'-1',"]"))) +
  xlab(bquote('Site')) +
  scale_y_continuous(limits = c(0,0.7), breaks = seq(from = 0, to = 0.7, by = 0.1)) +
  scale_x_discrete(labels = c('bci\nn = 5','luq\nn = 17','scbi\nn = 5','serc\nn = 5')) +
  adams_theme


RoANPP_obs_fig <- RoLandRoANPP_obs %>%
  filter(var == "R/ANPP") %>%
  ggplot(aes(site,value)) +
  geom_point(shape = 2, size = 5, stroke = 2) +
  geom_errorbar(aes(ymin = value-se,ymax = value+se, width = 0)) +
  ylab(expression(paste('R / ANPP', ' [g m'^'-2','yr'^'-1',"]"))) +
  xlab(bquote('Site')) +
  scale_y_continuous(limits = c(0,0.15)) +
  scale_x_discrete(labels = c('bci\nn = [5,1]','luq\nn = [1,1]')) +
  adams_theme


RoLandRoANPPfig <- plot_grid(RoL_obs_fig,RoANPP_obs_fig, rel_widths = c(1.3,1), labels = c("(a)","(b)"),
          label_fontface = "bold", label_x = 0, label_size = 18)

makePNG(RoLandRoANPPfig,'figures/',file_name = 'RoLandRoANPP_obs_fig',height = 5, width = 8.5)

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

# SI_Table_S5 <- rbind(serc_scbi, bci, luq) %>%
#   mutate(units = case_when(
#     var == "R/L" ~ "-",
#     var %in% c("R","L") ~ "g C m-2 yr-1",
#     TRUE ~ var
#   )) 
read_csv('data/ANPP_ests_BCI_Luquillo.csv')


# SI_Table_S5 <- RoLSCBIandSERC %>%
#   rbind(mean_sd_bci_luq) %>%
#   rbind(tibble(site = "bci", yr = "1969-2000", var = "ANPP", value = 1800)) %>%
#   rbind(tibble(site = "luq", yr = "1969-2000", var = "R", value = 1800)) %>%



SI_Table_S5 <- RoLSCBIandSERC %>% add_column(units = NA) %>%
  rename(year = yr) %>%
  rbind(read_csv('data/all_BCI_obs.csv') %>%
         rbind(read_csv('data/all_Luquillo_obs.csv')))

write_csv(x = SI_Table_S5,path = "data/Table_S5.csv")

# read_csv('data/all_BCI_obs.csv') %>%
#   rbind(read_csv('data/all_Luquillo_obs.csv')) %>%
#   rbind(RoLSCBIandSERC) %>%
#   write_csv('data/Table_S5.csv')


SI_Table_S4 <- vitalRates_allSites %>% left_join(spLatinMap,by = "sp") 
write_csv(SI_Table_S4,'data/Table_S4_v2.csv')



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

