
library(tidyverse)
library(lubridate)
library(ggpubr)

source('utils/figure_formatting_esm_tools.R')

########################### Correlation between R/L and R/NPP, ForC database
#Anderson-Teixeira, K.J., Wang, M.M.H., McGarvey, J.C., Herrmann, V., Tepley, A.J., Bond-Lamberty, B. and LeBauer, D.S. (2018), 
#ForC: a global database of forest carbon stocks and fluxes. Ecology, 99: 1507-1507. https://doi.org/10.1002/ecy.2229

#git hub: https://github.com/forc-db/ForC

##########################Script to format For C data

path <- "~/Desktop/NGEE-Tropics/Repro_metadata_synthesis/ForC_v2.0 2/"
`%notin%` <- Negate(`%in%`)

variables <- read_csv(paste0(path,"ForC_variables.csv"))
site_info<- read_csv(paste0(path, "ForC_sites.csv"))
plot_info <- read_csv(paste0(path, "ForC_plots.csv"))

site_id <- site_info %>% 
  select(site.ID,sites.sitename)

plot_id <- plot_info %>%
  select(plot.ID, sites.sitename, plot.name)


ForC_msmts <- read_csv(paste0(path,"ForC_measurements.csv")) %>% 
  left_join(plot_id, by = c("sites.sitename", "plot.name"))

###############select variables of interest
repro_vars <- c("ANPP_repro_OM", "ANPP_repro_C")
leaf_vars <- c("ANPP_foliage_OM", "ANPP_foliage_C") # leaf or needle

NPP_vars <- c("NPP_3_OM", "NPP_3_C", #these include repro structures
              "NPP_4_OM", "NPP_4_C", # + leaf herbivory (folivory)
              "NPP_5_OM", "NPP_5_C") # +  VOCs & exudation

#note ANPP methods include foliage + branch + stem components, not repro structures, not including here


############### function to extract data in useable format
get_var_df<- function(data, var_names){
  #data = ForC measurment data
  #var_names = character vector of variable names to include in df
  data %>% 
    filter(variable.name %in% var_names) %>% 
    select(plot.ID, measurement.ID, sites.sitename, plot.name, stand.age, variable.name, mean) %>%
    pivot_wider(names_from = variable.name, values_from = mean) 
}

############slight adjustments
repro_df <- get_var_df(ForC_msmts,repro_vars) %>% 
  pivot_longer(cols = c(contains("ANPP_repro")), names_to = "R_method", values_to = "R_flux") %>% 
  filter(!is.na(R_flux))

leaf_df <- get_var_df(ForC_msmts, leaf_vars) %>% 
  pivot_longer(cols = c(contains("ANPP_foliage")), names_to = "L_method", values_to = "L_flux") %>% 
  filter(!is.na(L_flux))

npp_df <- get_var_df(ForC_msmts, NPP_vars) %>% 
  pivot_longer(cols = c(contains("NPP_")), names_to = "NPP_method", values_to = "NPP_flux") %>% 
  filter(!is.na(NPP_flux)) 

########### specify variables to use to join dataframes
joining_vars<- c("plot.ID", "sites.sitename", "plot.name", "stand.age")

df <- repro_df %>% 
  left_join(leaf_df, by = joining_vars) %>% 
  #left_join(anpp_df, by = joining_vars) %>% 
  left_join(npp_df, by = joining_vars) %>% 
  select(- c(starts_with("Measur")))  %>% # remove measurment.id cols
  mutate(RL = R_flux / L_flux,
         #RANPP = R_flux / (ANPP_flux+ R_flux), # Repro cateogry not included in ANPP_3_ method (foliage, stem, branch)
         RNPP = R_flux / NPP_flux) %>% 
  #filter(RL > 1) cocoflux coconut plantation 
  filter(RL < 2) %>% #filtering out coconut plantations
  rowwise() %>% # the below code puts all leaf and npp fluxes in like terms (gC); these were used for exploratory work and 
  # are not currently used in anlaysis of this data
  mutate(L_flux_C = case_when(
    str_detect(L_method,"_C") ~ L_flux, 
    str_detect(L_method, "_OM") ~ (L_flux/2) # assuming .5*mass approx ~ massC
  )) %>% 
  mutate(NPP_flux_C = case_when(
    str_detect(NPP_method,"_C") ~ NPP_flux, 
    str_detect(NPP_method, "_OM") ~ (NPP_flux/2) # assuming .5*mass approx ~ massC
  )) 


write_csv(df, file = "~/Desktop/NGEE-Tropics/Earth-System-Model-Tools/data/ForCformatted.csv")

