#This script receives a set of netcdf files from one model simulation
#and plots the output for reproductive litter flux, leaf litter flux,
#NPP, R/L, and R/NPP.
#I use this for the default parameter set

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
#case_name <- 'bci-1pft-fates-main-from-inventory.Ccd4bc7aba-Fa0568083.2021-09-29_RAis0.079'
case_name <- 'bci-1pft-fates-main-from-inventory.Ccd4bc7aba-Fa0568083.2021-09-29_RAis0.149'

#case_alias <- "Ris0.079"
case_alias <- "Ris0.149"
files <- getFilePaths(data_path,case_name)



##############################################
#####Create dataframe of vdm data over time###
##############################################

#specify the variables of interest
siteLevelVars <- c('SEEDS_IN_LOCAL_ELEM','LEAF_LITTER_IN','NPP','AGB','RECRUITMENT','PFTbiomass','PFTnindivs','PFTnpp','NPP_FNRT_SCPF', 'NPP_BGSW_SCPF','NPP_BGDW_SCPF')


#put 1D site vars into a df where all units are in g m-2 d-1##
VDM_data <- standardize_units(getManySiteVarsOverTime(files,siteLevelVars),siteLevelVars,varUnitsMapFates)
VDM_data <- VDM_data %>% mutate(ANPP = NPP - (NPP_FNRT_SCPF + NPP_BGSW_SCPF + NPP_BGDW_SCPF))

#####################
#plot AGB over time##
#####################
AGB_plot <- plotVar(VDM_data,"date","AGB") +
  ylab(label = "g C m-2") +
  labs(title = "AGB")
makePNG(AGB_plot,'figures/',pasteo("AGB_",case_alias))

##############################################
#plot all site-level flux variables over time#
##############################################

#This plot is used mostly to check that fluxes reached dynamic equilibrium
y_axis_flux_units <- expression(paste("g C"," [m"^"-2"," day"^"-1","]"))

fluxVarsPlot <- plotGridofVars(VDM_data,
               c("ANPP","LEAF_LITTER_IN","SEEDS_IN_LOCAL_ELEM"),
               "date",
               y_axis_flux_units,prettyNames = FatesPrettyNamesMap)


makePNG(fluxVarsPlot,"figures/",paste0("fluxVarsPlot_",case_alias))

############################
#plot recruitment over time#
############################
rec_y_axis_units <- expression(paste("N recruits"," [ha"^"-1"," yr"^"-1","]"))

rec_plot <- plotVar(VDM_data %>%  #making recruitment in units of ha per yr for the plot
                      mutate_at(.vars = "RECRUITMENT",.funs = function(x){x*m2_per_ha*days_per_yr}),
                    "date","RECRUITMENT") +
            ylab(rec_y_axis_units) +
  adams_theme
makePNG(rec_plot,'figures/',paste0("RECRUITMENT_",case_alias))


#############################################################################
##calculate mean recruitment after simulation reaches dynamics equilibrium###
#############################################################################

#first we get the mean of each year  
meanRECperYrFATES <- getAnnualMeanAfterSpecificYr(VDM_data,"RECRUITMENT",start_date_str = "1905-01-01")



#and then plot a box plot which shows inter-annual variation
meanRECRUITMENTplot <- meanRECperYrFATES %>% 
  ggplot(aes("",value*m2_per_ha*days_per_yr)) +
  #facet_wrap(~var,scales = "free") +
  geom_boxplot() +
  geom_point(data = tibble(x = NA, value = 75 * has_per_m2 * yrs_per_day), shape = 17, size = 5) +
  ylab(rec_y_axis_units) +
  xlab("FATES Recruitment BCI \n triangle = obs") +
  theme_minimal() +
  adams_theme
makePNG(meanRECRUITMENTplot,'figures/',paste0("MEAN_RECRUITMENT_",case_alias))

############################################
#Calculate and plot R over L and R over NPP#
############################################
RoLandRoANPP_df <- VDM_data %>% filter(date > as.Date("1907-01-01")) %>%
  mutate(RoL = SEEDS_IN_LOCAL_ELEM / LEAF_LITTER_IN,
         RoANPP = SEEDS_IN_LOCAL_ELEM / ANPP) %>%
  select(simYr,date,RoL,RoANPP) %>%
  toLongForm(vars = c('RoL','RoANPP')) %>%
  group_by(simYr,var) %>%
  summarise(value = mean(value), n = length(value))


#plot the model predictions
RoL_and_RoANPP_plot <- RoLandRoANPP_df %>%
  ggplot(aes(var,value)) +
  facet_wrap(~var,scales = "free") +
  labs(title = case_alias) +
  geom_boxplot() +
  theme_minimal() +
  adams_theme

makePNG(RoL_and_RoANPP_plot,'figures/',paste0("RoL&RoNPP_",case_alias))


#SCRATCH
#these variables have multiple dimensions
#multidim_vars <- c('NPATCH_BY_AGE','DDBH_CANOPY_SCPF')








