#load packages
library(tidyverse)
library(ncdf4)
library(lubridate)

#load supporting functions
source('utils/supporting_funcs_esm_tools.R')
source('utils/figure_formatting_esm_tools.R')
source('generalFunctions.R')

############
#functions##
############
getSiteVar <- function(file,varName){ #input are netcdf files
  nc <- nc_open(file)
  r <- ncvar_get(nc,varName)
  nc_close(nc)
  return(r) #returns variable value of interest
}


getDateFromNC <- function(file){
  return(ymd(getSiteVar(file,"mcdate")))
}

getSiteVarOverTime <- function(files,varName,withDate=F){ #accepts a list of netcdf output files
  output <- c()
  date <- c()
  for(f in files){
    if(withDate == T){
      date <- append(date,getDateFromNC(f))
      }
    output <- append(output,getSiteVar(f,varName))
  }
  if(withDate == T){return(tibble(date,output))}else{
    return(output)
  }
}

getYrFromDate <- function(date){
  as.numeric(substr(date,1,4))
}

getSimYrFromDate <- function(date,start_date){
  getYrFromDate(date) - getYrFromDate(start_date) + 1
}




getManySiteVarsOverTime <- function(files,varNames){
  output <- matrix(nrow = length(files),ncol = length(varNames))
  c <- 0
  
  for(v in varNames){
    c <- c + 1
    output[,c] <- getSiteVarOverTime(files,v,F)
    print(paste("Done with variable",v,"of",length(varNames),"variables!"))
  }
  colnames(output) <- varNames
  
  dates <- c()
  for(f in files){
    datef <- getDateFromNC(f)
    dates <- append(dates,datef)
  }
  
  #add the date to the data
  output <- as_tibble(output) %>% add_column(date = dates)
  #add the simulation year to the data
  start_date <- output$date[1]
  output <- output %>% mutate(simYr = getSimYrFromDate(date,start_date))
  return(output)
}

plotVar <- function(d,x,y){
  p <- ggplot(data = d, mapping = aes_string(x,y)) +
    geom_point() +
    theme_minimal() +
    adams_theme
  return(p)
}


toLongForm <- function(d,vars){
    gather(data = d,as.name(head(vars,1)):as.name(tail(vars,1)),
           key = "var",value = "value")
}

toWideForm <- function(d,var_names_col,value_col){
  wide <- pivot_wider(data = d, 
                      names_from = as.name(var_names_col), 
                      values_from = as.name(value_col))
  return(wide)
}



plotGridofVars <- function(d,y_vars,x_var,y_axis_units,prettyNames){
  p <- d %>% select(date,simYr,y_vars) %>%
    toLongForm(y_vars) %>% 
    left_join(prettyNames,by = "var") %>%
    ggplot(aes_string(x_var,"value")) +
    geom_point() +
    theme_minimal() +
    facet_wrap(~prettyName,scales = "free") +
    ylab(y_axis_units) +
    xlab("SimYear") +
    scale_x_date(date_labels = "%y", date_breaks = "5 years") +
    adams_theme
  return(p)
  }

#this function receives a dataframe of variables (as columns)
#and standardizes their units based on a "variable unit map" (see script 'FATES_variable_units_map.R')
#It returns the dataframe with all variables in standard units (g m-2 day-1)

standardize_units <- function(d,vars,unit_map){
  
 
  
  output <- toLongForm(d,vars) %>%
    left_join(unit_map, by = "var") %>%
    mutate(std_mass = case_when(
      (value != 0 & grepl(' kg ',units)) ~ value*g_per_kg,
      TRUE ~ value
    )) %>%
    mutate(std_area = case_when(
      (value != 0 & grepl(' ha-1 ',units)) ~ std_mass*has_per_m2,
      TRUE ~ std_mass
    )) %>%
    mutate(std_time = case_when(
      (value != 0 & grepl(' yr-1 ',units)) ~ std_area*yrs_per_day,
      (value != 0 & grepl(' s-1 ',units)) ~ std_area*sec_per_day,
      TRUE ~ std_area
    )) %>%
    mutate(value = std_time) %>%
    select(simYr,date,var,value) 
  return(toWideForm(output,"var","value"))
}

#this function receives aata frame and returns the mean value of a variable 
#for each year of the data set
#NOTE: could add standard deviation to this
getAnnualMeanAfterSpecificYr <- function(df,var,start_date_str){
  df %>% filter(date > as.Date(start_date_str)) %>%
    select(simYr,date,as.name(var)) %>%
    toLongForm(vars = var) %>% 
    group_by(simYr) %>%
    summarise("value" = mean(value))
}



##############################################
#####Create dataframe of vdm data over time###
##############################################
data_path <- '~/cloud/gdrive/FATES/FATES_data'
case_name <- 'bci-1pft-fates-main-from-inventory.Ccd4bc7aba-Fa0568083.2021-09-29_RAis0.04680188'
files <- getFilePaths(data_path,case_name)

#variables of interest
#these vars have one dimension when there is one pft
siteLevelVars <- c('SEEDS_IN_LOCAL_ELEM','LEAF_LITTER_IN','NPP','AGB','RECRUITMENT','PFTbiomass','PFTnindivs','PFTnpp')

#create a map of the units of each variable
varUnitsMapFates <- tibble(
  var = c('SEEDS_IN_LOCAL_ELEM','LEAF_LITTER_IN','NPP','AGB','RECRUITMENT','PFTbiomass','PFTnindivs','PFTnpp'),
  units = c(' kg ha-1 d-1 ',' g m-2 s-1 ',' g m-2 s-1 ',' g m-2 ',' indiv ha-1 yr-1 ',' g m-2 ',' indiv m-2 ',' kg m-2 yr-1 ')
)

FatesPrettyNamesMap <- tibble(
  var = c('SEEDS_IN_LOCAL_ELEM','LEAF_LITTER_IN','NPP','AGB','RECRUITMENT','PFTbiomass','PFTnindivs','PFTnpp'),
  prettyName = c('ReproCFlux','LeafCFlux','NPP','AGB','Recruits','PFTbiomass','PFT_N','PFT_NPP')
)


#put 1D site vars into a df where all units are in g m-2 d-1##
VDM_data <- standardize_units(getManySiteVarsOverTime(files,siteLevelVars),siteLevelVars,varUnitsMapFates)




##########
#figures##
##########

###########
#plot AGB##
###########
AGB_plot <- plotVar(VDM_data,"date","AGB") +
  ylab(label = "g C m-2") +
  labs(title = "AGB")
makePNG(AGB_plot,'figures/',"AGB")

##############################################
#plot all site-level flux variables over time#
##############################################
#This plot is used mostly to check that fluxes reached dynamic equilibrium
y_axis_flux_units <- expression(paste("g C"," [m"^"-2"," day"^"-1","]"))

fluxVarsPlot <- plotGridofVars(VDM_data,
               c("NPP","LEAF_LITTER_IN","SEEDS_IN_LOCAL_ELEM"),
               "date",
               y_axis_flux_units,prettyNames = FatesPrettyNamesMap)


makePNG(fluxVarsPlot,"figures/","fluxVarsPlot")

############################
#plot recruitment over time#
############################
rec_y_axis_units <- expression(paste("N recruits"," [ha"^"-1"," yr"^"-1","]"))

rec_plot <- plotVar(VDM_data %>%  #making recruitment in units of ha per yr for the plot
                      mutate_at(.vars = "RECRUITMENT",.funs = function(x){x*m2_per_ha*days_per_yr}),
                    "date","RECRUITMENT") +
            ylab(rec_y_axis_units) +
            labs(title = "RECRUITMENT")
makePNG(rec_plot,'figures/',"RECRUITMENT")


#############################################################################
##calculate mean recruitment after simulation reaches dynamics equilibrium###
#############################################################################
#first we get the mean of the year and then plot a box plot which shows variation in inter-annual variation

meanRECperYrFATES <- getAnnualMeanAfterSpecificYr(VDM_data,"RECRUITMENT",start_date_str = "1905-01-01")
write_csv(meanRECperYrFATES,"data/FATES-recruitment-predictions")



meanRECRUITMENTplot <- meanRECperYrFATES %>% 
  ggplot(aes("",value*m2_per_ha*days_per_yr)) +
  #facet_wrap(~var,scales = "free") +
  geom_boxplot() +
  geom_point(data = tibble(x = NA, value = 75 * has_per_m2 * yrs_per_day), shape = 17, size = 5) +
  ylab(rec_y_axis_units) +
  xlab("FATES Recruitment BCI \n triangle = obs") +
  theme_minimal() +
  adams_theme
makePNG(meanRECRUITMENTplot,'figures/',"MEAN_RECRUITMENT")

############################################
#Calculate and plot R over L and R over NPP#
############################################

RoLandRoNPP_df <- VDM_data %>% filter(date > as.Date("1907-01-01")) %>%
  mutate(RoL = SEEDS_IN_LOCAL_ELEM / LEAF_LITTER_IN,
         RoNPP = SEEDS_IN_LOCAL_ELEM / NPP) %>%
  select(simYr,date,RoL,RoNPP) %>%
  toLongForm(vars = c('RoL','RoNPP')) %>%
  group_by(simYr,var) %>%
  summarise(value = mean(value), n = length(value))

RoL_and_RoNPP_plot <- RoLandRoNPP_df %>%
  ggplot(aes(var,value)) +
  facet_wrap(~var,scales = "free") +
  geom_boxplot() +
  theme_minimal() +
  adams_theme

makePNG(RoL_and_RoNPP_plot,'figures/',"RoL&RoNPP")


#SCRATCH
#these variables have multiple dimensions
multidim_vars <- c('NPATCH_BY_AGE','DDBH_CANOPY_SCPF')








