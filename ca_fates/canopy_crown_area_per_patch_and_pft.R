#install libraries
library(ncdf4)
library(tidyverse)
library(lubridate)

data_path <- '/home/adam/cloud/gdrive/postdoc/simulation_output/CZ2_I2000Clm51Fates_pine_cedar_fir_montane_shrub_102322'

getFilePaths <- function(base_data_path = data_path){
  path <- paste0(data_path,'/')
  fileNames <- list.files(path,pattern = '.nc')
  files <- paste0(path,fileNames)
  return(files)
}

files <- getFilePaths(data_path)

getMultiDimSiteVar <- function(file,varName){ #input are netcdf files
  nc <- nc_open(file)
  r <- ncvar_get(nc,varName)
  #if(length(r) > 1){r <- sum(tail(r,3))} #temporary change
  nc_close(nc)
  return(r) #returns variable value of interest
}

#create patch age and pft variables
pft <- c()
for(i in c("pine","cedar","fir","shrub")){
  pft <- c(pft,rep(i,7))
}
ages <- c(0,1,2,5,10,20,50)
age_dim <- rep(ages,4)


getDateFromNC <- function(file){
  return(ymd(getSiteVar(file,"mcdate")))
}


getCanopyCrownArea <- function(file){
  date <- getDateFromNC(file)
  mydata <- nc_open(file)
  ca <- ncvar_get(mydata, "FATES_CANOPYCROWNAREA_APPF")
  patch_area <- rep(ncvar_get(mydata, "FATES_PATCHAREA_AP"),4)
  data <- tibble(age = age_dim, pft = pft, cca = ca, date = date, patch_area = patch_area) %>%
    mutate(pft = factor(pft, levels = c("pine","cedar","fir","shrub"))) %>%
    mutate_at(.vars = "age", .funs = as.factor) %>%
    mutate(cca_adj = cca / patch_area)
  nc_close(mydata)
  return(data)
}


#dates <- map(files,getDateFromNC)


canopy_crown_area_all_dates <- c()
for(f in files[(50*12):(75:12)]){
  tmp  <- getCanopyCrownArea(f)
  canopy_crown_area_all_dates <- rbind(canopy_crown_area_all_dates,tmp)
}



mean_cca_df <- canopy_crown_area_all_dates %>%
  group_by(age, pft) %>%
  summarise(mean_cca = mean(cca_adj, na.rm = T), n = length(cca), sd = sd(cca_adj))


mean_cca_df %>%
  ggplot(aes(x = age, y = mean_cca, fill = pft)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") +
  #geom_point(size = 4, position = "dodge") +
  scale_fill_manual(values = c("yellow","darkorange1","chartreuse4","brown")) +
  xlab("Patch age [yrs]") +
  ylab("Canopy crown area [m2 m-2]") +
  labs(title = "Canopy crown area by pft and patch age \n (1950-1975)") +
  theme_minimal()


canopy_crown_area_all_dates %>%
  ggplot(aes(x = age, y = cca, fill = pft)) +
  geom_boxplot() +
  #geom_point(size = 4, position = "dodge") +
  scale_fill_manual(values = c("yellow","darkorange1","chartreuse4","brown")) +
  xlab("Patch age [yrs]") +
  ylab("Canopy crown area [m2 m-2]") +
  theme_minimal()

