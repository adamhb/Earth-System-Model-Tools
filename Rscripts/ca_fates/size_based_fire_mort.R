#install libraries
library(ncdf4)
library(tidyverse)
library(lubridate)

data_path <- '/home/adam/cloud/gdrive/postdoc/simulation_output/'

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



getDateFromNC <- function(file){
  return(ymd(getSiteVar(file,"mcdate")))
}

s <- c(0.,   5.,  10.,  15.,  20.,  30.,  40.,  50.,  60.,  70.,  80.,  90.,
  100)

scls <- rep(c(0.,   5.,  10.,  15.,  20.,  30.,  40.,  50.,  60.,  70.,  80.,  90.,
               100),4)

pfts <- c(rep("pine",length(scls)/4),rep("cedar",length(scls)/4),rep("fir",length(scls)/4),rep("shrub",length(scls)/4))


fire_mort_by_size <- tibble()
for(i in 1:length(files)){
  fire_mort_by_size <- rbind(fire_mort_by_size, tibble(fire_mort = getMultiDimSiteVar(file = files[i],varName = "FATES_MORTALITY_FIRE_SZPF"), scls = scls, pft = pfts))
}
fire_mort_by_size <- fire_mort_by_size %>% filter(fire_mort > 0)


mean_fire_mort <- fire_mort_by_size %>%
  group_by(scls, pft) %>%
  summarise(mean_f_mort = mean(fire_mort, na.rm = T), n = length(fire_mort), sd = sd(fire_mort)) %>%
  mutate_at(.vars = c("mean_f_mort","sd"), .funs = function(x){x * 1e4}) %>%
  mutate(pft = factor(pft, levels = c("pine","cedar","fir","shrub")))


mean_fire_mort %>%
  ggplot(aes(x = scls, y = mean_f_mort, fill = pft)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") +
  scale_fill_manual(values = c("yellow","darkorange1","chartreuse4","brown")) +
  xlab("size class lower bound [cm]") +
  ylab("Mortality rate [N ha-1 yr-1]") +
  labs(title = "Mean size-based mort rate from fire (1980-1989)") +
  scale_x_continuous(breaks = s, labels = s) +
  theme_minimal()


canopy_crown_area_all_dates %>%
  ggplot(aes(x = age, y = cca, fill = pft)) +
  geom_boxplot() +
  #geom_point(size = 4, position = "dodge") +
  scale_fill_manual(values = c("yellow","darkorange1","chartreuse4","brown")) +
  xlab("Patch age [yrs]") +
  ylab("Canopy crown area [m2 m-2]") +
  theme_minimal()

