#This script is used to visualize how fixed window running means for solar radiation and moisture deficit days behave relative to exponential moving averages of the variables. This is done to calibrate the parameters of the exponential moving average to find parameter values that best match the fixed window moving averages that were used by the Tree Recruitment Scheme.

#The data used to visualize the behavior of fixed window moving averages and exponential moving averages are from FATES (for par at the forest floor) and ED2 (for moisture deficit days)

library(ncdf4)
library(tidyverse)
library(lubridate)

source('generalFunctions.R')
source('utils/supporting_funcs_esm_tools.R')

#Fates data
data_path <- '~/cloud/gdrive/FATES/FATES_data'
case_name <- 'regenTest.Cfd1459da-F8a065b20.2021-07-28'
files <- getFilePaths(base_data_path = data_path, case_name = case_name)[1:12]

#ED2 data
ED2_data <- read_csv('data/ED_input_example.csv')


#get running means of smp and sr in ED2
ED2_rMeans <- getRunningMeans(ED2_data, "sr", 7, 14)

ggplot(data = ED2_rMeans, 
       mapping = aes(index,sr)) +
  geom_line() +
  geom_line(data = ED2_rMeans, mapping = aes(index,rMean), color = "red") +
  geom_line(data = ED2_rMeans, mapping = aes(index,ema), color = "blue") +
  scale_x_continuous(limits = c(0,100)) +
  theme_minimal()


#testing how moisture deficit days behaves as an ema
psi_crit <- -175912.9
window_vector <- 126
window_ema <- 126

ED2_mdds_ema <- mdd_ema(ED2_data$smp,window_ema,psi_crit,1)

ED2_data %>%
  mutate(index = 1:nrow(ED2_data)) %>%
  mutate(mdd_vect = def_func(soil_moist = smp,psi_crit.x = psi_crit,
                             window = 126)) %>%
  mutate(mdd_ema = ED2_mdds_ema) %>%
  ggplot(mapping = aes(index, smp)) +
  geom_line(color = "blue") +
  geom_line(mapping = aes(index,mdd_vect), color = "green") +
  geom_line(mapping = aes(index,mdd_ema)) +
  geom_hline(yintercept = psi_crit) +
  scale_x_continuous(limits = c(400,700)) +
  theme_minimal()


#get running means of par at forest floor for FATES data
FATES_parForestFloor <- tibble(par = map_dbl(files, getPARatLowestLeafLayer))
FATES_rMeans <- getRunningMeans(FATES_parForestFloor, "par", 3, 7)
#FATES_rMeans <- rbind(FATES_rMeans, FATES_rMeans) %>% rbind(FATES_rMeans)

ggplot(data = FATES_rMeans, 
       mapping = aes(index,par)) +
  geom_line() +
  geom_line(data = FATES_rMeans, mapping = aes(index,rMean), color = "red") +
  geom_line(data = FATES_rMeans, mapping = aes(index,ema), color = "blue") +
  theme_minimal()


#A factor of 3/7 seems to get the ema to match the regular moving average, but this needs more testing

