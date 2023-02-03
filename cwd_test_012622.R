#install libraries
library(ncdf4)
library(tidyverse)
library(lubridate)

#set paths
data_path_control <- '~/cloud/gdrive/postdoc/simulation_output/cwd_test_control/'
data_path_test <- '~/cloud/gdrive/postdoc/simulation_output/cwd_test/'



#case_name <- "oak_resprout_vs_not_012422"
#case_alias <- "oak_resprout_vs_not_012422"

vars <- c('FATES_NPP_PF','FATES_NPLANT_PF','FATES_MORTALITY_PF')
files <- getFilePaths(case_name = case_name)[200:250]

mort <- c()
for(f in files){
  file <- nc_open(f)
  mort <- c(mort,ncvar_get(file, vars[3])[2])
  nc_close(file)
}

n <- c()
for(f in files){
  file <- nc_open(f)
  n <- c(n,ncvar_get(file, vars[2])[2])
  nc_close(file)
}

tibble(mort = mort, n = n) %>% mutate(per_cap_mort = mort / n) %>%
  filter(per_cap_mort > 1)
