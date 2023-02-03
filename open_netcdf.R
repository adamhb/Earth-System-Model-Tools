
#install libraries
library(ncdf4)
library(tidyverse)
library(lubridate)

#set paths
data_path <- '~/cloud/gdrive/'
case_name <- "postdoc"
case_alias <- "-"
vars <- c('FATES_NPP_PF')

files <- getFilePaths(case_name = case_name)
myfile <- files[1]
mydata <- nc_open(myfile)

ncvar_get(mydata, "TSAI")




for(i in 1:3){
  myfile <- files[i]
  mydata <- nc_open(myfile)
  print(ncvar_get(mydata, 'FATES_NPP_PF'))
}



ncvar_get(mydata, 'RECRUITMENT')
ncvar_get(mydata, 'PATCH_AREA_BY_AGE')
ncvar_get(mydata, 'fates_lfmap_levcnlf')
ncvar_get(mydata, 'fates_levage')






