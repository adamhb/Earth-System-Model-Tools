
#install libraries
library(ncdf4)
library(tidyverse)
library(lubridate)

#set paths
data_path <- '~/cloud/gdrive/FATES/FATES_data'
case_name <- "repro_alloc_sensitivity.C63f3ea8a4-Fc1f6dddf.2021-08-03"
case_alias <- "-"
vars <- c('RECRUITMENT','DDBH_CANOPY_SCPF','DDBH_UNDERSTORY_SCPF','DDBH_SCPF','NPP_SCPF','PARPROF_DIR_CNLF','PARPROF_DIF_CNLF','PATCH_AREA_BY_AGE','NPATCH_BY_AGE','BA_SCPF', 'NPLANT_SCPF')


path <- paste0(data_path,'/',case_name,'/')
fileNames <- list.files(path,pattern = '.nc')
files <- paste0(path,fileNames)


myfile <- files[10]
mydata <- nc_open(myfile)
ncvar_get(mydata, 'DDBH_CANOPY_SCPF')





