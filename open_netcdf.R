
#install libraries
library(ncdf4)
library(tidyverse)
library(lubridate)

#set paths
data_path <- '~/cloud/gdrive/FATES/FATES_data'
case_name <- "bci-1pft-fates-main.Ccd4bc7aba-F70a3e133.2021-09-28"
case_alias <- "-"
vars <- c('RECRUITMENT','DDBH_CANOPY_SCPF','DDBH_UNDERSTORY_SCPF','DDBH_SCPF','NPP_SCPF','PARPROF_DIR_CNLF','PARPROF_DIF_CNLF','PATCH_AREA_BY_AGE','NPATCH_BY_AGE','BA_SCPF', 'NPLANT_SCPF')


path <- paste0(data_path,'/',case_name,'/')
fileNames <- list.files(path,pattern = '.nc')
files <- paste0(path,fileNames)


myfile <- files[3]
mydata <- nc_open(myfile)

ncvar_get(mydata, 'NPP')
ncvar_get(mydata, 'PARPROF_DIR_CNLF')
ncvar_get(mydata, 'PATCH_AREA_BY_AGE')
ncvar_get(mydata, 'fates_lfmap_levcnlf')
ncvar_get(mydata, 'fates_levage')






