
#install libraries
library(ncdf4)
library(tidyverse)
library(lubridate)

#set paths
data_path <- '~/cloud/gdrive/FATES/FATES_data'
case_name <- "bci-1pft-fates-main-from-inventory.Ccd4bc7aba-Fa0568083.2021-09-29_RAis0.079"
case_alias <- "-"
vars <- c('RECRUITMENT','DDBH_CANOPY_SCPF','DDBH_UNDERSTORY_SCPF','DDBH_SCPF','NPP_SCPF','PARPROF_DIR_CNLF','PARPROF_DIF_CNLF','PATCH_AREA_BY_AGE','NPATCH_BY_AGE','BA_SCPF', 'NPLANT_SCPF')

files <- getFilePaths(case_name = case_name)
myfile <- files[100]
mydata <- nc_open(myfile)


ncvar_get(mydata, 'fates_levscls')
ncvar_get(mydata, 'RECRUITMENT')
ncvar_get(mydata, 'PATCH_AREA_BY_AGE')
ncvar_get(mydata, 'fates_lfmap_levcnlf')
ncvar_get(mydata, 'fates_levage')






