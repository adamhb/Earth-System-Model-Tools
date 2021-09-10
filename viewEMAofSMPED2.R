library(ncdf4)
library(tidyverse)
library(lubridate)

source('generalFunctions.R')

data_path <- '~/cloud/gdrive/rec_submodel/ED2/'
case_name <- 'regenTest.Cfd1459da-F8a065b20.2021-07-28'

files <- getFilePaths(base_data_path = data_path, case_name = case_name)[1:12]