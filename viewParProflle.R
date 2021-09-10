library(ncdf4)
library(tidyverse)
library(lubridate)

source('generalFunctions.R')
data_path <- '~/cloud/gdrive/FATES/FATES_data/'
case_name <- 'regenTest.Cfd1459da-F8a065b20.2021-07-28'

files <- getFilePaths(base_data_path = data_path, case_name = case_name)[1:12]

ggplot(data = tibble(month = as.integer(1:length(files)), 
                     par = map_dbl(.x = files,.f = getPARatLowestLeafLayer)),
       mapping = aes(month,par)) +
         geom_line() +
         ylab(label = "par at seedling layer [MJ m2 per day]") +
         xlab(label = "month since bare ground") +
         geom_hline(yintercept = 0.334) + #mean understory shortwave radiation in ED2 simulations
         theme_minimal() 
