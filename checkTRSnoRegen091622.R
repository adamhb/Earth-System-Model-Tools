library(tidyverse)
library(ncdf4)
library(lubridate)

source('generalFunctions.R')
source('utils/figure_formatting_esm_tools.R')

data_path <- '~/cloud/gdrive/FATES/FATES_data'
case_name <- 'bci-12pfts-fates-TRS-API24-July-2022.Cba72f7272-F97ad792b.2022-07-21_checkForRachel091622'
case_alias <- "fates-trs-noRegen"
files <- getFilePaths(data_path,case_name)
FATES_trs_noRegen <- getManySiteVarsOverTime(files,varNames = c('FATES_VEGC_PF','FATES_RECRUITMENT_PF'))


FATES_trs_noRegen %>%
  ggplot(aes(date,FATES_RECRUITMENT_PF * 1e4)) +
  geom_point() +
  ylab("N recruits [ha-1 yr-1]") +
  adams_theme
