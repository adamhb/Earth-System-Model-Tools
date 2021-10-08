library(tidyverse)
library(lubridate)


##############  Ecosystem Level R/L in the Luquillo Experimental Forest

# Data from: 
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-luq.162.862754
# Citation: 
# RamIrez, A. 2020. Canopy Trimming Experiment (CTE) Litterfall ver 862754. 
# Environmental Data Initiative. https://doi.org/10.6073/pasta/725e8064673b05a53ef5bf5a45abb4db (Accessed 2021-10-08).



###################################cleaning raw data 
CTE <- read_csv("data/CTE-Litterfall-data_2002-2019.csv") %>% 
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"), 
         Year = year(Date), 
         Month = month(Date))

#sanity check
# CTE %>% arrange(Year) %>%
#   ggplot() +
#   geom_col(aes(x = Date,
#                y = Plot,
#                fill = Block))

## foramt for years 2017, 2019 needs to be fixed
CTE_YYYY <- CTE %>% 
  filter(Year >1900)


CTE_YY <- CTE %>% 
  filter(Year <1900) %>% 
  mutate(Date = as.Date(paste((year(Date) +2000),"/",month(Date),"/",day(Date), sep = ""))) 


CTE_datesfixed <- CTE_YYYY %>% 
  rbind(CTE_YY) %>% 
  mutate(Year = year(Date), 
         Month = month(Date))

#sanity check:
CTE_datesfixed %>% 
  arrange(Year) %>% 
  ggplot() +
  geom_col(aes(x = Date, 
               y = Plot, 
               fill = Block))


## double entry of data,  April 23 2015 data

CTE_datesfixed %>% arrange(Year) %>% 
  filter(Year >2014 & Year <2016) %>% 
  filter(Month >3& Month <7) %>% 
  ggplot() +
  geom_col(aes(x = Date, 
               y = Plot, 
               fill = Block))

CTE_datesfixed %>% 
  arrange(Year) %>% 
  filter(Date == "2015-04-23") %>% 
 # distinct() %>% 
  ggplot() +
  geom_col(aes(x = Plot, 
               y = Block))

April_cln <- CTE_datesfixed %>% 
  filter(Date == "2015-04-23") %>% 
  distinct()

CTE_clean <- CTE_datesfixed %>% 
  filter(Date != "2015-04-23") %>% 
  left_join(April_cln)

##sanity check: looks good 
CTE_clean %>% 
  ggplot() +
  geom_col (aes(x = Date, y = Plot, fill = Block))

#write_csv(CTE_clean, file = "data/CTE_clean.csv")


