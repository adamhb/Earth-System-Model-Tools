library(lubridate)
library(tidyverse)


data_path <- "~/Desktop/NGEE-Tropics/BCI_data/"
BCI_mo_RL <- read_csv(paste(data_path,"BCI_mo_RL.csv", sep = "")

BCI_mo_RL %>% 
  group_by(Year) %>% 
  summarise(Rm2yr = sum(Rm2), 
            Lm2yr = sum(Lm2),
            Sm2yr = sum(Sm2)) %>% 
  mutate(RL = Rm2yr/Lm2yr) %>% 
  ggplot() +
  geom_boxplot(aes(x = "", 
                   y = RL))+
  xlab("BCI R/L (all species)")
