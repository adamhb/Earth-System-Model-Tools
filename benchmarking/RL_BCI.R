library(tidyverse)
source('utils/figure_formatting_esm_tools.R')

data_path <- "data/"
BCI_mo_RL <- read_csv(paste(data_path,"BCI_mo_RL.csv", sep = ""))

BCI_mo_RL %>% 
  group_by(Year) %>% 
  summarise(Rm2yr = sum(Rm2), 
            Lm2yr = sum(Lm2),
            Sm2yr = sum(Sm2)) %>% 
  mutate(RL = Rm2yr/Lm2yr) %>% 
  ggplot() +
  geom_boxplot(aes(x = "", 
                   y = RL))+
  xlab("BCI obs. R/L (all species)") +
  adams_theme


#calculate leaf and repro litter fluxes in g C m-2 day-1
BCI_mo_RL %>% 
  group_by(Year) %>% 
  summarise(Rm2yr = sum(Rm2), 
            Lm2yr = sum(Lm2),
            Sm2yr = sum(Sm2)) %>% 
  mutate_at(.vars = c("Rm2yr","Lm2yr","Sm2yr"),.funs = function(x){x / 365}) 


