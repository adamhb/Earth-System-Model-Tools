
library(tidyverse)
library(lubridate)
source('utils/figure_formatting_esm_tools.R')

##############  Ecosystem Level R/L in the Luquillo Experimental Forest

# Data from: 
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-luq.162.862754
# Citation: 
# RamIrez, A. 2020. Canopy Trimming Experiment (CTE) Litterfall ver 862754. 
# Environmental Data Initiative. https://doi.org/10.6073/pasta/725e8064673b05a53ef5bf5a45abb4db (Accessed 2021-10-08).

nTraps_plot <-10 # 10 litterfall traps per plot
trap_size_m2 <- 1.75^2# 1.75 m X 1.75 m baskets


CTE <- read_csv("data/CTE_clean.csv")

# Adding treatment information, source: https://luq.lter.network/data/luqmetadata162

bA<- c("Control", "Trim&clear", "Trim+debris", "NoTrim+debris")
bB<- c("Control", "Trim+debris", "NoTrim+debris", "Trim&clear")
bC <- c("NoTrim+debris", "Trim+debris", "Trim&clear", "Control")

CTE_Block_desc <- tibble(Block = rep(c("A", "B", "C"), ea = 4),
                         Plot = rep(1:4, 3),
                         Treatment = c(bA, bB, bC))


CTE_Control <- CTE %>% 
  left_join(CTE_Block_desc, by = c("Block", "Plot")) %>% #add treatment information
  filter(Treatment == "Control") %>% #filter for control plots
  rename(L = `Leaves (in g)`,
         R = `Fruits;seeds;flower (g)`) %>% 
  mutate(Lgm2 = L / (nTraps_plot*trap_size_m2), 
         Rgm2 = R / (nTraps_plot*trap_size_m2)) %>% 
  select(Date, Year, Block, Plot, Lgm2, Rgm2) %>% 
  group_by(Year) %>% 
  summarise(Rgm2yr = sum(Rgm2)/3, # three unique block-plots
            Lgm2yr = sum(Lgm2)/3) %>% 
  filter(Year>2002) #full years


  
#observatsion of R g/m2yr to add to csv

Rdf <- CTE_control %>% 
  add_column(case = "Luquillo obs.", var = "Rgm2yr") %>% 
  rename(simYr = Year, value = Rgm2yr) %>% 
  select(case, simYr, var, value)


#write observations to a csv
CTE_control %>% 
  add_column(case = "Luquillo obs.", var = "RoL") %>% 
  rename(simYr = Year, value = RL) %>% 
  select(case,simYr,var,value) %>%
  ungroup() %>%
  bind_rows(Rdf) %>% 
  write_csv(file = "data/RoverL_Luquillo_obs.csv")

#plot 
  CTE_control %>% 
  ggplot() +
  geom_boxplot(aes(x = "", 
                   y = RL)) +
  xlab("Luquillo") +
  ggtitle("Control plots of CTE, El Verde, Luquillo (2003-2019)") +
  adams_theme
  
  
Rdf$value *10000
