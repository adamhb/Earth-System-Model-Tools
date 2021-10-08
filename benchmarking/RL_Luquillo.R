
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

# Adding treatment information

bA<- c("Control", "Trim&clear", "Trim+debris", "NoTrim+debris")
bB<- c("Control", "Trim+debris", "NoTrim+debris", "Trim&clear")
bC <- c("NoTrim+debris", "Trim+debris", "Trim&clear", "Control")

CTE_Block_desc <- tibble(Block = rep(c("A", "B", "C"), ea = 4),
                         Plot = rep(1:4, 3),
                         Treatment = c(bA, bB, bC))

CTE_control <- CTE %>% 
  left_join(CTE_Block_desc, by = c("Block", "Plot")) %>% #add treatment information
  filter(Treatment == "Control") %>% #filter for control plots
  rename(L = `Leaves (in g)`,
         R = `Fruits;seeds;flower (g)`) %>% 
  mutate(Lgm2 = L / (nTraps_plot*trap_size_m2),
         Rgm2 = R / (nTraps_plot*trap_size_m2)) %>%# observations from each plot are sum of material from 10 traps
  mutate(Yr_mo = as.Date(paste(month(Date), "01", year(Date), sep = "/"), format = "%m/%d/%Y")) %>% # create Yr-mo for easy grouping
  group_by(Yr_mo, Block, Plot) %>% # monthly g/m2month for each plot within each block
  summarise(Lgm2moplot = sum(Lgm2, na.rm = T),
            Rgm2moplot = sum(Rgm2, na.rm = T)) %>%
  group_by(Yr_mo) %>% # get monthly avg for all blocks 
  summarise(Lgm2mo = mean(Lgm2moplot, na.rm = T),
            Rgm2mo = mean(Rgm2moplot, na.rm = T)) %>% 
  group_by(Year = year(Yr_mo)) %>%  # sum over year
  summarise(Lgm2yr = sum(Lgm2mo, na.rm = T),
            Rgm2yr = sum(Rgm2mo, na.rm = T)) %>% 
  mutate(RL = Rgm2yr/Lgm2yr)
  


#write observations to a csv
CTE_control %>% 
  add_column(case = "Luquillo obs.", var = "RoL") %>% 
  rename(simYr = Year, value = RL) %>% 
  select(case,simYr,var,value) %>%
  ungroup() %>%
  write_csv(file = "data/RoverL_Luquillo_obs.csv")

  
  CTE_control %>% 
    filter(Year > 2002)%>% # full years 
  ggplot() +
  geom_boxplot(aes(x = "", 
                   y = RL)) +
  xlab("Luquillo") +
  ggtitle("Control plots of CTE, El Verde, Luquillo (2003-2019)") +
  adams_theme
  





