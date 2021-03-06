
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
trap_size_m2 <- 1.75^2 # 1.75 m X 1.75 m baskets


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



#Luquillo average R g/m2yr
LuquillomeanRgm2yr <- mean(CTE_control$Rgm2yr) 
# Not currently using this value


#write observations for R/ANPP to a csv
# source: https://daac.ornl.gov/NPP/guides/NPP_LQL.html
# file = "lql1_npp_r1.txt"
RANPP_Luquillo <- tibble(site = "Luquillo", 
                         year = NA, 
                         var = c("R", "ANPP", "R/ANPP"), 
                         value = c(51, (525*2), 51/(525*2)),
                         units = c(rep("g m-2 yr-1", 2), NA)) # units are g/m2yr

#write_csv(RANPP_Luquillo, path = "data/RANPP_Luquillo_obs.csv")


  
# write all Luquillo observatsion to .csv

Luquillodf <- CTE_control %>% 
  add_column(site = "Luquillo") %>% 
  rename(year = Year, 
         R = Rgm2yr, 
         L = Lgm2yr,
         `R/L` = RL) %>% 
  pivot_longer(cols = c(R, L, `R/L`), names_to = "var") %>% 
  mutate(units = ifelse((var == "R" | var == "L"), "g m-2 y-1", "-")) %>%
  arrange(var) %>% 
  select(site, year, var, value, units) %>% 
  bind_rows(RANPP_Luquillo) %>% 
  write_csv("data/all_Luquillo_obs.csv")
  



#write R/ L observations to a csv
CTE_control %>% 
  add_column(case = "Luquillo obs.", var = "R/L") %>% 
  rename(simYr = Year, value = RL) %>% 
  select(case,simYr,var,value) %>%
  ungroup() %>%
  #bind_rows(Rdf) %>% 
  write_csv(file = "data/RoverL_Luquillo_obs.csv")




#plot 
  CTE_control %>% 
  ggplot() +
  geom_boxplot(aes(x = "", 
                   y = RL)) +
  xlab("Luquillo") +
  ggtitle("Control plots of CTE, El Verde, Luquillo (2003-2019)") +
  adams_theme
  
  
  
  ######### raw values look anomalously low, by a factor of 10
  ######### when metadata says litter was pooled from the 10 baskets per site, do they mean averaged? 
  ######## *not* including the additional dividion by 10 makes these values make much more sense
  ###### I emailed Jess Zimmerman about this, hopefully he gets back to me quickly
  
CTE %>%
    left_join(CTE_Block_desc, by = c("Block", "Plot"))  %>%
  filter(Treatment == "Control") %>% #filter for control plots
  rename(L = `Leaves (in g)`,
         R = `Fruits;seeds;flower (g)`)  %>%
  select(Block, Plot, Treatment, Date, Year, R, L ) %>%
    group_by(Block, Plot) %>% 
    mutate(ID = cur_group_id()) %>% 
    mutate(ID = factor(ID)) %>% 
    ungroup() %>% 
    group_by(ID, Year) %>% 
  summarise(Ryr = sum(R, na.rm = T)/(10*1.75^2), #######################to divide, or not divide ...
            Lyr = sum(L, na.rm = T)/(10*1.75^2)) %>%
    pivot_longer(cols = c("Ryr", "Lyr"), names_to = "Component", values_to = "Value") %>% 
    ggplot(aes(x = Year, 
               y = Value,
               linetype = Component,
               color = ID)) +
    geom_line()





