
library(tidyverse)
library(lubridate)

source('utils/figure_formatting_esm_tools.R')
source('benchmarking/RL_functions.R')


##############################Ecosystem level R/L , all species###########################
data_path <- "~/cloud/gdrive/review_paper/quant_summary_data/"


BCI_mo_RL <- read_csv("data/BCI_mo_RL.csv")

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

############################## Canopy sp  R/L ##############################################

# load species list
canopy_sp <- read_csv("data/canopy_sp_bci.csv") %>%
  select(sp)
#load palm species
palms <- read_csv("data/palms_Arecaceae.csv")

#remove palm species from canopy species
canopy_sp <- canopy_sp %>% filter(!sp %in% palms$sp)

#load and clean the trap data
trap_data <- read_csv("data/Hanbury_Brown_50ha_rawdata_20190404.csv") %>% 
  mutate(part = as_factor(part)) %>% 
  mutate(sp = tolower(sp6)) %>% 
  mutate(Date = as.Date(fecha, format = "%Y%m%d")) %>% 
  mutate(Year = lubridate::year(Date)) %>% 
  filter(sp %in% canopy_sp$sp) %>% # filter for only canopy sp
  filter(Year>2013 & Year < 2019) ## full years of data

# save n traps, trap area to later scale to g/m2
ntraps = n_distinct(trap_data$trap)
trap_size = 0.5 #units = m2

# nest data by species, Year
Year_sp_dat <- trap_data %>% 
  mutate(sp2 = sp) %>% # create another sp var to use so that sp col stays in nested data
  group_by(Year, sp2) %>% 
  nest()

# map function to find repro, leaf, and corrected repro mass
Y2 <- Year_sp_dat %>% 
  mutate(sp_tib = map(data, get_RL_df, fruit_traits = fruit_traits, use_best_est = F)) 

# unnest
Y3 <- Y2 %>% 
  select(Year, sp_tib) %>% 
  unnest(cols = c(sp_tib)) %>%
  ungroup() 

# summarise over species, convert g/year --> g/m2year

Y4 <- Y3 %>%
  rowwise() %>%
  mutate(Repro_c = ifelse(is.na(all_R_corr), all_R, all_R_corr)) %>% #if no R corrrection, use "all_R"
  group_by(Year) %>% 
  summarise(Leaf = sum(all_L, na.rm =T)/(trap_size * n_traps),  ## 64 traps, each 0.5 m2
            Repro_corr = sum(Repro_c, na.rm = T)/(trap_size * n_traps),
            Repro = sum(all_R, na.rm = T)/(trap_size * n_traps)) %>% 
  mutate(L_corr =round(Leaf*1.097, 1)) %>% 
  ##leaf herbivory correction source: pers. communication with Joe Wright, Oct 13, 2021
  mutate(RL_corr = Repro_corr / L_corr,
         RL_no_corr = Repro/Leaf)


Y5 <- Y4 %>% 
  pivot_longer(cols = c(RL_corr, RL_no_corr), names_to = "Type", values_to = "Value")

Y5 %>%
  ggplot() + 
  geom_boxplot(aes(x = Type, 
                   y = Value, 
                   group = Type)) 


# BCI average R g/m2yr
BCImeanRgm2yr <- mean(Y5$Repro_corr[Y5$Type == "RL_corr"])


# ANPP estimate info from https://daac.ornl.gov/NPP/guides/NPP_BRR.html
RANPP_BCI <- tibble(site = "BCI", 
                    year = NA, 
                    var = c("R", "ANPP", "R/ANPP"), 
                    value = c(BCImeanRgm2yr, 1800, BCImeanRgm2yr/1800),
                    units = c(rep("g m-2 yr-1", 2), "-")) # units are g/m2yr

#write_csv(RANPP_BCI, path = "data/RANPP_BCI_obs.csv")

# join all BCI observations, write to .csv
BCIdf <- Y5 %>% 
  filter(Type == "RL_corr") %>% 
  add_column(site = "BCI") %>% 
  rename(year = Year, 
         L = L_corr, 
         R = Repro_corr, 
         "R/L" = Value) %>% 
  select(site, year, R, L, `R/L`) %>% 
  pivot_longer(cols = c(R, L, `R/L`), names_to = "var") %>% 
  mutate(units = ifelse((var == "R" | var == "L"), "g m-2 y-1", "-")) %>%
  arrange(var) %>% 
  bind_rows(RANPP_BCI) %>% 
  write_csv("data/all_BCI_obs.csv")



# write observations for R/L to a csv
Y5 %>%
  filter(Type == "RL_corr") %>%
  add_column(case = "BCI obs.", var = "R/L") %>%
  rename(simYr = Year, value = Value) %>% 
  select(case,simYr,var,value) %>%
  ungroup() %>%
  write_csv(path = "data/RoverL_BCI_obs.csv")









            