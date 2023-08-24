#This script is used to plot/analyze the recruitment data
library(ggplot2)
library(tidyverse)
source('utils/figure_formatting_esm_tools.R')
source('utils/supporting_funcs_esm_tools.R')
source('generalFunctions.R')

sdx <- 15

bci_all_obs <- read_csv('data/all_BCI_obs.csv')

bci_R <-  237.17687500000002 # g m-2 yr-1; source?
bci_Rvect <- c(242,202,280,289,172) #g m-2 yr-1
bci_ANPP <- 1800 # source:"Condit, R., D.M. Windsor, and S.P. Hubbell. 2013. NPP Tropical Forest: Barro Colorado, Panama, 1969-2000, R1. Data set. Available on-line [http://daac.ornl.gov] from Oak Ridge National Laboratory Distributed Active Archive Center, Oak Ridge, Tennessee, USA doi:10.3334/ORNLDAAC/157"
bci_RoANPP <- bci_Rvect / bci_ANPP

luq_R <- 51 #source? g m-2 yr-1
luq_Rvect <- rnorm(10,luq_R,sdx)
luq_ANPP <- 1050 #source? g m-2 yr-1
luq_RoANPP <- luq_Rvect / luq_ANPP


df1 <- tibble(RoANPP = bci_RoANPP) %>% add_column(site = "bci")
df2 <- tibble(RoANPP = luq_RoANPP) %>% add_column(site = "luq")

df <- rbind(df1,df2)

RoANPP_fig <- df %>%
  ggplot(aes(site,RoANPP)) +
  geom_boxplot() +
  ylab(expression(paste('R / ANPP', ' [g m'^'-2','yr'^'-1',"]"))) +
  xlab(bquote('Site')) +
  scale_x_discrete(labels = c('bci\nn = 5','luq\nn = 17')) +
  scale_y_continuous(limits = c(0,0.2)) +
  adams_theme
    
    







