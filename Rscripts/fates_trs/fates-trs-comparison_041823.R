library(tidyverse)
source("utils/supporting_funcs_esm_tools.R")
source("utils/figure_formatting_esm_tools.R")

path <- "~/cloud/gdrive/FATES/Earth-System-Model-Tools/bci_fates-trs_test"
files <- list.files(path,pattern = ".csv")
simulations <- c("fates-main","trs-off","trs-on","trs-on-060923","fates-main","trs-off","trs-on","trs-on-060923")
m2_per_ha <- 10000
s_per_yr <- 3600 * 24 * 365

df <- tibble()
j <- 0
for(f in files){
  j <- j + 1
  tmp <- read_csv(paste0(path,"/",f)) %>%
    mutate(simulation = simulations[j]) %>%
    gather(2, key = "var",value = "value")
  df <- rbind(df,tmp)
}

rec_fig <- df %>%
  filter(var == "FATES_RECRUITMENT_PF") %>%
  ggplot(aes(time,value * m2_per_ha,color = simulation, linetype = simulation)) +
  geom_line(size = 2) +
  scale_linetype_manual(values = c("solid","dashed","solid","dashed")) +
  scale_y_log10() +
  labs(title = "Recruitment") +
  geom_hline(yintercept = 70,color = "black") +
  ylab(label = "Recruitment rate [N ha-1 yr-1]") +
  xlab(label = "Simulation Year") +
  theme_minimal()

makePNG(fig = rec_fig,path_to_output.x = path,file_name = "rec")

npp_fig <- df %>%
  filter(var == "FATES_NPP_PF") %>%
  ggplot(aes(time,value * s_per_yr, color = simulation, linetype = simulation)) +
  geom_line(size = 2) +
  scale_linetype_manual(values = c("solid","dashed","solid","dashed")) +
  scale_y_log10() +
  labs(title = "NPP") +
  ylab(label = "NPP [Kg m-2 yr-1]") +
  xlab(label = "Simulation Year") +
  theme_minimal()

makePNG(fig = npp_fig,path_to_output.x = path,file_name = "npp")
