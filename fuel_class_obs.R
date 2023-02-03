library(tidyverse)
library(ggplot2)

fuel_obs <- read_csv('/home/adam/cloud/gdrive/postdoc/benchmarking/fuel_class_obs.csv')



fuel_obs %>%
  ggplot(aes(x = Fuel_class, y = Mean_kg_C_m2, fill = Site)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") +
  #geom_point(size = 4, position = "dodge") +
  scale_fill_manual(values = c("yellow","darkorange1","chartreuse4")) +
  xlab("Fuel Class") +
  ylab("Fuel Load [kg C m-2]") +
  #labs(title = "Canopy crown area by pft and patch age \n (1950-1975)") +
  theme_minimal()