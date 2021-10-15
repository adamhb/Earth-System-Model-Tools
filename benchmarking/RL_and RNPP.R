

library(tidyverse)
library(lubridate)
library(ggpubr)

source('utils/figure_formatting_esm_tools.R')

########################### Correlation between R/L and R/NPP, ForC database
#Anderson-Teixeira, K.J., Wang, M.M.H., McGarvey, J.C., Herrmann, V., Tepley, A.J., Bond-Lamberty, B. and LeBauer, D.S. (2018), 
#ForC: a global database of forest carbon stocks and fluxes. Ecology, 99: 1507-1507. https://doi.org/10.1002/ecy.2229

#git hub: https://github.com/forc-db/ForC


ForCdf <- read_csv("data/ForCformatted.csv")


ForCdf %>% 
  ggplot(aes(x = RNPP, 
             y = RL)) +
  geom_point() +
  geom_smooth(method = "lm", 
              color = "#007FFF") + 
  stat_regline_equation(aes(label = ..rr.label..))+
  adams_theme

# two sites were removed becase R/L was an outlier (> 2)
# these sites were from plot Id 539, site =  Cocoflux,  plot = coconut plantation

