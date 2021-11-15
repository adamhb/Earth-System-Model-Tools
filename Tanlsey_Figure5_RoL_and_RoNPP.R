

library(tidyverse)
library(ggpubr)
library(plotbiomes)
source('utils/figure_formatting_esm_tools.R')

############################################################
###########  Figures 5a and 5b in Tansley review ###########
############################################################


#load data
TableS2 <- read.csv("data/TableS2.csv")

# For Figure 5a: filter for only sites where NPP is reported
Fig5a_df <- TabelS2 %>% 
  filter(!is.na(NPP))

# For Figure 5b: filter to remove sites which fall beyond axis limits, non-forest biomes, duplicates
Fig5b_df <- TabelS2 %>% 
  filter(MAP < 5000) %>% # remove sites where precip is beyond y axis limit
  filter(!str_detect(`Whittaker biome`, "desert")) # remove sites in desert biomes
  distinct(MAT, MAP, `R/L`, .keep_all = T) # remove any duplicates



#####################################################################
######## Prepare Whittaker biome names and colors for legend ########
#####################################################################

# Assign biome names and colors
biomes_tbl <- data.frame(biome_id = 1:10,
                         biome = c("Tropical seasonal forest/savanna",
                                   "Subtropical desert",
                                   "Temperate rain forest",
                                   "Tropical rain forest",
                                   "Woodland/shrubland",
                                   "Tundra",
                                   "Boreal forest",
                                   "Temperate grassland/desert",
                                   "Temperate seasonal forest", 
                                   "NA"))

colors = c("#A09700", "#DCBB50", "#75A95E", "#317A22", "#D16E3F", "#C1E1DD", "#A5C790", "#FCD57A", "#97B669", "#989898")

# customizing colors, only diff is inclusion of an "NA" color
Ricklefs_colors <- colors
names(Ricklefs_colors) <- biomes_tbl$biome

# Sets a desired order for the names. This will affect the order in legend.
biomes_order <- c("Tundra",
                  "Boreal forest",
                  "Temperate seasonal forest",
                  "Temperate rain forest",
                  "Tropical rain forest",
                  "Tropical seasonal forest/savanna",
                  "Subtropical desert",
                  "Temperate grassland/desert",
                  "Woodland/shrubland",
                  "NA")
Ricklefs_colors <- Ricklefs_colors[order(factor(names(Ricklefs_colors), levels = biomes_order))]



#####################################################################
######## Figure 5 a: Plot R / NPP vs R / L ########
#####################################################################


# set axis size for all plots
axis_size = 15 


RNPP_plot<- Fig5a_df %>% 
  ggplot(aes(x = `R/NPP`, 
             y = `R/L`)) +
  geom_point(aes(x = `R/NPP`, 
                 y = `R/L`, 
                 color = `Whittaker biome`)) +
  scale_color_manual(name = "Whittaker biomes",
                     breaks = names(Ricklefs_colors),
                     labels = names(Ricklefs_colors),
                     values = Ricklefs_colors) +
  geom_smooth(method = "lm",
              color = "#666666") + 
  stat_regline_equation(aes(label = ..rr.label..),  label.x = .003, label.y = .78) +
  ylab("R/L") +
  xlab("R/NPP") +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        strip.text.x = element_text(size = axis_size),
        axis.title.x = element_text (size = axis_size), 
        axis.title.y = element_text (size = axis_size),
        axis.title.y.right = element_text (size = axis_size, color = "blue"),
        axis.text.x = element_text (size = axis_size, colour = "black"),
        axis.text.y = element_text (size = axis_size, colour = "black") )

#ggsave(RNPP_plot, device = "png", filename = "figures/Figure5_a.png", width = 7, height = 5 )

#####################################################################
######## Figure S1 : Plot R / ANPP vs R / L ########
#####################################################################

RANPP_plot<- Fig5a_df %>% 
  ggplot(aes(x = `R/ANPP`, 
             y = `R/L`)) +
  geom_point(aes(x = `R/ANPP`, 
                 y = `R/L`, 
                 color = `Whittaker biome`)) +
  scale_color_manual(name = "Whittaker biomes",
                     breaks = names(Ricklefs_colors),
                     labels = names(Ricklefs_colors),
                     values = Ricklefs_colors) +
  geom_smooth(method = "lm",
              color = "#666666") + 
  stat_regline_equation(aes(label = ..rr.label..), label.x = .003, label.y = .9) +
  ylab("R/L") +
  xlab("R/ANPP") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        strip.text.x = element_text(size = axis_size),
        axis.title.x = element_text (size = axis_size), 
        axis.title.y = element_text (size = axis_size),
        axis.title.y.right = element_text (size = axis_size, color = "blue"),
        axis.text.x = element_text (size = axis_size, colour = "black"),
        axis.text.y = element_text (size = axis_size, colour = "black") )


#ggsave(RANPP_plot, device = "png", filename = "figures/FigureS1.png", width = 7, height = 5 )


#####################################################################
######## Figure 5 b : Plot  R / L in Whittaker space ########
#####################################################################
base <- ggplot() +
  # add biome polygons
  geom_polygon(data = Whittaker_biomes,
               aes(x    = temp_c,
                   y    = precip_cm,
                   fill = biome),
               # adjust polygon borders
               colour = "gray98",
               size   = 1) + 
  xlab("Mean annual temperature (°C)") + 
  ylab("Mean annual precipitation (cm)") +
  theme_bw() +
  scale_fill_manual(name = "Whittaker biomes",
                    breaks = names(Ricklefs_colors),
                    labels = names(Ricklefs_colors),
                    values = Ricklefs_colors)



RL_plot<-  base + 
  geom_point(data =Fig5b_df,
             aes(x = MAT, 
                 y = (MAP / 10), #mm --> cm
                 size = `R/L`), 
             shape = 21, 
             color = "gray95",
             fill = "black", 
             stroke = .5, 
             alpha = .30) +
  theme_bw() + 
  scale_radius(name = "R/L", breaks = c(.05, .1, .5, 1), range = c(2, 15)) +
  
  guides(shape = guide_legend(override.aes = list(shape = c(1, 1),
                                                  size = c(2, 2),
                                                  stroke = c(1, 1),
                                                  color = c("black", "black"),
                                                  alpha= c(.4, .8))),
         size = guide_legend(override.aes = list(color = "gray80"))) +
  theme(aspect.ratio = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        strip.text.x = element_text(size = axis_size),
        axis.title.x = element_text (size = axis_size), 
        axis.title.y = element_text (size = axis_size),
        axis.title.y.right = element_text (size = axis_size, color = "blue"),
        axis.text.x = element_text (size = axis_size, colour = "black"),
        axis.text.y = element_text (size = axis_size, colour = "black") )


#ggsave(RL_plot, device = "png", filename = "figures/Figure5_b.png", width = 7, height = 5 )

#####################################################################
######## Create legend for overlapping sites, Figure 5 b  ########
#####################################################################

# This legend was added to Figure 5b in photoshop

# fake data to assess transparency at overlapping sites
testdf <- tibble (MAP = c(100, rep(150, 2), rep(200, 3), rep(250, 4)),
                  MAT = c(-10, rep(-10, 2), rep(-10, 3), rep(-10, 4)),
                  RL = c(rep(.3, 10)))

# test legend to match transparency levels at overlapping sites
test2<- tibble(MAP = c(100, 150, 200, 250), 
               MAT = c(rep(-15, 4)),
               RL = rep(.3, 4),
               alpha =  c("1", "2", "3", "\u2265 4"))

# convert to factor to set legend order
test2$alpha <- factor(test2$alpha, levels = c("1", "2", "3", "\u2265 4"))

Fig5b_legend<-  base + 
  geom_point(data = testdf,
             aes(x = MAT,
                 y = MAP),
             alpha = .3,  
             size = 3,
             shape = 21,
             color = "black",
             fill = "black",
             stroke = .5) +
  scale_radius(name = "R/L", breaks = c(.05, .1, .5, 1), range = c(2, 15)) +
  geom_point(data =test2,
             aes(x=MAT,
                 y=MAP, 
                 alpha = alpha), 
             size = 3,
             shape = 21, 
             color = "black",
             fill = "black", 
             stroke = .5) +
  labs(alpha = "N sites") +
  theme_classic()
