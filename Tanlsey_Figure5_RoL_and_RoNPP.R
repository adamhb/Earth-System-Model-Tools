
#Load libraries
library(tidyverse)
library(ggpubr)

library(plotbiomes)
#if plotbiomoes library won't install, try: 
#install.packages("devtools")
#devtools::install_github("valentinitnelav/plotbiomes")

library(Cairo) 
# needed to save as eps file and preserve transparency

library(cowplot) 
library(patchwork)
# needed to arrange legends for panel a and b of figure 5


source('utils/figure_formatting_esm_tools.R')

data_path <- "~/Desktop/Posters&Papers/Tansley_review/Table_S2_RAFluxData_submitted.csv"

############################################################
###########  Figures 5a and 5b in Tansley review ###########
############################################################


#load data
TableS2 <- read.csv(data_path)

# For Figure 5a: filter for only sites where NPP is reported
Fig5a_df <- TableS2 %>% 
  filter(!is.na(NPP))

# For Figure 5b: filter to remove sites which fall beyond axis limits, non-forest biomes, duplicates
Fig5b_df <- TableS2 %>% 
  mutate(Whittaker_biome = replace_na(Whittaker_biome, "NA")) %>% 
  filter(!str_detect(Whittaker_biome, "desert")) %>% 
  filter(MAP < 5000) #%>% # remove sites where precip is beyond y axis limit




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
  ggplot(aes(x = R.NPP, 
             y = R.L)) +
  geom_point(aes(x = R.NPP, 
                 y = R.L, 
                 color = Whittaker_biome)) +
  scale_color_manual(name = "Whittaker biomes",
                     breaks = names(Ricklefs_colors),
                     labels = names(Ricklefs_colors),
                     values = Ricklefs_colors) +
  geom_smooth(method = "lm",
              color = "#666666") + 
  stat_regline_equation(aes(label = ..rr.label..),  label.x = .003, label.y = .78, size = 5) +
  ylab("R/L") +
  xlab("R/NPP") +
  theme_bw() +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, size = title_size),
        strip.text.x = element_text(size = axis_size),
        axis.title.x = element_text (size = axis_size), 
        axis.title.y = element_text (size = axis_size),
        axis.title.y.right = element_text (size = axis_size, color = "blue"),
        axis.text.x = element_text (size = axis_size, colour = "black"),
        axis.text.y = element_text (size = axis_size, colour = "black") )

#ggsave(RNPP_plot, device = "png", filename = "figures/Figure5_a.png", width = 7, height = 5 )
#ggsave(RNPP_plot, device = cairo_ps, filename = "figures/Figure5_a_highres.eps", dpi = 600, width = 7, height = 5)


#####################################################################
######## Figure S1 : Plot R / ANPP vs R / L ########
#####################################################################

RANPP_plot<- Fig5a_df %>% 
  ggplot(aes(x = R.ANPP, 
             y = R.L)) +
  geom_point(aes(x = R.ANPP, 
                 y = R.L, 
                 color = Whittaker_biome)) +
  scale_color_manual(name = "Whittaker biomes",
                     breaks = names(Ricklefs_colors),
                     labels = names(Ricklefs_colors),
                     values = Ricklefs_colors) +
  geom_smooth(method = "lm",
              color = "#666666") + 
  stat_regline_equation(aes(label = ..rr.label..), label.x = .003, label.y = .9, size = 5) +
  ylab("R/L") +
  xlab("R/ANPP") +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        strip.text.x = element_text(size = axis_size),
        axis.title.x = element_text (size = axis_size), 
        axis.title.y = element_text (size = axis_size),
        axis.title.y.right = element_text (size = axis_size, color = "blue"),
        axis.text.x = element_text (size = axis_size, colour = "black"),
        axis.text.y = element_text (size = axis_size, colour = "black") )


#ggsave(RANPP_plot, device = "png", filename = "figures/FigureS1.png", width = 7, height = 5 )
#ggsave(RANPP_plot, device = cairo_ps, filename = "figures/FigureS1_highres.eps", dpi = 600, width = 7, height = 5)



#####################################################################
######## Figure S2 : Plot L vs NPP ########
#####################################################################

LNPP_plot<- Fig5a_df %>% 
  ggplot(aes(x = NPP, 
             y = L)) +
  geom_point(aes(x = NPP, 
                 y = L, 
                 color = Whittaker_biome)) +
  scale_color_manual(name = "Whittaker biomes",
                     breaks = names(Ricklefs_colors),
                     labels = names(Ricklefs_colors),
                     values = Ricklefs_colors) +
  geom_smooth(method = "lm",
              formula = y ~ 0 + x, 
              color = "#666666") + 

  stat_regline_equation(aes(label = ..rr.label..), label.x = 0, label.y = 600, size = 5) +
  ylab(expression(paste("L (gC ",m^{-2}, yr^{-1}, ")"))) +
  xlab(expression(paste("NPP (gC ",m^{-2}, yr^{-1}, ")"))) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  ylim(0, 600) +
  xlim(0, 2000) +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        strip.text.x = element_text(size = axis_size),
        axis.title.x = element_text (size = axis_size), 
        axis.title.y = element_text (size = axis_size),
        axis.title.y.right = element_text (size = axis_size, color = "blue"),
        axis.text.x = element_text (size = axis_size, colour = "black"),
        axis.text.y = element_text (size = axis_size, colour = "black") )

LNPP_plot
#ggsave(LNPP_plot, device = "png", filename = "figures/FigureS1.png", width = 7, height = 5 )
ggsave(LNPP_plot, device = cairo_ps, filename = "figures/FigureS2_highres.eps", dpi = 600, width = 7, height = 5)

#####################################################################
######## Figure 5 b : Plot  R / L in Whittaker space ########
#####################################################################
base <- ggplot() +
  # add biome polygons
  geom_polygon(data = Whittaker_biomes,
               aes(x = temp_c,
                   y = precp_cm,
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
                    values = Ricklefs_colors) +
  guides(fill = guide_legend(order = 1))



RL_plot<-  base + 
  geom_point(data =Fig5b_df,
             aes(x = MAT, 
                 y = (MAP / 10), #mm --> cm
                 size = R.L), 
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
         size = guide_legend(order = 2,
                             override.aes = list(color = "gray80"))) +

  theme(aspect.ratio = 1, 
        plot.title = element_text(hjust = 0.5, size = title_size),
        strip.text.x = element_text(size = axis_size),
        axis.title.x = element_text (size = axis_size), 
        axis.title.y = element_text (size = axis_size),
        axis.title.y.right = element_text (size = axis_size, color = "blue"),
        axis.text.x = element_text (size = axis_size, colour = "black"),
        axis.text.y = element_text (size = axis_size, colour = "black")) 

 


RL_plot
#ggsave(RL_plot, device = "png", filename = "figures/Figure5_b.png", width = 7, height = 5 )
#ggsave(RL_plot, device = cairo_ps, filename = "figures/Figure5_b_highres.eps", dpi = 600, width = 7, height = 5)

#####################################################################
######## Create legend for N_sites, Figure 5 b  ########
#####################################################################


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
  theme_classic() +
  guides(color = guide_legend(order = 1))


#ggsave(Fig5b_legend, device = cairo_ps, filename = "figures/Figure5_b_highres_legend.eps", dpi = 600, width = 7, height = 5)


#####################################################################
######## Combine Figure 5 a and b, arrange legend guides ########
#####################################################################

#Option one: combine using ggarrange, use photoshop to insert N_sites guide fromo Fig5b_legend
#Fig5_a_and_b <- ggarrange(RNPP_plot, RL_plot, labels = c("a)", "b)"), legend.grob = get_legend(RL_plot), legend = "right")
#ggsave(Fig5_a_and_b, device = cairo_ps, filename = "~/Desktop/Posters&Papers/Tansley_review/Figure5_a_and_b_highres.eps", dpi = 600, width = 11, height = 5)


#Option 2: use cowplot and patchwork libraries to extract legends and arrange using the function plot_grid()

#R/L ~ R/NPP plot with no legend
RNPP_plot_noleg <- RNPP_plot +
  guides(color = "none")


#R/L Whittaker biome plot with no legend
RL_plot_noleg<-  RL_plot +
  guides(fill = "none", 
         size = "none")

# Create plots where one legend can be extracted
P1 <- RL_plot +
  guides(size = "none")
P2 <- RL_plot +
  guides(fill = "none")
P3 <- Fig5b_legend +
  theme(aspect.ratio = 1)+
  guides(fill = "none", 
         size = "none")


#Extract legends from plots
whit <- get_legend(P1)
RLsize <- get_legend(P2)
Nsites <- get_legend(P3)

#Create blank space
blank_p <- plot_spacer() +theme_void()

#Legend row 1 = Whittaker biomes
pt1 <- plot_grid(whit, blank_p,blank_p, ncol = 3, rel_widths = c(1, .2, .2), axis = "tl") #, blank_p, ncol = 3, rel_widths = c(1, .2, .2)) # axis) #= "tl") #, ncol = 4)

#Legend row 2 = R/L size and N_sites
pt2 <- plot_grid(RLsize,  blank_p, Nsites,blank_p,  ncol = 4, axis = "l", rel_widths = c(.3, .3,1, .3)) #rel widths are wonky, no idea why this combo worked best

#Combine row 1 and row 2 to create legend
pt_1_and_2 <- plot_grid(pt1, blank_p, pt2, blank_p, nrow = 4, axis = "l", align = "v", rel_heights = c(1, .1, .5, .2)) # again, rel heights are wonky 

#add legend to two-panel figure
Fig5Final <- plot_grid(RNPP_plot_noleg, RL_plot_noleg, blank_p, pt_1_and_2, nrow = 1, ncol = 4, align = "h", axis = "t", labels = c("a)","b)", ""), rel_widths = c(1, 1, .2, .4))

#This will NOT look correct in R plots viewer, but does work when saved & viewed as .eps

ggsave(Fig5Final, device = cairo_ps, filename = "figures/Fig5_highres_w_legends.eps", dpi = 600, width = 13, height = 5)
