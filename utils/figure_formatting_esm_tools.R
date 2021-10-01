library(scales)
library(cowplot)


pft.cols <- c("darkolivegreen2","darkolivegreen4","lightskyblue", "midnightblue")

#this function adds n (number of samples) to a boxplot
give.n <- function(x){
  return(c(y = median(x), label = length(x))) 
  # experiment with the multiplier to find the perfect position
}



psize <- 5
axis_size <- 20
title_size <- 25
legend_symbol_size  <- 7


adams_guides <- guides(color = guide_legend(override.aes = list(size=6)),
                       shape = guide_legend(override.aes = list(size=6)),
       fill=guide_legend(title="PFT"))

color_guide <-  guides(color = guide_legend(override.aes = list(shape = 15)))

rec.y.axis <- ylab(expression(paste('N recruits [ha'^'-1','yr'^'-1',"]")))

#set theme for the plots
adams_theme <- theme(plot.title = element_text(hjust = 0.5, size = title_size),
                     strip.text.x = element_text(size = axis_size),
                     legend.title = element_blank (),
                     axis.title.x = element_text (size = axis_size), # change the axis title
                     axis.title.y = element_text (size = axis_size),
                     axis.title.y.right = element_text (size = axis_size, color = pft.cols[2]),
                     axis.text.x = element_text (size = axis_size, colour = "black"),
                     axis.text.y = element_text (size = axis_size, colour = "black"),
                     legend.text = element_text (size = axis_size),
                     legend.spacing.x = unit(0.3, 'cm'),
                     legend.spacing.y = unit(0.3, 'cm'), #this changes the spacing between groups of legend symbols
                     legend.key.size = unit(0.9, "cm"))
adams_theme <- theme_minimal() + adams_theme




multipanel_theme <- theme(plot.title = element_text(hjust = 0.5, size = 25),
                     strip.text.x = element_text(size = 20),
                     legend.title = element_blank (),
                     axis.title.x = element_text (size = 22), # change the axis title
                     axis.title.y = element_text (size = 22),
                     axis.title.y.right = element_text (size = 25, color = pft.cols[2]),
                     axis.text.x = element_text (size = 20, colour = "black"),
                     axis.text.y = element_text (size = 20, colour = "black"),
                     legend.position = "none")



long_term_sim_theme <- theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 22),
                             strip.text.x = element_text(size = 20),
                             legend.title = element_blank (),
                             axis.title.x = element_text (size = 18), # change the axis title
                             axis.title.y = element_text (size = 18),
                             axis.title.y.right = element_text (size = 18, color = pft.cols[2]),
                             axis.text.x = element_text (size = 16, colour = "black"),
                             axis.text.y = element_text (size = 16, colour = "black"),
                             legend.text = element_text (size = 18),
                             legend.spacing.x = unit(0.3, 'cm'),
                             legend.spacing.y = unit(0.3, 'cm'), #this changes the spacing between groups of legend symbols
                             legend.key.size = unit(0.9, "cm"),
                             panel.spacing = unit(1, "lines"))
x_axis <- scale_x_date(limits = c(as.Date("2005-01-01"),as.Date("2034-12-31")), 
                       breaks = c(as.Date("2005-01-01"), as.Date("2010-01-01"), as.Date("2015-01-01"),
                                  as.Date("2020-01-01"), as.Date("2025-01-01"), as.Date("2030-01-01"), as.Date("2035-01-01")),
                       #date_breaks = "5 year", 
                       labels = date_format("%y"))

point <- geom_point(size = 2.5, stroke = 1, alpha = 1)

smooth_line <- geom_smooth(size = 1.2, method = "loess", span = .01, se = F)
smoother_line <- geom_smooth(size = 1.8, method = "loess", span = .05, se = F)
smooth_line_black <- geom_smooth(size = 1.8, method = "loess", span = .01, se = F, color = "black")



if(exists("start_date")){
  plot_title <- labs(title = paste0(start_date,"--",end_date,"\n",
                                    "driver: ",basename(driver_data_path)))
}



#png options
PNGheight = 5
PNGwidth = 8 #usually 8
PNGunits = "in"
PNGres = 100
