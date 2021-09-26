#This script is used to plot/analyze the recruitment data
m2_per_ha <- 1e4

########################################################
#######Plotting recruitment data########################
########################################################
#preparing recruitment data for plotting
df_plot <- vitalRates_allSites %>% #turn this into a function?
  group_by(site, cen_interval) %>%
  summarise(R_ha_yr = sum(R_t1,na.rm = T) * m2_per_ha, #the sum of all species area-based recruitment rates
            M_ha_yr = sum(M_t1,na.rm = T) * m2_per_ha) #the sum of all species area-based mortality rates

#plotting the recruitment observations
rec_benchmarks_plot <- df_plot %>% 
  ggplot(mapping = aes(x = cen_interval, y = M_ha_yr, color = site)) + #get the census interval year to line up on this plot instead of the census number
  geom_point(size = 5) +
  geom_line() +
  ylab(label = "recruitment rate (# ind. per ha yr)") +
  xlab(label = "census interval") +
  theme_minimal()

rec_benchmarks_plot

#makePNG(fig = rec_benchmarks, path_to_output.x = path_to_benchmarking_output, file_name = "rec_benchmarks")