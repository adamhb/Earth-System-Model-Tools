#prep data for mean recruitment
df_plot <- vitalRates_allSites %>% #turn this into a function?
  group_by(site,cen_interval) %>%
  summarise(R_ha_yr = sum(R_t1,na.rm = T) * m2_per_ha, #the sum of all species area-based recruitment rates, for each census int at each site
            M_ha_yr = sum(M_t1,na.rm = T) * m2_per_ha) %>% #the sum of all species area-based mortality rates, for each census int at each site
  group_by(site) %>%
  summarise(R_ha_yr_mean_all_cen_ints = mean(R_ha_yr,na.rm = T), #the mean of the above across cen ints
            M_ha_yr_mean_all_cen_ints = mean(M_ha_yr,na.rm = T),
            sd_R_ha_yr_mean_all_cen_ints = sd(R_ha_yr,na.rm = T),
            sd_M_ha_yr_mean_all_cen_ints = sd(M_ha_yr,na.rm = T),
            n_cens = length(R_ha_yr)) #the mean of the above across cen ints

