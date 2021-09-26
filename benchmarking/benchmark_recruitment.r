#This script calculates annual mortality rates and recruitment rates
#in size class Z (dbh, mm)




#Function to calculate recruitment rates based on
#prior work by Kohyama et al., 2018 (Eqn 11) which accounts for
#recruits that may have died before the next census interval 

#Kohyama TS, Kohyama TI, Sheil D. 2018. Definition and estimation of vital rates from repeated 
#censuses: Choices, comparisons and bias corrections focusing on trees. Methods in Ecology and 
#Evolution 9: 809–821.


#This function returns the per area annual mortality rate (m2-1 yr-1) 
#Nt0 = the number of alive individuals during census i
#Nt1 = the number of alive individuals during census i+1
#S_t1 = the number of survivors from census i to census i+1
#M_t1 = the area-based (m2-1) annual mortality rate
#tCenInt = the census interval duration (years)
#plot_area = plot area in m2
getMortRate <- function(Nt0,S_t1,tCenInt,plot_area){
  M_t1 <- (Nt0 / plot_area) * (1 - (S_t1 / Nt0)^(1/tCenInt) )
  return(M_t1)
}


#This function returns the per area annual recruitment rate (m2-1 yr-1) 
getRecRate <- function(Nt0,Nt1,S_t1,M_t1){
  R <- M_t1 * (Nt1 - S_t1) / (Nt0 - S_t1)
  return(R) #area-based recruitment rate (m2-1 yr-1)
}

#function to calculate the number of alive trees during census i
#for species spp
#option to return the IDs of the alive trees
getAliveTrees <- function(df,i,spp,output = "N"){
  tmp <- df %>% filter(cen == i, status == "A", sp == spp)
  if(output == 'N'){
    return(nrow(tmp))
  }else{
    return(tmp %>% pull(treeID))
  }
}

#function returns species level Nt0, Nt1, S_t1, M_t1, and R over n census intervals for one species
#function requires cleaned census data as an input and outputs
getVitalRates_1spp <- function(cen_data,#needs to be 'clean' (see 'clean_census_data.R')
                                    spp, #species of interest
                                    plot_area) { #in m2
#number of censuses
n_cens <- length(unique(cen_data$cen))

vital_rates_spp <- tibble()  

for(i in 1:(n_cens-1) ){
      #calculate vital rates for species spp in census i
      cen_interval <- i
      tStart <- cen_data %>% filter(cen == i) %>% pull(date) %>% mean()
      tEnd <- cen_data %>% filter(cen == i+1) %>% pull(date) %>% mean()
      tCenInt <- time_length(interval(ymd(tStart), ymd(tEnd)),unit = "year")
      Nt0 <- getAliveTrees(cen_data,i,spp)
      Nt1 <- getAliveTrees(cen_data,(i+1),spp)
      TreesAlive_t0 <- getAliveTrees(cen_data,i,spp,output = "IDs")
      TreesAlive_t1 <- getAliveTrees(cen_data,(i+1),spp,output = "IDs")
      S_t1 = sum(TreesAlive_t1 %in% TreesAlive_t0)
      M_t1 <- getMortRate(Nt0,S_t1,tCenInt,plot_area)
      R_t1 <- getRecRate(Nt0,Nt1,S_t1,M_t1) #recruitment rate (yr-1 m-2)
      
    #put the vital rates calculated for census i into a df
    vital_rates_spp_cen_i <- tibble(sp = spp,
                                    cen_interval = cen_interval,
                                    tStart = tStart,
                                    tEnd = tEnd,
                                    tCenInt = tCenInt,
                                    Nt0 = Nt0,
                                    Nt1 = Nt1,
                                    S_t1 = S_t1,
                                    M_t1 = M_t1,
                                    R_t1 = R_t1)
    
    #add the vital rates for census i to the df of vital rates for all census
    vital_rates_spp <- rbind(vital_rates_spp,vital_rates_spp_cen_i)
}
return(vital_rates_spp)
}

#function to get vital rates of all species in a cleaned census dataset
getVitalRates <- function(cen_data,plot_area,site_name){
  vital_rates <- tibble()
  i <- 0
 for(sp in unique(cen_data$sp)){
  i <- i+1 
  tmp <- getVitalRates_1spp(cen_data,sp,plot_area) %>% add_column(site = site_name)
  vital_rates <- rbind(vital_rates,tmp)
  print(paste("done with",i,"of",length(unique(cen_data$sp)),"species"))
 }
return(vital_rates)  
}

#calculate recruitment and mortality rates from the census data
luq_size <- 500*320 #m2
bci_size <- 5e4 #m2

luqVitalRates <- getVitalRates(luqfull_clean,luq_size,"luq")
bciVitalRates <- getVitalRates(bcifull_clean,bci_size,"bci")




#plotting the recruitment observations
rec_benchmarks <- rec_data1 %>% ggplot(mapping = aes(x = int, y = R_Koh, color = pft)) +
  geom_point(size = 5) +
  geom_line() +
  scale_color_manual(values = pft.cols) +
  ylab(label = "recruitment rate (# ind. per ha yr)") +
  xlab(label = "census interval") +
  theme_minimal() +
  adams_theme

makePNG(fig = rec_benchmarks, path_to_output.x = path_to_benchmarking_output, file_name = "rec_benchmarks")



if(write_benchmark_csv == T){
  write_csv(rec_benchmarks_with_dates_long, path = "benchmarking/bci_rec_benchmarks_long.csv")
}

print("created benchmarks")

