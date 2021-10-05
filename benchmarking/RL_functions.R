

library(tidyverse)
#############################################Functions for trap data 
##get mass of designated litterfall part codes

get_part_mass_total <- function(data, parts, trap_size_m2 = 0.5, by_sp = F){
  #data is raw, in form sp, fecha, census, trap, part, quantity, mass
  #part codes is a vector of integerts corresponding to codes recorded in the BCI Seed Rain Metadata file
  #by_sp will return results grouped by sp if TRUE, else will return total mass
  n_traps <- n_distinct(data$trap)
  trap_size <- trap_size_m2
  
  if(by_sp == T){
    data %>% 
      filter(part %in% parts) %>% 
      group_by(sp, part) %>% 
      summarise(mass_sum_sp = sum(mass))%>% 
      mutate(mass_m2_sp = mass_sum_sp / (n_traps*trap_size))
  }else(
    data %>% 
      filter(part %in% parts) %>% 
      group_by(sp, part) %>% 
      summarise(mass_sum_sp = sum(mass))%>% 
      mutate(mass_m2_sp = mass_sum_sp / (n_traps*trap_size)) %>% 
      ungroup() %>%
      summarise(mass_m2_total = sum(mass_m2_sp)))
}


##get mature fruit equivalent mass from relevant part codes 
get_mf_eq_mass <- function (data, parts = c(3, 4, 7, 10), fruit_mass = fruit_mass_df) {
  
  #data cols are: sp, fecha, census, trap, part, quantity, mass
  #part codes are as in METADATA BCI SEED RAIN, 
  #fruit_mass_df cols are: "sp"  = sp6 code, "FRUIT_DRY" = dry mature fruit mass 
  #default parts that represent mature fruit equivalents:
  #3:capsules, 4:fragments,7:fruit w insect emergence holes, 10: damaged fruit
  # see email with Joe Wright April 3 2021; 
  # fragments should in fact be multiplied by with sp specific ratios, conservative estimate is 1:1 ratio
  
  mf_eq <- data %>% 
    filter(part %in% parts) %>% 
    group_by(sp, part) %>% 
    summarise(mf_equiv = sum(quantity)) %>% 
    ungroup() %>% 
    left_join(fruit_mass, by = "sp") %>%
    rowwise() %>% 
    mutate(mf_eq = mf_equiv * fruit_mass) %>% 
    dplyr::select(sp, part, mf_eq)
  return(mf_eq)
  
}

##get total reproductive mass total, including correction from herbivory 

get_repro_total_m2 <- function(data, fruit_mass = fruit_mass_df, trap_size_m2 = 0.5){
  
  #data cols are: sp, fecha, census, trap, part, quantity, mass
  #part codes are as in METADATA BCI SEED RAIN, 
  #fruit_mass_df cols are: "sp"  = sp6 code, "FRUIT_DRY" = dry mature fruit mass 
  
  repro_mass <- get_part_mass_total(data, parts = c(0, 1, 2, 5, 6, 8, 9), trap_size_m2 = 0.5, by_sp = T) %>% 
    #parts 3, 4, 7, 10 are excluded so as not to double count material 
    group_by(sp) %>% 
    summarise(repro_mass_sp = sum(mass_sum_sp), 
              repro_mass_sp_m2 = sum(mass_m2_sp))
  
  mfe_mass<- get_mf_eq_mass(data, fruit_mass = fruit_mass_df) %>% 
    #default parts for which to extrapolate mature fruit mass are 3, 4, 7, 10
    #per personal communication from Joe Wright: 
    #it's likely that part 4, fragments, should be multiplied by a sp specific ratio, 1:1 assumption is conservative
    group_by(sp) %>% 
    summarise(mfe_mass = sum(mf_eq, na.rm = T))
  
  repro_total <-  repro_mass %>% 
    select(sp, repro_mass_sp) %>% 
    left_join(mfe_mass, by = "sp") %>% 
    rowwise() %>%
    mutate(repro_total = sum(repro_mass_sp, mfe_mass, na.rm = T)) %>% 
    dplyr::select(sp, repro_total)
  
  n_traps <- n_distinct(data$trap)
  trap_size <- trap_size_m2
  
  repro_total_m2 <- repro_total %>% 
    mutate(repro_m2 = repro_total / (n_traps *trap_size)) %>% 
    ungroup() %>% 
    summarise(mass_m2_total = sum(repro_m2))
  
  return(as_tibble(repro_total_m2))
}


