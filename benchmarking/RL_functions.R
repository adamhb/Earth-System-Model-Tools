
library(tidyverse)
library(lubridate)

`%!in%` <- Negate(`%in%`)


# 
# fruit_traits <- bci_traits %>%
#   rename(sp = `SP$`) %>% ##might need to do this
#   mutate(sp = tolower(sp)) %>% 
#   select(sp, FRUIT_DRY, N_SEEDFULL, DSPR_DRY) %>%
#   mutate(DRY_MASS = ifelse(!is.na(FRUIT_DRY), FRUIT_DRY, DSPR_DRY)) %>%
#   write_csv("data/BCI_fruit_traits.csv")

fruit_traits <- read_csv("data/BCI_fruit_traits.csv")

# Helper functions: 

get_counts <- function(data, parts){
  # data cols = traps, part, quantity, mass
  # parts = correspond to litterfall categories
  data %>% 
    filter(part %in% parts) %>% 
    summarise(counts = sum(quantity, na.rm = T)) %>% 
    as.numeric(.)
}

get_mass <- function(data, parts){
  # data cols = traps, part, quantity, mass
  # parts = correspond to litterfall categories
  data %>% 
    filter(part %in% parts) %>% 
    summarise(mass = round(sum(mass, na.rm = T), 1)) %>% 
    as.numeric(.)
}


get_best_est <- function(data, parts, est){
  # data cols = traps, part, quantity, mass
  # parts = correspond to litterfall categories
  # est = fruit_dry or dspr_dry estimate from bci traits
  
  # for many species reproductive parts are small and light, 
  # average fruit or diaspore weight is < the resolution of the litterfall scale (0.1 g)
  # in this case, mass will only be recorded when a large enough quantity of a part are weighed together on the scale
  # one option = estsimate mass using bci trait data (multiply counts by the average dry mass of a fruit or diaspore)
  # if, however, the recorded mass is larger than the estimate mass, observed mass should be used
  
  op1 <- round(get_counts(data, parts) * est, 1)
  op2 <- round(get_mass(data, parts), 1)
  est <- max(op1, op2, na.rm = T)
  
}

get_fr_from_cap <- function(data, n_seeds, mf_mass){
  # capsule counts, fragment counts, may represent mature fruits, but seeds from capsules, fragments may be dispersed in other traps 
  # if we estimate more seeds (diaspores, part = 2) from capsules & fragments than are observed in the plot,
  # then we can use est fruit counts * fruit mass, and forget other seeds 
  # (these seeds could come from observed capsules, fragments)
  # if we observe more seeds than are estimated from capsule and framents 
  # then we add additional fruits represented by captured seeds 
  
  est_fr_count <- get_counts(data, parts = c(3, 4))  #estimated fruits from capsules, fragments
  est_seeds <- est_fr_count * n_seeds #estimated seeds given fruit count estimate from capsules, fragments               
  capt_seeds <- get_counts(data, parts = c(2)) # observed seeds
  
  if (est_seeds > capt_seeds){                    
    est_mfm_caps = round(est_fr_count * mf_mass, 1) # estimated fruit mass from capsules, forget other seeds
  }else{                                          
    est_mfm_caps = round((est_fr_count + (ceiling((capt_seeds-est_seeds)/n_seeds))) * mf_mass, 1) 
    #add additional fruits represented by captured seeds to estimated fruit mass from capsules
  }  
  
}


get_RL_df <- function(data, fruit_traits, use_best_est = T){
  # data = for a single species (nested or filtered), cols = sp, part, quantity, mass
  # fruit traits = from BCI traits, cols = FRUIT_DRY, DSP_DRY, N_SEEDFUL, DRY_MASS
  # use_best_est = logical, if T will use bci trait values and get_best_est() function described above,
  # if F will use get_mass() function rather than get_best_est()
  
  repro_parts = c(0:10) # corresponds to cateogories listed in metadata for BCI seed rain 
  sp <- data$sp[1]      
  if(length(unique(data$sp)) > 1){
    stop("data must be observations for 1 species only")
  } 
  
  all_R <- get_mass(data, repro_parts)
  all_L <- get_mass(data, parts = 11)
  
  
  ####################################get trait data if we have trait data for sp
  if(sp %in% fruit_traits$sp){
    
    mf_mass <- fruit_traits$FRUIT_DRY[fruit_traits$sp == sp] 
    dsp_mass <- fruit_traits$DSPR_DRY[fruit_traits$sp == sp]
    n_seeds <- fruit_traits$N_SEEDFULL[fruit_traits$sp == sp]
    avail_dry_mass <- fruit_traits$DRY_MASS[fruit_traits$sp == sp]
    # "avail_dry_mass" == dry fruit mass unless fruit mass is NA, in which case == dry diaspore mass
    # for several species, fruits are better described as diaspores but recorded as parts 1, 7, 10
    # avail_dry_mass is used in get_best_est() functions when the dsp_mass should be used if fruit_dry mass is not available
    
    
    ###########################if we have everything:fruit mass and n seed trait data                
    if(!is.na(mf_mass) & !is.na(n_seeds)){
      
      est_mfm_caps <- get_fr_from_cap(data, n_seeds = n_seeds, mf_mass = mf_mass)
      dsp <-  NA # dsp counts are incoporated into estimate from capsules, fragments 
      
      if(use_best_est == T){
        mfm <- get_best_est(data, parts = c(1, 7, 10), est = mf_mass)
      }else{
        mfm <- get_mass(data, parts = c(1, 7, 10))
      }
      other_part_mass = get_mass(data, parts = c(0, 5, 6, 8, 9)) # mass from buds, flowers, immature fruit, aborted fruit
      all_R_corr = sum(est_mfm_caps, dsp, mfm, other_part_mass, na.rm = T) # sum all components
      
      #############################if we have only fruit mass or dsp mass from trait data                       
    }else if((!is.na(mf_mass) | !is.na(dsp_mass)) & is.na(n_seeds)) {
      
      
      est_mfm_caps = NA  #can't make est of mature fruit from capsules without n_seeds
      
      if(use_best_est == T){
        dsp <- get_best_est(data, parts = 2, est = dsp_mass)
        mfm <- get_best_est(data, parts = c(1, 7, 10), est = avail_dry_mass)
        
        
      }else{
        dsp <- get_mass(data, parts = 2)
        mfm <- get_mass(data, parts = c(1, 7, 10))          
      }
      
      other_part_mass = get_mass(data, parts = c(0, 3, 4, 5, 6, 8, 9)) #other part mass summed, includes capsules, fragments
      all_R_corr = sum(est_mfm_caps, dsp, mfm, other_part_mass, na.rm = T)  # sum components
      
    }else if (is.na(mf_mass) & is.na(dsp_mass)){
      #############################if sp is present in trait data but we don't have fruit mass (only n_seeds)
      est_mfm_caps = NA
      dsp = NA
      mfm = NA
      other_part_mass = NA
      all_R_corr = NA
    }
    ####################################if sp not present in triat data       
  }else {
    est_mfm_caps = NA
    dsp = NA
    mfm = NA
    other_part_mass = NA
    all_R_corr = NA
  }
  ####################################create tibble  
  tibble(sp = sp,
         all_L = as.numeric(all_L),
         all_R = as.numeric(all_R),
         all_R_corr = as.numeric(all_R_corr),
         est_from_capsules= as.numeric(est_mfm_caps),
         est_dsp_mass = as.numeric(dsp),
         est_fr_mass = as.numeric(mfm),
         other_part_mass = as.numeric(other_part_mass)
  )
  
}

