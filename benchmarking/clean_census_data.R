source('settings.R')
print("Cleaning census data...")
#This script cleans census data from luquillo, BCI, SERC, and ?
#to make it ready to calculate vital demographic rates such as recruitment

#########################
#set path to census data#
#########################
path_to_data <- '~/cloud/gdrive/review_paper/quant_summary_data/'
path_to_bci_data <- '~/cloud/gdrive/rec_submodel/data/observations/'
#################
###functions#####
#################

#The below function selects and renames the data fields required to calculate demographic rates
#from census data. As input it takes:
#1) raw census data
#2) a vector of the raw data field names to be included*
#*see "field names" vector in the function for the fields required
#3) the site name
clean_cen_data <- function(df, from_names, site_name){
  
  field_names <- c("cen","date","treeID","latin","sp","status","dbh") #dbh in mm
  clean_df <- df %>% select(all_of(from_names)) %>% 
    rename_at(vars(all_of(from_names)), ~ field_names) %>%
    mutate_at(.vars = "sp", .funs = tolower) %>%
    mutate_at(.vars = "treeID", .funs = as.character) %>%
    mutate_at(.vars = "date",.funs = as.Date) %>%
    drop_na(date) %>%
    add_column(site_name)
  
  return(clean_df) #returned cleaned census data ready for demographic rate analysis
}

#The below function filters non-canopy trees out of census data
#by only including trees that can grow larger than size Z (dbh mm)
filter_for_canopy_trees <- function(df,Z = 200){
  canopy_sp <- df %>% 
    select(sp, dbh) %>% 
    filter(dbh > Z) %>% 
    pull(sp) %>% unique()
  output <- df %>% filter(sp %in% canopy_sp)
}


##########################################################
###prepping luquillo data for standard cleaning function##
##########################################################

#load luquillo data
files <- list.files(path = path_to_data,pattern = '*txt')
luqfull_raw <- tibble()
for(i in 1:length(files)){
  tmp <- read_delim(file = paste0(path_to_data,'luqfull',i,'.txt'), delim = "\t") 
  luqfull_raw <- rbind(luqfull_raw,tmp)
}


#Filter out sprouts (just focus on main stems)
luqfull_raw <- luqfull_raw %>% filter(grepl("MAIN",Codes))

#Anything with Status = "alive" or an "A" in the Codes
#is considered alive, and all everything else is considered dead
luqfull_raw <- luqfull_raw %>%
mutate(status_new = case_when(
  Status == "alive" ~ "A",
  Status == "dead" ~ "D",
  grepl(";A",Codes) ~ "A",
  grepl(";A;",Codes) ~ "A",
  grepl("A;",Codes) ~ "A",
  TRUE ~ "D"
))

#write out the field names of the luquillo census data that
#we want to keep for further analysis
luqFieldNames <- c("Census","Date","TreeID","Latin","Mnemonic","status_new","DBH")


##########################################################
###prepping BCI data for standard cleaning function#######
##########################################################

#loading bci data census 1 to 7
load(paste0(path_to_bci_data,"bcifull.RData")) #object is called bci.full
#loading bci data census 8
load(paste0(path_to_bci_data,"bci.full8.rdata")) #object is called bci.tree8

#joining BCI data censuses 1-7 with the 8th census
#selecting fields for joining to bci census data 1-7
bci.tree8_raw <- bci.tree8 %>% mutate(bid = 8) %>% 
  select(bid,sp,status,dbh,date,ExactDate,treeID) %>%
  mutate_at(.vars = "bid",.funs = as.numeric)
#prepping censuses 1-7 for join with 8th  
bci.full.ahb <- bci.full %>%
  select(bid, sp,status,dbh,date,ExactDate,treeID) %>%
  mutate_at(.vars = "bid",.funs = as.numeric)
#joining bci censuses together. This is bci censuses 1-8
bcifull_raw <- bci.full.ahb %>%
  rbind(bci.tree8_raw) %>% add_column(latin = NA) %>%
  mutate(date_new = as.Date(date,origin = "1960-01-01")) 

bci_from_names <- c("bid","date_new","treeID","latin","sp","status","dbh")

#########################################################
####Apply standard cleaning function to all data sets####
#########################################################

luqfull_clean <- filter_for_canopy_trees(clean_cen_data(df = luqfull_raw,from_names = luqFieldNames,site_name = "luq"))
bcifull_clean <- filter_for_canopy_trees(clean_cen_data(df = bcifull_raw,from_names = bci_from_names,site_name = "bci"))

print("Done cleaning census data!")




