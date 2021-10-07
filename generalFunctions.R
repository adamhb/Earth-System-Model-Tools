getFilePaths <- function(base_data_path = data_path, case_name){
  path <- paste0(data_path,'/',case_name,'/')
  fileNames <- list.files(path,pattern = '.nc')
  files <- paste0(path,fileNames)
  return(files)
}

filterByPFT <- function(file = files[1],pft_int = 1,
                        map = "fates_pftmap_levscpf",
                        var='DDBH_CANOPY_SCPF'){
  
  mapVector <- ncvar_get(nc_open(file),varid = map)
  output <-  ncvar_get(nc_open(file),varid = var)[mapVector == pft_int]
  return(output)
}


watts_to_MJ_per_day <- function(watts){
  return(watts * 3600 * 24 / 1e6)
}

#NEED TO UPDATE THIS TO INCLUDE PAR DIF AS WELL 
#add the indirect par (by getting the weighted average)
getPARatLowestLeafLayer <- function(file = files[1]){ #input are netcdf files with par in mean W m2-1 over history file interval
  r <- ncvar_get(nc_open(file),varid = 'PARPROF_DIR_CNLF')
  nonZero <- r > 0
  r <- watts_to_MJ_per_day(tail(r[nonZero],1))
  return(r) #return par at seedling layer in MJ m2 per day (mean over history file interval-- usually one month)
}


#update an exponential moving average
updateEMA <- function(currentAverage, newValue, window, alpha = 1, i){
  
  if ((i-window) < 0) {
    newValue <- currentAverage
  } else {
    wgtNew <- alpha / window
    wgtCurrent <- 1 - wgtNew
    newValue <- (currentAverage * wgtCurrent) + (newValue * wgtNew)
  }
  
  return(newValue)
}

calculateMovingWindowMean <- function(input, window){
  runningMean <- c()
  for(i in 1:length(input)){
    if (i-window >= 0) {
      values <- input[(i-window):i]
      tmp <- mean(values)
    }else{
      tmp <- NA
    }
    runningMean <- append(runningMean,tmp)
  }
  return(runningMean)
}

#create running means (exponential and moving window) of a variable
getRunningMeans <- function(df, var, ema_window, ma_window, alpha = 1){
  
  ema <- c()
  ema[1] <- as.numeric(df[1,var])
  for(i in 2:nrow(df)){
    ema[i] <- updateEMA(currentAverage = ema[i-1],
                        newValue = as.numeric(df[i,var]),
                        window = ema_window,
                        alpha = alpha,
                        i = i)
  }
  
  
  rMean = calculateMovingWindowMean(input = dplyr::pull(df,var), 
                                    window = ma_window)
  
  output <- df %>% mutate(ema = ema) %>% 
    mutate(rMean = rMean) %>%
    cbind(index = 1:nrow(df))
  
  return(output)
  
}

#moisture deficit days calculator from the regeneration submodel
def_func <- function(soil_moist, psi_crit.x = psi_crit[PFT], window){
  def <- (abs(psi_crit.x) - abs(soil_moist))*-1
  no_def <- def < 0 
  def[no_def] <- 0
  deficit_days <- c()
  for(i in 1:length(def)){
    deficit_days[i] <- ifelse(i < window, sum(def[1:i]), sum(def[(i-window):i]))
  }
  return(deficit_days)
}

mdd_ema <- function(smp, mdd_window, psi_crit, alpha = 1){
  
  ema <- c()
  ema[1] <- def_func(smp[1],psi_crit,mdd_window) * mdd_window
  for(i in 2:length(smp)){
    ema[i] <- updateEMA(currentAverage = ema[i-1],
                        newValue = def_func(smp[i],psi_crit,mdd_window) * mdd_window,
                        window = mdd_window,
                        alpha = alpha,
                        i = i)
  }
  return(ema)
}
  

getSiteVar <- function(file,varName){ #input are netcdf files
  nc <- nc_open(file)
  r <- ncvar_get(nc,varName)
  if(length(r) > 1){r <- sum(r)}
  nc_close(nc)
  return(r) #returns variable value of interest
}


getDateFromNC <- function(file){
  return(ymd(getSiteVar(file,"mcdate")))
}

getSiteVarOverTime <- function(files,varName,withDate=F){ #accepts a list of netcdf output files
  output <- c()
  date <- c()
  for(f in files){
    if(withDate == T){
      date <- append(date,getDateFromNC(f))
    }
    output <- append(output,getSiteVar(f,varName))
  }
  if(withDate == T){return(tibble(date,output))}else{
    return(output)
  }
}

getYrFromDate <- function(date){
  as.numeric(substr(date,1,4))
}

getSimYrFromDate <- function(date,start_date){
  getYrFromDate(date) - getYrFromDate(start_date) + 1
}




getManySiteVarsOverTime <- function(files,varNames){
  output <- matrix(nrow = length(files),ncol = length(varNames))
  c <- 0
  
  for(v in varNames){
    c <- c + 1
    output[,c] <- getSiteVarOverTime(files,v,F)
    print(paste("Done with variable",v,"of",length(varNames),"variables!"))
  }
  colnames(output) <- varNames
  
  dates <- c()
  for(f in files){
    datef <- getDateFromNC(f)
    dates <- append(dates,datef)
  }
  
  #add the date to the data
  output <- as_tibble(output) %>% add_column(date = dates)
  #add the simulation year to the data
  start_date <- output$date[1]
  output <- output %>% mutate(simYr = getSimYrFromDate(date,start_date))
  return(output)
}

plotVar <- function(d,x,y){
  p <- ggplot(data = d, mapping = aes_string(x,y)) +
    geom_point() +
    labs(title = case_alias)
    theme_minimal() +
    adams_theme
  return(p)
}


toLongForm <- function(d,vars){
  gather(data = d,as.name(head(vars,1)):as.name(tail(vars,1)),
         key = "var",value = "value")
}

toWideForm <- function(d,var_names_col,value_col){
  wide <- pivot_wider(data = d, 
                      names_from = as.name(var_names_col), 
                      values_from = as.name(value_col))
  return(wide)
}


plotGridofVars <- function(d,y_vars,x_var,y_axis_units,prettyNames){
  p <- d %>% select(date,simYr,y_vars) %>%
    toLongForm(y_vars) %>% 
    left_join(prettyNames,by = "var") %>%
    ggplot(aes_string(x_var,"value")) +
    geom_point() +
    theme_minimal() +
    facet_wrap(~prettyName,scales = "free") +
    ylab(y_axis_units) +
    xlab("SimYear") +
    labs(title = case_alias) +
    scale_x_date(date_labels = "%y", date_breaks = "5 years") +
    adams_theme
  return(p)
}

#this function receives a dataframe of variables (as columns)
#and standardizes their units based on a "variable unit map" (see script 'FATES_variable_units_map.R')
#It returns the dataframe with all variables in standard units (g m-2 day-1)

standardize_units <- function(d,vars,unit_map){
  
  
  
  output <- toLongForm(d,vars) %>%
    left_join(unit_map, by = "var") %>%
    mutate(std_mass = case_when(
      (value != 0 & grepl(' kg ',units)) ~ value*g_per_kg,
      TRUE ~ value
    )) %>%
    mutate(std_area = case_when(
      (value != 0 & grepl(' ha-1 ',units)) ~ std_mass*has_per_m2,
      TRUE ~ std_mass
    )) %>%
    mutate(std_time = case_when(
      (value != 0 & grepl(' yr-1 ',units)) ~ std_area*yrs_per_day,
      (value != 0 & grepl(' s-1 ',units)) ~ std_area*sec_per_day,
      TRUE ~ std_area
    )) %>%
    mutate(value = std_time) %>%
    select(simYr,date,var,value) 
  return(toWideForm(output,"var","value"))
}

#this function receives aata frame and returns the mean value of a variable 
#for each year of the data set
#NOTE: could add standard deviation to this
getAnnualMeanAfterSpecificYr <- function(df,var,start_date_str){
  df %>% filter(date > as.Date(start_date_str)) %>%
    select(simYr,date,as.name(var)) %>%
    toLongForm(vars = var) %>% 
    group_by(simYr) %>%
    summarise("value" = mean(value))
}






  
  # rMean = calculateMovingWindowMean(input = pull(df,var), 
  #                                   window = ma_window)
  # 
  # output <- df %>% mutate(ema = ema) %>% 
  #   mutate(rMean = rMean) %>%
  #   cbind(index = 1:nrow(df))
  # 
  # return(output)





