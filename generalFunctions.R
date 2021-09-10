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
  
  
  rMean = calculateMovingWindowMean(input = pull(df,var), 
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
  
  
  rMean = calculateMovingWindowMean(input = pull(df,var), 
                                    window = ma_window)
  
  output <- df %>% mutate(ema = ema) %>% 
    mutate(rMean = rMean) %>%
    cbind(index = 1:nrow(df))
  
  return(output)





