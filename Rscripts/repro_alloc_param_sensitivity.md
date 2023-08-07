Install Packages

    #install libraries
    library(ncdf4)
    library(tidyverse)

    ## ── Attaching packages ───────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.5
    ## ✓ tidyr   1.0.2     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.5.0

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(lubridate)

    ## 
    ## Attaching package: 'lubridate'

    ## The following object is masked from 'package:base':
    ## 
    ##     date

Set path to FATES output

    #set paths
    data_path <- 'data/'
    case_names <- c("repro_alloc_sensitivity.C63f3ea8a4-Fc1f6dddf.2021-08-03","repro_alloc_sensitivity_RA-0.4.C63f3ea8a4-Fc1f6dddf.2021-08-03")
    case_aliases <- c("RA:0.1","RA:0.4")

    getFilePaths <- function(base_data_path = data_path, case_name){
     path <- paste0(data_path,'/',case_name,'/')
     fileNames <- list.files(path,pattern = '.nc')
     files <- paste0(path,fileNames)
     return(files)
    }

Define variables

    vars <- c('RECRUITMENT','DDBH_CANOPY_SCPF','DDBH_UNDERSTORY_SCPF','DDBH_SCPF','NPP_SCPF','PARPROF_DIR_CNLF','PARPROF_DIF_CNLF','PATCH_AREA_BY_AGE','NPATCH_BY_AGE','BA_SCPF', 'NPLANT_SCPF','AGB_SCPF')

Filter data by pft

    filterByPFT <- function(file = files[1],pft_int = 1,
                            map = "fates_pftmap_levscpf",
                            var='DDBH_CANOPY_SCPF'){
      
      mapVector <- ncvar_get(nc_open(file),varid = map)
      output <-  ncvar_get(nc_open(file),varid = var)[mapVector == pft_int]
      return(output)
      }

Get growth rates (cm per ind.) for the tropical tree PFT in the final
timestep (at 10 yrs), dissaggregated by size class.

    n_size_classes = 13

    growth_rates <- matrix(nrow = n_size_classes,ncol = length(case_names)) 
    ind_per_ha <- matrix(nrow = n_size_classes,ncol = length(case_names))
    agb <- matrix(nrow = n_size_classes,ncol = length(case_names))

    for(c in 1:length(case_names)){
      files <- getFilePaths(case_name = case_names[c])
      lastTimeStep <- files[length(files)]
      
      growth_rates[,c] <- filterByPFT(file = lastTimeStep, var = 'DDBH_CANOPY_SCPF')
      ind_per_ha[,c] <- filterByPFT(file = lastTimeStep, var = 'NPLANT_SCPF')
      agb[,c] <- filterByPFT(file = lastTimeStep, var = 'AGB_SCPF')
      
    }

    size_class_lower_bounds_cm <- ncvar_get(nc_open(getFilePaths(case_name = case_names[1])),varid = 'fates_levscls')
    growth_rates_cm_perInd_perYr <- cbind(growth_rates / ind_per_ha, size_class_lower_bounds_cm)
    agb <- cbind(agb, size_class_lower_bounds_cm)

    colNames <- c(case_aliases,"size_class_lower_cm")

    colnames(x = growth_rates_cm_perInd_perYr) <- colNames
    colnames(x = agb) <- colNames

    print("growth rates (cm per ind. per year) by size class")

    ## [1] "growth rates (cm per ind. per year) by size class"

    print(growth_rates_cm_perInd_perYr)

    ##          RA:0.1     RA:0.4 size_class_lower_cm
    ##  [1,] 0.1247532 0.05488317                   0
    ##  [2,]       NaN        NaN                   5
    ##  [3,]       NaN        NaN                  10
    ##  [4,]       NaN        NaN                  15
    ##  [5,]       NaN        NaN                  20
    ##  [6,]       NaN        NaN                  30
    ##  [7,]       NaN        NaN                  40
    ##  [8,]       NaN        NaN                  50
    ##  [9,]       NaN        NaN                  60
    ## [10,]       NaN        NaN                  70
    ## [11,]       NaN        NaN                  80
    ## [12,]       NaN        NaN                  90
    ## [13,]       NaN        NaN                 100

    print("AGB (kg per m2) by size class")

    ## [1] "AGB (kg per m2) by size class"

    print(agb) #kg C m2-1

    ##           RA:0.1     RA:0.4 size_class_lower_cm
    ##  [1,] 0.01491347 0.01209539                   0
    ##  [2,] 0.00000000 0.00000000                   5
    ##  [3,] 0.00000000 0.00000000                  10
    ##  [4,] 0.00000000 0.00000000                  15
    ##  [5,] 0.00000000 0.00000000                  20
    ##  [6,] 0.00000000 0.00000000                  30
    ##  [7,] 0.00000000 0.00000000                  40
    ##  [8,] 0.00000000 0.00000000                  50
    ##  [9,] 0.00000000 0.00000000                  60
    ## [10,] 0.00000000 0.00000000                  70
    ## [11,] 0.00000000 0.00000000                  80
    ## [12,] 0.00000000 0.00000000                  90
    ## [13,] 0.00000000 0.00000000                 100

Summary stats

    percent_change <- function(start,end){
      pc <- (end - start) / start
      return(na.omit(pc))
    }

    growth_rate_change <- percent_change(growth_rates_cm_perInd_perYr[,1], growth_rates_cm_perInd_perYr[,2])
    agb_change <- percent_change(agb[,1], agb[,2])

    print(paste("increasing RA from 0.1 to 0.4 results in a reduction of", as.numeric(growth_rate_change), "in growth rates (cm of dbh per ind.) at year 10"))

    ## [1] "increasing RA from 0.1 to 0.4 results in a reduction of -0.560066137733111 in growth rates (cm of dbh per ind.) at year 10"

    print(paste("increasing RA from 0.1 to 0.4 results in a reduction of", as.numeric(agb_change), "in AGB (kg per m2) at year 10"))

    ## [1] "increasing RA from 0.1 to 0.4 results in a reduction of -0.188962319193911 in AGB (kg per m2) at year 10"
