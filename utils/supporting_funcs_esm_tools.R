#constants
g_per_kg <- 1000
yrs_per_day <- 1/365
days_per_yr <- 365
has_per_m2 <- 1e-4
m2_per_ha <- 1e4
sec_per_day <- 60*60*24   


#Units map of the units of each variable
varUnitsMapFates <- tibble(
  var = c('SEEDS_IN_LOCAL_ELEM','LEAF_LITTER_IN','NPP','AGB','RECRUITMENT','PFTbiomass','PFTnindivs','PFTnpp','NPP_FNRT_SCPF', 'NPP_BGSW_SCPF','NPP_BGDW_SCPF'),
  units = c(' kg ha-1 d-1 ',' g m-2 s-1 ',' g m-2 s-1 ',' g m-2 ',' indiv ha-1 yr-1 ',' g m-2 ',' indiv m-2 ',' kg m-2 yr-1 ', ' kg m-2 yr-1 ', ' kg m-2 yr-1 ', ' kg m-2 yr-1 ')
)

FatesPrettyNamesMap <- tibble(
  var = c('SEEDS_IN_LOCAL_ELEM','LEAF_LITTER_IN','NPP','AGB','RECRUITMENT','PFTbiomass','PFTnindivs','PFTnpp','ANPP'),
  prettyName = c('ReproCFlux','LeafCFlux','NPP','AGB','Recruits','PFTbiomass','PFT_N','PFT_NPP','ANPP')
)



pft_names <- c("LD_DI", "LD_DT", "ST_DI", "ST_DT")
path_to_figures <- "figures/"

options(dplyr.print_max = 1e5)
options(max.print=1e4)
options(tibble.print_max = 1e4)

#constants
megajoules_per_joule <- 1e-6
par_per_solar_rad <- 0.45 #Garcia-Rodriguez et al., 2020


dateFromFile <- function(filename){
  yr <- unlist(str_extract_all(filename, "(?<=-)[:digit:]{4}(?=-)"))
  month <- unlist(str_extract(filename, "(?<=-)[:digit:]{2}(?=-)"))
  date <- lubridate::ymd(paste(yr, month, "01" ,sep = "-"))
  return(date)
}


makePNG <- function(fig, path_to_output.x = path_to_output, file_name = "unamed_graph",
                    height=PNGheight,  width=PNGwidth, units=PNGunits, res = PNGres){
  
  #fig = SMP_fig
  #  file_name = "sMP_fig"
  model_run_time_stamp <- Sys.time() %>% 
    sub(pattern = ":", replacement = "-") %>%
    sub(pattern = ":", replacement = "-") %>%
    sub(pattern = " ", replacement = "-")
  
  png(paste0(path_to_output.x,file_name,"_",model_run_time_stamp,".png"), height=height, width=width, units=units, res = res)
  print(fig)
  dev.off()
}



#source only parts of a file
source2 <- function(file, start, end, ...) {
  file.lines <- scan(file, what=character(), skip=start-1, nlines=end-start+1, sep='\n')
  file.lines.collapsed <- paste(file.lines, collapse='\n')
  source(textConnection(file.lines.collapsed), ...)
}


#generate a log sequence
lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  round(exp(seq(log(from), log(to), length.out = length.out)))
}

