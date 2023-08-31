
from tabulate import tabulate
import sys
import getopt
import os
import xarray as xr
import functools
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import hist
import viz
import myfuncs
import warnings
from tabulate import tabulate
from myparams import *
from PyPDF2 import PdfMerger
from datetime import datetime

warnings.filterwarnings("ignore")

def arg_func(argv):
    arg_input = ""
    arg_output = ""
    arg_user = ""
    arg_help = "{0} -c <case-name> -s <start year> -e <end year>".format(argv[0])
    
    try:
        opts, args = getopt.getopt(argv[1:], "hc:s:e:d:ivp:", ["help", "case-name=","start=", 
        "end=","dbh-min=","verbose","pft-names="])
    except:
        print(arg_help)
        sys.exit(2)
    
   

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-c", "--case-name"):
            arg_case_name = arg
        elif opt in ("-s", "--start"):
            arg_start_year = arg
        elif opt in ("-e", "--end"):
            arg_end_year = arg
        elif opt in ("-d", "--dbh-min"):
            arg_dbh_min = arg
        elif opt in ("-p","--pft-names"):
            use_custom_pft_names = arg

   
    for opt, arg in opts:
        if opt in ("-v","--verbose"):
           
           print('case-name:', arg_case_name)
           print('start-year', arg_start_year)
           print('end-year:', arg_end_year)
    #print("\n".join(getFullFilePaths(arg_case_name)))

    arg_out = (arg_case_name,arg_start_year,arg_end_year,use_custom_pft_names,arg_dbh_min)
    return(arg_out)




if __name__ == "__main__":
    

    original_stdout = sys.stdout

    #parse the arguments
    arg_case_name,arg_start_year,arg_end_year,use_custom_pft_names,arg_dbh_min = arg_func(sys.argv)

    print("pft_names: ",use_custom_pft_names)

    #define paths to netcdf files
    full_paths = myfuncs.getFullFilePaths(arg_case_name,arg_start_year,arg_end_year)

    time = datetime.now().strftime("%Y%m%d-%H%M%S")

    #path to store output
    output_path = archive_path + "/" + "processed_output_" + arg_case_name + "_" +\
                  str(time)
    
    isExist = os.path.exists(output_path)
    if not isExist:
        os.makedirs(output_path)

    print("Importing data from netcdf....")
    #Import the data into an xarray
    ds = myfuncs.fix_time(xr.open_mfdataset(full_paths, decode_times=True,
                                    preprocess=functools.partial(myfuncs.preprocess, fields=hist.fields)))
   
    #get number of pfts
    n_pfts = len(ds.fates_levpft)

    if use_custom_pft_names == "True":
        pft_names = pd.Series(pft_names).iloc[list(ds.fates_levpft.values - 1)]
    else:
        pft_names = list(range(n_pfts))
    
    #Create output figures   
    
    #Environmental variables
    
    #Fire return interval
    fri = 1 / (ds.FATES_BURNFRAC.values.mean() * s_per_yr)
    
    #Burn frac and area-weighted fire intensity
    burn_frac = ds.FATES_BURNFRAC * s_per_day
    mean_burn_frac = burn_frac.values.mean()
    aw_fi = ds.FATES_FIRE_INTENSITY_BURNFRAC / burn_frac * kW_per_W
    mean_aw_fi = np.nanmean(aw_fi.values)
    viz.plot_array(aw_fi,ds,n_pfts,1,"01_FireLineIntensity","kW m-1",output_path)

    #get canopy crown area by patch age
    if 'FATES_CANOPYCROWNAREA_APPF' in hist.fields:
        canopy_crown_area_tab = myfuncs.get_rate_table(ds.FATES_CANOPYCROWNAREA_APPF,ds,"CCA",pft_names,"pft")

    #plot disturbance rate
    viz.plot_array(ds.FATES_DISTURBANCE_RATE_FIRE,ds,n_pfts,1,"01_Disturbance Rate Fire","m2 m-2 yr-1",output_path)
    mean_dist_rate_fire = np.nanmean(ds.FATES_DISTURBANCE_RATE_FIRE.values)

    #plot changes in patch age
    viz.plot_ap(ds.FATES_PATCHAREA_AP, ds, "01_Patch Area", "[m2 m-2]", output_path)
    patch_area_tab = myfuncs.get_rate_table(ds.FATES_PATCHAREA_AP, ds, 'Patch area',ds.fates_levage.values,"Patch age")

    #SMP
    viz.plot_smp(ds.SMP,ds,output_path,return_means=True)
    smp_tab = myfuncs.get_rate_table(ds.SMP,ds,"Mean SMP [MPa]",[4,7,8],"Depth [m]") 

    #Basal Area
    viz.plot_array(ds.FATES_BASALAREA_SZPF,ds,n_pfts,m2_per_ha,"02_Basal Area",'BA [m-2 ha-1]', output_path, dbh_min = arg_dbh_min)
    
    #Put basal area of last year into a table
    ba = viz.plot_array(ds.FATES_BASALAREA_SZPF,ds,n_pfts,m2_per_ha,"02_Basal Area",'BA [m-2 ha-1]', output_path, getData = True, dbh_min = arg_dbh_min)
    ba_last_year = ba.isel(time = slice(-13,-1))
    ba_last_year_tab = myfuncs.get_rate_table(ba_last_year,ds, "Mean Basal Area in last year", pft_names, "pft")

    #Stem density
    viz.plot_array(ds.FATES_NPLANT_PF,ds,n_pfts,m2_per_ha,"02_Stem Density",'N [ha-1]',output_path, subplots = True)
    sd = viz.plot_array(ds.FATES_NPLANT_PF,ds,n_pfts,m2_per_ha,"02_Stem Density",'N [ha-1]', output_path, subplots = True, getData = True)
    sd_last_year = sd.isel(time = slice(-13,-1))
    sd_last_year_tab = myfuncs.get_rate_table(sd_last_year,ds, "Stem density in last year", pft_names, "pft")

    #Canopy crown area
    viz.plot_array(ds.FATES_CANOPYCROWNAREA_PF,ds,n_pfts,1,"02_Canopy Crown Area",'[m2 m-2]', output_path)
    
    if 'FATES_CANOPYCROWNAREA_APPF' in hist.fields:

        viz.plot_appf(ds.FATES_CANOPYCROWNAREA_APPF, ds, n_pfts, "02_CCA by Age", "[m2 m2-1]", output_path)

    #Crown area
    viz.plot_array(ds.FATES_CROWNAREA_PF,ds,n_pfts,1,"02_Crown Area",'[m2 m-2]', output_path)

    #Plot total recruitment rates
    viz.plot_array(ds.FATES_RECRUITMENT_PF,ds,n_pfts,m2_per_ha,"03_Recruitment",'[N ha-1 yr-1]',
    output_path, subplots = True)

    #Calculate per capita recruitment
    rec_per_cap = myfuncs.per_capita_rate(ds.FATES_RECRUITMENT_PF, ds, unit_conversion = 1)
    viz.plot_array(rec_per_cap,ds,n_pfts,1,"03_Per Capita Recruitment",'[N per capita yr-1]',
    output_path, subplots = True)

    #Tabulate mean per capita recuitment rates
    rec_tab = myfuncs.get_rate_table(rec_per_cap,ds,"Per Capita Recruitment",pft_names,"pft")

    #Plot mortality rates
    mort_types = (ds.FATES_MORTALITY_PF, ds.FATES_MORTALITY_BACKGROUND_SZPF,
                  ds.FATES_MORTALITY_HYDRAULIC_SZPF, ds.FATES_MORTALITY_CSTARV_SZPF,
                  ds.FATES_MORTALITY_FIRE_SZPF)

    
    m_total, m_back, m_hydr, m_cstarv, m_fire = map(functools.partial(myfuncs.per_capita_rate,
                                               xds = ds, unit_conversion = yrs_per_month),mort_types)
    morts = (m_total, m_back, m_hydr, m_cstarv, m_fire)
    
    mort_titles = ["04_Total Mortality Rate","04_Background Mortality Rate","04_Hydraulic Failure Mortality",
    "04_Carbon Starvation Mortality","04_Fire Mortality"]

    for i in range(len(mort_titles)):
        viz.plot_array(morts[i],ds,n_pfts,1,mort_titles[i],'[N per capita month-1]', output_path)
  
    #Tabulate mean mortality rates
    mort_tabs = []
    for i,m in enumerate(morts):
        mort_tabs.append(myfuncs.get_rate_table(m,ds,mort_titles[i],pft_names,"pft"))

    #Plot par incident on each pft over time
    pft_par = myfuncs.incident_par(ds)
    viz.plot_array(pft_par,ds,n_pfts,1,"01_Mean PAR at Top of PFT","[W m-2 ground]", output_path)
    par_tab = myfuncs.get_rate_table(pft_par,ds,"PFT-level PAR [W m-2 ground]",pft_names,"pft")

    #Plot canopy fraction for each pft
    frac_in_canopy = myfuncs.frac_in_canopy(ds)
    viz.plot_array(frac_in_canopy,ds,n_pfts,1,"01_Fraction In Canopy","Fraction in canopy layer", output_path)
    frac_in_canopy_tab = myfuncs.get_rate_table(frac_in_canopy,ds,"Fraction In Canopy Layer",pft_names,"pft")
    

    #write text and tabular data to report
    with open(output_path + '/' + 'report.txt', 'w') as f:
        sys.stdout = f # Change the standard output to the file we created.
        
        #Environment
        print('Mean annual burn fraction:',mean_burn_frac)
        print('Fire return interval [yrs]: ',fri)
        print('Area-weighted burn intensity [kW m-1]: ',mean_aw_fi)
        print('Mean disturbance rate from fire [m2 m-2]: ',mean_dist_rate_fire)
        print("Patch area by age")
        print(patch_area_tab)
        print("Soil matric potential [MPa]")
        print(smp_tab)
        print(frac_in_canopy_tab)
        print(par_tab)
        
        if 'FATES_CANOPYCROWNAREA_APPF' in hist.fields:

            print("Canopy crown area by patch age")
            print(canopy_crown_area_tab)
        

        #Composition
        print("Mean basal area in last year [m2 ha-1]")
        print(ba_last_year_tab)
        print("Stem Density Last Year [N ha-1]")
        print(sd_last_year_tab)
        
        #Rates
        print('Recruitment [N per capita per year]')
        print(rec_tab)
        print('Mortality rates [N per capita per month]')
        for t in mort_tabs:
            print(t)
        
        sys.stdout = original_stdout
    print("finished making figs")
    #merge pdfs

    files = sorted(os.listdir(output_path))

    files = [os.path.join(output_path, f) for f in files]


    pdfs = []

    for f in files:
        if (f.endswith(".pdf")) & (f != "figs.pdf"):
                pdfs.append(f)


    merger = PdfMerger()

    for pdf in pdfs:
        merger.append(pdf)
    

    merger.write(output_path + '/' + "figs_" + str(time) + ".pdf")
    merger.close()


    print("Finished!")
    print("Output: ",output_path)


