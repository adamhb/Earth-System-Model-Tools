
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


warnings.filterwarnings("ignore")

def arg_func(argv):
    arg_input = ""
    arg_output = ""
    arg_user = ""
    arg_help = "{0} -c <case-name> -s <start year> -e <end year>".format(argv[0])
    
    try:
        opts, args = getopt.getopt(argv[1:], "hc:s:e:ivp", ["help", "case-name=","start=", 
        "end=","verbose","pft-names"])
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

    for opt, arg in opts:
        if opt in ("-p", "--pft-names"):
            use_custom_pft_names = True
        else:
            use_custom_pft_names = False
    
    for opt, arg in opts:
        if opt in ("-v","--verbose"):
           
           print('case-name:', arg_case_name)
           print('start-year', arg_start_year)
           print('end-year:', arg_end_year)
    #print("\n".join(getFullFilePaths(arg_case_name)))

    arg_out = (arg_case_name,arg_start_year,arg_end_year,use_custom_pft_names)
    return(arg_out)




if __name__ == "__main__":
    

    original_stdout = sys.stdout

    #parse the arguments
    arg_case_name,arg_start_year,arg_end_year,use_custom_pft_names = arg_func(sys.argv)
    
    #define paths to netcdf files
    full_paths = myfuncs.getFullFilePaths(arg_case_name,arg_start_year,arg_end_year)

    #path to store output
    output_path = archive_path + "/" + "processed_output_" + arg_case_name
    isExist = os.path.exists(output_path)
    if not isExist:
        os.makedirs(output_path)

    print("Importing data from netcdf....")
    #Import the data into an xarray
    ds = myfuncs.fix_time(xr.open_mfdataset(full_paths, decode_times=True,
                                    preprocess=functools.partial(myfuncs.preprocess, fields=hist.fields)))
   
    #get number of pfts
    n_pfts = len(ds.fates_levpft)

    if use_custom_pft_names == True:
        pft_names = pd.Series(pft_names).iloc[list(ds.fates_levpft.values - 1)]
    else:
        pft_names = list(range(n_pfts))
    

    #Create output figures   
    print("Making figures...")
    
    #Environmental variables
    
    #Fire return interval
    fri = 1 / (ds.FATES_BURNFRAC.values.mean() * s_per_yr)
    
    #Burn frac and area-weighted fire intensity
    burn_frac = ds.FATES_BURNFRAC * s_per_day
    mean_burn_frac = burn_frac.values.mean()
    aw_fi = ds.FATES_FIRE_INTENSITY_BURNFRAC / burn_frac * kW_per_W
    mean_aw_fi = np.nanmean(aw_fi.values)
    viz.plot_array(aw_fi,ds,n_pfts,1,"FireLineIntensity","kW m-1",output_path)

    #Basal Area
    viz.plot_array(ds.FATES_BASALAREA_SZPF,ds,n_pfts,m2_per_ha,"Basal Area",'BA [m-2 ha-1]', output_path)

    #Stem density
    viz.plot_array(ds.FATES_NPLANT_PF,ds,n_pfts,m2_per_ha,"Stem Density",'N [ha-1]',output_path, subplots = True)

    #Canopy crown area
    viz.plot_array(ds.FATES_CANOPYCROWNAREA_PF,ds,n_pfts,1,"Canopy Crown Area",'[m2 m-2]', output_path)

    #Crown area
    viz.plot_array(ds.FATES_CROWNAREA_PF,ds,n_pfts,1,"Crown Area",'[m2 m-2]', output_path)

    #Plot total recruitment rates
    viz.plot_array(ds.FATES_RECRUITMENT_PF,ds,n_pfts,m2_per_ha,"Recruitment",'[N ha-1 yr-1]',
    output_path, subplots = True)

    #Calculate per capita recruitment
    rec_per_cap = myfuncs.per_capita_rate(ds.FATES_RECRUITMENT_PF, ds, unit_conversion = 1)
    viz.plot_array(rec_per_cap,ds,n_pfts,1,"Per Capita Recruitment",'[N per capita yr-1]',
    output_path, subplots = True)

    #Tabulate mean per capita recuitment rates
    rec_tab = myfuncs.get_rate_table(rec_per_cap,"Per Capita Recruitment",pft_names)

    #Plot mortality rates
    mort_types = (ds.FATES_MORTALITY_PF, ds.FATES_MORTALITY_BACKGROUND_SZPF,
                  ds.FATES_MORTALITY_HYDRAULIC_SZPF, ds.FATES_MORTALITY_CSTARV_SZPF,
                  ds.FATES_MORTALITY_FIRE_SZPF)

    
    m_total, m_back, m_hydr, m_cstarv, m_fire = map(functools.partial(myfuncs.per_capita_rate,
                                               xds = ds, unit_conversion = yrs_per_month),mort_types)
    morts = (m_total, m_back, m_hydr, m_cstarv, m_fire)
    
    mort_titles = ["Total Mortality Rate","Background Mortality Rate","Hydraulic Failure Mortality",
    "Carbon Starvation Mortality","Fire Mortality"]

    for i in range(len(mort_titles)):
        viz.plot_array(morts[i],ds,n_pfts,1,mort_titles[i],'[N per capita month-1]', output_path)
  
    #Tabulate mean mortality rates
    mort_tabs = []
    for i,m in enumerate(morts):
        mort_tabs.append(myfuncs.get_rate_table(m,mort_titles[i],pft_names))

    #Plot par incident on each pft over time
    pft_par = myfuncs.incident_par(ds)
    viz.plot_array(pft_par,ds,n_pfts,1,"Mean PAR at Top of PFT","[W m-2 ground]", output_path)
    par_tab = myfuncs.get_rate_table(pft_par,"PFT-level PAR [W m-2 ground]",pft_names)

    #Plot canopy fraction for each pft
    frac_in_canopy = myfuncs.frac_in_canopy(ds)
    viz.plot_array(frac_in_canopy,ds,n_pfts,1,"Fraction In Canopy","Fraction in canopy layer", output_path)
    frac_in_canopy_tab = myfuncs.get_rate_table(frac_in_canopy,"Fraction In Canopy Layer",pft_names)
    

    #write text and tabular data to report
    with open(output_path + '/' + 'report.txt', 'w') as f:
        sys.stdout = f # Change the standard output to the file we created.
        print('Mean annual burn fraction:',mean_burn_frac)
        print('Fire return interval [yrs]: ',fri)
        print('Area-weighted burn intensity [kW m-1]: ',mean_aw_fi)
        print('Recruitment [N per capita per year]')
        print(rec_tab)
        print(par_tab)
        print(frac_in_canopy_tab)
        print('Mortality rates [N per capita per month]')
        
        for t in mort_tabs:
            print(t)
        
        sys.stdout = original_stdout

    #merge pdfs

    files = os.listdir(output_path)

    files = [os.path.join(output_path, f) for f in files]



    pdfs = []

    for f in files:
        if f.endswith(".pdf"):
                pdfs.append(f)


    merger = PdfMerger()

    for pdf in pdfs:
        merger.append(pdf)

    merger.write(output_path + '/' + "figs.pdf")
    merger.close()


    print("Finished!")
    print("Output: ",output_path)


