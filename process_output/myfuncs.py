#funcs
import sys
import getopt
import os
import xarray as xr
import functools
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import hist
from get_output import archive_path as archive_path
import viz
from myparams import *
from tabulate import tabulate


def arg_func(argv):
    arg_input = ""
    arg_output = ""
    arg_user = ""
    arg_help = "{0} -c <case-name> -start <start year> -end <end year>".format(argv[0])
    
    try:
        opts, args = getopt.getopt(argv[1:], "hc:s:e:iv", ["help", "case-name=","start=", 
        "end=","verbose"])
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
       if opt in ("-v","--verbose"):
           
           print('case-name:', arg_case_name)
           print('start-year', arg_start_year)
           print('end-year:', arg_end_year)
	   #print("\n".join(getFullFilePaths(arg_case_name)))

    arg_out = (arg_case_name,arg_start_year,arg_end_year)
    return(arg_out)

def getFullFilePaths(case,start_year,end_year):
    
    years = list(range(int(start_year), int(end_year))) 
    months = list(range(1, 13, 1)) 
    file_names = [f"{case}.clm2.h0.{str(year)}-{str(month).rjust(2, '0')}.nc"
              for year in years for month in months]

    full_paths = [os.path.join(archive_path, case, 'lnd/hist', fname) for fname in file_names]
    return full_paths


def preprocess(ds, fields):
    '''Selects the variables we want to read in 
       Drops lndgrid because we are on a single point'''

    return ds[fields].sel(lndgrid=0)


def fix_time(ds):
    '''Does a quick fix to adjust time vector for monthly data'''
    nmonths = len(ds.time)
    yr0 = ds['time.year'][0].values
    #ds['time'] = xr.cftime_range(str(yr0), periods=nmonths, freq='MS')
    ds['time'] = pd.date_range(start=str(yr0),periods=nmonths,freq="MS")
    return ds


def scpf_to_scls_by_pft(scpf_var, dataset):
    """function to reshape a fates multiplexed size and pft-indexed variable to one indexed by size class and pft
    first argument should be an xarray DataArray that has the FATES SCPF dimension
    second argument should be an xarray Dataset that has the FATES SCLS dimension 
    (possibly the dataset encompassing the dataarray being transformed)
    returns an Xarray DataArray with the size and pft dimensions disentangled"""
    n_scls = len(dataset.fates_levscls)
    ds_out = (scpf_var.rolling(fates_levscpf=n_scls, center=False)
            .construct("fates_levscls")
            .isel(fates_levscpf=slice(n_scls-1, None, n_scls))
            .rename({'fates_levscpf':'fates_levpft'})
            .assign_coords({'fates_levscls':dataset.fates_levscls})
            .assign_coords({'fates_levpft':dataset.fates_levpft}))
    #ds_out.attrs['long_name'] = scpf_var.attrs['long_name']
    #ds_out.attrs['units'] = scpf_var.attrs['units']
    return(ds_out)

def agefuel_to_age_by_fuel(agefuel_var, dataset):
    n_age = len(dataset.fates_levage)
    ds_out = (agefuel_var.rolling(fates_levagefuel = n_age, center=False).construct("fates_levage")
          .isel(fates_levagefuel=slice(n_age-1, None, n_age))
          .rename({'fates_levagefuel':'fates_levfuel'})
          .assign_coords({'fates_levage':dataset.fates_levage})
          .assign_coords({'fates_levfuel':np.array([1,2,3,4,5,6])}))
    return ds_out
    #ds_out.attrs['long_name'] = agefuel_var['long_name']
    #ds_out.attrs['units'] = agefuel_var['units']

def get_n_subplots(n_pfts): 
 
    if (n_pfts % 2 == 0) | (n_pfts == 1): 
        n_subplots = n_pfts 
    else: 
        n_subplots = n_pfts + 1

    if n_subplots == 1:
        ncol = 1
        nrow = 1

    else:
        ncol = 2
        nrow = n_subplots / ncol

    return (ncol,int(nrow))

def per_capita_rate(xarr,xds,unit_conversion):
    
    xarr = xarr * unit_conversion
    
    if xarr.dims == ('time', 'fates_levscpf'):
        xarr = scpf_to_scls_by_pft(xarr, xds)
        xarr = xarr.sum(axis=2) #sum across size classes
        
    xarr_per_cap = xarr / xds.FATES_NPLANT_PF
    
    return(xarr_per_cap)


def get_rate_table(arr_timeXpft, col_title, pft_names):
    
    series = pd.DataFrame(arr_timeXpft.mean(axis = 0).values,
                           index = pft_names).sort_values(by = 0, ascending=False).reset_index()

    my_dict = {'pft':list(series.iloc[:,0]), col_title:list(series.iloc[:,1])}
    my_df = pd.DataFrame.from_dict(my_dict).set_index('pft')
    mort_tab = tabulate(my_df, headers='keys', tablefmt='psql')
    return(mort_tab)

def weighted_avg_par(par_stream,frac_in_canopy):
    par_z = (par_stream.isel(fates_levcnlf = 0) * frac_in_canopy) +\
    (par_stream.isel(fates_levcnlf = 30) * (1 - frac_in_canopy))
    return(par_z)


def frac_in_canopy(xds):
    return(xds.FATES_CANOPYCROWNAREA_PF / xds.FATES_CROWNAREA_PF)

def incident_par(xds):
    f = frac_in_canopy(xds)

    par_z_dir = weighted_avg_par(xds.FATES_PARPROF_DIR_CLLL, f)
    par_z_dif = weighted_avg_par(xds.FATES_PARPROF_DIF_CLLL, f)
    par_total = par_z_dir + par_z_dif

    return(par_total.rolling(time=12, center=True).mean())

















