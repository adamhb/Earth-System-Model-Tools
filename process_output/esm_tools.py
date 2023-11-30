#funcs
import netCDF4 as nc4
import sys
import re
import getopt
import os
import xarray as xr
import functools
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import functools
from datetime import datetime
from datetime import date
#import hist
#from get_output import archive_path as archive_path
#import viz
#from myparams import *
from tabulate import tabulate
from scipy.io import netcdf as nc

# Constants
s_per_yr = 31536000
s_per_day = 3600 * 24
m2_per_ha = 1e4
m2_per_km2 = 1e6
g_per_kg = 1000
mm_per_m = 1000
months_per_yr = 12
s_per_month = 3600 * 24 * 30.4




#################################
# Directories, paths, and files #
#################################

def print_current_datetime_no_spaces():
    # Get the current date and time
    now = datetime.now()
    # Format the datetime as a string without spaces (e.g., YYYYMMDDHHMMSS)
    datetime_string = now.strftime('%Y%m%d%H%M%S')
    return datetime_string


def get_path_to_sim(case_name,case_output_root,suffix = 'lnd/hist'):
    '''case_name: name of cime case
       case_output_root: root directory where case output is store
       suffic: subdirectory where case output is stored within a case'''
    path_to_sim = os.path.join(case_output_root,case_name,suffix)
    return path_to_sim

def extract_digits(filename):
    '''Extracts the 4-digit instance tag from a land model output netcdf file'''
    
    match = re.search(r'_(\d{4})\.', filename)
    if match:
        return match.group(1)
    return None

def get_unique_inst_tags(full_case_path):
    '''Return unique list of instance tags for a case'''
    substring = ".h0."
    files = find_files_with_substring(full_case_path, substring)
    inst_tags = np.unique(np.array([extract_digits(f) for f in files]))
    return inst_tags

def find_files_with_substring(directory, substring):
    """
    Returns a list of filenames in the given directory that contain the given substring.
    
    :param directory: The path to the directory to search in.
    :param substring: The substring to search for in filenames.
    :return: A list of filenames containing the substring.
    """
    
    # List all files in the directory
    all_files = os.listdir(directory)
    
    # Filter the ones that contain the substring
    matching_files = [f for f in all_files if substring in f]
    
    return matching_files



def filter_files(file_list, before_year):
    # Regular expression to find YYYY-MM pattern
    date_pattern = re.compile(r'(\d{4})-(\d{2})')

    # Filtered list of files
    filtered_files = []

    # Iterate over files in the list
    for file_name in file_list:
        # Search for the date pattern in the filename
        match = date_pattern.search(file_name)
        if match:
            year, month = map(int, match.groups())
            file_date = date(year, month, 1)  # Create a date object

            # Check if the file date is before the specified year
            if file_date.year < before_year:
                filtered_files.append(file_name)

    return filtered_files

def get_files_of_inst(full_case_path,inst_tag,last_n_years,output_period = "monthly",
                      end_year = None):
    
    '''Returns files from the last n years (param: last_n_years)
    of a simulation (param: full_case_path) belonging to a specific ensemble member (param: inst_tag)'''
    
    substring = "clm2_" + inst_tag + ".h0"
    
    # Get the instance files
    inst_files = find_files_with_substring(full_case_path, substring)
    # Subset based on years
    if end_year is not None:
         files = filter_files(inst_files,end_year)
    
    full_files = [os.path.join(full_case_path,f) for f in files]
    
    # Get last n files from last n years
    # This assumes the output is monthly 
    if output_period == "monthly":
        last_n_files = int(last_n_years * 12)
    else:
        print("Ouput not monthly")
        return
              
    inst_files = sorted(full_files)[-last_n_files:]
    return inst_files



def create_directory(directory_path):
    try:
        os.mkdir(directory_path)
        print(f"Directory '{directory_path}' created successfully!")
    except FileExistsError:
        print(f"Directory '{directory_path}' already exists!")

def fancy_index(x,indices):
    return np.array(x)[np.array(indices)]

###################################
# Import netcdf files into xarray #
###################################

def preprocess(ds, fields):
    '''Selects the variables we want to read in 
       Drops lndgrid because we are on a single point'''
    
    return ds[fields].sel(lndgrid=0)


def fix_time(ds):
    '''Does a quick fix to adjust time vector for monthly data'''
    nmonths = len(ds.time)
    yr0 = ds['time.year'][0].values
    ds['time'] = xr.cftime_range(str(yr0), periods=nmonths, freq='MS')

    return ds

def multiple_netcdf_to_xarray(full_paths, fields):
    
    # open the dataset -- this may take a bit of time
    ds = fix_time(xr.open_mfdataset(full_paths, decode_times=True,
                                    preprocess=functools.partial(preprocess, fields=fields)))

    
    return(ds)


def load_fates_output_data(model_output_root, case_name, years, fields,
                           inst_tag = None, manual_path = None, save_processed_output = False):
    '''
    Load fates data from netcdf files

    '''
    if case_name == None:
        return None
    
    months = list(range(1, 13, 1))

    # build a list of file names based on the year and month
    
    if inst_tag == None:
    
        file_names = [f"{case_name}.clm2.h0.{str(year)}-{str(month).rjust(2, '0')}.nc"
                  for year in years for month in months]
    else:
        file_names = [f"{case_name}.clm2_{inst_tag}.h0.{str(year)}-{str(month).rjust(2, '0')}.nc"
                  for year in years for month in months]

    if manual_path != None:
        full_paths = [os.path.join(manual_path, fname) for fname in file_names]
    
    else:
    # create their full path
        full_paths = [os.path.join(model_output_root, case_name, 'lnd/hist', fname) for fname in file_names]

    # open the dataset -- this may take a bit of time
    ds = fix_time(xr.open_mfdataset(full_paths, decode_times=True,
                                    preprocess=functools.partial(preprocess, fields=fields)))
    
    print('-- your data have been read in -- ')
    
    return(ds)


def extract_variable_from_netcdf(file_path, variable_name,pft_index):
    """
    Extract a variable from a NetCDF file.

    Parameters:
    - file_path: The path to the NetCDF file.
    - variable_name: The name of the variable to extract.

    Returns:
    - The extracted variable data.
    """
    with nc4.Dataset(file_path, 'r') as dataset:
        # Check if the variable exists in the dataset
        if variable_name in dataset.variables:
            variable_data = dataset.variables[variable_name][:]
            if len(variable_data.shape) > 1:
                variable_data = variable_data[0,:]
            return variable_data.data[pft_index]
        else:
            raise ValueError(f"'{variable_name}' not found in the NetCDF file.")

def assign_var_to_nc(file_path,var_name,value,index):
    
    '''
    assigns a value to a netcdf for a particular pft and organ
    
    file = full path to netcdf
    var_name = fates parameter name
    value = parameter value to add to file
    '''
    
    # open nc file
    ncfile = nc.netcdf_file(file_path, 'a')
    
    # define param of interest
    var = ncfile.variables[var_name]
    
    # get number of dimensions
    ndim = len(var.dimensions)
    
    if var_name == "fates_stoich_nitr":
        organ_index = 0
        var[organ_index,index] = value
    
    elif "fates_leafage_class" in var.dimensions:
        var[:,index] = value
   
    elif ndim == 0:
        
        var[...] = value
    
    else:
        var[index] = value


def extract_variable_from_netcdf_specify_organ(file_path, variable_name,pft_index,organ_index):
    """
    Extract a variable from a NetCDF file.

    Parameters:
    - file_path: The path to the NetCDF file.
    - variable_name: The name of the variable to extract.
    - pft index (python index starting at 0)
    - organ index (python index starting at 0)

    Returns:
    - The extracted variable data.
    """
    with nc4.Dataset(file_path, 'r') as dataset:
        # Check if the variable exists in the dataset
        if variable_name in dataset.variables:
            variable_data = dataset.variables[variable_name][:]

            if organ_index == None:
                if len(variable_data.shape) > 1:
                    return variable_data[0,pft_index]
                elif len(variable_data.shape) ==  1:
                    return variable_data[pft_index]
                elif len(variable_data.shape) ==  0:
                    return variable_data
                else:
                    print("Dimension? Not returning anything")
                    
            else:
                if (pft_index == 0) & (organ_index is not None) & (variable_name != "fates_stoich_nitr"):
                    return variable_data[organ_index]
                else:
                    return variable_data[organ_index,pft_index]
        else:
            raise ValueError(f"'{variable_name}' not found in the NetCDF file.")



def get_parameter_file_of_inst(params_root,param_dir,inst):
    
    '''Inputs:
    1) root directory where param perturbation params are stored
    2) subdirectory for case of interest where instance-specific parameter files are stored
    3) the instance tag (e.g. 0001)
    
    Returns: full path to parameter file'''
    
    path_to_param_files = os.path.join(params_root,param_dir)
    
    substring = "_" + inst + ".nc"

    # Get the instance files
    file_oi = find_files_with_substring(path_to_param_files, substring)
    full_file_path = os.path.join(path_to_param_files,file_oi[0])
    
    return full_file_path

#####################################
# Unraveling multiplexed dimensions #
#####################################

# These functions were originally developed by
# Adriana Foster and Charlie Koven

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


def appf_to_ap_by_pft(appf_var, dataset):
    """function to reshape a fates multiplexed size and pft-indexed variable to one indexed by size class and pft
    first argument should be an xarray DataArray that has the FATES SCPF dimension
    second argument should be an xarray Dataset that has the FATES SCLS dimension 
    (possibly the dataset encompassing the dataarray being transformed)
    returns an Xarray DataArray with the size and pft dimensions disentangled"""
    n_ap = len(dataset.fates_levage)
    ds_out = (appf_var.rolling(fates_levagepft=n_ap, center=False)
            .construct("fates_levage")
            .isel(fates_levagepft=slice(n_ap-1, None, n_ap))
            .rename({'fates_levagepft':'fates_levpft'})
            .assign_coords({'fates_levage':dataset.fates_levage})
            .assign_coords({'fates_levpft':dataset.fates_levpft}))
    #ds_out.attrs['long_name'] = scpf_var.attrs['long_name']
    #ds_out.attrs['units'] = scpf_var.attrs['units']
    return(ds_out)



########
# Misc #
########

def per_capita_rate(xarr,xds,unit_conversion):
    
    xarr = xarr * unit_conversion
    
    if xarr.dims == ('time', 'fates_levscpf'):
        xarr = scpf_to_scls_by_pft(xarr, xds)
        xarr = xarr.sum(dim="fates_levscls") #sum across size classes
        
    xarr_per_cap = xarr / xds.FATES_NPLANT_PF
    
    return(xarr_per_cap)

####################
# Forest structure #
####################

def get_pft_level_basal_area(ds,dbh_min = None):
    '''Returns a numpy array of pft-specific basal area [m-2 ha-1]
       time-averaged over the timesteps in the dataset (ds)

       Input: xarray dataset containing FATES_BASALAREA_SZPF'''
    basal_area = scpf_to_scls_by_pft(ds.FATES_BASALAREA_SZPF, ds)
    basal_area = basal_area.sel(fates_levscls = slice(dbh_min,None))
    basal_area_pf = basal_area.sum(dim="fates_levscls").mean(dim = "time").values * m2_per_ha
    return basal_area_pf


def get_tree_basal_area_over_time(ds,dbh_min = None):
    '''Returns a numpy array of pft-specific basal area [m-2 ha-1]
      

       Input: xarray dataset containing FATES_BASALAREA_SZPF'''
    basal_area = scpf_to_scls_by_pft(ds.FATES_BASALAREA_SZPF, ds)
    basal_area = basal_area.sel(fates_levscls = slice(dbh_min,None))
    basal_area = basal_area.sel(fates_levpft = slice(0,3))
    basal_area_pf = basal_area.sum(dim="fates_levscls").sum(dim = "fates_levpft").values * m2_per_ha
    return basal_area_pf

def get_conifer_basal_area_over_time(ds,dbh_min = None):
    '''Returns a numpy array of pft-specific basal area [m-2 ha-1]
      

       Input: xarray dataset containing FATES_BASALAREA_SZPF'''
    basal_area = scpf_to_scls_by_pft(ds.FATES_BASALAREA_SZPF, ds)
    basal_area = basal_area.sel(fates_levscls = slice(dbh_min,None))
    basal_area = basal_area.sel(fates_levpft = slice(0,3))
    basal_area_pf = basal_area.sum(dim="fates_levscls").sum(dim = "fates_levpft").values * m2_per_ha
    return basal_area_pf

def get_oak_basal_area_over_time(ds,dbh_min = None):
    '''Returns a numpy array of pft-specific basal area [m-2 ha-1]
      

       Input: xarray dataset containing FATES_BASALAREA_SZPF'''
    basal_area = scpf_to_scls_by_pft(ds.FATES_BASALAREA_SZPF, ds)
    basal_area = basal_area.sel(fates_levscls = slice(dbh_min,None))
    basal_area = basal_area.isel(fates_levpft = 4)
    basal_area_pf = basal_area.sum(dim="fates_levscls").values * m2_per_ha
    return basal_area_pf


def get_pft_level_crown_area(ds,canopy_area_only = True, pft_index = None, over_time = False):
    '''Returns a numpy array of pft-specifi crown area [m2 m-2]
       time-averaged over the timesteps in the dataset (ds)

       Input: xarray dataset containing FATES_CANOPYCROWNAREA_PF
       pft_index starts at 0

    '''
    
    if canopy_area_only == True:
        crown_area = ds.FATES_CANOPYCROWNAREA_PF
    else:
        crown_area = ds.FATES_CROWNAREA_PF
    
    if over_time == False:
        crown_area = crown_area.mean(dim = "time")
    if over_time == True:
        crown_area = crown_area

    if pft_index is not None:
        crown_area = crown_area.isel(fates_levpft = pft_index)

    return crown_area.values    

def get_conifer_crown_area(ds,canopy_area_only = True, pft_index = None, over_time = False):
    '''Returns a numpy array of pft-specifi crown area [m2 m-2]
       time-averaged over the timesteps in the dataset (ds)

       Input: xarray dataset containing FATES_CANOPYCROWNAREA_PF
       pft_index starts at 0

    '''

    if canopy_area_only == True:
        crown_area = ds.FATES_CANOPYCROWNAREA_PF
    else:
        crown_area = ds.FATES_CROWNAREA_PF

    if over_time == False:
        crown_area = crown_area.mean(dim = "time")
    if over_time == True:
        crown_area = crown_area

    crown_area = crown_area.sel(fates_levpft = slice(0,3)).sum(dim = "fates_levpft")

    return crown_area.values


def shannon_equitability(arr):
    if isinstance(arr, np.ndarray) == True:
        total = arr.sum()
        p = arr / total
        h_i = []
        for i,a in enumerate(arr):
            h_i.append(p[i] * np.log(p[i]))
            
        return (np.array(h_i).sum() * -1) / np.log(len(arr))
    else:
        print("Error: non np array passed to shannon equatability")

def get_n_failed_pfts(arr,ba_thresh=0.1):
    '''Returns the number of pfts that had a mean basal area less than ba_thresh.
       ba_thresh should be set to a number below which we can assume they were
       competetively excluded'''
    failed = np.less(arr,ba_thresh)
    return failed.astype(int).sum()


def get_resprout_stem_den(ds,pft = None):
    
    if pft == None:
        den = ds.FATES_NPLANT_RESPROUT_PF.sum(dim = "fates_levpft").mean(dim = "time") * m2_per_ha
    else:
        den = ds.FATES_NPLANT_RESPROUT_PF.isel(fates_levpft = pft).mean(dim = "time") * m2_per_ha
    
    return den.values

def get_total_stem_den(ds,trees_only=True, dbh_min=None, resprout = False, over_time = False):
    
    '''This function returns a time-averaged value for
    stem density [N ha-1] with the option to exclude shrubs (pft 4)'''
    den = scpf_to_scls_by_pft(ds.FATES_NPLANT_SZPF, ds)

    if over_time == False:
        den = den.mean(dim = "time")
    if over_time == True:
        den = den

    if dbh_min != None:
        den = den.sel(fates_levscls = slice(dbh_min,None)).sum(dim = "fates_levscls")
    else:
        den = den.sum(dim = "fates_levscls")

    den_total = den.sum(dim="fates_levpft")
    
    if trees_only == False:
        return den_total.values * m2_per_ha
    else:
        den_shrub = den.isel(fates_levpft = 3)
        den_trees = den_total - den_shrub
        den_trees = den_trees.values * m2_per_ha
        return den_trees


def get_AGB(ds):
    '''Returns Total AGB [kg C m-2]'''
    agb_total = ds.FATES_VEGC_ABOVEGROUND.mean(dim = "time").values 
    return agb_total.item()





##################
## Productivity ##
##################

def get_total_npp(ds):
    '''Returns NPP [kg m-2 yr-1]'''
    npp_total = ds.FATES_NPP_PF.sum(dim="fates_levpft").mean(dim = "time").values * s_per_yr
    return npp_total


################
# Making plots #
################



#########
# Light #
#########

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
    
    
    
##########
# Output #
##########

def store_output_csv(case_name,file_name,case_output_df,processed_output_root):
    output_path_for_case = os.path.join(processed_output_root,case_name)
    create_directory(output_path_for_case)
    current_date_and_time = print_current_datetime_no_spaces()
    df_file_name = "ensemble_output_" + file_name + "_" + current_date_and_time + ".csv"
    case_output_df.to_csv(os.path.join(output_path_for_case,df_file_name))
       

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
     




def get_ignition_success(ds, ignition_density):
    successful_ignitions = ds.FATES_IGNITIONS.values.mean() * s_per_yr * m2_per_km2
    ignition_success_rate = successful_ignitions / ignition_density
    return np.round(ignition_success_rate,3)


def get_mean_annual_burn_frac(ds,start_date=None,end_date=None,over_time = False):
    if over_time == False:
        burnfrac = ds.FATES_BURNFRAC.sel(time = slice(start_date,end_date)).values.mean()  * s_per_yr
    else:
        burnfrac = ds.FATES_BURNFRAC.sel(time = slice(start_date,end_date)).values * s_per_yr
    return np.round(burnfrac,3)


    
def get_awfi(ds):
    '''
    Returns area-weighted fire intensity (kW m-1)
    '''
    aw_fi = ds.FATES_FIRE_INTENSITY_BURNFRAC / (ds.FATES_BURNFRAC * s_per_day) / 1000
    return aw_fi

def plot_area_weighted_fire_intensity(ds,case):
    aw_fi = get_awfi(ds)
    aw_fi.plot(marker = "o",linewidth = 0.5)
    plt.ylabel("Fire line intensity [kW m-1]")
    title = "Fire Intensity"
    plt.title(title)
    plt.savefig(output_path + "/" + case + "_" + title.replace(" ","-") + ".png")
    plt.clf()

def get_n_fire_months(ds):
    aw_fi = get_awfi(ds)
    n_months = len(aw_fi.values)
    aw_fi = aw_fi.where(~np.isnan(aw_fi), 0)
    n_fire_months_boolean = aw_fi > 0
    n_fire_months = np.sum(n_fire_months_boolean.values)
    return n_fire_months

def get_PHS_FLI_thresh(ds,FLI_thresh):
    '''
    Returns percent of fire months that had a mean severity above FLI_thresh 
    ''' 
    aw_fi = get_awfi(ds)
    n_months_greater_than_thresh_boolean = aw_fi > FLI_thresh
    n_months_greater_than_thresh = np.sum(n_months_greater_than_thresh_boolean.values)
    n_fire_months = get_n_fire_months(ds)
    PHS = n_months_greater_than_thresh / n_fire_months * 100
    return PHS

def get_PHS_FLI_thresh_isel(ds,i_start,i_end,FLI_thresh):
    '''
    Returns percent of fire months that had a mean severity above FLI_thresh 
    '''
    ds = ds.isel(time = slice(i_start,i_end))
    aw_fi = get_awfi(ds)
    n_months_greater_than_thresh_boolean = aw_fi > FLI_thresh
    n_months_greater_than_thresh = np.sum(n_months_greater_than_thresh_boolean.values)
    n_fire_months = get_n_fire_months(ds)
    PHS = n_months_greater_than_thresh / n_fire_months * 100
    return PHS





def get_combustible_fuel(ds):
    '''
    Returns the amount of combustible fuel on the landscape. Averages over the time dimension.
    '''

    age_by_fuel = agefuel_to_age_by_fuel(ds.FATES_FUEL_AMOUNT_APFC,ds)
    fates_fuel_amount_by_class = age_by_fuel.sum(dim = "fates_levage").mean(dim = "time")
    fates_trunk_fuel_amount = fates_fuel_amount_by_class.isel(fates_levfuel = 3)
    fates_combustible_fuel_amount = fates_fuel_amount_by_class.sum(dim = "fates_levfuel") - fates_trunk_fuel_amount
    return fates_combustible_fuel_amount.values

    
def get_PHS(ds,start_date,end_date):
    ds = ds.sel(time = slice(start_date,end_date))
    
    n_fire_months = get_n_fire_months(ds)
    
    #disentangle the multiplexed size class X pft dimension
    mort_fire_by_pft_and_scls = scpf_to_scls_by_pft(ds.FATES_MORTALITY_FIRE_SZPF, ds)

    #get the monthly burned area to calculate mortality rates just on the burned area
    monthly_burnfrac = ds.FATES_BURNFRAC  * s_per_month

    #sum across size classes to get pft-level mort from fire
    mort_fire_by_pft = mort_fire_by_pft_and_scls.sum(axis=2)

    #per capita mort per month per unit area that burned
    mort_fire_per_capita_per_month_per_burned_area = mort_fire_by_pft / ds.FATES_NPLANT_PF / months_per_yr / monthly_burnfrac
    greater_than_95_mort_bool = mort_fire_per_capita_per_month_per_burned_area.sel(fates_levpft = slice(1,3)).mean(axis = 1) > 0.95
    greater_than_95_mort = np.sum(greater_than_95_mort_bool.values)
    #print("conifer mort greater_than_95_mort",greater_than_95_mort)

    greater_than_95_mort_bool_all_pfts = mort_fire_per_capita_per_month_per_burned_area.mean(axis = 1) > 0.95
    greater_than_95_mort_all_pfts = np.sum(greater_than_95_mort_bool_all_pfts.values)

    frac_greater_than_95_mort = greater_than_95_mort / n_fire_months
    #print("Conifer PHS:",np.round(frac_greater_than_95_mort,3))

    frac_greater_than_95_mort_all_pfts = greater_than_95_mort_all_pfts / n_fire_months
    #print("PHS:",np.round(frac_greater_than_95_mort_all_pfts,3))
    
    return np.round(frac_greater_than_95_mort_all_pfts,3) * 100


def get_PHS_conifer(ds,start_date,end_date):
    ds = ds.sel(time = slice(start_date,end_date))
    
    n_fire_months = get_n_fire_months(ds)
    
    #disentangle the multiplexed size class X pft dimension
    mort_fire_by_pft_and_scls = scpf_to_scls_by_pft(ds.FATES_MORTALITY_FIRE_SZPF, ds)

    #get the monthly burned area to calculate mortality rates just on the burned area
    monthly_burnfrac = ds.FATES_BURNFRAC  * s_per_month

    #sum across size classes to get pft-level mort from fire
    mort_fire_by_pft = mort_fire_by_pft_and_scls.sum(axis=2)

    #per capita mort per month per unit area that burned
    mort_fire_per_capita_per_month_per_burned_area = mort_fire_by_pft / ds.FATES_NPLANT_PF / months_per_yr / monthly_burnfrac
    greater_than_95_mort_bool = mort_fire_per_capita_per_month_per_burned_area.sel(fates_levpft = slice(1,3)).mean(axis = 1) > 0.95
    greater_than_95_mort = np.sum(greater_than_95_mort_bool.values)
    #print("conifer mort greater_than_95_mort",greater_than_95_mort)

    greater_than_95_mort_bool_all_pfts = mort_fire_per_capita_per_month_per_burned_area.mean(axis = 1) > 0.95
    greater_than_95_mort_all_pfts = np.sum(greater_than_95_mort_bool_all_pfts.values)

    frac_greater_than_95_mort = greater_than_95_mort / n_fire_months
    #print("Conifer PHS:",np.round(frac_greater_than_95_mort,3))

    frac_greater_than_95_mort_all_pfts = greater_than_95_mort_all_pfts / n_fire_months
    #print("PHS:",np.round(frac_greater_than_95_mort_all_pfts,3))
    
    return np.round(frac_greater_than_95_mort,3) * 100


def get_frac_pft_level_basal_area(ds,pft_i,date,dbh_min = 0):
    basal_area = scpf_to_scls_by_pft(ds.FATES_BASALAREA_SZPF, ds) 
    basal_area = basal_area.sel(fates_levscls = slice(dbh_min,None)).sel(time = date)
    total_basal_area = basal_area.sum(axis=1).sum(axis = 0)
    basal_area_pf = basal_area.isel(fates_levpft = pft_i).sum(axis = 0)
    frac_ba = basal_area_pf.values / total_basal_area.values
    return frac_ba

def write_fire_report(ds,ignition_density,output_path,case):
    
    original_stdout = sys.stdout
    
    mean_burn_frac = get_mean_annual_burn_frac(ds)
    
    with open(output_path + '/' + 'fire_report.txt', 'w') as f:
        sys.stdout = f # Change the standard output to the file we created.
        
        print("case:",case)
        print("Mean annual burn frac:",mean_burn_frac)
        print("Mean FRI:",1 / mean_burn_frac)
        print('PHS (> 95% mort):',get_PHS(ds))
        print('PHS (> 3500 kW m-1):',get_PHS_FLI_thresh(ds,FLI_thresh))
        print('Ignition success:',get_ignition_success(ds, ignition_density))
        
        sys.stdout = original_stdout
        
        
def plus_minus_20_pct(number):
    plus_20 = number * 1.2
    minus_20 = number * 0.8
    return [plus_20,minus_20]



def check_fire_regime(ds):

    FLI_threshold = 3500

    burnfrac = ds.FATES_BURNFRAC  * s_per_yr
    print("Mean annual burn fraction",burnfrac.values.mean())
    fri = 1 / (ds.FATES_BURNFRAC.values.mean() * s_per_yr)
    print("Mean fire return interval (yrs):",fri)

    aw_fi = ds.FATES_FIRE_INTENSITY_BURNFRAC / (ds.FATES_BURNFRAC * s_per_day) / 1000
    n_months = len(aw_fi.values)
    print("n months",n_months)
    aw_fi = aw_fi.where(~np.isnan(aw_fi), 0)
    n_fire_months_boolean = aw_fi > 0
    n_fire_months = np.sum(n_fire_months_boolean.values)
    print("n fire months", n_fire_months)
    n_months_greater_than_thresh_boolean = aw_fi > FLI_threshold
    n_months_greater_than_thresh = np.sum(n_months_greater_than_thresh_boolean.values)
    print("n months > threshold",n_months_greater_than_thresh)

    print("Fraction of fire months that burned hotter than X",n_months_greater_than_thresh / n_fire_months)

    #disentangle the multiplexed size class X pft dimension
    mort_fire_by_pft_and_scls = scpf_to_scls_by_pft(ds.FATES_MORTALITY_FIRE_SZPF, ds)

    #get the monthly burned area to calculate mortality rates just on the burned area
    monthly_burnfrac = ds.FATES_BURNFRAC  * s_per_month

    #sum across size classes to get pft-level mort from fire
    mort_fire_by_pft = mort_fire_by_pft_and_scls.sum(axis=2)

    #per capita mort per month per unit area that burned
    mort_fire_per_capita_per_month_per_burned_area = mort_fire_by_pft / ds.FATES_NPLANT_PF / months_per_yr / monthly_burnfrac
    greater_than_95_mort_bool = mort_fire_per_capita_per_month_per_burned_area.sel(fates_levpft = slice(1,3)).mean(axis = 1) > 0.95
    greater_than_95_mort = np.sum(greater_than_95_mort_bool.values)
    print("conifer mort greater_than_95_mort",greater_than_95_mort)

    greater_than_95_mort_bool_all_pfts = mort_fire_per_capita_per_month_per_burned_area.mean(axis = 1) > 0.95
    greater_than_95_mort_all_pfts = np.sum(greater_than_95_mort_bool_all_pfts.values)


    frac_greater_than_95_mort = greater_than_95_mort / n_fire_months
    print("frac_greater_than_95_mort",frac_greater_than_95_mort)

    frac_greater_than_95_mort_all_pfts = greater_than_95_mort_all_pfts / n_fire_months
    print("frac_greater_than_95_mort_all_pfts",frac_greater_than_95_mort_all_pfts)
        
    
def getFullFilePaths(case,start_year,end_year):
    
    years = list(range(int(start_year), int(end_year))) 
    months = list(range(1, 13, 1)) 
    file_names = [f"{case}.clm2.h0.{str(year)}-{str(month).rjust(2, '0')}.nc"
              for year in years for month in months]

    full_paths = [os.path.join(archive_path, case, 'lnd/hist', fname) for fname in file_names]
    return full_paths

def get_rate_table(xarr,xds,var_title,indices,index_title):
    
    if xarr.dims == ('time', 'fates_levage'):
        xarr = xarr.isel(time = slice(-12,-1))
        series = pd.DataFrame(xarr.mean(axis = 0).values,
                     index=xds.fates_levage.values)

    if xarr.dims == ('time', 'fates_levagepft'):
        xarr = appf_to_ap_by_pft(xarr,xds)
        xarr = xarr / xds.FATES_PATCHAREA_AP
        series = pd.DataFrame(xarr.mean(axis = 0).values,
                     index = indices, columns=xds.fates_levage.values)
        series.loc["Total"] = series.sum()
        tab = tabulate(series, headers="keys", tablefmt="psql")
        return(tab)

    if xarr.dims == ('time', 'levgrnd'):
        grnd_depths = xds.levgrnd.values[indices]
        xarr = xarr.isel(levgrnd = indices).isel(time = slice(12,-1)) * MPa_per_mmh2o
        series = pd.DataFrame(xarr.mean(axis = 0).values,
                          index = grnd_depths).sort_values(by = 0, ascending=True).reset_index()
    else:
        series = pd.DataFrame(xarr.mean(axis = 0).values,
                          index = indices).sort_values(by = 0, ascending=False).reset_index()
    
    my_dict = {index_title:list(series.iloc[:,0]), var_title:list(series.iloc[:,1])}
    my_df = pd.DataFrame.from_dict(my_dict).set_index(index_title)
    tab = tabulate(my_df, headers='keys', tablefmt='psql')
    return(tab)




    
def is_xarray_dataset(obj):
    return isinstance(obj, xr.Dataset)

def filter_data(ds,start,stop):
    if is_xarray_dataset(ds):
        return ds.sel(time = slice(start, stop))
    else:  
        return None
    


def get_area_weighed_FLI(ds):    
    return (ds.FATES_FIRE_INTENSITY_BURNFRAC / (ds.FATES_BURNFRAC * s_per_day) / 1000)


def get_per_capita_fire_mort_by_scls(ds):
    mort_fire_by_pft_and_scls = scpf_to_scls_by_pft(ds.FATES_MORTALITY_FIRE_SZPF, ds)
    nplant_by_pft_scls = scpf_to_scls_by_pft(ds.FATES_NPLANT_SZPF, ds)
    return(mort_fire_by_pft_and_scls / nplant_by_pft_scls)
