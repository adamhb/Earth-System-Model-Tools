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
#import hist
#from get_output import archive_path as archive_path
#import viz
#from myparams import *
from tabulate import tabulate


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

def get_files_of_inst(full_case_path,inst_tag,last_n_years,output_period = "monthly"):
    
    '''Returns files from the last n years (param: last_n_years)
    of a simulation (param: full_case_path) belonging to a specific ensemble member (param: inst_tag)'''
    
    substring = "clm2_" + inst_tag + ".h0"
    
    # Get the instance files
    files = find_files_with_substring(full_case_path, substring)
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
    ds_out.attrs['long_name'] = scpf_var.attrs['long_name']
    ds_out.attrs['units'] = scpf_var.attrs['units']
    return(ds_out)





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


def get_total_stem_den(ds,trees_only=True, dbh_min=None):
    
    '''This function returns a time-averaged value for
    stem density [N ha-1] with the option to exclude shrubs (pft 4)'''
    den = scpf_to_scls_by_pft(ds.FATES_NPLANT_SZPF, ds)

    den = den.mean(dim = "time")

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

def plot_multi_panel(df, x_col, y_cols, figsize=(6, 8)):
    """
    Plots multiple y-columns against one x-column in a multi-panel figure.
    
    Args:
    - df (pd.DataFrame): DataFrame containing the data.
    - x_col (str): Name of the column for the x-axis.
    - y_cols (list of str): List of column names for the y-axes.
    - figsize (tuple, optional): Figure size. Defaults to (6, 8).
    """
    # Create subplots
    fig, axs = plt.subplots(nrows=len(y_cols), figsize=figsize, sharex=True)
    
    # If only one y_col is provided, axs is not a list; make it one for consistent indexing
    if len(y_cols) == 1:
        axs = [axs]
    
    # Plot each y_col
    for i, ycol in enumerate(y_cols):
        axs[i].scatter(df[x_col], df[ycol], label=ycol)
        axs[i].legend(loc='upper right')
        axs[i].set_ylabel(ycol)
    
    axs[-1].set_xlabel(x_col)
    plt.tight_layout()
    plt.show()


##########
# Output #
##########

def store_output(case_name,case_output_df,processed_output_root,makeFig = True):
    output_path_for_case = os.path.join(processed_output_root,case_name)
    create_directory(output_path_for_case)
    df_file_name = "ensemble_output_" + case_name + ".csv"
    case_output_df.to_csv(os.path.join(output_path_for_case,df_file_name))
    
    if makeFig == True:
        fig_file_name = "ensemble_fig_" + case_name + ".png"
        plt.savefig(os.path.join(output_path_for_case,fig_file_name))
