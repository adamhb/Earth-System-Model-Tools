#!/usr/bin/env python
# coding: utf-8

# # CZ2 Simualtion with a montane shrub 081122
# 

# In[1]:


import os

import xarray as xr
import functools
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

get_ipython().run_line_magic('matplotlib', 'inline')


# <div class="alert alert-block alert-info">
# 
# <b>NOTE:</b>  
#     Resources and tutorials
# </div>
# 
# - [NCAR python tutorial](https://ncar.github.io/python-tutorial/tutorials/yourfirst.html), which introduces python, the conda package manager, and more on github.
# - [NCAR ESDS tutorial series](https://ncar.github.io/esds/blog/tag/python-tutorial-series/), features several recorded tutorials on a wide variety of topics.
# - [GeoCAT examples](https://geocat-examples.readthedocs.io/en/latest/), with some nice plotting examples
# 
# ***
# 
# ## Reading and formatting data
# 
# ### Find the correct directories and file names

# In[2]:


n_pfts = 5
#pft_names = ["no_resprout","resprout"]
pft_names = ["pine","cedar","fir","shrub","oak"]
pft_colors = ['gold','darkorange','darkolivegreen','brown','springgreen']
#pft_colors = ["darkolivegreen","springgreen"]

#path = '/glade/u/home/adamhb/' 
#path = '/glade/u/home/adamhb/pollys_with_fire'

path = '/glade/scratch/adamhb/archive/'


# path to archived simulations
case = 'CZ2_I2000Clm51Fates_api25_main-branch_TEST.C17e2acb6a-Fe663a6e6.2023-01-19'       # case name
#subrun_name = 'run_14'

years = list(range(1901, 1999)) #python won't count the last year in the list)
months = list(range(1, 13, 1))     #python will go up to, but not including month 13)

# build a list of file names based on the year and month
file_names = [f"{case}.clm2.h0.{str(year)}-{str(month).rjust(2, '0')}.nc"
              for year in years for month in months]

# create their full path
full_paths = [os.path.join(path, case, 'lnd/hist', fname) for fname in file_names]
#full_paths = [os.path.join(path, fname) for fname in file_names]
#full_paths = [os.path.join(path, fname) for fname in file_names]

# print the last file in our list
print(file_names[0])
print(full_paths[-1])


# ### Choose variables to import

# In[3]:


# define the history variables to read in
fields = [
          #have on to import these dimensions
          'FATES_SEED_PROD_USTORY_SZ',
          'FATES_VEGC_AP',
          #patches and cohorts
          'FATES_NPATCHES',
          'FATES_PATCHAREA_AP','FATES_CANOPYAREA_AP',
          'FATES_NCOHORTS','FATES_NPATCH_AP',
          #structure
          #'FATES_LAI_AP',
          #density
          'FATES_NPLANT_PF',
          'FATES_NPLANT_SZAPPF',
          'FATES_NPLANT_SZPF',
          'FATES_NPLANT_ACPF',
          'FATES_NPLANT_CANOPY_SZPF',
          'FATES_NPLANT_USTORY_SZPF',
          #basal area
          'FATES_BASALAREA_SZPF',
          #crown_area
          'FATES_CANOPYCROWNAREA_PF',
          #'FATES_CANOPYAREA_HT',
          #'FATES_CROWNAREA_CLLL',
          'FATES_CROWNAREA_PF',
          #biomass
          'FATES_VEGC_PF','FATES_VEGC_AP','FATES_VEGC_ABOVEGROUND','FATES_VEGC_ABOVEGROUND_SZPF',
          #growth
          #'FATES_DDBH_SZPF',
          #'FATES_DDBH_CANOPY_SZAP','FATES_DDBH_USTORY_SZAP',
          #mortality
          'FATES_MORTALITY_PF',
          #'FATES_MORTALITY_CANOPY_SZAP','FATES_MORTALITY_USTORY_SZAP',
          'FATES_MORTALITY_BACKGROUND_SZPF','FATES_MORTALITY_HYDRAULIC_SZPF','FATES_MORTALITY_CSTARV_SZPF',
          #'FATES_MORTALITY_IMPACT_SZPF',
          'FATES_MORTALITY_FIRE_SZPF','FATES_MORTALITY_CROWNSCORCH_SZPF',
          #'FATES_MORTALITY_CANOPY_SZ','FATES_MORTALITY_USTORY_SZ',
          'FATES_MORTALITY_SENESCENCE_SZPF',
          #seed production and recruitment
          #'FATES_SEED_PROD_USTORY_SZ','FATES_SEED_PROD_CANOPY_SZ',
          #'FATES_SEEDS_IN',
          'FATES_SEED_BANK','FATES_SEED_ALLOC_SZPF',
          'FATES_RECRUITMENT_PF',
          #GPP and NPP
          #'FATES_GPP','FATES_GPP_SZPF',
          'FATES_NPP_PF','FATES_NPP_SZPF',
          'FATES_AUTORESP_SZPF','FATES_MAINTAR_SZPF',
          #physical environment
          #Light
          #'FATES_LAISUN_Z_CLLL','FATES_LAISHA_Z_CLLL',
          #'FATES_LAISUN_Z_CLLLPF','FATES_LAISHA_Z_CLLLPF',
          #'FATES_PARSUN_Z_CLLLPF','FATES_PARSHA_Z_CLLLPF',
          #'FATES_PARPROF_DIR_CLLL','FATES_PARPROF_DIF_CLLL',
          #'FATES_PARPROF_DIF_CLLLPF','FATES_PARPROF_DIR_CLLLPF',
          #Litter
          'FATES_CWD_ABOVEGROUND_DC',
          'FATES_FUEL_AMOUNT',
          #CLM
          #'QVEGT','QVEGE','QSOIL','TLAI','TBOT','RAIN','QBOT','Q2M',
          #'BTRAN',
          #H20
          #'SMP',
          #allocation
          #'FATES_STOREC_CANOPY_SZPF','FATES_STOREC_USTORY_SZPF',
          #fire
          'FATES_BURNFRAC','FATES_IGNITIONS','FATES_FIRE_INTENSITY_BURNFRAC',
          'FATES_FUEL_BULKD','FATES_FUEL_SAV',
          # 'FATES_DISTURBANCE_RATE_FIRE',
          # 'FATES_FUEL_AMOUNT_AP',
          # 'FATES_FIRE_INTENSITY_BURNFRAC_AP',
          # 'FATES_BURNFRAC_AP',
           'FATES_FUEL_AMOUNT_APFC',
           'FATES_FUEL_AMOUNT',
          # 'FATES_FDI',
           'FATES_FIRE_INTENSITY',
          # 'FATES_FUELCONSUMED',
          # 'FATES_NESTEROV_INDEX',
          # 'FATES_MORTALITY_CROWNSCORCH_SZPF',
          # 'FATES_SCORCH_HEIGHT_APPF','FATES_FUEL_MEF',
           'FATES_FUEL_EFF_MOIST','FATES_FUEL_MOISTURE_FC','FATES_ROS'
          # 'FATES_MORTALITY_FIRE_CFLUX_PF'
          ]


# ### Define constants

# In[4]:


#Define constants
s_per_yr = 31536000
s_per_day = 3600 * 24
m2_per_ha = 1e4
g_per_kg = 1000
mm_per_m = 1000
months_per_yr = 12

Jan1_final_year = "{}-01-01".format(years[-1])
Dec31_final_year = "{}-12-31".format(years[-1])


# <div class="alert alert-block alert-info">
# These are the raw history files that CTSM writes out.  
# By default, they include grid cell averaged monthly means for different state and flux variables.
# <br><br>
# Typically we post-process these data to generate single variable time series for an experiment. 
# This means that the full time series of model output for each variable, like rain, air temperature, or latent heat flux, are each in their own file.
# A post-processing tutorial will be available at a later date, but for now we'll just read in the monthly history files described above.
# 
# </div>

# ### Open files & load variables
# 
# This is done with the xarray function `open_mfdataset`
# 
# To make this go faster, we're going to preprocess the data so we're just reading the variables we want to look at.

# In[5]:


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


# open the dataset -- this may take a bit of time 
ds = fix_time(xr.open_mfdataset(full_paths, decode_times=True,
                                preprocess=functools.partial(preprocess, fields=fields)))

print('-- your data have been read in -- ')
print(fields)


# ### Define functions

# In[6]:


def pftLevelVarOverTime(n_pfts,var_name):
    """Function to plot a pft-level variable over time in facetted by pft"""
    fig, axes = plt.subplots(4,figsize=(10,10),sharey = True)
    for p in range(n_pfts):
        ds[var_name].isel(fates_levpft=p).plot(ax = axes[p], color = pft_colors[p])
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.subplots_adjust(left=1, bottom=1, right=2, top=2, wspace=1, hspace=2)
    
    
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


def SZPF_to_facetted_plot(var, axis_name, s_to_yr = False):
    
    var2 = scpf_to_scls_by_pft(var, ds)

    #sum across size classes to get pft-level mort from hyd. failure.
    var3 = var2.sum(axis=2)
    
    if s_to_yr == True:
        var3 = var3 * s_per_yr

    fig, axes = plt.subplots(ncols=2,nrows=2,figsize=(7,7), sharey = True)

    for p,ax in zip(range(n_pfts),axes.ravel()):
        var3.isel(fates_levpft=p).plot(x = "time", color = pft_colors[p],lw = 1, ax = ax)
        ax.set_ylabel(axis_name)
        ax.set_title(axis_name)

    plt.tight_layout()
    
    
    
#ds.FATES_PARSUN_Z_CLLLPF


#n_age = len(ds.fates_levage)

#ds.FATES_FUEL_AMOUNT_APFC.rolling(fates_levagefuel = n_age, center=False).construct("fates_levage").isel(fates_levagefuel=slice(n_age-1, None, n_age))
          # .rename({'fates_levagefuel':'fates_levfuel'})
          # .assign_coords({'fates_levage':dataset.fates_levage})
          # .assign_coords({'fates_levfuel':np.array([1,2,3,4,5,6])}))
    #return ds_out
    #ds_out.attrs['long_name'] = agefuel_var['long_name']
    #ds_out.attrs['units'] = agefuel_var['units']

# def clllpf_to_cl_ll_pf(clllpf_var, dataset):
#     n_cl = 2
#     cl_llpf = (clllpf_var.rolling(fates_levcnlfpf = n_cl, center=False).construct("fates_levcl")
#         .isel(fates_levcnlfpf=slice(n_cl-1, None, n_cl))
#         .rename({'fates_levcnlfpf':'fates_levllpf'})
#         .assign_coords({'fates_levcl':np.array([1,2])})
#         .assign_coords({'fates_levllpf':np.arange(0,120)}))
    
#     cl_ll_pf = (cl_llpf.rolling(fates_levllpf = n_pfts, center=False).construct("fates_levpft")
#                 .isel(fates_levllpf=slice(n_pfts-1, None, n_pfts))
#                 .rename({"fates_levllpf":"fates_levll"})
#                 .assign_coords({"fates_levll":np.arange(0,30)})
#                 .assign_coords({"fates_levpft":dataset.fates_levpft}))
    
#     return cl_ll_pf 

# test = clllpf_to_cl_ll_pf(ds.FATES_LAISUN_Z_CLLLPF,ds)
# test


# ## Patch dynamics

# In[ ]:


ds.FATES_PATCHAREA_AP.plot(x = "time")


# In[ ]:


ds.FATES_NCOHORTS.plot(x = "time")


# ## Biomass
# ### Total biomass over time

# In[ ]:


ds.FATES_VEGC_PF.values


# In[ ]:


for p in range(n_pfts):
    ds.FATES_VEGC_PF.isel(fates_levpft=p).plot(color = pft_colors[p], linewidth = 5)


# ## Basal area

# ### PFT-level basal area over time

# In[7]:


#disentangle the multiplexed size class X pft dimension
basal_area = scpf_to_scls_by_pft(ds.FATES_BASALAREA_SZPF, ds)

#sum across size classes to get pft-level ba
basal_area_pf = basal_area.sum(axis=2)

#plot pft-level basal area over time
for p in range(n_pfts):
    ba_per_ha = basal_area_pf.isel(fates_levpft=p) * m2_per_ha
    ba_per_ha.plot(x = "time", color = pft_colors[p],lw = 5, add_legend = True)

plt.ylabel('BA [m-2 ha-1]', fontsize=12)
plt.title('Basal Area')


# ### Mean size class distribution for time slice

# In[ ]:


#get basal area in final time step
basal_area_final_timestep = basal_area.sel(time = slice("1930-01-01","1931-01-01")).mean(axis = 0) 

fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7))

for p,ax in zip(range(n_pfts),axes.ravel()):
    ba_m2_ha = basal_area_final_timestep.isel(fates_levpft=p) * m2_per_ha
    ba_m2_ha.plot(ax=ax, x = "fates_levscls", marker = "o", lw = 0, 
                                                        color = pft_colors[p], markersize = 10)
    ax.title.set_text('{}'.format(pft_names[p]))
    ax.set_ylabel('BA [m-2 ha-1]')
    

plt.tight_layout()
plt.subplots_adjust(hspace=0.4,wspace=0.4)


# ## Crown area

# In[8]:


plt.rcParams.update({'axes.titlesize': 'large', 'axes.labelsize':'large'})

for p in range(n_pfts):
    ds.FATES_CANOPYCROWNAREA_PF.isel(fates_levpft=p).plot(color = pft_colors[p], linewidth = 5)
    
plt.title(" ")
plt.xlabel("year")


# ## Canopy area height

# In[ ]:


#ds.FATES_CANOPYAREA_HT.plot(x = "time")


# ## Stem density

# ### Change in stem density over time

# In[ ]:


#N per area
n_per_area_pft = ds.FATES_NPLANT_PF


# In[ ]:


ds.FATES_NPLANT_PF.values


# In[ ]:


#pftLevelVarOverTime(4,"FATES_NPLANT_PF")
stem_den = ds.FATES_NPLANT_PF * m2_per_ha


#plot pft-level basal area over time
for p in range(n_pfts):
    stem_den.isel(fates_levpft=p).plot(x = "time", color = pft_colors[p],lw = 5)

plt.ylabel('Stem density [N m-2]', fontsize=12)
plt.title('Stem density')
plt.show()


# ### On separate axes

# In[9]:


fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7))

for p,ax in zip(range(n_pfts),axes.ravel()):
    
    density = ds.FATES_NPLANT_PF * m2_per_ha
    
    density.isel(fates_levpft=p).plot(ax = ax, 
                                                 x = "time", 
                                                 color = pft_colors[p],
                                                 lw = 5)
    ax.set_title('{}'.format(pft_names[p]))
    ax.set_ylabel('Density [N ha-1]')
    
    if p == 3:
        ax.set_yscale("log")
    
plt.tight_layout()
plt.subplots_adjust(hspace=0.6,wspace=0.6)


# ### Stem density of smallest size class

# In[ ]:


density_by_pft_and_scls = scpf_to_scls_by_pft(ds.FATES_NPLANT_SZPF,ds)
density_in_smallest_size_class = density_by_pft_and_scls.isel(fates_levscls = 1)

fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7))

for p,ax in zip(range(n_pfts),axes.ravel()):
    
    density = density_in_smallest_size_class * m2_per_ha
    
    density.isel(fates_levpft=p).plot(ax = ax, 
                                                 x = "time", 
                                                 color = pft_colors[p],
                                                 lw = 2)
    ax.set_title('{}'.format(pft_names[p]))
    ax.set_ylabel('Density [N ha-1]')
    
    if p == 3:
        ax.set_yscale("log")
    
plt.tight_layout()
plt.subplots_adjust(hspace=0.6,wspace=0.6)


# ## NPP

# ### NPP over time

# In[ ]:


#ds.GPP.rolling(time=6, center=True).mean()
for p in range(n_pfts):
    mv_avg = ds.FATES_NPP_PF.isel(fates_levpft=p).rolling(time=12, center=True).mean()\
    * s_per_yr
    mv_avg.plot(x = "time", color = pft_colors[p],lw = 3)

plt.ylabel("NPP [kg m-2 yr-1]")
plt.title("NPP [kg m-2 yr-1]")
plt.show()


# ### Shrub NPP per individual in the smallest size class

# In[ ]:


size_class = 0

NPP_by_pft_and_scls = scpf_to_scls_by_pft(ds.FATES_NPP_SZPF, ds)

#Shrub NPP in smallest size class
shrub_NPP_smallest_scls = NPP_by_pft_and_scls[:,(n_pfts - 1),size_class]

#Shrub nunmber density in size class x
density_by_pft_and_scls = scpf_to_scls_by_pft(ds.FATES_NPLANT_SZPF,ds)
shrub_density_smallest_scls = density_by_pft_and_scls.isel(fates_levpft = 3).data[:,size_class]


#NPP per individual in smallest size class (0-5 cm)
shrub_NPP_per_ind_smallest_size_class = (shrub_NPP_smallest_scls / shrub_density_smallest_scls) * s_per_yr

plt.plot(shrub_NPP_per_ind_smallest_size_class, linewidth = 1, color = "black")
plt.title("Mean NPP per individual in smallest scls [kg C ind. yr-1]")
plt.show()


# ### Autotrophic respiration

# In[ ]:


respiration = scpf_to_scls_by_pft(ds.FATES_AUTORESP_SZPF, ds) * s_per_yr

respiration = respiration.rolling(time=12, center=True).mean()

for p in range(n_pfts):
    respiration.sum(axis = 2).isel(fates_levpft=p).plot(x = "time", color = pft_colors[p],lw = 2) 

plt.ylabel("Respiration [kg m-2 yr-1]")
plt.title("Respiration [kg m-2 yr-1]")


# ## Seed Production and Recruitment

# ### Allocation to Seed

# In[ ]:


alloc_to_seed = scpf_to_scls_by_pft(ds.FATES_SEED_ALLOC_SZPF,ds).sum(axis = 2)

for p in range(n_pfts):
    seed_alloc = alloc_to_seed.isel(fates_levpft=p) * s_per_yr
    seed_alloc.rolling(time = 12, center = True).mean().plot(x = "time", color = pft_colors[p],lw = 3)

plt.ylabel("Allocation to seed [Kg m-2 yr-1]")
plt.title("Allocation to seed")


# ### Seed production (all PFTs)

# This is local and external (i.e. seed rain)

# In[ ]:


seed_production_m2_yr = ds.FATES_SEEDS_IN * s_per_yr
seed_production_m2_yr.rolling(time = 12, center = True).mean().plot()
plt.ylabel("Seed Production [kg m-2 yr-1]")


# ### Seed Bank Size

# In[ ]:


seed_bank = ds.FATES_SEED_BANK * s_per_day #s per day fixes issue in fates code (issue # 900)
seed_bank.plot(linewidth = 3, color = "black")
plt.ylabel("Seed Bank [Kg m-2]")
plt.title("Seed Bank")


# ### Per area recruitment

# In[ ]:


fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7))

for p,ax in zip(range(n_pfts),axes.ravel()):
    rec = ds.FATES_RECRUITMENT_PF.isel(fates_levpft=p) * m2_per_ha
    rec.plot(x = "time", color = pft_colors[p], ax=ax)
    ax.set_ylabel("[N ha-1 yr-1]")
    ax.set_title('{} recruitment'.format(pft_names[p]))
    
plt.tight_layout()


# In[10]:


rec_per_cap = ds.FATES_RECRUITMENT_PF / ds.FATES_NPLANT_PF / months_per_yr

fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7))

for p,ax in zip(range(n_pfts),axes.ravel()):
    rec = rec_per_cap.isel(fates_levpft=p).rolling(time=24, center=True).mean()
    rec.plot(x = "time", color = pft_colors[p], ax=ax)
    ax.set_ylabel("[N ha-1 month-1]")
    ax.set_title('{} recruitment'.format(pft_names[p]))
    
plt.tight_layout()


# ### mean rec per capita

# In[13]:


for p in range(n_pfts):
    print(pft_names[p],":", rec_per_cap.isel(fates_levpft=p).values.mean())


# ## Growth

# ### Growth rates by size class (snapshot of timestep X)

# In[ ]:


ds.fates_levscls


# In[ ]:


growth = scpf_to_scls_by_pft(ds.FATES_DDBH_SZPF,ds)

fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7))

for p,ax in zip(range(n_pfts),axes.ravel()):
    rec = growth.isel(fates_levpft=p).sel(time = "1910-07-01")
    rec.plot(x = "fates_levscls", color = pft_colors[p], ax=ax, marker = "o")
    ax.set_ylabel("Growth [m m-2 yr-1]")
    ax.set_title(" ")

plt.tight_layout()


# ## Shrub growth over time in the smallest size class
# 
# This is an area-based measure?

# In[ ]:


#growth.isel(fates_levpft=3).isel(fates_levscls=0).plot()
#plt.title("Shrub relative growth rate")


# ## Mortality

# ### Total area-based mortality

# In[ ]:


fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7),sharey=True)

for p,ax in zip(range(n_pfts),axes.ravel()):
    mort = ds.FATES_MORTALITY_PF.isel(fates_levpft=p).rolling(time=12, center=True).mean() * m2_per_ha
    mort.plot(x = "time", color = pft_colors[p],lw = 3, ax = ax)
    ax.set_ylabel("Mortality [N ha-1 yr-1]")
    ax.set_title('{}'.format(pft_names[p]))

plt.tight_layout()


# Total per capita mortality

# In[11]:


mort_per_cap = ds.FATES_MORTALITY_PF / ds.FATES_NPLANT_PF / months_per_yr

fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7),sharey=False)

for p,ax in zip(range(n_pfts),axes.ravel()):
    mort_per_cap.isel(fates_levpft=p).plot(x = "time", color = pft_colors[p],lw = 1, ax = ax)
    ax.set_ylabel("Mortality [N per capita]")
    ax.set_title('{}'.format(pft_names[p]))

plt.tight_layout()


# In[12]:


for p in range(n_pfts):
    print(pft_names[p], ":", mort_per_cap.isel(fates_levpft=p).values.mean())


# ### Background mortality rate

# In[ ]:


#disentangle the multiplexed size class X pft dimension
mort_back = scpf_to_scls_by_pft(ds.FATES_MORTALITY_BACKGROUND_SZPF, ds)
#mort_hydr

#sum across size classes to get pft-level mort from hyd. failure.
mort_back_pf = mort_back.sum(axis=2) * m2_per_ha

fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7), sharey = False)

for p,ax in zip(range(n_pfts),axes.ravel()):
    mort_back_pf.isel(fates_levpft=p).plot(x = "time", color = pft_colors[p],lw = 1, ax = ax)
    ax.set_ylabel("[N ha-1 yr-1]")
    ax.set_title("Background Mortality")

plt.tight_layout()


# Background mort rate per capita

# In[ ]:


background_mort_per_capita = mort_back.sum(axis=2) / ds.FATES_NPLANT_PF / mont
background_mort_per_capita

fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7), sharey = True)

for p,ax in zip(range(n_pfts),axes.ravel()):
    background_mort_per_capita.isel(fates_levpft=p).plot(x = "time", color = pft_colors[p],lw = 1, ax = ax)
    ax.set_ylabel("[N per capita yr-1]")
    ax.set_title("Background Mortality [{}]".format(pft_names[p]))

plt.tight_layout()


# ### Mortality from hydraulic failure

# In[16]:


#disentangle the multiplexed size class X pft dimension
import matplotlib.dates as mdates

mort_hydr = scpf_to_scls_by_pft(ds.FATES_MORTALITY_HYDRAULIC_SZPF, ds)
#mort_hydr

#sum across size classes to get pft-level mort from hyd. failure.
mort_hydr_pf = mort_hydr.sum(axis=2) * m2_per_ha

fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7), sharey = False)

for p,ax in zip(range(n_pfts),axes.ravel()):
    mort_hydr_pf.isel(fates_levpft=p).plot(x = "time", color = pft_colors[p],lw = 1, ax = ax)
    ax.set_ylabel("[N ha-1 yr-1]")
    ax.title.set_text('{} Mort Hyd. Fail. '.format(pft_names[p]))
    #ax.xaxis.set_major_locator(mdates.MonthLocator(interval = 120))
    #ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

plt.tight_layout()


# ### Mortality from hydraulic failure (per capita)

# In[15]:


mort_hydr_per_capita = mort_hydr.sum(axis=2) / ds.FATES_NPLANT_PF / months_per_yr

fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7), sharey = True)

for p,ax in zip(range(n_pfts),axes.ravel()):
    mort_hydr_per_capita.isel(fates_levpft=p).plot(x = "time", color = pft_colors[p],lw = 1, ax = ax)
    ax.set_ylabel("[N per capita month-1]")
    ax.title.set_text('{} Mort Hyd. Fail. '.format(pft_names[p]))

plt.tight_layout()


# In[ ]:


for p in range(n_pfts):
    print(pft_names[p], ":", mort_hydr_per_capita.isel(fates_levpft=p).values.mean())


# ### Mortality from Carbon Starvation

# In[25]:


#disentangle the multiplexed size class X pft dimension
mort_c = scpf_to_scls_by_pft(ds.FATES_MORTALITY_CSTARV_SZPF, ds)

#sum across size classes to get pft-level mort from hyd. failure.
mort_c_pf = mort_c.sum(axis=2) * m2_per_ha

fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7), sharey = False)

for p,ax in zip(range(n_pfts),axes.ravel()):
    mort_c_pf.isel(fates_levpft=p).plot(x = "time", color = pft_colors[p],lw = 1, ax = ax)
    ax.set_ylabel("[N ha-1 yr-1]")
    ax.title.set_text('{} Mort Carb Starv. '.format(pft_names[p]))

plt.tight_layout()


# ### Mortality from Carbon Starvation (per capita)

# In[26]:


#disentangle the multiplexed size class X pft dimension
mort_c = scpf_to_scls_by_pft(ds.FATES_MORTALITY_CSTARV_SZPF, ds)

#sum across size classes to get pft-level mort from hyd. failure.
mort_c_pf = mort_c.sum(axis=2)



#Mort per capita
mort_c_pf = mort_c_pf / ds.FATES_NPLANT_PF / months_per_yr

fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7), sharey = True)

for p,ax in zip(range(n_pfts),axes.ravel()):
    mort_c_pf.isel(fates_levpft=p).plot(x = "time", color = pft_colors[p],lw = 1, ax = ax)
    ax.set_ylabel("[N per capita month-1]")
    ax.set_title("{} Carb. Starv.".format(pft_names[p]))

plt.tight_layout()


# #### Pft-level mean

# In[27]:


for p in range(n_pfts):
    print(pft_names[p], ":", mort_c_pf.isel(fates_levpft=p).values.mean())


# ### PFT-level mean carbon starvation normalized by hydraulic failure mort

# In[ ]:





# In[ ]:


output = {}
for i in range(n_pfts):
    output[pft_names[i]] = None

for p in range(n_pfts):
    cstarv = mort_c_pf.isel(fates_levpft=p).values.mean()
    hyrd_f = mort_hydr_per_capita.isel(fates_levpft=p).values.mean()
    output[pft_names[p]] = cstarv / hyrd_f

print(output)


# In[ ]:


pd.DataFrame.from_dict(output,orient='index',columns=['cstarv_per_hydr_f']).sort_values(by = 'cstarv_per_hydr_f',ascending=False)


# ### Mortality from fire

# In[ ]:


#disentangle the multiplexed size class X pft dimension
mort_fire_by_pft_and_scls = scpf_to_scls_by_pft(ds.FATES_MORTALITY_FIRE_SZPF, ds) * m2_per_ha

#sum across size classes to get pft-level mort from fire
mort_fire_by_pft = mort_fire_by_pft_and_scls.sum(axis=2)

fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7), sharey = False)

for p,ax in zip(range(n_pfts),axes.ravel()):
    mort_fire_pft_x = mort_fire_by_pft.isel(fates_levpft=p).plot(x = "time",lw = 1, ax = ax, color = pft_colors[p])
    ax.title.set_text('{} Mort Fire '.format(pft_names[p]))
    ax.set_ylabel("[N per ha-1 yr-1]")

plt.tight_layout()


# ### Mortality from fire (per capita)

# In[17]:


#disentangle the multiplexed size class X pft dimension
mort_fire_by_pft_and_scls = scpf_to_scls_by_pft(ds.FATES_MORTALITY_FIRE_SZPF, ds)

#sum across size classes to get pft-level mort from fire
mort_fire_by_pft = mort_fire_by_pft_and_scls.sum(axis=2)
                                                
mort_fire_per_capita = mort_fire_by_pft / ds.FATES_NPLANT_PF / months_per_yr

fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7), sharey = True)

for p,ax in zip(range(n_pfts),axes.ravel()):
    mort_fire_per_capita.isel(fates_levpft = p).plot(x = "time",lw = 1, ax = ax, color = pft_colors[p])
    ax.set_ylabel("[N per capita month -1]")
    ax.set_title("{} fire mort.".format(pft_names[p]))

plt.tight_layout()


# #### PFT-level mean

# In[23]:


for p in range(n_pfts):
    print(pft_names[p], ":", mort_fire_per_capita.isel(fates_levpft=p).values.mean())


# ### Fire mort from crown scorch for pft x in size class x

# In[ ]:


#fire mortality from crown scorch by pft/size in number of plants per m2 per year
crown_scorch_pf_sc = scpf_to_scls_by_pft(ds.FATES_MORTALITY_CROWNSCORCH_SZPF, ds)
crown_scorch_pf_sc.isel(fates_levpft = 0).isel(fates_levscls = 7).plot(x = "time")


# ### Mortality from senescence 

# In[ ]:


#SZPF_to_facetted_plot(ds.FATES_MORTALITY_SENESCENCE_SZPF, " ")
#'FATES_AUTORESP_SZPF','FATES_MAINTAR_SZPF','SMP','FATES_CROWNAREA_PF'

#disentangle the multiplexed size class X pft dimension
mort_senescence = scpf_to_scls_by_pft(ds.FATES_MORTALITY_SENESCENCE_SZPF, ds)
#mort_fire

#sum across size classes to get pft-level mort from hyd. failure.
mort_senescence_pf = mort_senescence.sum(axis=2) * m2_per_ha

fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7,7), sharey = True)

for p,ax in zip(range(n_pfts),axes.ravel()):
    mort_senescence_pf.isel(fates_levpft=p).plot(x = "time", color = pft_colors[p],lw = 1, ax = ax)
    ax.set_ylabel("Mort. Senescence [N m-2 yr-1]")
    ax.set_title("Mort. Senescence")

plt.tight_layout()


# ### Shrub mortality from senescence (per capita)

# In[ ]:


n_per_area = ds.FATES_NPLANT_PF.isel(fates_levpft = 3)

mort_sen_per_capita = mort_senescence_pf.isel(fates_levpft=3) / n_per_area
mort_sen_per_capita.plot()

plt.title("Shrub Mortality From Senescence")
plt.ylabel("Mortality rate [per capita]")
plt.xlabel("Year")


# ## Disturbance

# ### Fire

# #### Number of ignitions

# In[ ]:


ignitions_per_ha_per_yr = ds.FATES_IGNITIONS * s_per_yr * m2_per_ha
ignitions_per_ha_per_yr.plot()
plt.ylabel("N ignitions [ha-1 yr-1]")

plt.show()


# #### Burned fraction

# In[18]:


plt.rcParams.update({'axes.titlesize': 'large', 'axes.labelsize':'large'})

burnfrac = ds.FATES_BURNFRAC  * s_per_yr
burnfrac.plot()
plt.ylabel("Annual burned fraction")
plt.xlabel("Year")
plt.show()


# In[19]:


print("mean annual burn frac: ", ds.FATES_BURNFRAC.values.mean() * s_per_yr)


# In[20]:


fri = 1 / (ds.FATES_BURNFRAC.values.mean() * s_per_yr)
print("Mean fire return interval (yrs):",fri)
#','FATES_IGNITIONS','FATES_FIRE_INTENSITY_BURNFRAC',


# ### Area-weighted fire intensity

# In[ ]:


# fireline_intensity = ds.FATES_FIRE_INTENSITY / 1000
# fireline_intensity.plot()
# plt.ylabel("Fire line intensity [kW m-1]")


# In[ ]:


#ds.FATES_FIRE_INTENSITY_BURNFRAC.plot()


# In[21]:


aw_fi = ds.FATES_FIRE_INTENSITY_BURNFRAC / (ds.FATES_BURNFRAC * s_per_day) / 1000
aw_fi.plot()
plt.ylabel("Fire line intensity [kW m-1]")


# In[22]:


print(np.nanmean(aw_fi.values))


# In[ ]:


#area_weighted_mean_fireline_intensity = (ds.FATES_FIRE_INTENSITY_BURNFRAC / ds.FATES_BURNFRAC) / 1000
area_weighted_mean_fireline_intensity.dims


# #### Fire return interval

# #### Nesterov Index

# In[ ]:


#ds.FATES_NESTEROV_INDEX.plot()


# ## Physical Environment

# ### Fuel Load

# In[ ]:


age_by_fuel = agefuel_to_age_by_fuel(ds.FATES_FUEL_AMOUNT_APFC,ds)
fates_fuel_amount_by_class = age_by_fuel.sum(axis = 2)

plt.rcParams.update({'axes.titlesize': 'large', 'axes.labelsize':'large'})

fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(10,10), sharey = False)

fuel_class_names = ['twig','small branch','large branch','trunk','dead leaves','live grass']

for p,ax in zip(range(6),axes.ravel()):
    fates_fuel_amount_by_class.isel(fates_levfuel=p).sel(time = slice("1901-01-01","1999-01-01")).plot(x = "time",
                                                                                                       lw = 1, ax = ax)
    ax.set_ylabel("Fuel amount [Kg C m-2]")
    ax.set_title(fuel_class_names[p])
    
plt.tight_layout()


# ### Fuel Bulk Density

# In[ ]:


ds.FATES_FUEL_BULKD.plot()


# ### Fuel SAV

# In[ ]:


ds.FATES_FUEL_SAV.plot()


# In[ ]:


x = fates_fuel_amount_by_class.isel(fates_levfuel = 3).values
plt.plot(range(len(x)),x)


# ### Soil moisture

# In[ ]:


smp_at_half_meter_mpa = ds.SMP.isel(levgrnd = 7) / 1e5
smp_at_half_meter_mpa.plot()
plt.show()


# ### Light

# #### Direct PAR profile for timestep x

# In[ ]:


ds.FATES_PARPROF_DIR_CLLL.sel(time="1910-07-01").plot(marker = "o")


# #### Diffuse PAR profile for timestep x

# In[ ]:


ds.FATES_PARPROF_DIF_CLLL.sel(time="1910-07-01").plot(marker = "o")


# #### Shrub PAR absorbtion per unit of LAI

# In[ ]:


shrubLAI_shade = ds.FATES_LAISHA_Z_CLLLPF.sel(fates_levcnlfpf = slice(180,240)).sum(axis = 1)
shrubPARZ_shade = ds.FATES_PARSHA_Z_CLLLPF.sel(fates_levcnlfpf = slice(180,240)).sum(axis = 1)
shrubPARZ_shade_per_m2_leaf_area = shrubPARZ_shade / shrubLAI_shade

shrubLAI_sun = ds.FATES_LAISUN_Z_CLLLPF.sel(fates_levcnlfpf = slice(180,240)).sum(axis = 1)
shrubPARZ_sun = ds.FATES_PARSUN_Z_CLLLPF.sel(fates_levcnlfpf = slice(180,240)).sum(axis = 1)
shrubPARZ_sun_per_m2_leaf_area = shrubPARZ_sun / shrubLAI_sun

shrubPARZ_per_m2_leaf_area = shrubPARZ_shade_per_m2_leaf_area + shrubPARZ_sun_per_m2_leaf_area

shrubPARZ_per_m2_leaf_area.plot(marker = "o")

plt.ylabel("PAR absorbed by shrubs [W m-2 leaf area]")


# #### Total PAR absorbed by shrubs

# In[ ]:


PAR_Z_total = shrubPARZ_shade + shrubPARZ_sun
PAR_Z_total.plot()
plt.ylabel("Shrub PAR absorbed [W m-2]")


# In[ ]:




