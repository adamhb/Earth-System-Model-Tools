import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import pandas as pd
import xarray as xr
import numpy as np
import datetime as datetime
import sys
import myfuncs
from myparams import *

#funcs
def plot_array(xarr,xds,n_pfts,conversion,title,ylabel,output_path,subplots = False):
    
    xarr = xarr * conversion
    nyears = len(np.unique(pd.to_datetime(xarr.time).year))
    nbase = max(nyears // n_ticks, 1)

    original_stdout = sys.stdout

    with open('log.txt', 'w') as f:
        sys.stdout = f # Change the standard output to the file we created.
        print('nyears: ',nyears)
        print('nbase: ',nbase)
        sys.stdout = original_stdout


    #prep data if it has multiplexed dimensions
    if xarr.dims == ('time', 'fates_levscpf'):
        xarr = myfuncs.scpf_to_scls_by_pft(xarr, xds)
        xarr = xarr.sum(axis=2)
    
    if subplots == False:

        fig, ax = plt.subplots(figsize = (7,7))
        
        if xarr.dims == ('time', 'fates_levpft'):

            for p in range(n_pfts):
                xarr.isel(fates_levpft=p).plot(x = "time",
                color = pft_colors[p],lw = line_width,add_legend = True, label = pft_names[p])
                    
            plt.legend()

        if xarr.dims == ('time',):

            xarr.plot(x = "time", lw = line_width)
        

        ax.xaxis.set_major_formatter(DateFormatter('%Y'))
        ax.xaxis.set_major_locator(mdates.YearLocator(base=nbase))
        ax.set_ylabel(ylabel,fontsize = fontsize)
        ax.title.set_text(title)


    if subplots == True:
        
        ncol,nrow = myfuncs.get_n_subplots(n_pfts)
        fig, axes = plt.subplots(ncols=ncol,nrows=nrow,figsize=(7,7))
        
        for p,ax in zip(range(n_pfts),axes.ravel()):
            xarr.isel(fates_levpft = p).plot(ax = ax, x = "time", color = pft_colors[p],
                                             lw = line_width)
            ax.set_title('{}'.format(pft_names[p]))
            ax.set_ylabel(ylabel,fontsize = fontsize)
            ax.xaxis.set_major_formatter(DateFormatter('%Y'))
            ax.xaxis.set_major_locator(mdates.YearLocator(base=nbase))


            if (pft_names[p] == "shrub") & (title == "Stem Density"):
                ax.set_yscale("log")

        plt.tight_layout()
        plt.subplots_adjust(hspace=0.6,wspace=0.6)

    plt.savefig(output_path + "/" + title.replace(" ","") + ".png")

