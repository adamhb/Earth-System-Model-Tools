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


def getNBase(xarr):
    nyears = len(np.unique(pd.to_datetime(xarr.time).year))
    nbase = max(nyears // n_ticks, 1)
    return(nbase)

#funcs
def plot_array(xarr,xds,n_pfts,conversion,title,ylabel,output_path,subplots = False, getData = False, dbh_min = None):
    
    xarr = xarr * conversion
    nbase = getNBase(xarr)

    #prep data if it has multiplexed dimensions
    if xarr.dims == ('time', 'fates_levscpf'):
        xarr = myfuncs.scpf_to_scls_by_pft(xarr, xds)
        if dbh_min != None:
            xarr = xarr.sel(fates_levscls = slice(dbh_min,None))
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
            ax.xaxis.set_major_locator(mdates.YearLocator(base=nbase*2))


            if (pft_names[p] == "shrub") & (title == "Stem Density"):
                ax.set_yscale("log")

        plt.tight_layout()
        plt.subplots_adjust(hspace=0.6,wspace=0.6)
        fig.suptitle(title, fontsize=fontsize)

    plt.savefig(output_path + "/" + title.replace(" ","") + ".pdf")
    plt.savefig(output_path + "/" + title.replace(" ","") + ".png")
    
    if getData == True:
        return(xarr)
    
def plot_smp(xarr,xds,output_path,depths = [4,7,8],return_means = False):
    title = "SMP"
    nbase = getNBase(xarr)
 
    fig, ax = plt.subplots(figsize = (7,7))

    mean_smp = []
    for d in depths:
        smp = xarr.isel(levgrnd = d).isel(time = slice(12,-1)) * MPa_per_mmh2o
        smp.plot(label = xds.levgrnd.values[d])
        plt.legend()

    ax.xaxis.set_major_formatter(DateFormatter('%Y'))
    ax.xaxis.set_major_locator(mdates.YearLocator(base=nbase))
    ax.set_ylabel("SMP [MPa]",fontsize = fontsize)
    ax.title.set_text("SMP [MPa]")

    plt.savefig(output_path + "/" + title.replace(" ","") + ".pdf")
    plt.savefig(output_path + "/" + title.replace(" ","") + ".png")

def plot_appf(xarr, xds, n_pfts, sup_title, ylabel, output_path):
    
    xarr = myfuncs.appf_to_ap_by_pft(xarr, xds)

    n_age = len(xds.fates_levage)
    
    ncol,nrow = myfuncs.get_n_subplots(n_age)
    
    nbase = getNBase(xarr) * 2
    
    fig, axes = plt.subplots(ncols=ncol,nrows=nrow,figsize=(12,10))

    for age,ax in zip(range(n_age),axes.ravel()):
    
         cca = xarr.isel(fates_levage = age) / xds.FATES_PATCHAREA_AP.isel(fates_levage = age)
    
         for p in range(n_pfts):
             cca.isel(fates_levpft=p).plot(x = "time",
                      color = pft_colors[p],lw = int(line_width/0.75),add_legend = True,
                      label = pft_names[p], ax = ax)
    
             #plt.legend()
         ax.set_title('{} yr old patches'.format(xds.fates_levage.values[age]))
         ax.set_ylabel(ylabel,fontsize = int(fontsize * 0.75))
         ax.xaxis.set_major_formatter(DateFormatter('%Y'))
         ax.xaxis.set_major_locator(mdates.YearLocator(base=nbase))
    
    plt.tight_layout()
    plt.subplots_adjust(hspace=1,wspace=0.2)
    fig.suptitle(sup_title, fontsize=fontsize,y=0.99)

    plt.savefig(output_path + "/" + sup_title.replace(" ","") + ".pdf")
    plt.savefig(output_path + "/" + sup_title.replace(" ","") + ".png")

def plot_ap(xarr, xds, sup_title, ylabel, output_path):
    nbase = getNBase(xarr) * 2
    n_age = len(xds.fates_levage)
    ncol,nrow = myfuncs.get_n_subplots(n_age)
    fig, axes = plt.subplots(ncols=ncol,nrows=nrow,figsize=(10,10))

    for age,ax in zip(range(n_age),axes.ravel()):
        xds.FATES_PATCHAREA_AP.isel(fates_levage = age).plot(ax = ax)
        ax.set_title('{} yr old patches'.format(xds.fates_levage.values[age]))
        ax.set_ylabel(ylabel,fontsize = int(fontsize * 0.75))
        ax.xaxis.set_major_formatter(DateFormatter('%Y'))
        ax.xaxis.set_major_locator(mdates.YearLocator(base=nbase))

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.6,wspace=0.6)
    fig.suptitle(sup_title, fontsize=fontsize)

    plt.savefig(output_path + "/" + sup_title.replace(" ","") + ".pdf")
    plt.savefig(output_path + "/" + sup_title.replace(" ","") + ".png")
