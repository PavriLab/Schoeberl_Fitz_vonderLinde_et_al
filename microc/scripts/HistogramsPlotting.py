import numpy as np
from .NonNumToNan import NonNumToNan
import matplotlib.pyplot as plt
import pandas as pd

def Plot2DHist(
    x, 
    y, 
    ax=None, 
    bins=100, 
    xlabel="", 
    ylabel="", 
    colorbar=False, 
    hist2dKwargs={}
):
    """
        Plot 2D histogram (ignoring Nans and infs)
        
        PARAMETERS
        ----------
        x: 1D array
        y: 1D array
        ax: matplotlib axis object / None
        bins: int, default 100
            bins for histogram
        xlabel: str, default: ""
            label for x-axis
        ylabel: str, default: ""
            label for y-axis
        colorbar: bool, default False
            plot colorbar
        hist2dKwargs: dict, default {}
            keyword arguments for 2d histogram
            
        OUTPUT
        ------        
        matplotlib axis object
    """
    if ax is None:
        ax = plt.gca()
    xmin = np.nanmin(NonNumToNan(x))
    xmax = np.nanmax(NonNumToNan(x))
    ymin = np.nanmin(NonNumToNan(y))
    ymax = np.nanmax(NonNumToNan(y))
    c = ax.hist2d(x,y,bins=bins,range=((xmin,xmax),(ymin,ymax)),**hist2dKwargs)
    # c = ax.hexbin(x,y,bins=bins,**hist2dKwargs)
    if colorbar:
        plt.colorbar(c[3],ax=ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return ax