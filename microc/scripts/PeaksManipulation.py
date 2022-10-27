import numpy as np
import pandas as pd
import bioframe as bf

def GetIntervalsFromPeaks(
    peaksDF, 
    flank=200,
    cols=["chrom","start","end"]
):
    """
        Get irregular peaks into regular sized interval centered at the peaks 
        centers
        
        PARAMETERS
        ----------
            peaksDF: pandas DataFrame
                peaks dataframe
            flank: int, default 200
                flank for each side of the peak center
            cols: list, default ["chrom","start","end"]
                column names
        OUTPUT
        ------
    """
    peaksCopy = peaksDF.copy()
    peakCenter = (peaksCopy[cols[1]] + peaksCopy[cols[2]])//2
    peaksCopy.start = peakCenter - flank
    peaksCopy.end = peakCenter + flank
    return peaksCopy

def GetTSS(df, cols=["chrom","start","end","strand"], dropDups=True):
    """
        Get TSS from a dataframe with strand column
        
        PARAMETERS
        ----------
        df: pandas dataframe
            dataframe with bedlike format
        cols: list, default, ["chrom","start","end","strand"]
            column names of the dataset

        OUTPUT
        ------
        pandas dataframe
    """
    pos = []
    for i in np.array(df[cols]):
        if i[3] == '-':
            pos.append(i[2])
        elif i[3] == '+':
            pos.append(i[1])
        else:
            pos.append(-1)
    pos = np.array(pos)
    tssDF = df.copy()
    tssDF["start"] = pos - 50
    tssDF["end"] = pos + 50
    if dropDups:
        return tssDF.drop_duplicates(subset=cols[:3]).reset_index(drop=True)
    else:
        return tssDF

def GetTTS(df, cols=["chrom","start","end","strand"], dropDups=True):
    """
        Get TTS from a dataframe with strand column

        PARAMETERS
        ----------
        df: pandas dataframe
            dataframe with bedlike format
        cols: list, default, ["chrom","start","end","strand"]
            column names of the dataset
        
        OUTPUT
        ------
        pandas dataframe

    """
    pos = []
    for i in np.array(df[cols]):
        if i[3] == '-':
            pos.append(i[1])
        elif i[3] == '+':
            pos.append(i[2])
        else:
            pos.append(-1)
    pos = np.array(pos)
    ttsDF = df.copy()
    ttsDF["start"] = pos - 50
    ttsDF["end"] = pos + 50
    if dropDups:
        return ttsDF.drop_duplicates(subset=cols[:3]).reset_index(drop=True)
    else:
        return ttsDF