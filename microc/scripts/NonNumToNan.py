import numpy as np

def NonNumToNan(x, changeTo='nan'):
    """
        Convert all non-numbers to Nans
        
        PARAMETERS
        ----------
        x: numpy array
        changeTo: int/'nan', default 'nan'
        
        OUTPUT
        ------
        numpy array
    """
    if changeTo == 'nan':
        changeTo = np.nan
    return np.nan_to_num(x, 
                         posinf=changeTo, 
                         neginf=changeTo, 
                         nan=changeTo
                        )