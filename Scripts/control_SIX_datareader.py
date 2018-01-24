"""
Functions Read in NetCDF4 data files for CCSM Climate Model
"""
from netCDF4 import Dataset

### Changing working directory
directory = '/volumes/eas-shared/ault/ecrl/spring-indices/data/'

def ccsm(years):
    """
    Function reads CESM Climate Model
    
    
    Parameters
    ----------
    years : 100 year time slice from CESM control (as string)
    
    Returns
    ----------
    lon : array of longitudes
    lat : array of latitudes 
    leaf_index : array of first leaf indices
    bloom_index : array of bloom indices
    lstfrz_index : array of last freezes
    """
    filename = directory + years
    values = Dataset(filename)
    lon = values.variables['lon'][:]
    lat = values.variables['lat'][:]
    leaf_index = values.variables['leaf_index'][:,:,:]
    bloom_index = values.variables['bloom_index'][:,:,:]
    lstfrz_index = values.variables['lstfrz_index'][:,:,:]
    values.close()
    
    return lon, lat, leaf_index, bloom_index, lstfrz_index