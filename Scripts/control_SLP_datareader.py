"""
Functions Read in NetCDF4 data files for CESM Climate Model
"""

from netCDF4 import Dataset

def ocean(oceanfilename):
    """
    Function reads CESM data for the oceans
    
    
    Parameters
    ----------
    oceanfilename : file (string)
    
    Returns
    ----------
    lon : array of longitudes
    lat : array of latitudes
    nino34 : time series array of nino3.4 anoms
    nino12 : time series array of nino1+2 anoms
    nino3 : time series array of nino3 anoms
    nino4 : time series array of nino4 anoms
    pdo : time series array of PDO indices
    """
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Data/'
    filename = directory + oceanfilename
    values = Dataset(filename)
    lon = values.variables['lon'][:]
    lat = values.variables['lat'][:]
    nino34 = values.variables['nino34'][:]
    nino12 = values.variables['nino12'][:]
    nino3 = values.variables['nino3'][:]
    nino4 = values.variables['nino4'][:]
    pdo = values.variables['pdo_timeseries_mon'][:]
    values.close()

    return lon, lat, nino34, nino12, nino3, nino4, pdo   

#SLPfilename = 'b.e11.B1850C5CN.f09_g16.005.cam.h1.PSL.04020101-04991231.nc'
    
def SLP(SLPfilename):
    """
    Function reads CESM data for SLP (CESM-LE control)
    
    Parameters
    ----------
    filename : file (string)
    
    Returns
    ----------
    lon : array of longitudes
    lat : array of latitudes
    date : list of available dates
    slp : sea level pressure in Pa
    """
    directory = '/volumes/zml5/scripts/'
    filename = directory + SLPfilename
    values = Dataset(filename)
    date = values.variables['date'][:]
    lon = values.variables['lon'][185:244]
    lat = values.variables['lat'][123:155]
    SLP = values.variables['PSL'][:,123:155,185:244] 
    values.close()
    
    return lon, lat, date, SLP