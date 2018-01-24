"""
Functions Read NCEP Reanalysis II Data in NetCDF4 Files
"""
from netCDF4 import Dataset
import numpy as np

def wind(cmpt,yr):
    """
    Function reads in NCEP Reanalysis wind data
    
    
    Parameters
    ----------
    cmpt : 'u' or 'v' (wind components)
    yr : year as integer
    
    Returns
    ----------
    wind : u or v component of wind in array
    levelw : list of vertical available levels
    latitudew : array of latitudes
    longitudew : array of longitudes
    
    """
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/data/ncep/'
    filename = directory + '%swnd.%d.nc' % (cmpt,yr)
    values = Dataset(filename)
    wind = values.variables['%swnd' % cmpt][:]
    levelw = values.variables['level'][:]
    latitudew = values.variables['lat'][:]
    longitudew = values.variables['lon'][:]
    
    wind = np.squeeze(wind)  
    
    return wind,levelw,latitudew,longitudew

def hgt(yr):
    """
    Function reads in NCEP Reanalysis height data
    
    Parameters
    ----------
    yr : year as integer
    
    Returns
    ----------
    hgts: varying heights array
    levelh : list of vertical available levels
    latitudeh : array of latitudes
    longitudeh : array of longitudes
    """
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/data/ncep/'
    filename = directory + 'hgt.%d.nc' % yr
    values = Dataset(filename)
    hgts = values.variables['hgt'][:]
    levelh = values.variables['level'][:]
    latitudeh = values.variables['lat'][:]
    longitudeh = values.variables['lon'][:]   
    
    return hgts,levelh,latitudeh,longitudeh
    
def temp(yr):
    """
    Function reads in NCEP Reanalysis temperature data
    
    
    Parameters
    ----------
    yr : year as integer
    
    Returns
    ----------
    temp : surface temperature array
    levelt : list of vertical available levels
    latitudet : array of latitudes
    longitudet : array of longitudes
    """
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/data/ncep/'
    filename = directory + 'air.%d.nc' % yr
    values = Dataset(filename)
    temp = values.variables['air'][:]
    levelt = values.variables['level'][:]
    latitudet = values.variables['lat'][:]
    longitudet = values.variables['lon'][:]   
    
    return temp,levelt,latitudet,longitudet
    
def climo(var):
    """
    Function reads in NCEP climate data for daily heights
        
    
    Parameters
    ----------
    var : variable as string for climatologies
    mghts : array of height fields
    levelmh : available levels for analysis
    latitudemh : array of latitudes
    longitudemh : array of longitudes
    """
    directory = '/volumes/data/zml5/NCEP_21st_Reanalysis/'
    filename = directory + '%s.day.1981-2010.ltm.nc' % var
    values = Dataset(filename)
    mhgts = values.variables['hgt'][:]
    levelmh = values.variables['hgt'][:]
    latitudemh = values.variables['lat'][:]
    longitudemh = values.variables['lon'][:]
    
    return mhgts,levelmh,latitudemh,longitudemh
    
    
def lftx(yr):
    """
    Function reads in NCEP Reanalysis Lifted Index Surface data
    
    
    Parameters
    ----------
    yr : year as integer
    
    Returns
    ----------
    lftx : array of lifted indices
    latitudel : array of latitudes
    longitudel : array of longitudes
    """
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/data/ncep/'
    filename = directory + 'lftx.sfc.%d.nc' %yr
    values = Dataset(filename)
    lftx = values.variables['lftx'][:]
    latitudeL = values.variables['lat'][:]
    longitudeL = values.variables['lon'][:]
    
    return lftx,latitudeL,longitudeL    
    
def MSLP(yr):
    """
    Function reads in NCEP Sea Level Pressure
    
    
        
    Parameters
    ----------
    yr : year as integer
    
    Returns
    ----------
    slp12 : array of sea level pressures
    latitudep : array of latitudes
    longitudep : array of longitudes
    """
    directory = '/volumes/data/zml5/NCEP_21st_Reanalysis/'
    filename = directory + 'slp.%d.nc' %yr
    values = Dataset(filename)
    slp12 = values.variables['slp'][:]
    latitudep = values.variables['lat'][:]
    longitudep = values.variables['lon'][:]
    
    return slp12,latitudep,longitudep 
    

def wind20(cmpt,yr):
    """
    Function reads in 20th century NCEP Reanalysis wind data
    
    
    Parameters
    ----------
    cmpt : 'u' or 'v' (wind components)
    yr : year as integer
    
    Returns
    ----------
    wind : u or v component of wind in array
    levelw : list of vertical available levels
    latitudew : array of latitudes
    longitudew : array of longitudes
    """
    directory = '/volumes/data/zml5/NCEP_20th_Reanalysis/'
    filename = directory + '%swnd.%d.nc' % (cmpt,yr)
    values = Dataset(filename)
    wind = values.variables['%swnd' % cmpt][:]
    levelw = values.variables['level'][:]
    latitudew = values.variables['lat'][:]
    longitudew = values.variables['lon'][:]
    
    wind = np.squeeze(wind)  
    
    return wind,levelw,latitudew,longitudew

def hgt20(yr):
    """
    Function reads in 20th century NCEP Reanalysis height data
    
    
    Parameters
    ----------
    yr : year as integer
    
    Returns
    ----------
    hgts: varying heights array
    levelh : list of vertical available levels
    latitudeh : array of latitudes
    longitudeh : array of longitudes
    """

    directory ='/volumes/data/zml5/NCEP_20th_Reanalysis/'
    filename = directory + 'hgt.%d.nc' % yr
    values = Dataset(filename)
    hgts = values.variables['hgt'][:]
    levelh = values.variables['level'][:]
    latitudeh = values.variables['lat'][:]
    longitudeh = values.variables['lon'][:]   
    
    return hgts,levelh,latitudeh,longitudeh
    
def temp20(yr):
    """
    Function reads in 20th century NCEP Reanalysis temperature data
        
    
    Parameters
    ----------
    yr : year as integer
    
    Returns
    ----------
    temp : surface temperature array
    levelt : list of vertical available levels
    latitudet : array of latitudes
    longitudet : array of longitudes
    """
    directory ='/volumes/data/zml5/NCEP_20th_Reanalysis/'
    filename = directory + 'air.%d.nc' % yr
    values = Dataset(filename)
    temp = values.variables['air'][:]
    levelt = values.variables['level'][:]
    latitudet = values.variables['lat'][:]
    longitudet = values.variables['lon'][:]   
    
    return temp,levelt,latitudet,longitudet
    
