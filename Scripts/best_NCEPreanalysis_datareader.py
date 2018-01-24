"""
Data reader for daily and decadal netCDF4 files
"""
from netCDF4 import Dataset

### Change Working Directory
directory = '/volumes/eas-shared/ault/ecrl/spring-indices/data/best/'

class Reader(object):
    """
    Functions read varying daily/decadal temperature data from NCEP Reanalysis II Files
    
    
    Parameters
    ----------
    yearsSD : year as integer
    decade : decade as integer
    
    Returns
    ----------
    latitude : array of latitudes
    longitude : array of longitudes
    month : list of months in data file
    tmin : minimum temperature array
    tmax : maximum temperature array
    """    
    
    def __init__(self,yrsSD,decade):
        self.yrsSD = yrsSD
        self.decade = decade

    def dlytmin(self):
        """
        Import Daily tmin NetCDF4 Datasets
        """
        tminfile = directory + 'BEST_Daily_TMIN_%d.25N_to_85N.nc' % self.yrsSD
        data = Dataset(tminfile)
        latitude = data.variables['latitude'][:]
        longitude = data.variables['longitude'][:]
        month = data.variables['month'][:]
        tmin = data.variables['tmin'][:]
        data.close()
        return latitude,longitude,month,tmin
    
    def dlytmax(self):
        """
        Import Daily tmax NetCDF4 Datasets
        """        
        tmaxfile = directory + 'BEST_Daily_TMAX_%d.25N_to_85N.nc' % self.yrsSD
        data = Dataset(tmaxfile)
        tmax = data.variables['tmax'][:]
        data.close()
        return tmax
        
    def yrtmin(self):
        """
        Import Decade tmin NetCDF4 Datasets
        """
        tminfile = directory + 'Complete_TMIN_Daily_LatLong1_%d.25N_to_85N.nc' % self.decade
        data = Dataset(tminfile)
        latitude = data.variables['latitude'][:]
        longitude = data.variables['longitude'][:]
        days = data.variables['day_of_year'][:]
        units = data.variables['temperature'].units
        dimensions = data.variables['temperature'].dimensions
        tmin = data.variables['temperature'][:]
        data.close()
        return latitude,longitude,days,units,dimensions,tmin 
    
    def yrtmax(self):
        """
        Import Decade tmax NetCDF4 Datasets
        """     
        tmaxfile = directory + 'Complete_TMAX_Daily_LatLong1_%d.25N_to_85N.nc' % self.decade
        data = Dataset(tmaxfile)
        tmax = data.variables['temperature'][:]
        data.close()
        return tmax