"""
*Calculates average temperature per 100 years in CESM control*
"""
import numpy as np
from netCDF4 import Dataset
from scipy.stats import nanmean

years = '13000101-13991231'

### Import Tmax
def avetemps(years):
    """
    Calculates average 2m maxT for a 100 year period
    
    
    Parameters
    ----------
    years : 100 year period from CESM control
    
    Returns
    ----------
    tq : average 2m max temperature (time x lat x lon)
    lat : array of latitudes
    lon : array of longitudes
    """
    directory = '/volumes/data/gcm/cesm-lens/B1850C5CN/Aday/trefhtmx/cesm1-cam5/005/'
    filename = 'b.e11.B1850C5CN.f09_g16.005.cam.h1.TREFHTMN.%s.NH.nc' % years
    data = directory + filename
    value = Dataset(data)
    lat = value.variables['lat'][:44]
    lon = value.variables['lon'][57:240]
    tmax1 = value.variables['TREFHTMN'][:,:44,57:240]
    value.close()
    
    # Convert K to F
    tmax2 = (9/5.*(tmax1-273.))+32.
    
    ### Calculate climatologiesv - very slow
    doy = list(xrange(61))
    
    avetemp = []
    for j in xrange(59,36500,365):
        start = j
        end = j+61
        q = np.array(xrange(start,end))
        temps = tmax2[q,:,:]
        avetemp.append(temps)
    century_averages = np.asarray(avetemp)
    
    tq = []
    for m in xrange(len(doy)):
        a = century_averages[:,m,:,:]
        b = nanmean(a)
        tq.append(b)
    tq = np.asarray(tq)
    
    return tq, lat, lon

def netcdf(tq, lat, lon, years):
    """
    Creates netcdf file of temperature climatologies
    
    
    Parameters
    ----------
    tq : average 500mb height (time x lat x lon)
    lat : array of latitudes
    lon : array of longitudes
    years : 100 year period from CESM control
    """
    directory = '/volumes/data/zml5/cesmclimo/'
    name = 'marchaprilclimotempstmin.%s.nc' % years
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = 'March and April Climatology from CESM Control'
    
    #dimensions
    ncfile.createDimension('time',None)
    ncfile.createDimension('lat',lat.shape[0])
    ncfile.createDimension('lon',lon.shape[0])
    
    #variables
    times = ncfile.createVariable('time','f4',('time'))
    latitude = ncfile.createVariable('latitude','f4',('lat'))
    longitude = ncfile.createVariable('longitude','f4',('lon'))
    temps = ncfile.createVariable('temps','f4',('time','lat','lon',))
    
    #data
    times[:] = list(xrange(tq.shape[0]))
    latitude[:] = lat
    longitude[:] = lon
    temps[:] = tq
    
    ncfile.close()
    
def climoMarch():
    """
    Function reads in netcdf files for temperature climatologies
    
    Returns
    ----------
    temps : average 500mb height (time x lat x lon)
    lat : array of latitudes
    lon : array of longitudes
    """
    years = '13000101-13991231'
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/'
    name = 'marchaprilclimotempstmin.%s.nc' % years
    filename = directory + name
    data = Dataset(filename,'r')
    temps = data.variables['temps'][:]
    lat = data.variables['latitude'][:]
    lon = data.variables['longitude'][:]
    
    return temps, lat, lon
    
