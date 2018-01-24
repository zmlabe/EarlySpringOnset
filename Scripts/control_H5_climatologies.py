"""
*Calculates average 500mb height per 100 years in CESM control*
"""
import numpy as np
from netCDF4 import Dataset
from scipy.stats import nanmean

years = '13000101-13991231' #arbitrary period selected

### Import Tmax
def aveH5(years):
    """
    Calculates average 500mb height for a 100 year period
    
    
    Parameters
    ----------
    years : 100 year period from CESM control
    
    Returns
    ----------
    tq : average 500mb height (time x lat x lon)
    lat : array of latitudes
    lon : array of longitudes
    """
    directory = '/volumes/data/gcm/cesm-lens/B1850C5CN/Aday/Z500/CESM1-CAM5/005/'
    filename = 'b.e11.B1850C5CN.f09_g16.005.cam.h1.Z500.%s.nc' % years
    data = directory + filename
    value = Dataset(data)
    lat = value.variables['lat'][112:167]
    lon = value.variables['lon'][58:240]
    z500 = value.variables['Z500'][:,112:167,58:240]
    value.close()
    
    ### Calculate climatologies - very slow
    doy = list(xrange(30))
    
    aveheights = []
    for j in xrange(59,36500,365):
        start = j
        end = j+30
        q = np.array(xrange(start,end))
        Z500 = z500[q,:,:]
        aveheights.append(Z500)
    century_averages = np.asarray(aveheights)
    
    tq = []
    for m in xrange(len(doy)):
        a = century_averages[:,m,:,:]
        b = nanmean(a)
        tq.append(b)
    tq = np.asarray(tq)
    
    return tq, lat, lon
tq, lat, lon = aveH5(years)

def netcdf(tq, lat, lon, years):
    """
    Creates netcdf file of H5 climatologies
    
    
    Parameters
    ----------
    tq : average 500mb height (time x lat x lon)
    lat : array of latitudes
    lon : array of longitudes
    years : 100 year period from CESM control
    """
    
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/'
    name = 'marchclimoZ500.%s.nc' % years
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = 'March Climatology from CESM Control'
    
    #dimensions
    ncfile.createDimension('time',None)
    ncfile.createDimension('lat',lat.shape[0])
    ncfile.createDimension('lon',lon.shape[0])
    
    #variables
    times = ncfile.createVariable('time','f4',('time'))
    latitude = ncfile.createVariable('latitude','f4',('lat'))
    longitude = ncfile.createVariable('longitude','f4',('lon'))
    heights = ncfile.createVariable('H5','f4',('time','lat','lon',))
    
    #data
    times[:] = list(xrange(tq.shape[0]))
    latitude[:] = lat
    longitude[:] = lon
    heights[:] = tq
    
    ncfile.close()
netcdf(tq,lat,lon,years)
    
def climoMarch():
    """
    Function reads in netcdf files for H5 climatologies
    
    Returns
    ----------
    heights : average 500mb height (time x lat x lon)
    lat : array of latitudes
    lon : array of longitudes
    """
    years = '13000101-13991231'
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/'
    name = 'marchclimoZ500.%s.nc' % years
    filename = directory + name
    data = Dataset(filename,'r')
    heights = data.variables['H5'][:]
    lat = data.variables['latitude'][:]
    lon = data.variables['longitude'][:]
    
    return heights, lat, lon
    