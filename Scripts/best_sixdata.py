"""
*Script creates line plot for BEST SI-x and calculates z-scores/trends*
"""
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import scipy.stats as sts
from scipy.stats import zscore

directory= '/volumes/eas-shared/ault/ecrl/spring-indices/data/SI-x/'

### 2012
filename = directory + 'BEST-Cat_SI-X_Daily_LatLong1_1880-to-2013.25N_to_85N.nc'
values = Dataset(filename)
longitude = values.variables['lon'][78:105]
latitude = values.variables['lat'][12:26]
lfs = values.variables['leaf_index'][40:,12:26,78:105]
bls = values.variables['bloom_index'][40:,12:26,78:105]
lstfrz = values.variables['lstfrz_index'][40:,12:26,78:105]
values.close()

### Restrict domain to United States (Great Lakes)
lonq = np.where((longitude > -102) & (longitude < -75))
lonq = np.squeeze(lonq)
lon = longitude[lonq]
latq = np.where((latitude > 37) & (latitude < 51))
latq = np.squeeze(latq)
lat = latitude[latq]

### Calculate trend in first leaf and last freeze
trendlf = []
trendlst = []
timeq = list(xrange(lfs.shape[0]))
for i in xrange(lfs.shape[1]):
    for j in xrange(lfs.shape[2]):
        lfslope, lfintercept, r_value, p_value, std_err = sts.linregress(timeq,lfs[:,i,j])
        latslope, latintercept, r_value, p_value, std_err = sts.linregress(timeq,lstfrz[:,i,j])
        trendlf.append(lfslope*10.)
        trendlst.append(latslope*10.)
trendlf = np.reshape(np.asarray(trendlf),(lfs.shape[1],lfs.shape[2]))
trendlst = np.reshape(np.asarray(trendlst),(lfs.shape[1],lfs.shape[2]))

lstfrz_ave = []
for i in xrange(len(lstfrz)):
    ave = np.nanmean(lstfrz[i,:,:])
    lstfrz_ave.append(ave)
lstfrz_ave = np.asarray(lstfrz_ave)

lon = lon+360.
lons,lats = np.meshgrid(lon,lat)

### Calculate mean leaf and bloom 
mean_lfs = []
mean_bls = []
for i in xrange(len(lfs)):
    lfsq = np.nanmean(lfs[i,:,:])
    blsq = np.nanmean(bls[i,:,:])
    mean_lfs.append(lfsq)
    mean_bls.append(blsq)
mean_lfs = np.asarray(mean_lfs)
mean_bls = np.asarray(mean_bls)
    
totalmean_lf = np.nanmean(mean_lfs)
totalmean_bl = np.nanmean(mean_bls)    
    
timeq = np.asarray(list(xrange(len(lfs))))
lfslope, lfintercept, r_value, p_value, std_err = sts.linregress(timeq,mean_lfs)
lf_line = lfslope*timeq+lfintercept
blslope, blintercept, r_value, p_value, std_err = sts.linregress(timeq,lstfrz_ave)
bl_line = blslope*timeq+blintercept

### Plot time series of first leaf and last freeze
fig = plt.figure()
plt.title('First Leaf and Last Freeze Index',fontsize=12)
plt.plot(mean_lfs,color='DarkOrange',linewidth=2,label='First Leaf')
plt.plot(lstfrz_ave,color='darkblue',linewidth=2,label='Last Freeze')
plt.xlim([0,len(bl_line)])
x = list(xrange(0,len(bl_line),20))
labels = [1920,1940,1960,1980,2000]
plt.xticks(x,labels)
plt.xlabel('Years')
fig.suptitle('Gridded BEST 1920-2013 SI-x Indices',fontsize=20)
fig.text(0.04,0.5,'DOY',va='center',rotation='vertical')

### Subtract regression
lf_diff = lf_line-totalmean_lf
bl_diff = bl_line-totalmean_bl

lf_notrend = mean_lfs-lf_diff
bl_notrend = mean_bls-bl_diff

plt.plot(bl_line,linewidth=3,color='black')
plt.plot(lf_line,linewidth=3,color='black',label='Trends')
plt.legend(loc=1,prop={'size':9},shadow=True)
plt.grid(True)
#
#### New z scores
#lf_zscore = zscore(mean_lfs)
#bl_zscore = zscore(mean_bls)
#yrlf = np.where(lf_zscore<=-2)
#yrbl = np.where(bl_zscore<=-2)
#
#damage = lfs - bls
#
#damages = np.empty((damage.shape[0]))
#for i in xrange(damages.shape[0]):
#    damages[i] = np.nanmean(damage[i,:,:])
#damage_zscore = zscore(damages)
#
#lf_zscore2 = zscore(lf_notrend)
#bl_zscore2 = zscore(bl_notrend)
#yrlf2 = np.where(lf_zscore2<=-2)
#yrbl2 = np.where(bl_zscore2<=-2)
#
#lf12z = lf_zscore[-2]
#bl12z = bl_zscore[-2]
#
#lf12z_adj = lf_zscore2[-2]
#bl12z_adj = bl_zscore2[-2]
#
#print (lf12z,lf12z_adj), 'Leaf Z Score 2012'
#print (bl12z,bl12z_adj), 'Bloom Z Score 2012'
plt.savefig('/Users/zlabe/documents/CESMspring/Fig1.eps',dpi=400,format='eps')
#
#### Output to text file
#np.savetxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/bestlf.txt',(mean_lfs),delimiter=',')
#np.savetxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/bestbl.txt',(mean_bls),delimiter=',')
#np.savetxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/bestlstfrz.txt',(lstfrz_ave),delimiter=',')

###############################################################################
### Create netcdf files of trends
def netcdf(trend, lat, lon):
    """
    Creates netcdf files for BEST SI-x trends
    
    
    Parameters
    ----------
    trend : SI-x trend array
    lat : array of latitudes
    lon : array of longitudes 

    """
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/'
    name = 'lstfrztrends_BEST.nc'
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = 'Lstfrz Index Trends BEST 1880-2013'
    
    #dimensions
    ncfile.createDimension('lat',lat.shape[0])
    ncfile.createDimension('lon',lon.shape[0])
    
    #variables
    latitude = ncfile.createVariable('latitude','f4',('lat'))
    longitude = ncfile.createVariable('longitude','f4',('lon'))
    trends = ncfile.createVariable('trend','f4',('lat','lon',))
    
    #data
    latitude[:] = lat
    longitude[:] = lon
    trends[:] = trend
    
    ncfile.close()
#netcdf(trendlst,lat,lon)
###############################################################################
