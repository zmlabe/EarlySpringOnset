"""
*Script analyzes BEST data and produces z-scores, anomalies, and trends*
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.stats import zscore, linregress
from mpl_toolkits.basemap import Basemap

directory= '/volumes/eas-shared/ault/ecrl/spring-indices/data/SI-x/'

### Extract BEST SIx Data
filename = directory + 'BEST-Cat_SI-X_Daily_LatLong1_1880-to-2013.25N_to_85N.nc'
values = Dataset(filename)
lon = values.variables['lon'][:]
lat = values.variables['lat'][:]
lfs = values.variables['leaf_index'][:,:,:]
bls = values.variables['bloom_index'][:,:,:]
lstfrz = values.variables['lstfrz_index'][:,:,:]
values.close()

lon = lon+360.
lons,lats = np.meshgrid(lon,lat)

### Average SI per year
lfmean = []
blmean = []
for i in xrange(len(lfs)):
    leafs = np.nanmean(lfs[i,:,:])
    blooms = np.nanmean(bls[i,:,:])
    lfmean.append(leafs)
    blmean.append(blooms)
    
### Climatology for SIx
lf_totalmean = np.nanmean(lfmean)
bl_totalmean = np.nanmean(blmean) 
   
### Anom per year
lf_anom = []
bl_anom = []
for j in xrange(len(lfmean)):
    anomlf = lf_totalmean - lfmean[j]
    anombl = bl_totalmean - blmean[j]
    lf_anom.append(anomlf)
    bl_anom.append(anombl)
    
### Damge Index trends
damage = lfs - lstfrz    
timeq = np.asarray(list(xrange(lstfrz.shape[0])))
damagetrends = np.empty((lstfrz.shape[1],lstfrz.shape[2]))
for i in xrange(lstfrz.shape[1]):
    for j in xrange(lstfrz.shape[2]):
        damageslope, damageintercept, r_value, p_value, std_err = linregress(timeq,damage[:,i,j])
        damagetrends[i,j] = damageslope*10.

### z-score per year
lf_zscore = zscore(lfmean)
bl_zscore = zscore(blmean)

yrlf = np.where(lf_zscore<=-2)
yrbl = np.where(bl_zscore<=-2)

damagetrends[np.where(damagetrends > 2)] = 2
damagetrends[np.where(damagetrends < -2)] = -2

### Plot Damage Trends
m = Basemap(projection='merc',llcrnrlon=235.5,llcrnrlat=26,urcrnrlon=300,
            urcrnrlat=54,resolution='l')           
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color='black',linewidth=0.5)
m.drawlsmask(land_color='grey',ocean_color='w')
x,y = m(lons,lats)
cs = m.contourf(x,y,damagetrends,np.arange(-2,2.1,.1))
cs.set_cmap('bwr_r')
cbar = m.colorbar(cs,location='bottom',pad='5%')
cbar.set_label('days/decade')
cbar.set_ticks(np.arange(-2,3,1))
plt.title('Historical Damage Trends')