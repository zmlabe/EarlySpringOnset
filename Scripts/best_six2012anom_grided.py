"""
*Script Reads and Plots Leaf,Bloom,Freeze Data from 1880-2013*
"""
from netCDF4 import Dataset
import numpy as np
from scipy.stats import nanmean
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as c

directory = '/volumes/eas-shared/ault/ecrl/spring-indices/data/SI-x/'

### Read in Netcdf file
filename = directory + 'BEST-Cat_SI-X_Daily_LatLong1_1880-to-2013.25N_to_85N.nc'
values = Dataset(filename)
longitude = values.variables['lon'][:]
latitude = values.variables['lat'][:]
leaf = values.variables['leaf_index'][:]
bloom = values.variables['bloom_index'][:]
lstfrz = values.variables['lstfrz_index'][:]
years = values.variables['time'][:]
values.close()

### Restrict Domain for contiguous United States
lonq = np.where((longitude > -125) & (longitude < -55))
lonq = np.squeeze(lonq)
lon = longitude[lonq]
latq = np.where((latitude > 25) & (latitude < 55))
latq = np.squeeze(latq)
lat = latitude[latq]

lon,lat = np.meshgrid(lon,lat)

leaf = leaf[:,latq,55:125]
bloom = bloom[:,latq,55:125]
lstfrz = lstfrz[:,latq,55:125]
lats = nanmean(lstfrz)

### Extract Particular Years
leaf1 = leaf[-2,:,:]
bloom1 = bloom[-2,:,:]
lstfrz1 = lstfrz[-2,:,:]
damage = leaf1-lstfrz1

### Draw Polygon
def plot_rec(bmap, lonmin,lonmax,latmin,latmax):
    xs = [lonmin,lonmax,lonmax,lonmin,lonmin]
    ys = [latmin,latmin,latmax,latmax,latmin]
    bmap.plot(xs, ys, latlon = True, color='k',linewidth=2.5,linestyle='-')
lonmin = -101.5
lonmax =  -75.5
latmin =  37.5
latmax =  50.5

### Calculate Anomaly
leafq = []
for yr in xrange(len(leaf)):
    leafn = leaf[yr,:,:]
    leafq.append(leafn)
leafq = np.asarray(leafq)

aveleaf = nanmean(leafq)
anomleaf = leaf1 - aveleaf

bloomq = []
for yr in xrange(len(bloom)):
    bloomn = bloom[yr,:,:]
    bloomq.append(bloomn)
bloomq = np.asarray(bloomq)

avebloom = nanmean(bloomq)
anombloom = bloom1 - avebloom

### Basemap Plot Indices
#m = Basemap(projection='merc',llcrnrlon=-124.7,llcrnrlat=25.5,urcrnrlon=-60,urcrnrlat=54,resolution='l') 
#m.drawstates()
#m.drawcountries()
#m.drawmapboundary(fill_color = 'white')
#m.drawcoastlines(color='black',linewidth=0.5)
#m.drawlsmask(land_color='grey',ocean_color='w')
#x,y = m(lon,lat)
##anomleaf[np.where(anomleaf < -40)] = -40
##anomleaf[np.where(anomleaf > 40)] = 40
#cs = m.contourf(x,y,bloom1,range(0,226,5))
#cs.set_cmap('jet_r')
#cbar = m.colorbar(cs,location='bottom',pad='5%')
#cbar.set_label('Day of Year')
#plt.title('March 2012 Bloom Index')
#
#directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/'
#plt.savefig(directory + '2012.bloom.png',dpi=1000)
#plt.show()

### Make Colormap
cmap = plt.get_cmap('spring')
cmap2 = plt.get_cmap('summer_r')
cmaplist = [cmap(i) for i in range(cmap.N)]
cmaplist2 = [cmap2(i) for i in range(cmap2.N)]

cms=c.ListedColormap(cmaplist+cmaplist2)

fig = plt.figure()
## Panel 1
ax1 = fig.add_subplot(2,1,1)   
m = Basemap(projection='merc',llcrnrlon=-124.7,llcrnrlat=25.5,urcrnrlon=-60,
            urcrnrlat=54,resolution='l')           
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color='black',linewidth=0.5)
m.drawlsmask(land_color='grey',ocean_color='w')
x,y = m(lon,lat)
plt.title('Average First Leaf Index',fontsize=16)
cs = m.contourf(x,y,aveleaf,range(0,233,8))
cbar = m.colorbar(cs,location='right',pad='5%',ticks=xrange(0,226,45))
cbar.set_label('DOY') 
cs.set_cmap(cms)

#### Panel 2
#ax2 = fig.add_subplot(212)   
#m1 = Basemap(projection='merc',llcrnrlon=-124.7,llcrnrlat=25.5,urcrnrlon=-60,
#            urcrnrlat=54,resolution='l')           
#m1.drawstates()
#m1.drawcountries()
#m1.drawmapboundary(fill_color = 'white')
#m1.drawcoastlines(color='black',linewidth=0.5)
#m1.drawlsmask(land_color='grey',ocean_color='w')
#x,y = m1(lon,lat)
#cs1 = m1.contourf(x,y,avebloom,range(0,233,8))
#cs1.set_cmap(cms)

#fig.subplots_adjust(bottom=0.15)
#cbar_ax = fig.add_axes([0.264,0.08, 0.5, 0.03])
#cbar = fig.colorbar(cs, cax=cbar_ax, orientation = 'horizontal',
#                    extend='both',extendfrac='auto',
#                    ticks=list(xrange(0,226,45)))
#cbar.set_label('DOY') 

### Calculate Leaf and Bloom Anomalies 2012
#fig = plt.figure()
### Panel 1
ax1 = fig.add_subplot(2,1,2)   
m = Basemap(projection='merc',llcrnrlon=-124.7,llcrnrlat=25.5,urcrnrlon=-60,
            urcrnrlat=54,resolution='l')           
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color='black',linewidth=0.5)
m.drawlsmask(land_color='grey',ocean_color='w')
x,y = m(lon,lat)
cs2 = m.contourf(x,y,anomleaf,range(-40,41,5))
cbar = m.colorbar(cs2,location='right',pad='5%',ticks=list(xrange(-40,41,10)))
cbar.set_label('Difference (Days)') 
plot_rec(m,lonmin,lonmax,latmin,latmax)
plt.title('March 2012, First Leaf Index Anomalies',fontsize=16)
cs2.set_cmap('RdYlGn')
#
#### Panel 2
#ax2 = fig.add_subplot(212)   
#m1 = Basemap(projection='merc',llcrnrlon=-124.7,llcrnrlat=25.5,urcrnrlon=-60,
#            urcrnrlat=54,resolution='l')           
#m1.drawstates()
#m1.drawcountries()
#m1.drawmapboundary(fill_color = 'white')
#m1.drawcoastlines(color='black',linewidth=0.5)
#m1.drawlsmask(land_color='grey',ocean_color='w')
#x,y = m1(lon,lat)
#cs1 = m1.contourf(x,y,anombloom,range(-40,41,5))
#plot_rec(m,lonmin,lonmax,latmin,latmax)
#cs1.set_cmap('RdYlGn')
#
#fig.subplots_adjust(bottom=0.15)
#cbar_ax = fig.add_axes([0.264,0.08, 0.5, 0.03])
#cbar = fig.colorbar(cs, cax=cbar_ax, orientation = 'horizontal',
#                    extend='both',extendfrac='auto',
#                    ticks=list(xrange(-40,41,10)))
#cbar.set_label('Difference (Days)') 
#
#plt.subplots_adjust(wspace=0.1)
#directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/'
##plt.savefig(directory + '2012_leaf_bloom_anom.png',dpi=300)
##plt.show()

plt.subplots_adjust(wspace=0.1)
directory = '/Users/zlabe/documents/CESMspring/'
plt.savefig(directory + 'Fig2.eps',format='eps',dpi=400)