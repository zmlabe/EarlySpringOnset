"""
*Script Reads and Plots Leaf,Bloom,Freeze Data from 1880-2013*
"""
from netCDF4 import Dataset
import numpy as np
from scipy.stats import nanmean
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as c

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

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
    bmap.plot(xs, ys, latlon = True, color='darkblue',linewidth=3.5,linestyle='-')
lonmin = -101.5
lonmax =  -75.5
latmin =  37.5
latmax =  50

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

### Make Colormap
cmap = plt.get_cmap('spring')
cmap2 = plt.get_cmap('summer_r')
cmaplist = [cmap(i) for i in range(cmap.N)]
cmaplist2 = [cmap2(i) for i in range(cmap2.N)]

cms=c.ListedColormap(cmaplist+cmaplist2)

fig = plt.figure()
## Panel 1
ax2 = fig.add_subplot(2,1,1)   
m = Basemap(projection='merc',llcrnrlon=-124.7,llcrnrlat=25.5,urcrnrlon=-60,
            urcrnrlat=54,resolution='l')           
m.drawstates()
m.drawcountries()
parallels = np.arange(0,61,5)
m.drawparallels(parallels,labels=[True,False,True,True],linewidth=0.35)
m.drawmapboundary(fill_color = 'whitesmoke')
m.drawcoastlines(color='black',linewidth=0.5)
m.drawlsmask(land_color='grey',ocean_color='whitesmoke')
x,y = m(lon,lat)
plt.title('Average First Leaf Index',fontsize=16)
cs = m.contourf(x,y,aveleaf,range(0,211,8),extend='max')
cbar = m.colorbar(cs,location='right',pad='5%',ticks=xrange(0,211,30),extend='max')
cbar.set_ticklabels(['January','','March','','May','','July',''])
cs.set_cmap(cms)
ax2.text(0.176,0.015,'Climatological',size='15',horizontalalignment= 'center',
                backgroundcolor='white',verticalalignment= 'center',
                bbox=dict(facecolor='white',edgecolor='black',alpha=1),
                transform=ax2.transAxes) 

### Panel 2
ax1 = fig.add_subplot(2,1,2)   
m = Basemap(projection='merc',llcrnrlon=-124.7,llcrnrlat=25.5,urcrnrlon=-60,
            urcrnrlat=54,resolution='l')           
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color = 'whitesmoke')
m.drawcoastlines(color='black',linewidth=0.5)
m.drawlsmask(land_color='grey',ocean_color='whitesmoke')
x,y = m(lon,lat)
parallels = np.arange(0,61,5)
m.drawparallels(parallels,labels=[True,False,True,True],linewidth=0.35)
cs2 = m.contourf(x,y,anomleaf,range(-40,41,5),extend='both')
cbar = m.colorbar(cs2,location='right',pad='5%',
                  ticks=list(xrange(-40,41,10)))
cbar.set_label('Difference (Days)') 
plot_rec(m,lonmin,lonmax,latmin,latmax)
plt.title('First Leaf Index Anomalies',fontsize=16)
cs2.set_cmap('RdYlGn')
ax1.text(0.14,0.015,'March 2012',size='15',horizontalalignment= 'center',
                backgroundcolor='white',verticalalignment= 'center',
                bbox=dict(facecolor='white',edgecolor='black',alpha=1),
                transform=ax1.transAxes)   

plt.subplots_adjust(wspace=0.1)
directory = '/Users/zlabe/desktop/'
plt.savefig(directory + 'spring2012.png',dpi=400)