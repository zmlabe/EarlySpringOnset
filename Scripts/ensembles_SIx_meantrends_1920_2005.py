"""
*Reads in Trends for Mean Leaf and Mean LSTFRZ for 1920-2005*
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import nanmean
from mpl_toolkits.basemap import Basemap

directory = '/volumes/data/zml5/cesmclimo/'
filename1 = 'leaftrends_2005.nc'
filename2 = 'lstfrztrends_2005.nc'

files1 = directory + filename1
data1 = Dataset(files1)
trendlf = data1.variables['trend'][:]
lat = data1.variables['latitude'][:]
lon = data1.variables['longitude'][:]
data1.close()

files2 = directory + filename2
data2 = Dataset(files2)
trendlstfrz = data2.variables['trend'][:]
data2.close()

trendlf = trendlf*10.
trendlstfrz = trendlstfrz*10.
lons,lats=np.meshgrid(lon,lat)
meantrendlf = nanmean(trendlf)
meantrendlstfrz = nanmean(trendlstfrz)

### Plot 2x1
def plot_rec(bmap, lonmin,lonmax,latmin,latmax):
    xs = [lonmin,lonmax,lonmax,lonmin,lonmin]
    ys = [latmin,latmin,latmax,latmax,latmin]
    bmap.plot(xs, ys, latlon = True, color='k',linewidth=1.5,linestyle='solid')
lonmin = -101.5
lonmax =  -75.5
latmin =  37.5
latmax =  50.5

fig = plt.figure()   
figure_title = 'LENS 1920-2005, Mean SI-x Trends'
fig.text(0.5, .973, figure_title,
     horizontalalignment='center',
     fontsize=14)
ax1 = plt.subplot(2,1,1)
m = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=31,urcrnrlon=298,
            urcrnrlat=54,resolution='l')           
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color='black',linewidth=0.5)
m.drawlsmask(land_color='grey',ocean_color='w')
x,y = m(lons,lats)

meantrendlf[np.where(meantrendlf > 3)] = 3
meantrendlf[np.where(meantrendlf < -3)] = -3
ax1.text(0.1565,.08,'First Leaf Index',size='12',horizontalalignment= 'center',
        backgroundcolor='white',verticalalignment= 'center',
        bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
        transform=ax1.transAxes) 

cs = m.contourf(x,y,meantrendlf,np.arange(-3.,3.1,.1))
plot_rec(m,lonmin,lonmax,latmin,latmax)
cs.set_cmap('bwr_r')

### 
ax2 = plt.subplot(2,1,2)
m = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=31,urcrnrlon=298,
            urcrnrlat=54,resolution='l')           
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color='black',linewidth=0.5)
m.drawlsmask(land_color='grey',ocean_color='w')
x,y = m(lons,lats)

meantrendlstfrz[np.where(meantrendlstfrz > 3)] = 3
meantrendlstfrz[np.where(meantrendlstfrz < -3)] = -3

cs1 = m.contourf(x,y,meantrendlstfrz,np.arange(-3.,3.1,.1))
plot_rec(m,lonmin,lonmax,latmin,latmax)
cs1.set_cmap('bwr_r')
ax2.text(0.143,.08,'LSTFRZ Index',size='12',horizontalalignment= 'center',
        backgroundcolor='white',verticalalignment= 'center',
        bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
        transform=ax2.transAxes) 

plt.tight_layout()
fig.subplots_adjust(bottom=0.11)
cbar_ax = fig.add_axes([0.15, 0.08, 0.7, 0.01])
cbar = fig.colorbar(cs, cax=cbar_ax, orientation = 'horizontal',
                    extend='both',extendfrac='auto',ticks=np.arange(-3.,4,1))
cbar.set_label('days/decade') 
plt.savefig('/Users/zlabe/documents/CESMspring/Fig8.eps',dpi=400,format='eps')
 

