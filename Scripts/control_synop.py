"""
*Script reads CESM-control meteorological parameters for case studies*
"""
from netCDF4 import Dataset
import numpy as np
from scipy.stats import nanmean
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

### Import Windu at 200 hPa
### Year 653
directory1 = '/volumes/data/zml5/'
filename1 = 'b.e11.B1850C5CN.f09_g16.005.cam.h1.U200.06000101-06991231.nc'
data1 = directory1 + filename1
values1 = Dataset(data1)
lon = values1.variables['lon'][58:240]
lat = values1.variables['lat'][112:167]
windu = values1.variables['U200'][19344:19709,112:167,58:240]
values1.close()

lons,lats = np.meshgrid(lon,lat)

### Import Windv at 500 hPa and Tmax
filename2 = 'b.e11.B1850C5CN.f09_g16.005.cam.h1.V200.06000101-06991231.nc'
data2 = directory1 +filename2
values2 = Dataset(data2)
windv = values2.variables['V200'][19344:19709,112:167,58:240]
values2.close()

### Dates
#doy = list(xrange(65,74))
#fig = plt.figure()    
#fig.suptitle('Year 653, CESM Control 200mb Winds',fontsize=16)
#for i in xrange(len(doy)):
#    ax = plt.subplot(3,3,i+1)
#    windu1 = windu[doy[i],:,:]
#    windv1 = windv[doy[i],:,:]
#    
#    ### Wind Speed Plot
#    speed = np.sqrt(windu1**2+windv1**2)
#    speed[np.where(speed <25)] = np.nan
#    #  
#    m = Basemap(projection='merc',llcrnrlon=183,llcrnrlat=25,urcrnrlon=297,
#                urcrnrlat=61,resolution='l')           
#    m.drawstates()
#    m.drawcountries()
#    m.drawmapboundary(fill_color = 'white')
#    m.drawcoastlines(color='black',linewidth=0.5)
#    m.drawlsmask(land_color='grey',ocean_color='w')
#    x,y = m(lons,lats)
#    cs1 = m.contourf(x,y,speed,50,cmap='jet')
#    cs =    m.quiver(x[::3,::3],y[::3,::3],windu1[::3,::3],windv1[::3,::3],scale=450) 
#    
#    ax.text(0.1,0.1,'doy=%d' % doy[i],size='8',horizontalalignment= 'center',
#            backgroundcolor='white',verticalalignment= 'center',
#            bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
#            transform=ax.transAxes)
#plt.tight_layout()
#fig.subplots_adjust(bottom=0.15)
#cbar_ax = fig.add_axes([0.15, 0.08, 0.7, 0.05])
#cbar = fig.colorbar(cs1, cax=cbar_ax, orientation = 'horizontal',
#                    extend='both',extendfrac='auto',ticks=(
#                    xrange(25,61,5)))
#cbar.set_label('Knots') 
#plt.savefig('/volumes/zml5/research/ccsm/results/CESM_control_T/year653_wind200.png',dpi=300)
