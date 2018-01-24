"""
*Script takes cases of early spring synoptic variables and creates figures from LENS future members*
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import control_H5_climatologies as C

version = '002'

directory = '/volumes/data/gcm/cesm-lens/BRCP85C5CNBDRD/Aday/Z500/CESM1-CAM5/%s/' % version
path = 'b.e11.BRCP85C5CNBDRD.f09_g16.%s.cam.h1.Z500.20060101-20801231.nc' % version
filename = directory + path
values = Dataset(filename)
longitude = values.variables['lon'][58:240]
latitude = values.variables['lat'][112:167]
heights = values.variables['Z500'][:,112:167,58:240]
values.close()

lons,lats = np.meshgrid(longitude,latitude)

### Import Climo
aveh,lat,lon = C.climoMarch()

heights = np.reshape(heights,(75,365,55,182))

month = []
anom = []
for i in xrange(len(heights)):
    yrsh = heights[i,59:89,:,:]
    anoms = yrsh - aveh
    month.append(yrsh)
    anom.append(anoms)
month = np.asarray(month)
anom = np.asarray(anom)

doy = list(xrange(21,30))
anom = anom[66,:,:,:]
heights = heights[66,:,:,:]

fig = plt.figure()    
fig.suptitle('CESM-LE Member 2 March 500mb Heights and Height Anomalies',fontsize=16)
for i in xrange(len(doy)):
    ax = plt.subplot(3,3,i+1)
 
    m = Basemap(projection='merc',llcrnrlon=183,llcrnrlat=25,urcrnrlon=297,
                urcrnrlat=61,resolution='l')           
    m.drawstates()
    m.drawcountries()
    m.drawmapboundary(fill_color = 'white')
    m.drawcoastlines(color='black',linewidth=0.5)
    m.drawlsmask(land_color='grey',ocean_color='w')
    x,y = m(lons,lats)
    
    heightq = heights[doy[i],:,:]  
    
    anomq = anom[doy[i],:,:]
    anomq[np.where(anomq < -400)] = -400
    anomq[np.where(anomq > 400)] = 400  

    cs1 = m.contour(x,y,heightq,12,colors='k')
    cs = m.contourf(x,y,anomq,xrange(-400,401,5))
    
    cs.set_cmap('RdBu_r')
    ax.text(0.1,0.1,'doy=%d' % (doy[i]+1),size='8',horizontalalignment= 'center',
            backgroundcolor='white',verticalalignment= 'center',
            bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
            transform=ax.transAxes)
plt.tight_layout()
fig.subplots_adjust(bottom=0.15)
cbar_ax = fig.add_axes([0.15, 0.08, 0.7, 0.05])
cbar = fig.colorbar(cs, cax=cbar_ax, orientation = 'horizontal',
                    extend='both',extendfrac='auto',ticks=xrange(-400,401,100))
cbar.set_label('Geopotential Meters (gpm)') 
plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/anommember2.png',dpi=300)