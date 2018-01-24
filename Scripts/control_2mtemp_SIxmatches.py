"""
*Script reads CESM-control meteorological parameters.*
"""
from netCDF4 import Dataset
import numpy as np
from scipy.stats import nanmean
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

### Import Tmax
directory = '/volumes/data/zml5/'
filename = 'b.e11.B1850C5CN.f09_g16.005.cam.h1.TREFHTMX.06000101-06991231.nc'
data = directory + filename
value = Dataset(data)
lon = value.variables['lon'][58:240]
lat = value.variables['lat'][112:167]
tmax1 = value.variables['TREFHTMX'][:,112:167,58:240]
value.close()

# Convert K to F
tmax2 = (9/5.*(tmax1-273.))+32.
tmax = tmax2[19344:19709,:,:]

### Calculate climatologiesv - very slow
doy = list(xrange(65,74))

avetemp = []
for j in xrange(65,36565,365):
    start = j
    end = j+9
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

### Import Windu at 500 hPa
filename1 = 'b.e11.B1850C5CN.f09_g16.005.cam.h1.U500.06000101-06991231.nc'
data1 = directory + filename1
values1 = Dataset(data1)
windu = values1.variables['U500'][19344:19709,112:167,58:240]
values1.close()

### Import Windv at 500 hPa
filename2 = 'b.e11.B1850C5CN.f09_g16.005.cam.h1.V500.06000101-06991231.nc'
data2 = directory +filename2
values2 = Dataset(data2)
windv = values2.variables['V500'][19344:19709,112:167,58:240]
values2.close()

lons,lats = np.meshgrid(lon,lat)

### Plot Temperature    
fig = plt.figure()   
fig.suptitle('Year 653, CESM Control 2m Temperature Anomaly and 500mb Winds',fontsize=16)
for i in xrange(len(doy)):
    ax = plt.subplot(3,3,i+1)
    m = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=26,urcrnrlon=298,
                urcrnrlat=54,resolution='l')           
    m.drawstates()
    m.drawcountries()
    m.drawmapboundary(fill_color = 'white')
    m.drawcoastlines(color='black',linewidth=0.5)
    m.drawlsmask(land_color='grey',ocean_color='w')
    x,y = m(lons,lats)
    
    T = tmax[doy[i],:,:]
    Tave = T - tq[0]
    Tave[np.where(Tave > 30)] = 30
    Tave[np.where(Tave < -30)] = -30
    
    winduq = windu[doy[i],:,:]
    windvq = windv[doy[i],:,:]
    speed = np.sqrt((winduq**2) + (windvq**2))
    
    cs = m.contourf(x,y,Tave,range(-30,31,1))
    cs.set_cmap('RdBu_r')
    cs2 = m.quiver(x[::2,::2],y[::2,::2],winduq[::2,::2],windvq[::2,::2],scale=400)

    ax.text(0.1,0.1,'doy=%d' % doy[i],size='8',horizontalalignment= 'center',
            backgroundcolor='white',verticalalignment= 'center',
            bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
            transform=ax.transAxes)    
plt.tight_layout()
fig.subplots_adjust(bottom=0.15)
cbar_ax = fig.add_axes([0.15, 0.08, 0.7, 0.05])
cbar = fig.colorbar(cs, cax=cbar_ax, orientation = 'horizontal',
                    extend='both',extendfrac='auto',ticks=(xrange(-30,31,5)))
cbar.set_label('degrees Fahrenheit') 
plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/year653_anoms_temp.png',dpi=300)
