"""
*Script gathers minimum temperatures and plots anomalies over USA for selected dates in future LENS*
"""
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.stats import nanmean
import control_temps_climatologies as C
from mpl_toolkits.basemap import Basemap

### Import Climo
#temps,lat,lon = C.climoMarch()
#lons,lats = np.meshgrid(lon,lat)
#
#
#versions=['002']
#
#damagevalues = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/damagevalues_2006-2080.txt',delimiter=',')
#lstfrzvalues = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/lstfrzvalues_2006-2080.txt',delimiter=',')
#damagetimes = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/damagetimes_2006-2080.txt',delimiter=',')
#
#ensemble, yr = damagetimes[:,:12]
#lstfrzvalues = lstfrzvalues[:12]
#damagevalues = damagevalues[:12]

tmaxq = []
for version in versions:
    directory = '/volumes/data/gcm/cesm-lens/BRCP85C5CNBDRD/aday/trefhtmn/cesm1-cam5/%s/' % version
    years = 'b.e11.BRCP85C5CNBDRD.f09_g16.%s.cam.h1.TREFHTMN.20060101-20801231.NH.nc' % version
    filename = directory + years
    values = Dataset(filename)
    longitude = values.variables['lon'][57:240]
    latitude = values.variables['lat'][:44]
    tmax = values.variables['TREFHTMN'][:,:44,57:240]
    tmaxq.append(tmax)
    values.close()
tmaxq = np.asarray(tmaxq)
tmaxq = np.squeeze(tmaxq)

### Conversion K to F
tmaxf = ((9./5.)*(tmaxq-273.))+32.
tmaxf = np.reshape(tmaxf,(75,365,44,183)) 

### Gather March Anomalous Temps
tmarch = []
for i in xrange(tmaxf.shape[0]):
    ts = tmaxf[i,59:120,:,:]
    anom = ts - temps
    tmarch.append(anom)
tmarch_anom = np.asarray(tmarch)

day = 29
doy = day+60
year = 69
Tave = tmarch_anom[year,day,:,:]
tmaxf = tmaxf[year,doy,:,:]

tmaxf[np.where(tmaxf<0)]=0
tmaxf[np.where(tmaxf>70)]=70
Tave[np.where(Tave<-20)]=-20
Tave[np.where(Tave>20)]=20
#Thigh = np.ones((Tave.shape[1],Tave.shape[2]))
fig = plt.figure()
fig.add_subplot(121)
m = Basemap(projection='merc',llcrnrlon=235.5,llcrnrlat=26,urcrnrlon=298,
            urcrnrlat=54,resolution='l')           
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color='black',linewidth=0.5)
m.drawlsmask(land_color='grey',ocean_color='w')
x,y = m(lons,lats)
cs = m.contourf(x,y,tmaxf,xrange(0,72,2))
cs1 = m.contour(x,y,tmaxf,xrange(32,33,1),colors='b',linestyles='dashed',linewidths=2.3)
cs.set_cmap('rainbow')
cbar = m.colorbar(cs,location='bottom',pad='5%')
cbar.set_label('degrees Fahrenheit')
cbar.set_ticks(np.arange(0,100,10))

fig.add_subplot(122)
m = Basemap(projection='merc',llcrnrlon=235.5,llcrnrlat=26,urcrnrlon=298,
            urcrnrlat=54,resolution='l')           
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color='black',linewidth=0.5)
m.drawlsmask(land_color='grey',ocean_color='w')
x,y = m(lons,lats)
cs = m.contourf(x,y,Tave,xrange(-20,21,1))
cs1 = m.contour(x,y,tmaxf,xrange(32,33,1),colors='b',linestyles='dashed',linewidths=2.3)
cs.set_cmap('bwr')
cbar = m.colorbar(cs,location='bottom',pad='5%')
cbar.set_label('degrees Fahrenheit')
cbar.set_ticks(np.arange(-20,21,5))
fig.suptitle('LENS Year %s Day %s, 2m Temperatures and Anomalies' % (year+1,doy+1),fontsize=18)

#doy = list(xrange(21,30))
#### Plot Temperature Anomalies  
#fig = plt.figure()   
#fig.suptitle('LENS, 2m Temperature Anomaly',fontsize=16)
#for i in xrange(len(doy)):
#    ax = plt.subplot(3,3,i+1)
#    m = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=26,urcrnrlon=298,
#                urcrnrlat=54,resolution='l')           
#    m.drawstates()
#    m.drawcountries()
#    m.drawmapboundary(fill_color = 'white')
#    m.drawcoastlines(color='black',linewidth=0.5)
#    m.drawlsmask(land_color='grey',ocean_color='w')
#    x,y = m(lons,lats)
#    
#    Taves = Tave[doy[i],:,:]
#    highs = np.where(Taves > 25)
#    for j in xrange(len(highs[0])):
#        values = Taves[highs[0][j],highs[1][j]]
#        Thigh[highs[0][j],highs[1][j]] = values
#
#    Taves[np.where(Taves > 30)] = 30
#    Taves[np.where(Taves < -30)] = -30
#    
##    if np.nanmax(Thigh) != 1.0:
##        cs1 = m.contour(x,y,Thigh,colors='r',linewidth=1.5)
#        
#    cs = m.contourf(x,y,Taves,xrange(-30,31,1))
#    cs.set_cmap('RdBu_r')
#
#    ax.text(0.1,0.1,'doy=%i' % (doy[i]+1),size='8',horizontalalignment= 'center',
#            backgroundcolor='white',verticalalignment= 'center',
#            bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
#            transform=ax.transAxes)    
#plt.tight_layout()
#fig.subplots_adjust(bottom=0.15)
#cbar_ax = fig.add_axes([0.15, 0.08, 0.7, 0.05])
#cbar = fig.colorbar(cs, cax=cbar_ax, orientation = 'horizontal',
#                    extend='both',extendfrac='auto',ticks=xrange(-30,31,5))
#cbar.set_label('degrees Fahrenheit') 

### Temperature Trends
#trendq = tmaxf[:,59,:,:]
#trendn = []
#for i in xrange(trendq.shape[0]):
#    trends = np.mean(trendq[i,:,:])
#    trendn.append(trends)
#trendn = np.asarray(trendn)
    
 
