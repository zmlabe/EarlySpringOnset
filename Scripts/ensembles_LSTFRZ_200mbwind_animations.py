"""
*Script creates animations for early LENS LSTFRZ cases and produces 200mb wind plots*
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

### Read in temperature 
def LENStemps(yr, ens, LENS):
    """
    Produces plots for 200mb winds (animation)
    
    
    Parameters
    ----------
    yr : 20060101-20801231 or 19200101-20051231
    ens : ensemble member(s) (list format)
    LENS : fut or hist
    
    Returns
    ----------
    U200 : array of 200mb winds (yr x doy x lat x lon)
    V200 : array of 200mb winds (yr x doy x lat x lon)
    Latitude : array of latitudes
    Longitude : array of longitudes 
    """    
    if LENS == 'fut':
        U200 = []
        V200 = []
        for i in ens:
            directory = '/volumes/data/gcm/cesm-lens/BRCP85C5CNBDRD/aday/U200/cesm1-cam5/%s/' % i          
            years = 'b.e11.BRCP85C5CNBDRD.f09_g16.%s.cam.h1.U200.%s.nc' % (i,yr)
            filename = directory + years
            values = Dataset(filename)
            longitude = values.variables['lon'][57:240]
            latitude = values.variables['lat'][112:167]
            u200 = values.variables['U200'][:14965,112:167,57:240]
            values.close()
            U200.append(u200)

            directory2 = '/volumes/data/gcm/cesm-lens/BRCP85C5CNBDRD/aday/V200/cesm1-cam5/%s/' % i              
            years2 = 'b.e11.BRCP85C5CNBDRD.f09_g16.%s.cam.h1.V200.%s.nc' % (i,yr)
            filename2 = directory2 + years2
            values2 = Dataset(filename2)
            v200 = values2.variables['V200'][:14965,112:167,57:240]
            values2.close()
            V200.append(v200)
            
        U200 = np.asarray(U200)
        V200 = np.asarray(V200)
        U200 = np.reshape(U200, (2,41,365,55,183))
        V200 = np.reshape(V200, (2,41,365,55,183))
        
    elif LENS == 'hist':
        U200 = []
        V200 = []
        for i in ens:
            directory = '/volumes/data/gcm/cesm-lens/B20TRC5CNBDRD/aday/U200/cesm1-cam5/%s/' % i 
            years = 'b.e11.B20TRC5CNBDRD.f09_g16.%s.cam.h1.U200.%s.nc' % (i,yr)
            filename = directory + years
            values = Dataset(filename)
            longitude = values.variables['lon'][57:240]
            latitude = values.variables['lat'][112:167]
            u200 = values.variables['U200'][:23725,112:167,57:240]
            values.close()
            U200.append(u200)
            
            directory2 = '/volumes/data/gcm/cesm-lens/B20TRC5CNBDRD/aday/V200/cesm1-cam5/%s/' % i 
            years2 = 'b.e11.B20TRC5CNBDRD.f09_g16.%s.cam.h1.V200.%s.nc' % (i,yr)
            filename2 = directory2 + years2
            values2 = Dataset(filename2)
            v200 = values2.variables['V200'][:23725,112:167,57:240]
            values2.close()
            V200.append(v200)
            
        U200 = np.asarray(U200)
        V200 = np.asarray(V200)
        U200 = np.reshape(U200, (65,365,55,183))
        V200 = np.reshape(V200, (65,365,55,183))
    
    return U200, V200, latitude, longitude
#U200f, V200f, latitude, longitude = LENStemps('20060101-20801231',['002','003'],'fut')
#U200h, V200h, latitude, longitude = LENStemps('19200101-20051231',['002'],'hist')

def plotWind(U200,V200,latitude, longitude, start, stop, yr, year,ens):
    """
    Produces LENS 200mb wind plots that can be used for animations
    
    
    Parameters
    ----------
    temps : historical or future temperature arrays
    latitude : latitude array 
    longitude : longitude array
    start : doy for starting plot
    stop : doy for ending plot
    yr : yr in LENS for plotting
    year : actual year for figure caption
    ens : ensemble member for figure caption
    
    Attributes
    ----------
    climoMarch : function to read in climatologies for temperatures from CESM
    """   
    lons, lats = np.meshgrid(longitude,latitude)
    doy = np.arange(start,stop,1)
    time = ['1','2','3','4','5','6','7','8','9','10']

    meanU = np.mean(U200)
    
    U200 = U200[yr,doy,:,:]    
    V200 = V200[yr,doy,:,:]
    
    speed = np.sqrt(U200**2 + V200**2)
    
    speeds = speed - meanU
    
    for i in xrange(len(doy)):
        plt.figure()
        plt.title('LENS Historic Year %s, Days %s' % (year,doy[i]))
        m = Basemap(projection='merc',llcrnrlon=235.5,llcrnrlat=26,urcrnrlon=298,
            urcrnrlat=54,resolution='l')           
        m.drawstates()
        m.drawcountries()
        m.drawmapboundary(fill_color = 'white')
        m.drawcoastlines(color='black',linewidth=0.5)
        m.drawlsmask(land_color='grey',ocean_color='w')
        x,y = m(lons,lats)
        cs = m.contourf(x,y,speeds[i,:,:],xrange(0,71,1))
        U200s = U200[i,:,:]
        V200s = V200[i,:,:]
        cs1 = m.quiver(x[::3,::3],y[::3,::3],U200s[::3,::3],V200s[::3,::3],scale=450) 
        cbar = m.colorbar(cs,location='bottom',pad='5%')
        cs.set_cmap('gist_ncar')
        cbar.set_label('knots')
        cbar.set_ticks(np.arange(0,80,10))
        plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/lens_wind_%s.png' % (time[i]), dpi=300)

plotWind(U200h,V200h,latitude,longitude,105,115,64,'65',0)