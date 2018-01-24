"""
*Script creates animations for early LENS LSTFRZ cases and produces 2m temperature plots*
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

### Read in LSTFRZ Historical
hdamagevalues = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/damagevalues_1920-2005.txt',delimiter=',')
hlstfrzvalues = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/lstfrzvalues_1920-2005.txt',delimiter=',')
hdamagetimes = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/damagetimes_1920-2005.txt',delimiter=',')

### Read in LSTFRZ Future
fdamagevalues = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/damagevalues_2006-2080.txt',delimiter=',')
flstfrzvalues = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/lstfrzvalues_2006-2080.txt',delimiter=',')
fdamagetimes = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/damagetimes_2006-2080.txt',delimiter=',')

### Slice only 5 cases
hist_damage = hdamagevalues[:5]
hist_lstfrz = hlstfrzvalues[:5]
hist_ens, hist_yr = hdamagetimes[:5,:5]
fut_damage = fdamagevalues[:5]
fut_lstfrz = flstfrzvalues[:5]
fut_ens, fut_yr = fdamagetimes[:5,:5]

### Read in temperature 
def LENStemps(yr, ens, LENS):
    """
    Produces plots for 2m temperatures (animation)
    
    
    Parameters
    ----------
    yr : 20060101-20801231 or 19200101-20051231
    ens : ensemble member(s) (list format)
    LENS : fut or hist
    
    Returns
    ----------
    tmaxf : array of temperatures in Fahrenheit (yr x doy x lat x lon)
    Latitude : array of latitudes
    Longitude : array of longitudes 
    """
    
    if LENS == 'fut':
        tmaxq = []
        for i in ens:
            directory = '/volumes/data/gcm/cesm-lens/BRCP85C5CNBDRD/aday/trefhtmn/cesm1-cam5/%s/' % i
            years = 'b.e11.BRCP85C5CNBDRD.f09_g16.%s.cam.h1.TREFHTMN.%s.NH.nc' % (i,yr)
            filename = directory + years
            values = Dataset(filename)
            longitude = values.variables['lon'][57:240]
            latitude = values.variables['lat'][:44]
            tmax = values.variables['TREFHTMN'][:14965,:44,57:240]
            values.close()
            tmaxq.append(tmax)
        fut_temp = np.asarray(tmaxq)
        tmaxf = np.reshape(fut_temp, (2,41,365,44,183))
        
    elif LENS == 'hist':
        tmaxq = []
        for i in ens:
            directory = '/volumes/data/gcm/cesm-lens/B20TRC5CNBDRD/aday/trefhtmn/cesm1-cam5/%s/' % i 
            years = 'b.e11.B20TRC5CNBDRD.f09_g16.%s.cam.h1.TREFHTMN.%s.NH.nc' % (i,yr)
            filename = directory + years
            values = Dataset(filename)
            longitude = values.variables['lon'][57:240]
            latitude = values.variables['lat'][:44]
            tmax = values.variables['TREFHTMN'][:23725,:44,57:240]
            values.close()
            tmaxq.append(tmax)
        tmax = np.asarray(tmaxq)
        hist_temp = np.squeeze(tmax)
        tmaxf = np.reshape(hist_temp, (65,365,44,183))

    ### Conversion K to F
    tmaxf = ((9./5.)*(tmaxf-273.))+32.
    
    return tmaxf, latitude, longitude

#fut_temp, latitude, longitude = LENStemps('20060101-20801231',['002','003'],'fut')
#hist_temp, latitude, longitude = LENStemps('19200101-20051231',['002'],'hist')

def plotTemp(temps,latitude, longitude, start, stop, yr, year,ens):
    """
    Produces LENS 2m temperature plots that can be used for animations
    
    
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
    
    import cesmcontrol_avet as C
    
    lons, lats = np.meshgrid(longitude,latitude)
    doy = np.arange(start,stop,1)
    time = ['1','2','3','4','5','6','7','8','9','10']
    
    temps = temps[ens,yr,doy,:,:]    
    
    tempclimo,lat,lon = C.climoMarch()  
    slice_anom = doy - 60
    tempclimo = tempclimo[slice_anom,:,:]
    
    anom = temps - tempclimo
    
    anom[np.where(anom<-20)]=-20
    anom[np.where(anom>20)]=20
    
    for i in xrange(len(doy)):
        plt.figure()
        plt.title('LENS Future Year %s, Days %s' % (year,doy[i]))
        m = Basemap(projection='merc',llcrnrlon=235.5,llcrnrlat=26,urcrnrlon=298,
            urcrnrlat=54,resolution='l')           
        m.drawstates()
        m.drawcountries()
        m.drawmapboundary(fill_color = 'white')
        m.drawcoastlines(color='black',linewidth=0.5)
        m.drawlsmask(land_color='grey',ocean_color='w')
        x,y = m(lons,lats)
        cs = m.contourf(x,y,anom[i,:,:],xrange(-20,21,1))
        cs1 = m.contour(x,y,temps[i,:,:],xrange(32,33,1),colors='b',linestyles='dashed',linewidths=2.3)
        cbar = m.colorbar(cs,location='bottom',pad='5%')
        cs.set_cmap('bwr')
        cbar.set_label('degrees Fahrenheit')
        cbar.set_ticks(np.arange(-20,21,5))
        plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/lens_temps_%s.png' % (time[i]), dpi=300)

plotTemp(fut_temp,latitude,longitude,108,118,22,'23',1)
        
        