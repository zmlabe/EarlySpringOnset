"""
*Script plots temperatures for March 2012 from NCEP Reanalysis*
"""

import numpy as np
import best_NCEPreanalysis_datareader as R #Reads in temperature data from NCEP Reanalysis II
#import best_NCEPreanalysis_synop_datareader as N
import temp_converter as T
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

### Create String of Years
yrsq = list(xrange(1880,2014))
yrsSD = 2012 # year with >1+SD

### Call Functions
a = R.Reader(yrsSD,None)
latitude,longitude,month,tmin = a.dlytmin()
tmax = a.dlytmax()

### Restrict Domain for Contiguous United States
lonq = np.where((longitude > -125) & (longitude < -55))
lonq = np.squeeze(lonq)
lon = longitude[lonq]
latq = np.where((latitude > 25) & (latitude < 55))
latq = np.squeeze(latq)
lat = latitude[latq]

lon,lat = np.meshgrid(lon,lat)

tmax = tmax[:,latq,55:125]
tmin = tmin[:,latq,55:125]

### Convert Temperature Units
tmin = T.tempunits(tmin,'C','F')
tmax = T.tempunits(tmax,'C','F')

### Extract Particular Dates
tmax = tmax[62:94,:,:]

doy = list(xrange(62,94))

### Basemap Plot Temperature
for i in xrange(len(doy)):
    
    fig = plt.figure()   
    m = Basemap(projection='merc',llcrnrlon=-124,llcrnrlat=26,urcrnrlon=-60,
                urcrnrlat=54,resolution='l')           
    m.drawstates()
    m.drawcountries()
    m.drawmapboundary(fill_color = 'white')
    m.drawcoastlines(color='black',linewidth=1.0)
    m.drawlsmask(land_color='grey',ocean_color='w')
    x,y = m(lon,lat)
    cs = m.contourf(x,y,tmax[i],range(-5,100,2))
    cbar = m.colorbar(cs,location='bottom',pad='5%')
    cbar.set_label('degrees Fahrenheit')
    plt.title('Maximum Surface Temperature 2012 (doy=%d)' % doy[i])
    
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/'
#    if i < 10:
#        plt.savefig(directory + '2012.marchtemp.00%d.png' % i,dpi=100)
#    else:
#        plt.savefig(directory + '2012.marchtemp.0%d.png' % i,dpi=100)
    plt.show()    
    fig.clear()
