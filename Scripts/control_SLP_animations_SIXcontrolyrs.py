"""
*Function reads SLP and ENSO data from CESM control. 
Years listed below are for leaf/bloom March 2012-CESM correlations. 
Plots can create animations for SLP data*
"""

from control_SLP_datareader import SLP
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from scipy.stats import nanmean

### Read in SLP files
years = ['04020101-04991231.nc','05000101-05991231.nc','06000101-06991231.nc',
         '07000101-07991231.nc','08000101-08991231.nc','09000101-09991231.nc',
         '10000101-10991231.nc','11000101-11991231.nc','12000101-12991231.nc',
         '13000101-13991231.nc','13000101-13991231.nc','15000101-15991231.nc',
         '16000101-16991231.nc','17000101-17991231.nc','18000101-18991231.nc',
         '19000101-19991231.nc','20000101-20991231.nc','21000101-21991231.nc',
         '21000101-22001231.nc']

rankyrs = years[2]
SLPfilename = 'b.e11.B1850C5CN.f09_g16.005.cam.h1.PSL.' + rankyrs
directory = '/volumes/zml5/scripts/'
filename = directory + SLPfilename

values = Dataset(filename)
date = values.variables['date'][:]
lon = values.variables['lon'][58:240]
lat = values.variables['lat'][112:167]
SLP = values.variables['PSL'][:,112:167,58:240] 
values.close()

### Years for CESM (402-999)
blrank9 = [444,251,16,4,204,556,214,156,364]
lfrankn = [4,251,333,156,64,401,12,44,183]

#### BL Year 1 (yr 846 at slice 444) 
#date1 = 846-800
#start1 = (date1 * 365) - 1
##time = [16789:17154] 
#
#### BL Year 2 (yr 653 at slice 251)
#date2 = 653-600
#start2 = (date2 * 365) - 1
##time = [19344:19709] 
#
#### BL Year 3 (yr 418 at slice 16)
#date3 = 418-400
#start3 = (date3 * 365) - 1
##time = [6569:6934] 
#
#### BL Year 4 (yr 406 at slice 4)
#date4 = 406-400
#start4 = (date4 * 365) - 1
##time = [2189:2554] 
#
#### BL Year 5 (yr 606 at slice 204)
#date5 = 606-600
#start5 = (date5 * 365) - 1
##time = [2189:2554] 
#
#### BL Year 6 (yr 958 at slice 556)
#date6 = 958-900
#start6 = (date6 * 365) - 1
##time = [21169:21534] 
#
#### BL Year 7 (yr 616 at slice 214)
#date7 = 616-600
#start7 = (date7 * 365) - 1
##time = [5839:6204] 
#
#### BL Year 8 (yr 558 at slice 156)
#date8 = 558-500
#start8 = (date8 * 365) - 1
##time = [21169:21534] 
#
#### BL Year 9 (yr 766 at slice 364)
#date9 = 766-700
#start9 = (date9 * 365) - 1
##time = [24089:24454] 

### LF Year 1 (yr 406 at slice 4) 
date1 = 406-400
start1 = (date1 * 365) - 1
#time = [2189:2554] 

### LF Year 2 (yr 653 at slice 251)
date2 = 653-600
start2 = (date2 * 365) - 1
#time = [19344:19709] 

### LF Year 3 (yr 735 at slice 333)
date3 = 735-700
start3 = (date3 * 365) - 1
#time = [12774:13139] 

### LF Year 4 (yr 558 at slice 156)
date4 = 558-500
start4 = (date4 * 365) - 1
#time = [21169:21534] 

### LF Year 5 (yr 466 at slice 64)
date5 = 466-400
start5 = (date5 * 365) - 1
#time = [24089:24454] 

### LF Year 6 (yr 803 at slice 401)
date6 = 803-800
start6 = (date6 * 365) - 1
#time = [1094:1459] 

### LF Year 7 (yr 414 at slice 12)
date7 = 414-400
start7 = (date7 * 365) - 1
#time = [5109:5474] 

### LF Year 8 (yr 446 at slice 44)
date8 = 446-400
start8 = (date8 * 365) - 1
#time = [16789:17154] 

### LF Year 9 (yr 585 at slice 183)
date9 = 585-500
start9 = (date9 * 365) - 1
#time = [31024:31389] 

### Conversion Pa to hPa
PSL = SLP/100.
PSL[np.where(PSL < 970.)] = 970.
PSL[np.where(PSL > 1040.)] = 1040.

### Plot Sea Level Pressure
doy = list(xrange(60,76))
lon,lat = np.meshgrid(lon,lat)

fig = plt.figure()

### Plot Mean Sea Level Pressure
slpq = []
for doy in xrange(74,89):
    slpn = PSL[doy,:,:]
    slpq.append(slpn)
slpq = np.asarray(slpq)
slp_mean = nanmean(slpq)

m = Basemap(projection='merc',llcrnrlon=183,llcrnrlat=25,urcrnrlon=297,
            urcrnrlat=61,resolution='l')
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color='black',linewidth=0.5)
m.drawlsmask(land_color='grey',ocean_color='w')
x,y = m(lon,lat)
cs = m.contour(x,y,slp_mean,11,colors='k')
cs1 = m.contourf(x,y,slp_mean,np.arange(970,1041,1))
cbar = m.colorbar(cs1,location='bottom',pad='5%',ticks=[np.arange(970,1041,5)])
cbar.set_label('Pressure (hPa)')
plt.title('CESM Year 418, Days 75-90, Sea Level Pressure',fontsize=20)
directory = '/volumes/zml5/research/CCSM/results/'
plt.savefig(directory + 'bl.rank3.meanslp.png',dpi=300)
plt.show()

#for K in xrange(len(doy)):
#    fig.suptitle('CESM Year 418, Sea Level Pressure',fontsize=20)
#    ax = plt.subplot(4,4,K+1)
#    
#    m = Basemap(projection='merc',llcrnrlon=130,llcrnrlat=16,urcrnrlon=299,
#            urcrnrlat=64,resolution='l')        
#    m.drawmapboundary(fill_color = 'white')
#    m.drawcoastlines(color='black',linewidth=0.5)
#    m.drawlsmask(land_color='grey',ocean_color='w')
#    x,y = m(lon,lat) 
#    
#    P = PSL[doy[K],:,:]    
#    cs1 = m.contour(x,y,P,11,colors='k')
#    cs = m.contourf(x,y,P,np.arange(970,1041,1))
#    
#    ax.text(0.1,0.1,'doy=%d' % doy[K],size='5',horizontalalignment= 'center',
#            backgroundcolor='white',verticalalignment= 'center',
#            bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
#            transform=ax.transAxes)
#
#plt.tight_layout()
#fig.subplots_adjust(bottom=0.15)
#cbar_ax = fig.add_axes([0.15, 0.08, 0.7, 0.05])
#cbar = fig.colorbar(cs, cax=cbar_ax, orientation = 'horizontal',
#                    extend='both',extendfrac='auto',
#                    ticks=list(xrange(970,1041,10)))
#cbar.set_label('Pressure (hPa)') 
#directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/'
#plt.savefig(directory + 'rank3.bl.png',dpi=300)

#for i in xrange(len(doy)):
#    
#    fig = plt.figure()   
#    m = Basemap(projection='merc',llcrnrlon=130,llcrnrlat=16,urcrnrlon=299,
#            urcrnrlat=64,resolution='l')         
#    m.drawstates()
#    m.drawcountries()
#    m.drawmapboundary(fill_color = 'white')
#    m.drawcoastlines(color='black',linewidth=0.5)
#    m.drawlsmask(land_color='grey',ocean_color='w')
#    x,y = m(lon,lat)
#    
#    P = PSL[doy[i],:,:]
#    
#    cs1 = m.contour(x,y,P,30,colors='k')
#    cs = m.contourf(x,y,P,np.arange(980,1041,1))
#    cbar = m.colorbar(cs,location='bottom',pad='5%',extend='both',extendfrac='auto')
#    cbar.set_label('Pressure (hPa)')
#    plt.title('CESM Year 653, Sea Level Pressure (doy=%d)' % doy[i])
    
#    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/'
#    if i < 10:
#        plt.savefig(directory + 'rank2.bl.00%d.png' % i,dpi=300)
#    else:
#        plt.savefig(directory + 'rank2.bl.0%d.png' % i,dpi=300)
#    fig.clear()    