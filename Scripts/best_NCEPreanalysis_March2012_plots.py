"""
*Script plots March 2012 200mb winds and height fields.
Data for March 2010 is also available for plotting*
"""
import best_NCEPreanalysis_synop_datareader as N #function reads in data from NCEP
import numpy as np
from scipy.stats import nanmean
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, interp

### Call NCEP Functions
u12,level,latitude,longitude = N.wind('u',2012)
v12,level,latitude,longitude = N.wind('v',2012)
hgts12,level,latitude,longitude = N.hgt(2012)
#temp12,level,latitude,longitude = N.temp(2012)
#u10,level1,latitude1,longitude1 = N.wind20('u',1910)
#v10,level1,latitude1,longitude1 = N.wind20('v',1910)
#hgts10,level,latitude1,longitude1 = N.hgt20(1910)
#lftx,latitude,longitude = N.lftx(2012)
#mhgts,levelmh,latitudemh,longitudemh = N.climo('hgt')
#slp12,latitude,longitude = N.MSLP(2012)

lonq = np.where((longitude > 180) & (longitude < 305))
lonq = np.squeeze(lonq)
lon = longitude[lonq]
latq = np.where((latitude > 20) & (latitude < 65))
latq = np.squeeze(latq)
lat = latitude[latq]

#lonq1 = np.where((longitude1 > 180) & (longitude1 < 305))
#lonq1 = np.squeeze(lonq1)
#lon1 = longitude1[lonq1]
#latq1 = np.where((latitude1 > 20) & (latitude1 < 65))
#latq1 = np.squeeze(latq1)
#lat1 = latitude1[latq1]

lons,lats = np.meshgrid(lon,lat)
#lon1,lat1 = np.meshgrid(lon1,lat1)

### Restrict Domain Over United States
u12 = u12[:,9,latq,73:122]
v12 = v12[:,9,latq,73:122]
hgt12 = hgts12[:,9,latq,73:122]
#temp12 = temp12[:,0,latq,73:122]
#u10 = u10[:,16,latq1,91:153]
#v10 = v10[:,16,latq1,91:153]
#hgt10 = hgts10[:,16,latq1,91:153]
#lftx = lftx[:,latq,93:122]
#mhgts12 = mhgts[84,9,latq,73:122]
#mhgts10 = mhgts[84,9,latq,73:122]
#slp12 = slp12[:,latq,73:122]

### Calculate Mean SLP Proceeding 20 days
#slpq = []
#for doy in xrange(71,86):
#    slpn = slp12[doy,:,:]
#    slpq.append(slpn)
#slpq = np.asarray(slpq)
#aveslp = nanmean(slpq)
#slp_mean = aveslp/100.
#
#slp_mean[np.where(slp_mean < 970.)] = 970.
#slp_mean[np.where(slp_mean > 1041.)] = 1041.

###Calculate Mean Winds
uq12 = []
for doy in xrange(83,86):
    un12 = u12[doy,:,:]
    uq12.append(un12)
uq12 = np.asarray(uq12)
aveu12 = nanmean(uq12)

#uq10 = []
#for doy in xrange(76,88):
#    un10 = u10[doy,:,:]
#    uq10.append(un10)
#uq10 = np.asarray(uq10)
#aveu10 = nanmean(uq10)

vq12 = []
for doy in xrange(83,86):
    vn12 = v12[doy,:,:]
    vq12.append(vn12)
vq12 = np.asarray(vq12)
avev12 = nanmean(vq12)

#vq10 = []
#for doy in xrange(76,88):
#    vn10 = v10[doy,:,:]
#    vq10.append(vn10)
#vq10 = np.asarray(vq10)
#avev10 = nanmean(vq10)

### Calculate Mean Geopotential Heights Proceeding 20 Days
hgtq12 = []
for doy in xrange(83,86):
    hgtn12 = hgt12[doy,:,:]
    hgtq12.append(hgtn12)
hgtq12 = np.asarray(hgtq12)
avehgts12 = nanmean(hgtq12)

#hgtq10 = []
#for doy in xrange(76,88):
#    hgtn10 = hgt10[doy,:,:]
#    hgtq10.append(hgtn10)
#hgtq10 = np.asarray(hgtq10)
#avehgts10 = nanmean(hgtq10)

### Calculate Geopotential Height Anomaly 
#hgt12 = hgt12[75,:,:]
#ahgt12 = hgt12 - avehgts12
#
#hgt10 = hgt10[84,:,:]
#ahgt10 = hgt10 - avehgts10

#### Daily Values for a particular level
#u = u[67,0,:,:]
#v = v[67,0,:,:]
#lftx = lftx[67,:,:]

### Basemap Plot SLP for March 1910 and 2012
#m = Basemap(projection='merc',llcrnrlon=183,llcrnrlat=25,urcrnrlon=297,
#            urcrnrlat=61,resolution='l')
#m.drawstates()
#m.drawcountries()
#m.drawmapboundary(fill_color = 'white')
#m.drawcoastlines(color='black',linewidth=0.5)
#m.drawlsmask(land_color='grey',ocean_color='w')
#x,y = m(lon,lat)
#cs = m.contour(x,y,slp_mean,11,colors='k')
#cs1 = m.contourf(x,y,slp_mean,np.arange(970,1041,1))
#cbar = m.colorbar(cs1,location='bottom',pad='5%',ticks=[np.arange(970,1041,5)])
#cbar.set_label('Pressure (hPa)')
#plt.title('10-25 March 2012, Sea Level Pressure',fontsize=20)
#directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/'
#plt.savefig(directory + 'meanslp.2012.png',dpi=200)
#plt.show()

#### Basemap Plot Heights
#m = Basemap(projection='merc',llcrnrlon=235,llcrnrlat=25,urcrnrlon=300,
#            urcrnrlat=54,resolution='l')
#m.drawstates()
#m.drawcountries()
#m.drawmapboundary(fill_color = 'white')
#m.drawcoastlines(color='black',linewidth=0.5)
#m.drawlsmask(land_color='grey',ocean_color='w')
#x,y = m(lon,lat)
#cs = m.contour(x,y,ahgt,15,colors='k')
#cs1 = m.contourf(x,y,ahgt,range(-450,600,2))
#cs = m.barbs(x,y,u,v,15)
#cbar = m.colorbar(cs1,location='bottom',pad='5%')
#cbar.set_label('Meters')
#plt.title('Geopotential Height (250mb) Trend (20 days) March 13, 2012')
#
#directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/'
#plt.savefig(directory + '2012.hgttrend.march.007.png',dpi=300)
#plt.show()

speed12 = np.sqrt(aveu12**2+avev12**2)
#speed10 = np.sqrt(aveu10**2+avev10**2)

speed12[np.where(speed12 <25)] = np.nan
#speed10[np.where(speed10 <25)] = np.nan

time = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20']
days = list(xrange(60,81))

### Create Figure  
#for i in xrange(len(time)):
#    fig = plt.figure() 
#    us12 = u12[days[i]]
#    vs12 = v12[days[i]]
#    speeds12 = speed12[days[i]]
#    speeds12[np.where(speeds12<25)]=25
#    speeds12[np.where(speeds12>55)]=55
#    hgtss12 = hgt12[days[i]]
#    #fig.suptitle('200 mb Daily Mean Winds and Heights',fontsize=16)
#    ### Panel 1
#    #ax1 = fig.add_subplot(211) 
#    m = Basemap(projection='merc',llcrnrlon=183,llcrnrlat=25,urcrnrlon=297,
#                urcrnrlat=61,resolution='l')           
#    m.drawstates()
#    m.drawcountries()
#    m.drawmapboundary(fill_color = 'white')
#    m.drawcoastlines(color='black',linewidth=0.5)
#    m.drawlsmask(land_color='grey',ocean_color='w')
#    x,y = m(lons,lats)
#    cs2 = m.contourf(x,y,speeds12,range(25,56,1))
#    cs1 = m.contour(x,y,hgtss12,20,colors='r',linewidth=1,linestyles='dashed')
#    cs = m.quiver(x[::2,::2],y[::2,::2],us12[::2,::2],vs12[::2,::2],scale=450)
#    cbar = m.colorbar(cs2,location='bottom',pad='5%',ticks=(xrange(25,61,5)))
#    cbar.set_label('Knots')
#    plt.title('March %s, 2012, 200 mb Daily Mean Winds and Heights' % time[i],fontsize=16) 
#    plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/2012winds.%d.png' % i,dpi=300) 
    
### Panels for March 2012
fig=plt.figure()
fig.suptitle('200mb Zonal Mean Wind and Geopotential Height',fontsize=16)
ax1 = fig.add_subplot(211) 
m = Basemap(projection='merc',llcrnrlon=183,llcrnrlat=25,urcrnrlon=297,
            urcrnrlat=61,resolution='l')           
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color='black',linewidth=0.5)
m.drawlsmask(land_color='grey',ocean_color='w')
x,y = m(lons,lats)
cs2 = m.contourf(x,y,speed12,range(25,61,1))
cs2.set_cmap('jet')
cs = m.quiver(x[::2,::2],y[::2,::2],aveu12[::2,::2],avev12[::2,::2],scale=450,color='darkred')
cbar = m.colorbar(cs2,location='right',pad='5%',ticks=list(xrange(25,61,5)))
cbar.set_label('Knots')
plt.title('March 23-26, 2012')  
    
ax1 = fig.add_subplot(212) 
m = Basemap(projection='merc',llcrnrlon=183,llcrnrlat=25,urcrnrlon=297,
            urcrnrlat=61,resolution='l')           
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color='black',linewidth=0.5)
m.drawlsmask(land_color='grey',ocean_color='w')
cs2 = m.contour(x,y,avehgts12,range(11000,12500,100),linestyles='dashed',linewidth=1,colors='k')
cs1 = m.contourf(x,y,avehgts12,range(11000,12500,50))
cs1.set_cmap('jet')
cbar1 = m.colorbar(cs1,location='right',pad='5%',ticks=range(11000,12600,200))
cbar1.set_label('Meters')

cbar.set_label('Knots') 
plt.subplots_adjust(wspace=0.1)

plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/march2012_1910_ncep.eps',dpi=400,format='eps') 