"""
*Calculates plots for SIx from historical LENS*
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import scipy.stats as sts
from mpl_toolkits.basemap import Basemap

directory = '/volumes/eas-shared/ault/ecrl/spring-indices/data/'

def SIxHistorical():
    """
    Reads in historical LENS SI-x data

    
    Returns
    ----------
    leafmean : array leaf indices (ens x year x lat x lon)
    latmean : array last freeze indices (ens x year x lat x lon)
    lat : array of latitudes
    lon : array of longitudes 
    lstfrz : list last freeze indices 
    """
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/data/'
    
    versions=['002','003','004','005','006','007','008','009','010','011','012','013','014','015','016','017','018','019','020','021','022','023','024','025','026','027','028','029','030']
    
    lstfrz = []
    leaf=[]
    for version in versions:
        years = 'b.e11.B20TRC5CNBDRD.f09_g16.%s.cam.h1.SI-x.1920-2005.nc' % version
        filename = directory + years
        values = Dataset(filename)
        lon = values.variables['lon'][189:240]
        lat = values.variables['lat'][:32]
        lstfrz_index = values.variables['lstfrz_index'][:,:32,189:240]
        leaf_index = values.variables['leaf_index'][:,:32,189:240]
        values.close()
        
        lstfrz.append(lstfrz_index)
        leaf.append(leaf_index)
        
    latmean = np.asarray(lstfrz)
    leafmean = np.asarray(leaf)
    
    return leafmean,latmean,lat,lon,lstfrz

leafmean,latmean,lat,lon,lstfrz = SIxHistorical()

def bestData():
    """
    Reads in BEST SI-x data

    
    Returns
    ----------
    bestlf : array leaf indices (year x lat x lon)
    bestlst : array last freeze indices (year x lat x lon)
    lats1 : array of latitudes
    lons1 : array of longitudes 
    """
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/'
    filename1 = 'lstfrztrends_BEST.nc'
    filename2 = 'leaftrends_BEST.nc'
    files1 = directory + filename1
    files2 = directory + filename2
    
    data1 = Dataset(files1)
    lat1 = data1.variables['latitude'][:]
    lon1 = data1.variables['longitude'][:]
    bestlst = data1.variables['trend'][:]
    data1.close()
    
    data2 = Dataset(files2)
    bestlf = data2.variables['trend'][:]
    data2.close()
    
    lons1,lats1 = np.meshgrid(lon1,lat1)
    
    return bestlst, bestlf, lons1, lats1
    
bestlst,bestlf,lons1,lats1 = bestData()
    
#lonq = lon - 180.  
#    
#lonqn = np.where((lonq < 120) & (lonq > 55))
#lonqn = np.squeeze(lonqn)
#lonr = lonq[lonqn]
#latqn = np.where((lat > 25) & (lat < 56))
#latqn = np.squeeze(latqn)
#latr = lat[latqn]    

lons,lats = np.meshgrid(lon,lat)

lstave = sts.nanmean(lstfrz)
lstfrz_totalmean = np.nanmean(latmean)

### Mean for each lat/lon
#latmean = []
#for i in xrange(lstfrz.shape[0]):
#    latmeans = sts.nanmean(lstfrz[i,:,:,:])
#    latmean.append(latmeans)
#latmean = np.asarray(latmean)

trend = []
lines = []
latdif = []
for i in xrange(latmean.shape[0]):
    timeq = np.asarray(list(xrange(latmean.shape[1])))
    for j in xrange(latmean.shape[2]):
        for k in xrange(latmean.shape[3]):
            lstslope, lstintercept, r_value, p_value, std_err = sts.linregress(timeq,latmean[i,:,j,k])
            lst_line = lstslope*timeq+lstintercept
            trend.append(lstslope)
            lines.append(lst_line)
            
            latdifq = lst_line-lstfrz_totalmean
            latdif.append(latdifq)
#latdif = np.reshape(np.asarray(latdif),(latmean.shape[0],latmean.shape[2],latmean.shape[3]))
trend = np.reshape(np.asarray(trend),(latmean.shape[0],latmean.shape[2],latmean.shape[3]))
lines = np.asarray(lines)

### Mean for each member
#meanlstfrz = []
#for i in xrange(len(lstave)):
#    yrslstfrz = np.nanmean(lstave[i,:,:])
#    meanlstfrz.append(yrslstfrz)
#meanlstfrz = np.asarray(meanlstfrz)
#
### Save Files and Import Leaf
##np.savetxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/1920-2005lstfrz.txt',(meanlstfrz),delimiter=',')
##
#hist_lf = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/1920-2005lf.txt',delimiter=',')
#ave_lstfrz = sts.nanmean(hist_lf)
#
### Calculate Damage Index
#avedamage = meanlstfrz - ave_lstfrz

### Standard Deviation and Mean of LSTFRZ
std1920 = np.nanstd(trend)
mean1920 = np.nanmean(trend)

### Damage Index 
damage_members = leafmean - latmean
stddamage1920 = np.nanstd(damage_members)
meandamage1920 = np.nanmean(damage_members)

damagez = sts.zscore(damage_members)

damagezscores = np.empty(trend.shape)
for i in xrange(latmean.shape[0]):
    for j in xrange(latmean.shape[2]):
        for k in xrange(latmean.shape[3]):
            damagezscores[i,j,k] = np.mean(damagez[i,:,j,k])
            
trenddamage = []
linesdamage = []
for i in xrange(latmean.shape[0]):
    timeq = np.asarray(list(xrange(latmean.shape[1])))
    for j in xrange(latmean.shape[2]):
        for k in xrange(latmean.shape[3]):
            damageslope, damageintercept, r_value, p_value, std_err = sts.linregress(timeq,damage_members[i,:,j,k])
            damage_line = damageslope*timeq+damageintercept
            trenddamage.append(damageslope)
            linesdamage.append(damage_line)
damagetrend = np.reshape(np.asarray(trenddamage),(latmean.shape[0],latmean.shape[2],latmean.shape[3]))
trenddamage = damagetrend*10.
linesdamage = np.asarray(linesdamage)

meandamagetrend=sts.nanmean(damagetrend)
meantrend = sts.nanmean(trend)
    
### Draw Polygon
def plot_rec(bmap, lonmin,lonmax,latmin,latmax):
    xs = [lonmin,lonmax,lonmax,lonmin,lonmin]
    ys = [latmin,latmin,latmax,latmax,latmin]
    bmap.plot(xs, ys, latlon = True, color='k',linewidth=1.5,linestyle='solid')
lonmin = -101.5
lonmax =  -75.5
latmin =  37.5
latmax =  50.5

member = list(xrange(1,30))
### Plot Trends
fig = plt.figure()   
ax1 = plt.subplot(6,5,1)
m = Basemap(projection='merc',llcrnrlon=-124.7,llcrnrlat=30,urcrnrlon=-62,
            urcrnrlat=54,resolution='l')           
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color='black',linewidth=0.5)
m.drawlsmask(land_color='grey',ocean_color='w')
x,y = m(lons1,lats1)

bestlst[np.where(bestlst > 3)] = 3
bestlst[np.where(bestlst< -3)] = -3

cs = m.contourf(x,y,bestlst,np.arange(-3.,3.1,.1))
plot_rec(m,lonmin,lonmax,latmin,latmax)
cs.set_cmap('bwr_r')

ax1.spines['top'].set_linewidth(3)
ax1.spines['right'].set_linewidth(3)
ax1.spines['bottom'].set_linewidth(3)
ax1.spines['left'].set_linewidth(3)

ax1.text(0.155,0.015,'BEST Data',size='8',horizontalalignment= 'center',
        backgroundcolor='white',verticalalignment= 'center',
        bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
        transform=ax1.transAxes) 
for i in xrange(len(trend)):
    ax = plt.subplot(6,5,i+2)
    m = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=31,urcrnrlon=298,
                urcrnrlat=54,resolution='l')           
    m.drawstates()
    m.drawcountries()
    m.drawmapboundary(fill_color = 'white')
    m.drawcoastlines(color='black',linewidth=0.5)
    m.drawlsmask(land_color='grey',ocean_color='w')
    x,y = m(lons,lats)

    trend[np.where(trend > 3)] = 3
    trend[np.where(trend < -3)] = -3
        
    cs = m.contourf(x,y,trend[i,:,:]*10.,np.arange(-3.,3.1,.1))
    cs.set_cmap('bwr_r')

    ax.text(0.16,0.015,'Member %i' % (member[i]+1),size='8',horizontalalignment= 'center',
            backgroundcolor='white',verticalalignment= 'center',
            bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
            transform=ax.transAxes)    
plt.tight_layout()
fig.subplots_adjust(bottom=0.098)
cbar_ax = fig.add_axes([0.15, 0.08, 0.7, 0.01])
cbar = fig.colorbar(cs, cax=cbar_ax, orientation = 'horizontal',
                    extend='both',extendfrac='auto',ticks=np.arange(-6.,7,1))
cbar.set_label('days/decade') 
figure_title = 'LENS 1920-2005, LSTFRZ Index Trends'
fig.text(0.5, .97, figure_title,
     horizontalalignment='center',
     fontsize=14)
plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/ensemble_lstfrztrend.eps',dpi=400,format='eps')

#meandamagetrend = sts.nanmean(damagetrend)
#
#meantrend[np.where(meantrend > 2)] = 2
#meantrend[np.where(meantrend < -2)] = -2
##
#m = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=31,urcrnrlon=298,
#                urcrnrlat=54,resolution='l')           
#m.drawstates()
#m.drawcountries()
#m.drawmapboundary(fill_color = 'white')
#m.drawcoastlines(color='black')
#m.drawlsmask(land_color='grey',ocean_color='w')
#x,y = m(lons,lats)
#cs = m.contourf(x,y,meantrend,np.arange(-2,2.1,.1))
#cs.set_cmap('bwr_r')
#cbar = m.colorbar(cs,location='bottom',pad='5%')
#cbar.set_label('days/decade')
#cbar.set_ticks(np.arange(-2.,3,1))
#plt.title('LENS, 1920-2005 Mean Leaf Index Trends')
#plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/meandamage_2005.png',dpi=300)

###############################################################################
### Calculate largest damage index
#damage = leafmean - latmean
#
#damage_index = np.empty((damage.shape[0],damage.shape[1]))
#lstfrz_index = np.empty((damage.shape[0],damage.shape[1]))
#for i in xrange(damage.shape[0]):
#    for j in xrange(damage.shape[1]):
#        damage_index[i,j] = np.nanmean(damage[i,j,:,:])
#        lstfrz_index[i,j] = np.nanmean(latmean[i,j,:,:])
#        
#hist_damage = -13.864080 # mean damage over historical LENS
#best_2std = -20. # 2std for BEST
#
#damagetimes = np.where(damage_index <= -20)
#damagevalues = np.empty((len(damagetimes[0])))
#lstfrzvalues = np.empty((len(damagetimes[0])))
#for i in xrange(len(damagetimes[0])):
#        damagevalues[i] = damage_index[damagetimes[0][i],damagetimes[1][i]]
#        lstfrzvalues[i] = lstfrz_index[damagetimes[0][i],damagetimes[1][i]]
#        
#np.savetxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/damagevalues_1920-2005.txt',(damagevalues),delimiter=',')
#np.savetxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/lstfrzvalues_1920-2005.txt',(lstfrzvalues),delimiter=',')
#np.savetxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/damagetimes_1920-2005.txt',(damagetimes),delimiter=',')


###############################################################################
### Create netcdf files of trends
#
#def netcdf(trend, lat, lon):
#    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/'
#    name = 'leaftrends_2005.nc'
#    filename = directory + name
#    ncfile = Dataset(filename,'w',format='NETCDF4')
#    ncfile.description = 'Leaf Index Trends 1970-2005'
#    
#    #dimensions
#    ncfile.createDimension('time',trend.shape[0])
#    ncfile.createDimension('lat',lat.shape[0])
#    ncfile.createDimension('lon',lon.shape[0])
#    
#    #variables
#    times = ncfile.createVariable('time','f4',('time'))
#    latitude = ncfile.createVariable('latitude','f4',('lat'))
#    longitude = ncfile.createVariable('longitude','f4',('lon'))
#    trends = ncfile.createVariable('trend','f4',('time','lat','lon',))
#    
#    #data
#    times[:] = list(xrange(trend.shape[0]))
#    latitude[:] = lat
#    longitude[:] = lon
#    trends[:] = trend
#    
#    ncfile.close()
###############################################################################
