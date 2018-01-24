"""
*Script reads historical and future LENS 500mb heights. Makes thresholds for plotting
based on json date files*
"""

from netCDF4 import Dataset
import numpy as np
import json as J
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import scipy.stats as sts

### Read in Chunks
histstd = open('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/json_hists.txt')
hist_std = histstd.readline()
histstd.close()
futstd = open('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/json_futs.txt')
fut_std = futstd.readline()
futstd.close()
stdf = J.loads(fut_std)
stdh = J.loads(hist_std)

histyr = stdh['year']
histdoy = stdh['doy']
futyr = stdf['year']
futdoy = stdf['doy']

def files(year,doy,period):
    """
    Reads in LENS historical and future files for 500mb heights
    
    
    Parameters
    ----------
    year : futyr or histyr
    doy : futdoy or histdoy
    period : 1920 or 2006 (integer arguments)
    
    Returns
    ----------
    z500 : 500mb heights (ens x doy x lat x lon)
    lat : array of longitudes
    lon : array of latitudes
    """
    ensembles = ['002','003','004','005','006','007','008','009','010','011','012','013','014','015','016','017','018','019','020','021','022','023','024','025','026','027','028','029','030']
    
    z500 = []
    for i in xrange(len(year)):
        end = int((year[i])*365. + (doy[i]+1))
        initial = end-3
        
        if period == 1920:
            directory = '/volumes/data/gcm/cesm-lens/B20TRC5CNBDRD/Aday/Z500/CESM1-CAM5/%s/' % ensembles[i]
            path = 'b.e11.B20TRC5CNBDRD.f09_g16.%s.cam.h1.Z500.19200101-20051231.nc' % ensembles[i]
            filepath = directory + path
            data = Dataset(filepath)
            Z500 = data.variables['Z500'][initial:end,:,:]
            lat = data.variables['lat'][:]
            lon = data.variables['lon'][:]
            data.close()
            
            z500.append(Z500)
        elif period == 2006:
            directory = '/volumes/data/gcm/cesm-lens/BRCP85C5CNBDRD/Aday/Z500/CESM1-CAM5/%s/' % ensembles[i]
            path = 'b.e11.BRCP85C5CNBDRD.f09_g16.%s.cam.h1.Z500.20060101-20801231.nc' % ensembles[i]
            filepath = directory + path
            data = Dataset(filepath)
            Z500 = data.variables['Z500'][initial:end,:,:]
            lat = data.variables['lat'][:]
            lon = data.variables['lon'][:]
            data.close()
            
            z500.append(Z500)
            
    z500 = np.asarray(z500)
            
    return z500, lat, lon
    
z500f, lat, lon = files(futyr,futdoy,2006)    
z500h, lat, lon = files(histyr,histdoy,1920)

def zonal(z500):
    """
    Calculates 500mb zonal height anomaly
    
    
    Parameters
    ----------
    z500 : 500mb heights (ens x doy x lat x lon)
    
    Returns
    ----------
    zonal : 500mb zonal anomaly (ens x lat x lon)
    meandomain : 500mb zonal anomaly mean across Great lakes (ens x lat x lon)
    """
    mean = np.empty((z500.shape[0],z500.shape[1],z500.shape[2]))
    for i in xrange(z500.shape[0]):
        for j in xrange(z500.shape[1]):
            for k in xrange(z500.shape[2]):
                mean[i,j,k] = np.nanmean(z500[i,j,k,:])   
    zonals = np.empty(z500.shape)
    for i in xrange(z500.shape[3]):
        zonals[:,:,:,i] = mean.copy()       
    zonalheights = z500 - zonals
    zonal = np.nanmean(zonalheights,axis=1)
    
    meandomain = np.empty((zonal.shape[0]))
    for i in xrange(zonal.shape[0]):
        meandomain[i] = np.nanmean(zonal[i,123:154,185:244])
                    
    return zonal, meandomain
zonalh, meandomainh = zonal(z500h)
zonalf, meandomainf = zonal(z500f)
  
def plot(z500,zonal,lat,lon):
    """
    Plots height anomalies from previous functions
    
    
    Parameters
    ----------
    z500 : 500mb heights (ens x doy x lat x lon)
    zonal : 500mb zonal anomaly (ens x lat x lon)
    lat : array of latitudes
    lon : array of longitudes
    
    Returns
    ----------
    zonal : 500mb zonal anomaly (ens x lat x lon)
    meandomain : 500mb zonal anomaly mean across Great lakes (ens x lat x lon)
    """
    
    lons, lats = np.meshgrid(lon,lat)    
    
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
    m = Basemap(projection='merc',llcrnrlon=183,llcrnrlat=25,urcrnrlon=297,
            urcrnrlat=61,resolution='l')           
    m.drawstates()
    m.drawcountries()
    m.drawmapboundary(fill_color = 'white')
    m.drawcoastlines(color='black',linewidth=0.5)
    m.drawlsmask(land_color='grey',ocean_color='w')
    x,y = m(lons,lats)
    
#    cs = m.contourf(x,y,sts.nanmean(z500[0][0]))
    plot_rec(m,lonmin,lonmax,latmin,latmax)
#    cs.set_cmap('jet')
    
    ax1.spines['top'].set_linewidth(3)
    ax1.spines['right'].set_linewidth(3)
    ax1.spines['bottom'].set_linewidth(3)
    ax1.spines['left'].set_linewidth(3)
    
    ax1.text(0.18,0.015,'Average LENS',size='8',horizontalalignment= 'center',
            backgroundcolor='white',verticalalignment= 'center',
            bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
            transform=ax1.transAxes) 
            
    for i in xrange(len(zonal)):
        ax = plt.subplot(6,5,i+2)
        m = Basemap(projection='merc',llcrnrlon=183,llcrnrlat=25,urcrnrlon=297,
            urcrnrlat=61,resolution='l')           
        m.drawstates()
        m.drawcountries()
        m.drawmapboundary(fill_color = 'white')
        m.drawcoastlines(color='black',linewidth=0.5)
        m.drawlsmask(land_color='grey',ocean_color='w')
        x,y = m(lons,lats)
        
        z500m = zonal[i,:,:]
        
        z500m[np.where(z500m)<-500]=-500
        z500m[np.where(z500m)>500]=500
        cs = m.contour(x,y,z500m,range(-500,600,100),colors='k')
        cs = m.contourf(x,y,z500m,range(-500,520,10))
        cs.set_cmap('RdYlBu_r')
    
        ax.text(0.16,0.015,'Member %i' % (member[i]+1),size='8',horizontalalignment= 'center',
                backgroundcolor='white',verticalalignment= 'center',
                bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
                transform=ax.transAxes)    
    plt.tight_layout()
    fig.subplots_adjust(bottom=0.098)
    cbar_ax = fig.add_axes([0.15, 0.08, 0.7, 0.01])
    cbar = fig.colorbar(cs, cax=cbar_ax, orientation = 'horizontal',
                        extend='both',extendfrac='auto',ticks=np.arange(-500,600,100))
    cbar.set_label('Geopotential Heights (m)') 
    figure_title = 'LENS 1920-2005, 500mb Zonal Height Anomaly'
    fig.text(0.5, .97, figure_title,
         horizontalalignment='center',
         fontsize=14)
#    plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/ensemble_histheights.eps',dpi=400,format='eps')    
#plot(z500h,zonal,lat,lon)
#plot(z500f,zonal,lat,lon)
         
### Plot Data in Box/Whisker Plot
dataq = [meandomainh,meandomainf]
         
fig=plt.figure()
ax = fig.add_subplot(111)
bp = plt.boxplot(dataq, patch_artist=True)
for box in bp['boxes']:
    # change outline color
    box.set( color='k', linewidth=2)
    # change fill color
    box.set( facecolor = 'lightgray',alpha=0.7)
for median in bp['medians']:
    median.set(color='k', linewidth=3,linestyle='solid')
for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=2)
for cap in bp['caps']:
    cap.set(color='k', linewidth=1)
for flier in bp['fliers']:
    flier.set(marker='o', color='k', alpha=1)
ax.set_xticklabels(['1920-2005','2006-2080'])
plt.xlabel('Years',fontsize=13)
plt.ylabel('Anomaly',fontsize=13)
fig.suptitle('LENS Early Springs March 15-30',fontsize=18)
plt.title('500mb Geopotential Height Anomaly',fontsize=13)
plt.grid(True)
plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/boxplot.png',dpi=400)