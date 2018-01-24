"""
*Calculates plots for SIx from future LENS*
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import scipy.stats as sts
from mpl_toolkits.basemap import Basemap

directory = '/volumes/eas-shared/ault/ecrl/spring-indices/data/'

def SIx():
    """
    Reads in future LENS SI-x data

    
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
    
    leaf=[]
    lstfrz = []
    for version in versions:
        years = 'b.e11.BRCP85C5CNBDRD.f09_g16.%s.cam.h.SI-x.2006-2080.nc' % version
        filename = directory + years
        values = Dataset(filename)
        lon = values.variables['lon'][189:240]
        lat = values.variables['lat'][:32]
        lstfrz_index = values.variables['lstfrz_index'][:,:32,189:240]
        leaf_index = values.variables['leaf_index'][:,:32,189:240]
        values.close()
        
        leaf.append(leaf_index)
        lstfrz.append(lstfrz_index)
    latmean = np.asarray(lstfrz)
    leafmean = np.asarray(leaf)
    print 'Done! 1'
    return leafmean, latmean, lstfrz, lat, lon
    
leafmean, latmean, lstfrz, lat, lon = SIx()

def Trends(leafmean, latmean, lstfrz, lat, lon):
    """
    Reads in BEST SI-x data

    
    Parameters
    ----------
    leafmean : array leaf indices (ens x year x lat x lon)
    latmean : array last freeze indices (ens x year x lat x lon)
    lat : array of latitudes
    lon : array of longitudes 
    lstfrz : list last freeze indices 
    
    Returns
    ----------
    trend : array of SIx trends (time x lat x lon)
    lat : array of latitudes
    lon : array of longitudes
    """
    
    lons,lats = np.meshgrid(lon,lat)
    lstave = sts.nanmean(lstfrz)
    lstfrz_totalmean = np.nanmean(latmean)
    
    #meanlstfrz = []
    #for i in xrange(len(lstave)):
    #    yrslstfrz = np.nanmean(lstave[i,:,:])
    #    meanlstfrz.append(yrslstfrz)
    #meanlstfrz = np.asarray(meanlstfrz)
    
    ### Save Files and Import Leaf
    #np.savetxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/2006_2008lstfrz.txt',(meanlstfrz),delimiter=',')
    #
    fut_lf = np.genfromtxt('/volumes/zml5/research/ccsm/text/2006_2008lf.txt',delimiter=',')
    ave_lstfrz = sts.nanmean(fut_lf)
    #
    #### Calculate Damage Index
    #avedamage = meanlstfrz - ave_lstfrz
    
    trend = []
    lines = []
    latdif =[]
    for i in xrange(latmean.shape[0]):
        timeq = np.asarray(list(xrange(latmean.shape[1])))
        for j in xrange(latmean.shape[2]):
            for k in xrange(latmean.shape[3]):
                lstslope, lstintercept, r_value, p_value, std_err = sts.linregress(timeq,leafmean[i,:,j,k])
                lst_line = lstslope*timeq+lstintercept
                trend.append(lstslope)
                lines.append(lst_line)
                
                latdifq = lst_line-lstfrz_totalmean
                latdif.append(latdifq)
    #latdif = np.reshape(np.asarray(latdif),(latmean.shape[0],latmean.shape[2],latmean.shape[3]))
    trend = np.reshape(np.asarray(trend),(latmean.shape[0],latmean.shape[2],latmean.shape[3]))
    trend = trend*10.
    lines = np.asarray(lines)
    
    ### Standard Deviation and Mean
    std2006= np.nanstd(trend)
    mean2080 = np.nanmean(trend)
    
    ### Damage Index
    damage_members = leafmean - latmean
    stddamage2006 = np.nanstd(damage_members)
    meandamage2006 = np.nanmean(damage_members)
    
    std1920 = 27.444
    mean1920 = -13.8640801
    
    damagez = (damage_members-mean1920)/std1920
    
    damagetrend = []
    linesdamage = []
    for i in xrange(latmean.shape[0]):
        timeq = np.asarray(list(xrange(latmean.shape[1])))
        for j in xrange(latmean.shape[2]):
            for k in xrange(latmean.shape[3]):
                damageslope, damageintercept, r_value, p_value, std_err = sts.linregress(timeq,damage_members[i,:,j,k])
                damage_line = damageslope*timeq+damageintercept
                damagetrend.append(damageslope)
                linesdamage.append(damage_line)
    damagetrend = np.reshape(np.asarray(damagetrend),(latmean.shape[0],latmean.shape[2],latmean.shape[3]))
    damagetrend = damagetrend*10.
    linesdamage = np.asarray(linesdamage)
    meantrends = sts.nanmean(trend)
    
    damagezscores = np.empty(trend.shape)
    for i in xrange(latmean.shape[0]):
        for j in xrange(latmean.shape[2]):
            for k in xrange(latmean.shape[3]):
                damagezscores[i,j,k] = np.mean(damagez[i,:,j,k])
                
    ### Draw Polygon
    def plot_rec(bmap, lonmin,lonmax,latmin,latmax):
        xs = [lonmin,lonmax,lonmax,lonmin,lonmin]
        ys = [latmin,latmin,latmax,latmax,latmin]
        bmap.plot(xs, ys, latlon = True, color='k',linewidth=1.5,linestyle='solid')
    lonmin = -101.5
    lonmax =  -75.5
    latmin =  37.5
    latmax =  50.5
    
    meandamagetrends = sts.nanmean(damagetrend)
    meandamagetrends[np.where(meandamagetrends > 6)] = 6
    meandamagetrends[np.where(meandamagetrends < -6)] = -6
    
    member = list(xrange(1,30))
    ### Plot Trends
    fig = plt.figure()   
#    fig.suptitle('LENS 2006-2040, LSTFRZ Index Trends',fontsize=10)
    ax1 = plt.subplot(6,5,1)
    m = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=31,urcrnrlon=298,
                urcrnrlat=54,resolution='l')           
    m.drawstates()
    m.drawcountries()
    m.drawmapboundary(fill_color = 'white')
    m.drawcoastlines(color='black',linewidth=0.5)
    m.drawlsmask(land_color='grey',ocean_color='w')
    x,y = m(lons,lats)
    
#    meantrends[np.where(meantrends > 6)] = 6
#    meantrends[np.where(meantrends< -6)] = -6
    
    cs = m.contourf(x,y,meandamagetrends,np.arange(-6.,6.1,.1))
    plot_rec(m,lonmin,lonmax,latmin,latmax)
    cs.set_cmap('bwr_r')
    
    ax1.spines['top'].set_linewidth(3)
    ax1.spines['right'].set_linewidth(3)
    ax1.spines['bottom'].set_linewidth(3)
    ax1.spines['left'].set_linewidth(3)
    
    ax1.text(0.18,0.015,'Average LENS',size='8',horizontalalignment= 'center',
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
    
        damagetrend[np.where(damagetrend > 6)] = 6
        damagetrend[np.where(damagetrend < -6)] = -6
        
        cs = m.contourf(x,y,damagetrend[i,:,:],np.arange(-6.,6.1,.1))
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
    figure_title = 'LENS 2006-2080, Damage Index Trends'
    fig.text(0.5, .97, figure_title,
         horizontalalignment='center',
         fontsize=14)
    plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/ensemble_damagetrend_0680.eps',dpi=400,format='eps')
    
#    m = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=31,urcrnrlon=298,
#                    urcrnrlat=54,resolution='l')           
#    m.drawstates()
#    m.drawcountries()
#    m.drawmapboundary(fill_color = 'white')
#    m.drawcoastlines(color='black')
#    m.drawlsmask(land_color='grey',ocean_color='w')
#    x,y = m(lons,lats)
#    cs = m.contourf(x,y,meandamagetrend,np.arange(-6,6.1,.1))
#    cs.set_cmap('bwr_r')
#    cbar = m.colorbar(cs,location='bottom',pad='5%')
#    cbar.set_label('days/decade')
#    cbar.set_ticks(np.arange(-6.,6.5,1))
#    plt.title('LENS, 2006-2080 Mean Damage Index Trends')
#    plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/meandamage_0680.png',dpi=300)
#    
#    d = sts.nanmean(damage_members)
#    e = np.empty((d.shape[0]))
#    for i in xrange(d.shape[0]):
#        e[i] = np.nanmean(d[i,:-5,-5])
    
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
    #np.savetxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/damagevalues_2006-2080.txt',(damagevalues),delimiter=',')
    #np.savetxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/lstfrzvalues_2006-2080.txt',(lstfrzvalues),delimiter=',')
    #np.savetxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/damagetimes_2006-2080.txt',(damagetimes),delimiter=',')
    return trend, lat, lon
    
trend,lat,lon = Trends(leafmean, latmean, lstfrz, lat, lon)

###############################################################################
### Create netcdf files of trends

#def netcdf(trend, lat, lon):
#    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/'
#    name = 'LSTFRZ_2040.nc'
#    filename = directory + name
#    ncfile = Dataset(filename,'w',format='NETCDF4')
#    ncfile.description = 'LSTFRZ Trends 2006-2040'
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
#netcdf(trend,lat,lon)