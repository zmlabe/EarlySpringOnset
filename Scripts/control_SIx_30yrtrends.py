"""
*Script creates 30 year SI-x trends for CESM Control (per 100 year chunks)*
"""

from control_SIX_datareader import ccsm
import numpy as np
from scipy import stats as sts
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

### Call CCSM Data
years = 'b.e11.B1850C5CN.f09_g16.005.cam.h1.TRE.SI-x.500-599.nc'
longitude, latitude, leaf_index, bloom_index, lstfrz_index = ccsm(years)
print 'DONE!'
yrs = np.array(xrange(0,30))

### Slice Grid for contiguous United States
lonq = np.where((longitude > 235) & (longitude < 305))
lonq = np.squeeze(lonq)
lon = longitude[lonq]
latq = np.where((latitude > 25) & (latitude < 55))
latq = np.squeeze(latq)
lat = latitude[latq]

lons,lats = np.meshgrid(lon,lat)
lf = leaf_index[0:90,latq,189:244]
bl = bloom_index[0:90,latq,189:244]
lstfrz = lstfrz_index[0:90,latq,189:244]

def overallTrend(var):
    """
    Calculates 100 Year Trends for SI-x    
    
    
    Parameters
    ----------
    var : Leaf, bloom , or last freeze index

    Returns            
    ----------
    meantrend : Overal SI-x trend
    """
    yrs = list(xrange(var.shape[0]))
    meantrend = np.empty((var.shape[1],var.shape[2]))
    for i in xrange(var.shape[1]):
        for j in xrange(var.shape[2]):
            SI = var[:,i,j]
            slope, intercept, r_value, p_value, std_err = sts.linregress(yrs,SI)
            meantrend[i,j] = slope*10.
    return meantrend

meanleaf = overallTrend(leaf_index[:,latq,189:244])
meanlstfrz = overallTrend(lstfrz_index[:,latq,189:244])
   
def trend(var,pd):  
    """
    Calculates 30 Year Trends for SI-x
    
    
    Parameters
    ----------
    var : Leaf, bloom, or last freeze index
    pd : 1, 2, 3 (first 30 years, second 30 years, third 30 years)

    Returns            
    ----------
    trendn : 30 year chunk for SI-x
    """
    varn = []
    for i in xrange(0,90,30):
        start = i
        end = i+30
        q = np.array(xrange(start,end))
        varq = var[q]
        varn.append(varq)
    
    if pd == '1':
        varq = varn[0]
    elif pd == '2':  
        varq = varn[1]
    elif pd == '3':
        varq = varn[2]    
    
    r = []
    for j in xrange(len(lat)):
        for k in xrange(len(lon)):
            SI = varq[:,j,k]
            slope, intercept, r_value, p_value, std_err = sts.linregress(yrs,SI)
            r.append(slope)
    
    trend = np.array(r)
    trendn = trend.reshape(31,55)  

    return trendn  

### Draw Polygon for Great Lakes
def plot_rec(bmap, lonmin,lonmax,latmin,latmax):
    xs = [lonmin,lonmax,lonmax,lonmin,lonmin]
    ys = [latmin,latmin,latmax,latmax,latmin]
    bmap.plot(xs, ys, latlon = True, color='k',linewidth=1.5,linestyle='solid')
lonmin = -101.5
lonmax =  -75.5
latmin =  37.5
latmax =  50.5         
   
leaf1 = trend(lf,'1')
leaf2 = trend(lf,'2')
leaf3 = trend(lf,'3')
#bloom1  = trend(bl,'1')
#bloom2 = trend(bl,'2')
#bloom3 = trend(bl,'3')
lstfrz1 = trend(lstfrz,'1')
lstfrz2 = trend(lstfrz,'2')  
lstfrz3 = trend(lstfrz,'3') 

leaf1[np.where(leaf1 < -.6)] = -.6
leaf2[np.where(leaf2 < -.6)] = -.6
leaf3[np.where(leaf3 < -.6)] = -.6
leaf1[np.where(leaf1 > .6)] = .6
leaf2[np.where(leaf2 > .6)] = .6
leaf3[np.where(leaf3 > .46)] = .6
#bloom1[np.where(bloom1 < -.4)] = -.4
#bloom2[np.where(bloom2 < -.4)] = -.4
#bloom3[np.where(bloom3 < -.4)] = -.4
#bloom1[np.where(bloom1 > .4)] = .4
#bloom2[np.where(bloom2 > .4)] = .4
#bloom3[np.where(bloom3 > .4)] = .4
lstfrz1[np.where(lstfrz1 < -.6)] = -.6
lstfrz2[np.where(lstfrz2 < -.6)] = -.6
lstfrz3[np.where(lstfrz3 < -.6)] = -.6
lstfrz1[np.where(lstfrz1 >= .6)] = .6
lstfrz2[np.where(lstfrz2 >= .6)] = .6
lstfrz3[np.where(lstfrz3 >= .6)] = .6

#   
### Plot Trends  
fig = plt.figure() 
ax9 = plt.subplot(421)
m = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=31,urcrnrlon=298,
            urcrnrlat=54,resolution='l')           
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color='black',linewidth=0.5)
m.drawlsmask(land_color='grey',ocean_color='w')
x,y = m(lons,lats)

meanleaf[np.where(meanleaf > 6)] = 6
meanleaf[np.where(meanleaf< -6)] = -6

cs9 = m.contourf(x,y,meanleaf,np.arange(-6.,6.1,.1))
plot_rec(m,lonmin,lonmax,latmin,latmax)
cs9.set_cmap('bwr_r')

ax9.spines['top'].set_linewidth(3)
ax9.spines['right'].set_linewidth(3)
ax9.spines['bottom'].set_linewidth(3)
ax9.spines['left'].set_linewidth(3)

plt.title('First Leaf Index',fontsize=9)

ax9.text(0.18,0.015,'Century Trend',size='8',horizontalalignment= 'center',
        backgroundcolor='white',verticalalignment= 'center',
        bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
        transform=ax9.transAxes) 
 
### Panel 1
ax1 = fig.add_subplot(423)   
m = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=31,urcrnrlon=298,
            urcrnrlat=54,resolution='l')           
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color='black',linewidth=0.5)
m.drawlsmask(land_color='grey',ocean_color='w')
x,y = m(lons,lats)
cs = m.contourf(x,y,leaf1*10.,np.arange(-6,6.1,0.1))
cs.set_cmap('bwr_r')

ax1.text(0.11,0.015,'Period 1',size='8',horizontalalignment= 'center',
        backgroundcolor='white',verticalalignment= 'center',
        bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
        transform=ax1.transAxes) 

### Panel 2
ax2 = fig.add_subplot(425)   
m1 = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=31,urcrnrlon=298,
            urcrnrlat=54,resolution='l')           
m1.drawstates()
m1.drawcountries()
m1.drawmapboundary(fill_color = 'white')
m1.drawcoastlines(color='black',linewidth=0.5)
m1.drawlsmask(land_color='grey',ocean_color='w')
x,y = m1(lons,lats)
cs1 = m1.contourf(x,y,leaf2*10.,np.arange(-6,6.1,0.1))
cs1.set_cmap('bwr_r')

ax2.text(0.11,0.015,'Period 2',size='8',horizontalalignment= 'center',
        backgroundcolor='white',verticalalignment= 'center',
        bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
        transform=ax2.transAxes) 

#### Panel 3
ax3 = fig.add_subplot(427)   
m3 = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=31,urcrnrlon=298,
            urcrnrlat=54,resolution='l')           
m3.drawstates()
m3.drawcountries()
m3.drawmapboundary(fill_color = 'white')
m3.drawcoastlines(color='black',linewidth=0.5)
m3.drawlsmask(land_color='grey',ocean_color='w')
x,y = m3(lons,lats)
cs3 = m3.contourf(x,y,leaf3*10.,np.arange(-6,6.1,0.1))
cs3.set_cmap('bwr_r')

ax3.text(0.11,0.015,'Period 3',size='8',horizontalalignment= 'center',
        backgroundcolor='white',verticalalignment= 'center',
        bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
        transform=ax3.transAxes) 

#### Panel 4
#ax4 = fig.add_subplot(332)   
#m4 = Basemap(projection='merc',llcrnrlon=235.5,llcrnrlat=26,urcrnrlon=300,
#            urcrnrlat=54,resolution='l')           
#m4.drawstates()
#m4.drawcountries()
#m4.drawmapboundary(fill_color = 'white')
#m4.drawcoastlines(color='black',linewidth=0.5)
#m4.drawlsmask(land_color='grey',ocean_color='w')
#x,y = m4(lons,lats)
#cs4 = m4.contourf(x,y,bloom1*10.,np.arange(-4,4.1,0.1))
#cs4.set_cmap('jet')
#cbar4 = m4.colorbar(cs4,location='bottom',pad='5%',ticks=list(xrange(-4,5,1)))
#cbar4.set_label('trend per decade')
#plt.title('Bloom Index Trends')
#
#### Panel 5
#ax5 = fig.add_subplot(335)   
#m5 = Basemap(projection='merc',llcrnrlon=235.5,llcrnrlat=26,urcrnrlon=300,
#            urcrnrlat=54,resolution='l')           
#m5.drawstates()
#m5.drawcountries()
#m5.drawmapboundary(fill_color = 'white')
#m5.drawcoastlines(color='black',linewidth=0.5)
#m5.drawlsmask(land_color='grey',ocean_color='w')
#x,y = m5(lons,lats)
#cs5 = m5.contourf(x,y,bloom2*10.,np.arange(-4,4.1,0.1))
#cs5.set_cmap('jet')
#cbar5 = m4.colorbar(cs5,location='bottom',pad='5%',ticks=list(xrange(-4,5,1)))
#
#### Panel 6
#ax6 = fig.add_subplot(338)   
#m6 = Basemap(projection='merc',llcrnrlon=235.5,llcrnrlat=26,urcrnrlon=300,
#            urcrnrlat=54,resolution='l')           
#m6.drawstates()
#m6.drawcountries()
#m6.drawmapboundary(fill_color = 'white')
#m6.drawcoastlines(color='black',linewidth=0.5)
#m6.drawlsmask(land_color='grey',ocean_color='w')
#x,y = m6(lons,lats)
#cs6 = m6.contourf(x,y,bloom3*10.,np.arange(-4,4.1,0.1))
#cs6.set_cmap('jet')
#cbar6 = m4.colorbar(cs6,location='bottom',pad='5%',ticks=list(xrange(-4,5,1)))

### Panel 7
ax10 = plt.subplot(422)
m = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=31,urcrnrlon=298,
            urcrnrlat=54,resolution='l')           
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color='black',linewidth=0.5)
m.drawlsmask(land_color='grey',ocean_color='w')
x,y = m(lons,lats)

meanlstfrz[np.where(meanlstfrz >= 6)] = 6
meanlstfrz[np.where(meanlstfrz< -6)] = -6

cs10 = m.contourf(x,y,meanlstfrz,np.arange(-6.,6.1,.1))
plot_rec(m,lonmin,lonmax,latmin,latmax)
cs10.set_cmap('bwr_r')

ax10.spines['top'].set_linewidth(3)
ax10.spines['right'].set_linewidth(3)
ax10.spines['bottom'].set_linewidth(3)
ax10.spines['left'].set_linewidth(3)
plt.title('LSTFRZ Index',fontsize=9)

ax10.text(0.18,0.015,'Century Trend',size='8',horizontalalignment= 'center',
        backgroundcolor='white',verticalalignment= 'center',
        bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
        transform=ax10.transAxes) 

ax7 = fig.add_subplot(424)   
m7 = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=31,urcrnrlon=298,
            urcrnrlat=54,resolution='l')           
m7.drawstates()
m7.drawcountries()
m7.drawmapboundary(fill_color = 'white')
m7.drawcoastlines(color='black',linewidth=0.5)
m7.drawlsmask(land_color='grey',ocean_color='w')
x,y = m7(lons,lats)
cs7 = m7.contourf(x,y,lstfrz1*10.,np.arange(-6,6.1,0.1))
cs7.set_cmap('bwr_r')

ax7.text(0.11,0.015,'Period 1',size='8',horizontalalignment= 'center',
        backgroundcolor='white',verticalalignment= 'center',
        bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
        transform=ax7.transAxes) 

### Panel 8
ax8 = fig.add_subplot(426)   
m8 = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=31,urcrnrlon=298,
            urcrnrlat=54,resolution='l')           
m8.drawstates()
m8.drawcountries()
m8.drawmapboundary(fill_color = 'white')
m8.drawcoastlines(color='black',linewidth=0.5)
m8.drawlsmask(land_color='grey',ocean_color='w')
x,y = m8(lons,lats)
cs8 = m8.contourf(x,y,lstfrz2*10.,np.arange(-6,6.1,0.1))
cs8.set_cmap('bwr_r')

ax8.text(0.11,0.015,'Period 2',size='8',horizontalalignment= 'center',
        backgroundcolor='white',verticalalignment= 'center',
        bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
        transform=ax8.transAxes) 

### Panel 9
ax9 = fig.add_subplot(428)   
m9 = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=31,urcrnrlon=298,
            urcrnrlat=54,resolution='l')           
m9.drawstates()
m9.drawcountries()
m9.drawmapboundary(fill_color = 'white')
m9.drawcoastlines(color='black',linewidth=0.5)
m9.drawlsmask(land_color='grey',ocean_color='w')
x,y = m9(lons,lats)
cs9 = m9.contourf(x,y,lstfrz3*10.,np.arange(-6,6.1,0.1))
cs9.set_cmap('bwr_r')

ax9.text(0.11,0.015,'Period 3',size='8',horizontalalignment= 'center',
        backgroundcolor='white',verticalalignment= 'center',
        bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
        transform=ax9.transAxes) 

plt.tight_layout()
fig.subplots_adjust(bottom=0.12)
cbar_ax = fig.add_axes([0.15, 0.08, 0.7, 0.01])
cbar = fig.colorbar(cs, cax=cbar_ax, orientation = 'horizontal',
                    extend='both',extendfrac='auto',ticks=np.arange(-6.,7,1))
cbar.set_label('days/decade')
figure_title = 'CESM 500-599, SI-x Trends'
fig.text(0.5, .965, figure_title,
     horizontalalignment='center',
     fontsize=14)
plt.savefig('/Users/zlabe/documents/CESMspring/Fig4.eps',dpi=400,format='eps')