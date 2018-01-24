"""
*Script calculates SI-x trends over the entire CESM Control run*
"""

from control_SIX_datareader import ccsm
import numpy as np
from scipy import stats as sts
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

### Call CESM Control Data
years = 'b.e11.B1850C5CN.f09_g16.005.cam.h1.TRE.SI-x.402-999.nc'
longitude, latitude, leaf_index, bloom_index, lstfrz_index = ccsm(years)
yrs = np.array(xrange(0,30))

### Slice Grid for contiguous United States
lonq = np.where((longitude > 180) & (longitude < 304))
lonq = np.squeeze(lonq)
lon = longitude[lonq]
latq = np.where((latitude > 23) & (latitude < 77))
latq = np.squeeze(latq)
lat = latitude[latq]

lons,lats = np.meshgrid(lon,lat)
lf = leaf_index[0:570,latq,:]
bl = bloom_index[0:570,latq,:]
lstfrz = lstfrz_index[0:570,latq,:]

def trend(var):
    """
    Calculates 100 Year Trends for SI-x    
    
    
    Parameters
    ----------
    var : Leaf, bloom , or last freeze index (lf, bl, or lstfrz (as string))

    Returns            
    ----------
    meantrend : Overal SI-x trend
    """
    varn = []
    for i in xrange(0,570,30):
        start = i
        end = i+30
        q = np.array(xrange(start,end))
        SI = var[q]
        varn.append(SI)
    
    r=[]
    for ys in xrange(len(varn)):
        varq = varn[ys]
        for j in xrange(len(lat)):
            for k in xrange(len(lon)):
                SI = varq[:,j,k]
                slope, intercept, r_value, p_value, std_err = sts.linregress(yrs,SI)
                r.append(slope)
    
    trend = np.array(r)
    trendn = trend.reshape(19,53,99)  
    
    return trendn

trendn = trend(bl)

fig = plt.figure() 
fig.suptitle('Bloom index trends from 401-999 (30-yr periods) \n',fontsize=15)

### Plot SIx trends from CESM Control
for years in xrange(len(trendn)):
    ax = plt.subplot(5,4,years+1)

    m = Basemap(projection='cyl',llcrnrlon=180,llcrnrlat=25,urcrnrlon=304,
            urcrnrlat=74,resolution='l')           
    m.drawmapboundary(fill_color = 'white')
    m.drawcoastlines(color='black',linewidth=0.5)
    m.drawlsmask(land_color='grey',ocean_color='w')
    x,y = m(lons,lats)
    
    trendq = trendn[years]*10.
    trendq[np.where(trendq > 4.0)] = 4.0
    trendq[np.where(trendq < -4.0)] = -4.0
    
    cs = m.contourf(x,y,trendq,np.arange(-4,4.1,0.1))
    
    ax.text(0.1,0.2,years+1,size='small',horizontalalignment= 'center',
            backgroundcolor='white',verticalalignment= 'center',
            bbox=dict(facecolor='white',edgecolor='black',alpha=1.0),
            transform=ax.transAxes)

cs.set_cmap('jet')
plt.tight_layout()
fig.subplots_adjust(bottom=0.15)
cbar_ax = fig.add_axes([0.15, 0.08, 0.7, 0.05])
cbar = fig.colorbar(cs, cax=cbar_ax, orientation = 'horizontal',ticks=[-4,-2,0,2,4])
cbar.set_label('Trend (days/decade)') 
plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/cesm_control_bl_402_999.eps',
            format='eps',dpi=300)      