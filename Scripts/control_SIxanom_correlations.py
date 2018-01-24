"""
*Script compares March 2012 to CESM control for spatial pattern correlations.
A variety of plots are produced*
"""

from control_SIX_datareader import ccsm
import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, interp
import scipy.stats as sts
from scipy.stats import nanmean
import matplotlib.pyplot as plt
import matplotlib.colors as c

### Make Colormap for SI-x
cmap = plt.get_cmap('spring')
cmap2 = plt.get_cmap('summer_r')
cmaplist = [cmap(i) for i in range(cmap.N)]
cmaplist2 = [cmap2(i) for i in range(cmap2.N)]

cms=c.ListedColormap(cmaplist+cmaplist2)

directory = '/volumes/eas-shared/ault/ecrl/spring-indices/data/SI-x/'

### Read in Netcdf file for SIx Averages
filename = directory + 'BEST-Cat_SI-X_Daily_LatLong1_1880-to-2013.25N_to_85N.nc'
values = Dataset(filename)
longituden = values.variables['lon'][:]
latituden = values.variables['lat'][:]
lfs = values.variables['leaf_index'][:]
bls = values.variables['bloom_index'][:]
values.close()

print 'DONE 1'

lonqn = np.where((longituden > -120) & (longituden < -55))
lonqn = np.squeeze(lonqn)
lonr = longituden[lonqn]
latqn = np.where((latituden > 25) & (latituden < 56))
latqn = np.squeeze(latqn)
latr = latituden[latqn]

new_lon = (lonr + 360.)

lon_12, lat_12 = np.meshgrid(new_lon,latr)

leafn = lfs[:,latqn,55:125]
bloomn = bls[:,latqn,55:125]

leafqn = []
for yr in xrange(101,131):
    leafnn = leafn[yr,:,:]
    leafqn.append(leafnn)
leafqn = np.asarray(leafqn)
aveleaf = nanmean(leafqn)

bloomqn = []
for yrs in xrange(len(bloomn)):
    bloomnn = bloomn[yrs,:,:]
    bloomqn.append(bloomnn)
bloomqn = np.asarray(bloomqn)  
avebloom = nanmean(bloomqn)  

### Call CESM Data
years = 'b.e11.B1850C5CN.f09_g16.005.cam.h1.TRE.SI-x.402-999.nc'
longitude, latitude, leaf_index, bloom_index,lstfrz_index = ccsm(years)

print 'DONE 2'

lonq = np.where((longitude > 230) & (longitude < 305))
lonq = np.squeeze(lonq)
lon = longitude[lonq]
latq = np.where((latitude > 25) & (latitude < 55))
latq = np.squeeze(latq)
lat = latitude[latq]

lons,lats = np.meshgrid(lon,lat)
lf = leaf_index[:,latq,185:244] #0:570
bl = bloom_index[:,latq,185:244] #0:570
lstfrz = lstfrz_index[:,latq,185:244] 

### Call BEST SI-x Files for 2012
leaf_12 = lfs[-2,latqn,55:125]
bloom_12 = bls[-2,latqn,55:125]

### Interp BEST Grid onto CCSM Grid
new_leaf12 = interp(np.transpose(leaf_12),latr,new_lon,lats,lons)
new_bloom12 = interp(np.transpose(bloom_12),latr,new_lon,lats,lons)

#new_aveleaf = interp(np.transpose(aveleaf),latr,new_lon,lats,lons)
#new_avebloom = interp(np.transpose(avebloom),latr,new_lon,lats,lons)

### BEST Anomalies
leaf12_anom = leaf_12 - avebloom
bloom12_anom = bloom_12 - avebloom

new_leaf12_anom = interp(np.transpose(leaf12_anom),latr,new_lon,lats,lons)
new_bloom12_anom = interp(np.transpose(bloom12_anom),latr,new_lon,lats,lons)

### CCSM Anomalies
lfave = []
for ys in xrange(len(lf)):
    lf_ave = lf[ys,:,:]
    lfave.append(lf_ave)
lfave = np.asarray(lfave)
ccsm_aveleaf = nanmean(lfave)

blave = []
for yr in xrange(len(bl)):
    bl_ave = bl[yr,:,:]
    blave.append(bl_ave)
blave = np.asarray(blave)  
ccsm_avebloom = nanmean(blave)  

lf_anom = []
bl_anom = []
for k in xrange(len(lf)):
    lf_anomn = lf[k,:,:] - ccsm_aveleaf
    bl_anomn = bl[k,:,:] - ccsm_avebloom
    lf_anom.append(lf_anomn)
    bl_anom.append(bl_anomn)
lf_anom = np.asarray(lf_anom)
bl_anom = np.asarray(bl_anom)    

### Make into 1D arrays
lf12_anom = np.ravel(new_leaf12_anom)
bl12_anom = np.ravel(new_bloom12_anom)
lf12_anom[np.where(np.isnan(lf12_anom))] = 0.0
bl12_anom[np.where(np.isnan(bl12_anom))] = 0.0

lfq = []
blq = []
for i in xrange(len(lf)):
    lfn = np.ravel(lf_anom[i,:,:])
    lfn[np.where(np.isnan(lfn))] = 0.0
    lfq.append(lfn)
    
    bln = np.ravel(bl_anom[i,:,:])
    bln[np.where(np.isnan(bln))] = 0.0
    blq.append(bln)
  
### Correlation each year
lfcorr = []
blcorr = []

for j in xrange(len(lf)):
    lfcorrn = sts.pearsonr(lfq[j],lf12_anom)
    lfcorr.append(lfcorrn)
    blcorrn = sts.pearsonr(blq[j],bl12_anom)
    blcorr.append(blcorrn)

lfcorrs = np.asarray(lfcorr)  
blcorrs = np.asarray(blcorr)

lfcorrs = lfcorrs[:,0]
blcorrs = blcorrs[:,0]  
    
### Correlation Extremes        
lfmax = np.nanmax(lfcorrs)
lfmin = np.nanmin(lfcorrs)
blmax = np.nanmax(blcorrs)
blmin = np.nanmin(blcorrs)

### Plot Correlations
timeq = np.asarray(list(xrange(0,598)))
lfslope, lfintercept, r_value, p_value, std_err = sts.linregress(timeq,lfcorrs)
lf_line = lfslope*timeq+lfintercept
blslope, blintercept, r_value, p_value, std_err = sts.linregress(timeq,blcorrs)
bl_line = blslope*timeq+blintercept

lfcorrsq = np.sort(lfcorrs)
blcorrsq = np.sort(blcorrs)
#
#fig = plt.figure()
##
#ax1 = plt.subplot(2,1,1)
#plt.title('Leaf Index',fontsize=12)
#plt.plot(lfcorrs,color='DarkOrange')
##line1 = plt.plot(lf_line,linewidth=3,color='black')
##plt.legend(line1,['Regression'],loc=1,prop={'size':7})
#x = list(xrange(0,len(lfcorrs)+3,100))
#labels = [400,500,600,700,800,900,1000]
#plt.xticks(x,labels)
#plt.grid(True)
#ax1=plt.subplot(2,1,2,sharey=ax1)
#plt.title('Bloom Index',fontsize=12)
#plt.plot(blcorrs,color='darkolivegreen')
##line2 = plt.plot(bl_line,linewidth=3,color='black')
#plt.xlabel('Years',fontsize=11)
##plt.legend(line2,['Regression'],loc=1,prop={'size':7})
#labels = [400,500,600,700,800,900,1000]
#plt.xticks(x,labels)
#plt.grid(True)
#fig.suptitle('SI-x Pattern Correlations for CESM and March 2012',fontsize=20)
#fig.text(0.04,0.5,'Pearson Correlation Coefficient',va='center',rotation='vertical')

#plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/correlations_CESM.png',dpi=1000)

#### Rank top 9 correlations for leaf
lfrank = lfcorrsq[-9:]
blranks = np.where(lfcorrs >= 0.60040)[0]
lfranks = np.where(lfcorrs >= 0.44186)[0]
blrankn = [444,251,16,4,204,556,214,156,364]
lfrankn = [4,251,333,156,64,401,12,44,183]

blyear = [846,653,418,406,606,958,616,558,766]
lfyear = [406,653,735,558,466,803,414,446,585]

lfplot = []
for yrcorr in lfrankn:
    lfyr = lf[yrcorr,:,:]
    lfplot.append(lfyr)
lfplot = np.asarray(lfplot)  

lf_rank9 = []
for k in xrange(len(lfplot)):
    lf9 = lfplot[k,:,:] - ccsm_aveleaf
    lf_rank9.append(lf9)
lf_rank9 = np.asarray(lf_rank9) 

def plot_rec(bmap, lonmin,lonmax,latmin,latmax):
    xs = [lonmin,lonmax,lonmax,lonmin,lonmin]
    ys = [latmin,latmin,latmax,latmax,latmin]
    bmap.plot(xs, ys, latlon = True, color='k',linewidth=1.5,linestyle='solid')
lonmin = -101.5
lonmax =  -75.5
latmin =  37.5
latmax =  50.5

fig=plt.figure()
for K in xrange(len(lfplot)):
    ax = plt.subplot(3,3,K+1)

    m = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=26,urcrnrlon=300,
        urcrnrlat=54,resolution='l')           
    m.drawmapboundary(fill_color = 'white')
    m.drawcoastlines(color='black',linewidth=0.5)
    m.drawlsmask(land_color='grey',ocean_color='w')
    x,y = m(lons,lats) 
    
    lf_ranks = lf_rank9[K,:,:]
    lf_ranks[np.where(lf_ranks < -40)] = -40
    lf_ranks[np.where(lf_ranks > 40)] = 40
    
    cs = m.contourf(x,y,lf_ranks,range(-40,45,5))
    cs.set_cmap('RdYlGn')
    plot_rec(m,lonmin,lonmax,latmin,latmax)
    ax.text(0.1,0.125,lfyear[K],size='medium',horizontalalignment= 'center',
            backgroundcolor='white',verticalalignment= 'center',
            bbox=dict(facecolor='white',edgecolor='black',alpha=1.0),
            transform=ax.transAxes)

plt.tight_layout()
fig.subplots_adjust(bottom=0.13)
cbar_ax = fig.add_axes([0.15, 0.08, 0.7, 0.05])
cbar = fig.colorbar(cs, cax=cbar_ax, orientation = 'horizontal',
                    extend='both',extendfrac='auto',
                    ticks=list(xrange(-40,50,10)))
cbar.set_label('Difference (Days)') 
figure_title = 'CESM Control, First Leaf Index Anomalies'
fig.text(0.5, .955, figure_title,
         horizontalalignment='center',
         fontsize=20)
plt.savefig('/Users/zlabe/documents/CESMspring/Fig5.eps',dpi=400,format='eps')

#### Plot Leaf/Bloom Index for Ranks
#fig=plt.figure()
#for K in xrange(len(lfplot)):
##    fig.suptitle('CESM Leaf Index',fontsize=20)
#    ax = plt.subplot(3,3,K+1)
#
#    m = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=26,urcrnrlon=300,
#        urcrnrlat=54,resolution='l')           
#    m.drawmapboundary(fill_color = 'white')
#    m.drawstates()
#    m.drawcountries()
#    m.drawcoastlines(color='black',linewidth=0.5)
#    m.drawlsmask(land_color='grey',ocean_color='w')
#    x,y = m(lons,lats) 
#    
#    lf_ranks = lfplot[K,:,:]
#    lf_ranks[np.where(lf_ranks > 226)] = 226
#    
#    cs = m.contourf(x,y,lf_ranks,range(0,233,8))
#    cs.set_cmap(cms)
#    
#    ax.text(0.1,0.1,blyear[K],size='medium',horizontalalignment='center',
#            backgroundcolor='white',verticalalignment= 'center',
#            bbox=dict(facecolor='white',edgecolor='black',alpha=1.0),
#            transform=ax.transAxes)
#
#plt.tight_layout()
#fig.subplots_adjust(bottom=0.13)
#cbar_ax = fig.add_axes([0.15, 0.08, 0.7, 0.05])
#cbar = fig.colorbar(cs, cax=cbar_ax, orientation = 'horizontal',
#                    extend='both',extendfrac='auto',
#                    ticks=list(xrange(0,226,45)))
#figure_title = 'CESM Bloom Index'
#fig.text(0.5, .955, figure_title,
#         horizontalalignment='center',
#         fontsize=20)
#cbar.set_label('DOY') 
#plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/bl_rank9.png',dpi=300)  

### Calculate Covariance
lf_cov = []
for m in xrange(len(lfq)):
    si = np.vstack((lfq[m],lf12_anom))
    var = np.cov(si)
    lf_cov.append(var)
lf_cov = np.asarray(lf_cov)
lf_cov = lf_cov[:,0,1]

bl_cov = []
for m in xrange(len(blq)):
    si = np.vstack((blq[m],bl12_anom))
    var = np.cov(si)
    bl_cov.append(var)
bl_cov = np.asarray(bl_cov)
bl_cov = bl_cov[:,0,1]

#lfslope, lfintercept, r_value, p_value, std_err = sts.linregress(timeq,lf_cov)
#blslope, blintercept, r_value, p_value, std_err = sts.linregress(timeq,bl_cov)
#lf_line = lfslope*timeq+lfintercept
#bl_line = blslope*timeq+blintercept
#
#ax1 = plt.subplot(2,1,1)
#plt.title('CESM Leaf Index Covariance',fontsize=12)
#plt.plot(lf_cov,color='g')
#line3 = plt.plot(lf_line,linewidth=3,color='black')
#plt.legend(line3,['Regression'],loc=1,prop={'size':7})
#
#plt.subplot(2,1,2,sharey=ax1)
#plt.title('CESM Bloom Index Covariance',fontsize=12)
#plt.plot(bl_cov,color='g')
#line4 = plt.plot(bl_line,linewidth=3,color='black')
#plt.xlabel('Years',fontsize=11)
#plt.legend(line4,['Regression'],loc=1,prop={'size':7})
#
#fig.suptitle('SI-x Covariance for CESM and March 2012',fontsize=20)
#fig.text(0.04,0.5,'Covariance',va='center',rotation='vertical')
#
#plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/covariance.CESM.png',dpi=1000)

### Output to text files
lfss=[]
blss=[]
lst=[]
for i in xrange(len(lf)):
    out = np.nanmean(lf[i,:,:])
    out2 = np.nanmean(bl[i,:,:])
    out3 = np.nanmean(lstfrz[i,:,:])
    lst.append(out3)
    lfss.append(out)
    blss.append(out2)
    
#np.savetxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/controltrends.txt',(lfslope,blslope),delimiter=',')
#np.savetxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/controllf.txt',(lfss),delimiter=',')
#np.savetxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/controlbl.txt',(blss),delimiter=',')
#np.savetxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/controllstfrz.txt',(lst),delimiter=',')