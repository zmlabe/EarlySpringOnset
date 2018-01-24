"""
*Script creates figure for CESM Control LF values*
"""
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.stats import zscore
from mpl_toolkits.basemap import Basemap
import scipy.stats as sts

directory= '/volumes/eas-shared/ault/ecrl/spring-indices/data/'

### CESM control for Great Lakes region
years = 'b.e11.B1850C5CN.f09_g16.005.cam.h1.TRE.SI-x.402-999.nc'
filename = directory + years
values = Dataset(filename)
lon = values.variables['lon'][207:229]
lat = values.variables['lat'][12:27]
leaf_index = values.variables['leaf_index'][:,12:27,207:229]
bloom_index = values.variables['bloom_index'][:,12:27,207:229]
values.close()

lons,lats = np.meshgrid(lon,lat)

### Average SI per year
lfmean = []
blmean = []
for i in xrange(len(leaf_index)):
    leafs = np.nanmean(leaf_index[i,:,:])
    blooms = np.nanmean(bloom_index[i,:,:])
    lfmean.append(leafs)
    blmean.append(blooms)
    
### Climatology
lf_totalmean = np.nanmean(lfmean)
bl_totalmean = np.nanmean(blmean)

### Anom per year
lf_anom = []
bl_anom = []
for j in xrange(len(lfmean)):
    anomlf = lf_totalmean - lfmean[j]
    anombl = bl_totalmean - blmean[j]
    lf_anom.append(anomlf)
    bl_anom.append(anombl)

### z-score per year
lf_zscorez = zscore(lfmean)
bl_zscorez = zscore(blmean)
yrlfz2 = np.where(lf_zscorez<=-2)
yrblz2 = np.where(bl_zscorez<=-2)
yrlfz3 = np.where(lf_zscorez<=-3)
yrblz3 = np.where(bl_zscorez<=-3)

### Adjust z-score
timeq = np.asarray(list(xrange(0,598)))
lfslope, lfintercept, r_value, p_value, std_err = sts.linregress(timeq,lfmean)
lf_line = lfslope*timeq+lfintercept
blslope, blintercept, r_value, p_value, std_err = sts.linregress(timeq,blmean)
bl_line = blslope*timeq+blintercept

lf_diff = lf_line-lf_totalmean
bl_diff = bl_line-bl_totalmean

lf_notrend = lfmean-lf_diff
bl_notrend = blmean-bl_diff

lf_zscore = zscore(lf_notrend)
bl_zscore = zscore(bl_notrend)
yrlf2 = np.where(lf_zscore<=-2)
yrbl2= np.where(bl_zscore<=-2)
yrlf3 = np.where(lf_zscore<=-3)
yrbl3 = np.where(bl_zscore<=-3)

### Plot Time Series
fig = plt.figure()
threshold3 = [81.5]*len(lfmean)
threshold2 = [86.0]*len(lfmean)
plt.plot(lfmean,color='DarkOrange',linewidth=1,label=' Leaf Dates')
line2 = plt.plot(threshold2,'darkred',linewidth=3,label='-2$\sigma$ Threshold',
                 linestyle='dashed')
line3 = plt.plot(threshold3,'r',linewidth=3,label='-3$\sigma$ Threshold',
                 linestyle='dashed')
plt.legend(loc=1,prop={'size':9},shadow=True)
x = list(xrange(0,len(lfmean)+3,100))
labels = [400,500,600,700,800,900,1000]
plt.xticks(x,labels)
plt.xlabel('Years')
plt.grid(True)
fig.suptitle('CESM Control, First Leaf Index',fontsize=20)
fig.text(0.04,0.5,'DOY',va='center',rotation='vertical')
plt.savefig('/Users/zlabe/documents/CESMspring/Fig6.eps',dpi=400,format='eps')

### Output to text file
#np.savetxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/controllf.txt',(lfmean),delimiter=',')