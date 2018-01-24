import numpy as np
from netCDF4 import Dataset
from scipy.stats import zscore
import matplotlib.pyplot as plt
import scipy.stats as sts
from mpl_toolkits.basemap import Basemap

directory = '/volumes/eas-shared/ault/ecrl/spring-indices/data/'

versions = ['002','003','004','005','006','007','008','009','010','011','012',
'013','014','015','016','017','018','019','020','021','022','023','024','025',
'026','027','028','029','030']

lf_z = []
bl_z = []
lfs = []
bls = []
years_lf =[]
years_lf3 = []
years_bl = []
years_anomlf = []
years_anombl = []
for version in versions:
    years = 'b.e11.B20TRC5CNBDRD.f09_g16.%s.cam.h1.SI-x.1920-2005.nc' % version
    filename = directory + years
    values = Dataset(filename)
    lon = values.variables['lon'][207:229]
    lat = values.variables['lat'][12:27]
    leaf_index = values.variables['leaf_index'][:,12:27,207:229] 
    bloom_index = values.variables['bloom_index'][:,12:27,207:229]
    values.close()
    
    lfmean = []
    blmean = []
    for i in xrange(len(leaf_index)):
        leafs = np.nanmean(leaf_index[i,:,:])
        blooms = np.nanmean(bloom_index[i,:,:])
        lfmean.append(leafs)
        blmean.append(blooms)
    lfs.append(lfmean)
    bls.append(blmean)
    
    lf_totalmean = np.nanmean(lfmean)
    bl_totalmean = np.nanmean(blmean)
    lf_anom = []
    bl_anom = []
    for j in xrange(len(lfmean)):
        anomlf = lf_totalmean - lfmean[j]
        anombl = bl_totalmean - blmean[j]
        lf_anom.append(anomlf)
        bl_anom.append(anombl)
    years_anomlf.append(lf_anom)
    years_anombl.append(bl_anom)
    
    lf_zscore = zscore(lfmean)
    bl_zscore = zscore(blmean)
    lf_z.append(lf_zscore)
    bl_z.append(bl_zscore)
    
    yrlf = np.where(lf_zscore<=-2)
    yrbl = np.where(bl_zscore<=-2)
    yrlf3 = np.where(lf_zscore<=-2.8)
    years_lf.append(yrlf)
    years_bl.append(yrbl)
    years_lf3.append(yrlf3)

### Calculate total mean
lf_totalmean = np.nanmean(lfs)
bl_totalmean = np.nanmean(bls) 

### Detrend data
lf_slope=[]
bl_slope=[]
lfs_notrend = []
bls_notrend = []
lfz_notrend = []
blz_notrend = []
years_lf2_notrend = []
years_bl2_notrend = []
years_lf3_notrend = []
years_bl3_notrend = []    
trendleaf = [] 
for i in xrange(len(lfs)):
    timeq = np.asarray(list(xrange(len(lfs[i]))))
    lfslope, lfintercept, r_value, p_value, std_err = sts.linregress(timeq,lfs[i])
    lf_line = lfslope*timeq+lfintercept
    blslope, blintercept, r_value, p_value, std_err = sts.linregress(timeq,bls[i])
    bl_line = blslope*timeq+blintercept
    
    lf_slope.append(lfslope)
    bl_slope.append(blslope)
       
    lf_diff = lf_line-lf_totalmean
    bl_diff = bl_line-bl_totalmean

    lfs_no = lfs[i]-lf_diff
    bls_no = bls[i]-bl_diff
    lfs_notrend.append(lfs_no)
    bls_notrend.append(bls_no)
    
    #lf_zscore_no = zscore(lfs_no)
    lfave = 96.856182177083653
    lfstd = 5.4175636772359557
    lf_zscore_no = (lfs_no-lfave)/lfstd
    bl_zscore_no = zscore(bls_no)
    lfz_notrend.append(lf_zscore_no)
    blz_notrend.append(bl_zscore_no)
    
    yrlf_no = np.where(lf_zscore_no<=-2)
    yrbl_no = np.where(bl_zscore_no<=-2)
    yrlf3_no = np.where(lf_zscore_no<=-2.8)
    yrbl3_no = np.where(bl_zscore_no<=-2.8)
    years_lf2_notrend.append(yrlf_no)
    years_bl2_notrend.append(yrbl_no)
    years_lf3_notrend.append(yrlf3_no)
    years_bl3_notrend.append(yrbl3_no)

    trendleaf.append(lfslope*10.)
### Create plot of leaf indices
#fig = plt.figure()
#threshold1 = [81.5]*len(lfmean)
#threshold2 = [88.]*len(lfmean)
#threshold3 = [lf_totalmean]*len(lfmean)
#plt.plot(lfs[-2],color='DarkOrange',linewidth=1,label='Leaf Dates')
#plt.plot(lfs_notrend[4],color='b',linewidth=1,label='Leaf No trend')
#plt.plot(threshold1,'r--',linewidth=3,label='-3SD Threshold')
#plt.plot(threshold2,'r--',linewidth=3,label='-2SD Threshold')
#plt.plot(threshold3,'k--',linewidth=3,label='Average Leaf')
#plt.legend(loc=1,prop={'size':9},shadow=True)
#plt.xlabel('Years')
#plt.grid(True)
#fig.suptitle('Leaf Index for CESM-LE (1920-2005)',fontsize=20)
#fig.text(0.04,0.5,'Day of the Year (doy)',va='center',rotation='vertical')
#plt.savefig('/volumes/zml5/research/CCSM/results/CESMlE_LFdates.png',dpi=300)

### Create Bar Graph No trend
#n_groups = len(lfs)
#
#thresholdlf = []
#thresholdlf3 = []
#for i in xrange(len(lfs)):
#    sizelf = np.size(years_lf2_notrend[i])
#    sizelf3 = np.size(years_lf3_notrend[i])
#    thresholdlf.append(sizelf)
#    thresholdlf3.append(sizelf3)
#    
#fig, ax = plt.subplots()
#index = np.arange(n_groups)
#bar_width = 0.35
#
#rect1 = plt.bar(index,thresholdlf,bar_width,color='darkgreen',
#                label='Leaf -2SD',alpha=0.85)
#rect2 = plt.bar(index+bar_width,thresholdlf3,bar_width,color='indianred',
#                label='Leaf -3SD',alpha=0.85)
#plt.xlabel('CESM Individual Ensemble Projections')
#plt.ylabel('Occurences')
#fig.suptitle('Detrended Number of Anomalously Early Springs on CESM-LE (2006-2080)',fontsize=15)
#plt.yticks(np.arange(0,5,1.0),('0','1','2','3','4'))
#plt.xticks(index+bar_width,('02','','04','','06','','08'))
## ('02','','04','','06','','08','','10','','12','','14','','16','','18','','20','','22','','24','','26','','28','','30'))
#plt.legend(loc=1,prop={'size':9},shadow=True)
#plt.grid(True)
#plt.savefig('/volumes/zml5/research/CCSM/results/CESMlE_future_earlyspring_thresholds.png',dpi=300)

### Create Bar Graph with Trend
n_groups = len(lfs)

thresholdlf = []
thresholdbl = []
thresholdlf3 = []
for i in xrange(len(lfs)):
    sizelf = np.size(years_lf2_notrend[i])
    sizelf3 = np.size(years_lf3_notrend[i])
    sizebl = np.size(years_bl[i])
    thresholdlf.append(sizelf)
    thresholdbl.append(sizebl)
    thresholdlf3.append(sizelf3)
    
fig, ax = plt.subplots()
index = np.arange(n_groups)
bar_width = 0.35

rect1 = plt.bar(index,thresholdlf,bar_width,color='darkgreen',
                label='Leaf -2SD',alpha=0.85)
rect2 = plt.bar(index+bar_width,thresholdlf3,bar_width,color='indianred',
                label='Leaf -3SD',alpha=0.85)
plt.xlabel('CESM Individual Ensemble Projections')
plt.ylabel('Occurences')
fig.suptitle('Number of Anomalously Early Springs on CESM-LE (1920-2005)',fontsize=15)
plt.yticks(np.arange(0,4,1.0),('0','1','2','3'))
plt.xticks(index+bar_width,('02','','04','','06','','08','','10','','12','','14','','16','','18','','20','','22','','24','','26','','28','','30'))
plt.legend(loc=1,prop={'size':9},shadow=True)
plt.grid(True)
#plt.savefig('/volumes/zml5/research/CCSM/results/CESMlE_earlyspring_thresholds.png',dpi=300)

np.savetxt('/volumes/zml5/research/ccsm/text/1920-2005trends.txt',(lf_slope,bl_slope),delimiter=',')
np.savetxt('/volumes/zml5/research/ccsm/text/1920-2005lf.txt',(lfs),delimiter=',')
np.savetxt('/volumes/zml5/research/ccsm/text/1920-2005bl.txt',(bls),delimiter=',')