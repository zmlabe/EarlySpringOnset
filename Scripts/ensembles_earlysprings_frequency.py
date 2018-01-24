"""
*Functions calculate frequency of early springs in LENS and Control*
"""

import numpy as np
from netCDF4 import Dataset
from scipy.stats import zscore
import matplotlib.pyplot as plt
import scipy.stats as sts

def Future():
    """
    Function reads in future LENS SI-x data

    
    Returns
    ----------
    years_lf2_notrend : years with early springs above 2std
    years_lf3_notrend : years with early springs above 3std
    lfs : actual SIx values (yr x lat x lon)
    """
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/data/'
    
    versions=['002','003','004','005','006','007','008','009','010','011','012','013','014','015','016','017','018','019','020','021','022','023','024','025','026','027','028','029','030']
    
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
        years = 'b.e11.BRCP85C5CNBDRD.f09_g16.%s.cam.h.SI-x.2006-2080.nc' % version
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
    lfslopes = []
    blslopes = []
    lfs_notrend = []
    bls_notrend = []
    lfz_notrend = []
    blz_notrend = []
    years_lf2_notrend = []
    years_bl2_notrend = []
    years_lf3_notrend = []
    years_bl3_notrend = []    
    for i in xrange(len(lfs)):
        timeq = np.asarray(list(xrange(len(lfs[i]))))
        lfslope, lfintercept, r_value, p_value, std_err = sts.linregress(timeq,lfs[i])
        lf_line = lfslope*timeq+lfintercept
        blslope, blintercept, r_value, p_value, std_err = sts.linregress(timeq,bls[i])
        bl_line = blslope*timeq+blintercept
           
        lf_diff = lf_line-lf_totalmean
        bl_diff = bl_line-bl_totalmean
        
        lfslopes.append(lfslope*10.)
        blslopes.append(blslope*10.)
    
        lfs_no = lfs[i]-lf_diff
        bls_no = bls[i]-bl_diff
        lfs_notrend.append(lfs_no)
        bls_notrend.append(bls_no)
        
    #    lf_zscore_no = zscore(lfs_no)
        bl_zscore_no = zscore(bls_no)
        lfave = 98.562333
        lfstd = 5.483669
        lf_zscore_no = (lfs_no-lfave)/lfstd
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
        
    return years_lf2_notrend, years_lf3_notrend, lfs

def Historical():
    """
    Function reads in historical LENS SI-x data

    
    Returns
    ----------
    years_lf2_notrend : years with early springs above 2std
    years_lf3_notrend : years with early springs above 3std
    lfs : actual SIx values (yr x lat x lon)
    """
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
    return years_lf2_notrend, years_lf3_notrend, lfs

def Control():
    """
    Function reads in CESM control SI-x data

    
    Returns
    ----------
    totalz2 : years with early springs above 2std
    totalz3 : years with early springs above 3std
    """
    directory= '/volumes/eas-shared/ault/ecrl/spring-indices/data/'
    
    ### CESM control for Great Lakes region
    years = 'b.e11.B1850C5CN.f09_g16.005.cam.h1.TRE.SI-x.402-999.nc'
    filename = directory + years
    values = Dataset(filename)
    lon = values.variables['lon'][207:229]
    lat = values.variables['lat'][12:27]
    leaf_index = values.variables['leaf_index'][:,12:27,207:229]
    values.close()
    
    ### Average SI per year
    lfmean = []
    for i in xrange(len(leaf_index)):
        leafs = np.nanmean(leaf_index[i,:,:])
        lfmean.append(leafs)
        
    ### Climatology
    lf_totalmean = np.nanmean(lfmean)
    
    ### Anom per year
    lf_anom = []
    for j in xrange(len(lfmean)):
        anomlf = lf_totalmean - lfmean[j]
        lf_anom.append(anomlf)
    
    ### z-score per year
    lf_zscorez = zscore(lfmean)
    yrlfz2_1 = np.where(lf_zscorez[:98]<=-2)[0]
    yrlfz3_1 = np.where(lf_zscorez[:98]<=-3)[0]
    
    yrs2_1 = np.size(yrlfz2_1)
    yrs3_1 = np.size(yrlfz3_1)
     
    lfslicez = lf_zscorez[98:] 
     
    yrlfz2 = []
    yrlfz3 = []
    for i in xrange(0,len(lf_zscorez[98:]),100):
        start = i
        end = i + 100
        q = np.arange(start,end)
        yr2 = np.where(lfslicez[q]<=-2)[0]
        yr2f = np.size(yr2)
        yr3 = np.where(lfslicez[q]<=-3)[0]
        yr3f = np.size(yr3)
        yrlfz2.append(yr2f)
        yrlfz3.append(yr3f)
    
    totalz2 = np.append(yrs2_1,np.asarray(yrlfz2))
    totalz3 = np.append(yrs3_1,np.asarray(yrlfz3))
    return totalz2, totalz3    
    
years_lf2_notrend, years_lf3_notrend, lfs = Future()
years_lf2_notrendH, years_lf3_notrendH, lfsH = Historical()
totalz2, totalz3 = Control()

# Create Bar Graph 
fig = plt.figure()
n_groups = len(lfs)
n_groups1 = len(totalz2)

thresholdlf = []
thresholdlf3 = []
for i in xrange(len(lfs)):
    size2 = np.size(years_lf2_notrend[i])
    size3 = np.size(years_lf3_notrend[i])
    
    sizelf = (float(size2) / len(lfs[0])) * 100 
    sizelf3 = (float(size3) / len(lfs[0])) * 100    
    
    thresholdlf.append(sizelf)
    thresholdlf3.append(sizelf3)

thresholdlfH = []
thresholdlf3H = []
for i in xrange(len(lfsH)):
    size2H = np.size(years_lf2_notrendH[i])
    size3H = np.size(years_lf3_notrendH[i])
    
    sizelfH = (float(size2H) / len(lfsH[0])) * 100 
    sizelf3H = (float(size3H) / len(lfsH[0])) * 100    
    
    thresholdlfH.append(sizelfH)
    thresholdlf3H.append(sizelf3H)

fig = plt.figure()
fig.suptitle('CESM-LE Frequency of Anomalously Early Springs',fontsize=16)
index = np.arange(n_groups)
index1 = np.arange(n_groups1)
bar_width = 0.35

ax2 = plt.subplot(3,1,1)
rect1 = plt.bar(index1,totalz2,bar_width,color='darkgreen',
                label='Leaf -2$\sigma$',alpha=0.85)
rect2 = plt.bar(index1+bar_width,totalz3,bar_width,color='indianred',
                label='Leaf -3$\sigma$',alpha=0.85)
#plt.xlabel('Century',fontsize=8)
ax2.text(6.04,41.56,'Control \nCenturies',fontsize=10,rotation='horizontal')
plt.ylabel('Frequency')
plt.xticks(index+bar_width,('400','500','600','700','800','900'))
plt.ylim([0,75])
plt.xlim([0,6])
plt.yticks(np.arange(0,76,25))

plt.legend(loc=1,prop={'size':10},shadow=True)
plt.grid(True)

ax1 = plt.subplot(3,1,2)
rect1 = plt.bar(index,thresholdlfH,bar_width,color='darkgreen',
                label='Leaf -2$\sigma$',alpha=0.85)
rect2 = plt.bar(index+bar_width,thresholdlf3H,bar_width,color='indianred',
                label='Leaf -3$\sigma$',alpha=0.85)
#plt.xlabel('LENS Members',fontsize=8)
ax1.text(29.2,41,'Historical \nMembers \n1920-2005',fontsize=10,rotation='horizontal')
plt.ylabel('Frequency')
plt.ylim([0,75])
plt.xlim([0,29])
plt.yticks(np.arange(0,76,25))

#plt.legend(loc=1,prop={'size':9},shadow=True)
plt.grid(True)

ax4 =plt.subplot(3,1,3,sharex=ax1)
rect1 = plt.bar(index,thresholdlf,bar_width,color='darkgreen',
                label='Leaf -2SD',alpha=0.85)
rect2 = plt.bar(index+bar_width,thresholdlf3,bar_width,color='indianred',
                label='Leaf -3SD',alpha=0.85)
#plt.xlabel('LENS Members',fontsize=8)
ax4.text(29.2,41,'Future \nMembers \n2006-2080',fontsize=10,rotation='horizontal')
plt.ylabel('Frequency')
plt.xticks(index+bar_width,('02','','04','','06','','08','','10','','12','','14','','16','','18','','20','','22','','24','','26','','28','','30'))
plt.ylim([0,75])
plt.yticks(np.arange(0,76,25))
plt.xlim([0,29])
#plt.legend(loc=1,prop={'size':8},shadow=True)
plt.grid(True)
plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/total_CESMstdsprings.eps',dpi=400,format='eps')
