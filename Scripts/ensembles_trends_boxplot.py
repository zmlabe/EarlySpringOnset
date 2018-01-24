"""
*Script creates box and whisker plot for SIx LENS/CESM control trends*
"""

import ensembles_SIxtrend_datareader as T
import numpy as np
import scipy.stats as sts
import matplotlib.pyplot as plt

### Import LENS Trends
lstfrz_trendh,damage_trendh,leaf_trendh, lat, lon = T.histtrends()
lstfrz_trendh2,leaf_trendh2 = T.histtrends2()
lstfrz_trendf,damage_trendf,leaf_trendf, lat, lon = T.futtrends()
lstfrz_trendf2,leaf_trendf2 = T.futtrends2()

print 'FINISHED!'

### Calculate Mean Trends
lstfrz_meanh = sts.nanmean(lstfrz_trendh)
damage_meanh = sts.nanmean(damage_trendh)
leaf_meanh = sts.nanmean(leaf_trendh)
lstfrz_meanh2 = sts.nanmean(lstfrz_trendh2)
leaf_meanh2 = sts.nanmean(leaf_trendh2)

lstfrz_meanf = sts.nanmean(lstfrz_trendf)
damage_meanf = sts.nanmean(damage_trendf)
leaf_meanf = sts.nanmean(leaf_trendf)
lstfrz_meanf2 = sts.nanmean(lstfrz_trendf2)
leaf_meanf2 = sts.nanmean(leaf_trendf2)

### Covariance Function
def Cov(a,b):
    """
    Calculates covariance

    
    Parameters
    ----------
    a : dependent variable
    b : independent variable
    
    Returns
    ----------
    cov : covariance array
    """
    a_mean = np.nanmean(a)
    b_mean = np.nanmean(b)
    
    e = np.empty((len(a)))
    for i in xrange(0, len(a)):
        c = (a[i] - a_mean)
        d = (b[i] - b_mean)
        e[i] = c*d
        
    cov = np.nansum(e)/(len(a)-1)
    return cov

### Calculate Correlation Coefficents
def ensCoeff(gridded,mean):
    """
    Calculates correlation coefficients

    
    Parameters
    ----------
    gridded : dependent variable
    mean : independent variable
    
    Returns
    ----------
    corr : correlation coefficent
    """
    corr = np.empty((gridded.shape[0]))
    for i in xrange(gridded.shape[0]):
        griddedn = np.ravel(gridded[i,:,:])
        meann = np.ravel(mean)
        
        cov = Cov(griddedn,meann)
        
        std1 = np.nanstd(griddedn)
        std2 = np.nanstd(meann)
        
        corr[i] = cov/(std1*std2)
    return corr

### Coefficents
leafh_coef1 = ensCoeff(leaf_trendh,leaf_meanh)
leafh_coef2 = ensCoeff(leaf_trendh2,leaf_meanh2)

leaff_coef1 = ensCoeff(leaf_trendf,leaf_meanf)
leaff_coef2 = ensCoeff(leaf_trendf2,leaf_meanf2)

dataq=[leafh_coef2,leafh_coef1,leaff_coef2,leaff_coef1]

### Plot correlation coefficient box plot for each LENS time period
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
ax.set_xticklabels(['1970-2005','1920-2005','2006-2040','2006-2080'])
plt.xlabel('Years',fontsize=13)
plt.ylabel('Correlation Coefficient',fontsize=13)
fig.suptitle('First Leaf Trend Correlation Coefficients',fontsize=18)
plt.title('LENS and Ensemble Mean',fontsize=13)
plt.grid(True)
plt.savefig('/Users/zlabe/documents/CESMspring/Fig17.eps',dpi=400,format='eps')