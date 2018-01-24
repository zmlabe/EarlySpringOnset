"""
*Functions gather the earliest springs from each member of the historical and future LENS*
"""

import numpy as np
from netCDF4 import Dataset
import scipy.stats as sts
import json

def SIxHistorical():
    """
    Calcules earliest springs from historical LENS

    
    Returns
    ----------
    earlyh : doy of earliest first leafs
    yrsh : years of earliest first leafs from each ensemble member
    pathsh : historical LENS path
    """
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/data/'
    
    versions=['002','003','004','005','006','007','008','009','010','011','012','013','014','015','016','017','018','019','020','021','022','023','024','025','026','027','028','029','030']
    leaf=[]
    pathsh=[]
    for version in versions:
        years = 'b.e11.B20TRC5CNBDRD.f09_g16.%s.cam.h1.SI-x.1920-2005.nc' % version
        path = 'B20TRC5CNBDRD.f09_g16.%s' % version
        filename = directory + years
        values = Dataset(filename)
        lon = values.variables['lon'][207:229]
        lat = values.variables['lat'][12:27]
        leaf_index = values.variables['leaf_index'][:,12:27,207:229]
        values.close()
        
        leaf.append(leaf_index)
        pathsh.append(path)
    leafs = np.asarray(leaf)
    
    lfmean = np.empty((leafs.shape[0],leafs.shape[1]))
    for i in xrange(leafs.shape[0]):
        for j in xrange(leafs.shape[1]):
            lfmean[i,j] = np.nanmean(leafs[i,j,:,:])
            
    yrsh = np.empty((leafs.shape[0]))
    earlyh = np.empty((leafs.shape[0]))
    for i in xrange(leafs.shape[0]):       
        yrshh = np.where((lfmean[i] >= 74) & (lfmean[i] <= 90))[0][-1]
        yrsh[i] = yrshh
        earlyh[i] = lfmean[i,yrshh]
    
#    earlyh = np.empty((leafs.shape[0]))
#    yrsh = np.empty((leafs.shape[0]))
#    for i in xrange(leafs.shape [0]):
#        early = np.nanmin(lfmean[i,74:88])
#        earlyh[i] = early
#        yrsh[i] = np.where(lfmean[i,74:88] == np.nanmin(lfmean[i,74:88]))[0][-1]
        
    return earlyh,yrsh,pathsh
earlyh,yrsh,pathsh = SIxHistorical()

def SIxFuture():
    """
    Calcules earliest springs from future LENS

    
    Returns
    ----------
    earlyf : doy of earliest first leafs
    yrsf : years of earliest first leafs from each ensemble member
    pathsf : future LENS path
    """
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/data/'
    versions=['002','003','004','005','006','007','008','009','010','011','012','013','014','015','016','017','018','019','020','021','022','023','024','025','026','027','028','029','030']
    
    leaf=[]
    pathsf=[]
    for version in versions:
        years = 'b.e11.BRCP85C5CNBDRD.f09_g16.%s.cam.h.SI-x.2006-2080.nc' % version
        filename = directory + years
        path = 'B20TRC5CNBDRD.f09_g16.%s' % version
        values = Dataset(filename)
        leaf_index = values.variables['leaf_index'][:,12:27,207:229]
        values.close()
        
        pathsf.append(path)
        leaf.append(leaf_index)
    leafs = np.asarray(leaf)
    
    lfmean = np.empty((leafs.shape[0],leafs.shape[1]))
    for i in xrange(leafs.shape[0]):
        for j in xrange(leafs.shape[1]):
            lfmean[i,j] = np.nanmean(leafs[i,j,:,:])
            
    yrsf = np.empty((leafs.shape[0]))
    earlyf = np.empty((leafs.shape[0]))
    for i in xrange(leafs.shape[0]):       
        yrsff = np.where((lfmean[i] >= 74) & (lfmean[i] <= 88))[0][-1]
        yrsf[i] = yrsff
        earlyf[i] = lfmean[i,yrsff]
            
    #earlyf = np.empty((leafs.shape[0]))
    #yrsf = np.empty((leafs.shape[0]))
    #for i in xrange(leafs.shape [0]):
    #    early = np.nanmin(lfmean[i,74:88]) 
    #    earlyf[i] = early
    #    yrsf[i] = np.where(lfmean[i,74:88] == np.nanmin(lfmean[i,74:88]))[0][0]
    return earlyf,yrsf,pathsf
#earlyf,yrsf,pathsf = SIxFuture()

### Round Values to DOY
#earlyh = np.round(earlyh)
#earlyf = np.round(earlyf)
#dh = {'id':pathsh,'year':list(yrsh),'doy':list(earlyh)}
#json_hist = json.dumps(dh)
#df ={'id':pathsf,'year':list(yrsf),'doy':list(earlyf)}
#json_fut = json.dumps(df)

### Create Json files
#obj1 = open('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/json_hists.txt','w')
#obj1.write(json_hist)
#obj1.close()
#
#obj2 = open('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/json_futs.txt','w')
#obj2.write(json_fut)
#obj2.close()