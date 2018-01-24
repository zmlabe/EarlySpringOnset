"""
*Script calculates and plots correlation coefficients between SIx and climate modes for historical/future LENS*
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import scipy.stats as sts
from CESM_2006_2080_lstfrz import SIx
from CESM_1920_2005_lstfrz import SIxHistorical
from mpl_toolkits.basemap import Basemap

def readData(period):
    """
    Reads in climate mode data for future and historical LENS
    
    
    Parameters
    ----------
    period : future or historical (string)
    
    Returns
    ----------
    PDOyr : PDO array for each year for each ensemble
    PDOave : average PDO for each ensemble
    NAOyr : NAO array for each year for each ensemble
    NAOave : average NAO for each ensemble
    NINOyr : NINO array for each year for each ensemble
    NINOave : average NINO for each ensemble
    leafmean : average leaf for each ensemble
    latmean : average last freeze for each ensemble
    lat : array of latitudes
    lon : array of longitudes
    """
    if period == 'future':
        directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/cesm1.lens.1920-2005.cvdp_data/'
        NAO = []
        PDO = []
        NINO = []
        ens = list(xrange(2,31))
        for i in xrange(len(ens)):
            files = 'CESM1-CAM5-BGC-LE_%s.cvdp_data.2013-2100.nc' % ens[i]
            filename = directory + files
            values = Dataset(filename)
            time = values.variables['time'][:]
            pdo = values.variables['pdo_timeseries_mon'][:]
            nao = values.variables['nao_pc_mon'][:]
            nino = values.variables['nino34'][:]
            values.close()
            
            NAO.append(nao)
            PDO.append(pdo)
            NINO.append(nino)
        time = np.asarray(time)
        PDO = np.asarray(PDO)
        NINO = np.asarray(NINO)
        NAO = np.asarray(NAO)
        PDOyr = np.reshape(PDO,(PDO.shape[0],PDO.shape[1]/12.,12.))
        PDOave = np.nanmean(PDOyr,axis=2)
        NAOyr = np.reshape(NAO,(NAO.shape[0],NAO.shape[1]/12.,12.))
        NAOave = np.nanmean(NAOyr,axis=2)
        NINOyr = np.reshape(NINO,(NINO.shape[0],NINO.shape[1]/12.,12.))
        NINOave = np.nanmean(NINOyr,axis=2)
        
        leafmean, latmean, lstfrz, lat, lon = SIx()        
        leafmean = leafmean[:,7:,:,:]
        latmean = latmean[:,7:,:,:]
        PDOave = PDOave[:,:-20]
        NAOave = NAOave[:,:-20]
        NINOave = NINOave[:,:-20]
        return PDOyr,PDOave,NAOyr,NAOave,NINOyr,NINOave,leafmean,latmean,lat,lon
    elif period == 'historical':
        directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/cesm1.lens.1920-2005.cvdp_data/'
        NAO = []
        PDO = []
        NINO = []
        ens = list(xrange(2,31))
        for i in xrange(len(ens)):
            files = 'CESM1-CAM5-BGC-LE_%s.cvdp_data.1920-2005.nc' % ens[i]
            filename = directory + files
            values = Dataset(filename)
            time = values.variables['time'][:]
            pdo = values.variables['pdo_timeseries_mon'][:]
            nao = values.variables['nao_pc_mon'][:]
            nino = values.variables['nino34'][:]
            values.close()
            
            NAO.append(nao)
            PDO.append(pdo)
            NINO.append(nino)
        time = np.asarray(time)
        PDO = np.asarray(PDO)
        NINO = np.asarray(NINO)
        NAO = np.asarray(NAO)
        PDOyr = np.reshape(PDO,(PDO.shape[0],PDO.shape[1]/12.,12.))
        PDOave = np.nanmean(PDOyr,axis=2)
        NAOyr = np.reshape(NAO,(NAO.shape[0],NAO.shape[1]/12.,12.))
        NAOave = np.nanmean(NAOyr,axis=2)
        NINOyr = np.reshape(NINO,(NINO.shape[0],NINO.shape[1]/12.,12.))
        NINOave = np.nanmean(NINOyr,axis=2)
        
        leafmean, latmean, lat, lon = SIxHistorical()
        return PDOyr,PDOave,NAOyr,NAOave,NINOyr,NINOave,leafmean,latmean,lat,lon

#PDOyr,PDOave,NAOyr,NAOave,NINOyr,NINOave,leafmean,latmean,lat,lon = readData('future')
#PDOyr,PDOave,NAOyr,NAOave,NINOyr,NINOave,leafmean,latmean,lat,lon = readData('historical')

def Corr(x,y):
    """
    Calculates pearson correlation coefficent between two variables
    
    
    Parameters
    ----------
    x : dependent variable
    y : independent variable
    
    Returns
    ----------
    cocoeff1 : pearson correlation coefficient 
    cocoeff2 : not important correlation coefficient 
    """
    
    cocoeff1 = np.empty((x.shape[0],y.shape[2],y.shape[3]))
    cocoeff2 = np.empty((x.shape[0],y.shape[2],y.shape[3]))
    for ens in xrange(x.shape[0]):
        for i in xrange(y.shape[2]):
            for j in xrange(y.shape[3]):
                cocoeff1[ens,i,j],cocoeff2[ens,i,j] = sts.pearsonr(x[ens,:],y[ens,:,i,j])
                
    return cocoeff1, cocoeff2

cocoeff1, cocoeff2 = Corr(PDOave,leafmean)

def PlotCorr(cocoeff,lat,lon):
    """
    Plots correlation coefficients on a grid
    
    
    Parameters
    ----------
    cocoeff1 : 2d array of correlation coefficients
    lat : array of latitudes
    lon : array of longitudes
    """

    lons,lats = np.meshgrid(lon,lat)
    
    def plot_rec(bmap, lonmin,lonmax,latmin,latmax):
        xs = [lonmin,lonmax,lonmax,lonmin,lonmin]
        ys = [latmin,latmin,latmax,latmax,latmin]
        bmap.plot(xs, ys, latlon = True, color='k',linewidth=1.5,linestyle='solid')
    lonmin = -101.5
    lonmax =  -75.5
    latmin =  37.5
    latmax =  50.5
    
    meancoef = sts.nanmean(cocoeff1)
    
    member = list(xrange(1,30))
    ### Plot Coefficients
    fig = plt.figure()   
    plt.title('Mean LENS 2013-2080, Leaf Index PDO Correlations',fontsize=16)
    ax1 = fig.add_subplot(1,1,1)
    m = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=31,urcrnrlon=298,
                urcrnrlat=54,resolution='l')           
    m.drawstates()
    m.drawcountries()
    m.drawmapboundary(fill_color = 'white')
    m.drawcoastlines(color='black',linewidth=0.5)
    m.drawlsmask(land_color='grey',ocean_color='w')
    x,y = m(lons,lats)
    
    meancoef[np.where(meancoef > 1)] = 1
    meancoef[np.where(meancoef < -1)] = -1
    
    cs = m.contourf(x,y,meancoef,np.arange(-1.,1.1,.05))
    plot_rec(m,lonmin,lonmax,latmin,latmax)
    cs.set_cmap('RdYlBu')
    
    cbar = m.colorbar(cs,location='bottom',pad='5%',ticks=np.arange(-1.,1.2,0.2))
    cbar.set_label('Correlation Coefficient') 
    
    plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/meanLF_PDOcorr.eps',dpi=400,format='eps')
            
PlotCorr(cocoeff1,lat,lon)    
