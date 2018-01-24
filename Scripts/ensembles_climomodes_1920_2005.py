"""
*Script plots climate modes and trends for historical LENS*
"""
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
#import lens_SIxtrends_19202080 as L
import scipy.stats as sts
from mpl_toolkits.basemap import Basemap

def SIxHistorical():
    """
    Reads in SIx for historical LENS


    Returns
    ----------
    leafmean : average leaf for each ensemble
    latmean : average last freeze for each ensemble
    lat : array of latitudes
    lon : array of longitudes
    """
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/data/'
    
    versions=['002','003','004','005','006','007','008','009','010','011','012','013','014','015','016','017','018','019','020','021','022','023','024','025','026','027','028','029','030']
    
    lstfrz = []
    leaf=[]
    for version in versions:
        years = 'b.e11.B20TRC5CNBDRD.f09_g16.%s.cam.h1.SI-x.1920-2005.nc' % version
        filename = directory + years
        values = Dataset(filename)
        lon = values.variables['lon'][189:240]
        lat = values.variables['lat'][:32]
        lstfrz_index = values.variables['lstfrz_index'][:,:32,189:240]
        leaf_index = values.variables['leaf_index'][:,:32,189:240]
        values.close()
        print version
        lstfrz.append(lstfrz_index)
        leaf.append(leaf_index)
        
    latmean = np.asarray(lstfrz)
    leafmean = np.asarray(leaf)
    print 'DONE 2'
    return leafmean,latmean,lat,lon

directory2 = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/cesm1.lens.1920-2005.cvdp_data/'

NAO = []
PDO = []
NINO = []
ens = list(xrange(2,31))
for i in xrange(len(ens)):
    files = 'CESM1-CAM5-BGC-LE_#%s.cvdp_data.1920-2005.nc' % ens[i]
#    files = 'CESM1-CAM5-BGC-LE_%s.cvdp_data.2013-2100.nc' % ens[i]
    filename = directory2 + files
    values = Dataset(filename)
    time = values.variables['time'][:]
    pdo = values.variables['pdo_timeseries_mon'][:]
    nao = values.variables['nao_pc_mon'][:]
    nino = values.variables['nino34'][:]
    values.close()
    
    NAO.append(nao)
    PDO.append(pdo)
    NINO.append(nino)
PDO = np.asarray(PDO)
NAO = np.asarray(NAO)
time = np.asarray(time)
NINO = np.asarray(NINO)

print 'DONE 1'

#fig = plt.figure()
#fig.suptitle('LENS 1920-2005 Monthly PDO and NAO',fontsize=16)
#for i in xrange(len(ens)):
#    ax = plt.subplot(6,5,i+1)
#    plt.plot(NAO[i,:],color='b',linewidth=0.5,label='NAO')
#    plt.plot(PDO[i,:],color='r',linewidth=1.5,label='PDO',linestyle='solid')
#    plt.ylim([-3,3])
#    plt.yticks([-3,-2,-1,0,1,2,3],fontsize=6)
#    plt.grid(True)
#    x=list(xrange(0,PDO.shape[1]+1,258))
#    labels = ['1920','1941','1962','1983','2004']
#    plt.xticks(x,labels,fontsize=6)
#    plt.xlim([0,PDO.shape[1]])
##    plt.legend(loc=1,prop={'size':3},shadow=True)
#plt.savefig('/volumes/zml5/research/ccsm/results/CESM_ensembles/LE_historical/NAOPDO.png',dpi=300)
##fig.clear()
#
#fig = plt.figure()
#fig.suptitle('LENS 1920-2005 Monthly Nino3.4',fontsize=16)
#for i in xrange(len(ens)):
#    ax = plt.subplot(6,5,i+1)
#    plt.plot(NINO[i,:],color='darkorange',linewidth=1.2,label='ENSO')
#    plt.ylim([-5,5])
#    plt.yticks([-5,-3,0,3,5],fontsize=6)
#    plt.grid(True)
#    x=list(xrange(0,PDO.shape[1]+1,258))
#    labels = ['1920','1941','1962','1983','2004']
#    plt.xticks(x,labels,fontsize=6)
#    plt.xlim([0,PDO.shape[1]])
##   plt.legend(loc=1,prop={'size':3},shadow=True)
#plt.savefig('/volumes/zml5/research/ccsm/results/CESM_ensembles/LE_historical/Nino34.png',dpi=300)
#fig.clear()

### Create Yearly Averages
PDOyr = np.reshape(PDO,(PDO.shape[0],PDO.shape[1]/12.,12.))
PDOave = np.nanmean(PDOyr,axis=2)

NAOyr = np.reshape(NAO,(NAO.shape[0],NAO.shape[1]/12.,12.))
NAOave = np.nanmean(NAOyr,axis=2)

NINOyr = np.reshape(NINO,(NINO.shape[0],NINO.shape[1]/12.,12.))
NINOave = np.nanmean(NINOyr,axis=2)

### Import SIx
#lstfrz_trendh,damage_trendh,leaf_trendh, lat, lon = L.histtrends()
#lstfrz_trendf,damage_trendf,leaf_trendf, lat, lon = L.futtrends()

#leafmean, latmean, lstfrz, lat, lon = SIx()
#leafmean, latmean, lat, lon = SIxHistorical()

leafmean = leafmean[:,:,:,:]
latmean = latmean[:,:,:,:]
PDOave = PDOave[:,:]
NAOave = NAOave[:,:]
NINOave = NINOave[:,:]

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
    
#cocoeff1,cocoeff2 = Corr(PDOave,leafmean)

def PlotCorr(cocoeff1,lat,lon):
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
#    fig.suptitle('LENS 1920-2005, Leaf Index PDO Correlations',fontsize=10)
    ax1 = plt.subplot(6,5,1)
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
    
    ax1.spines['top'].set_linewidth(3)
    ax1.spines['right'].set_linewidth(3)
    ax1.spines['bottom'].set_linewidth(3)
    ax1.spines['left'].set_linewidth(3)
    
    ax1.text(0.29,0.015,'Average Correlation',size='8',horizontalalignment= 'center',
            backgroundcolor='white',verticalalignment= 'center',
            bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
            transform=ax1.transAxes) 
            
    for i in xrange(len(cocoeff1)):
        ax = plt.subplot(6,5,i+2)
        m = Basemap(projection='merc',llcrnrlon=236,llcrnrlat=31,urcrnrlon=298,
                    urcrnrlat=54,resolution='l')           
        m.drawstates()
        m.drawcountries()
        m.drawmapboundary(fill_color = 'white')
        m.drawcoastlines(color='black',linewidth=0.5)
        m.drawlsmask(land_color='grey',ocean_color='w')
        x,y = m(lons,lats)
    
        cocoeff1[np.where(cocoeff1 > 1)] = 1
        cocoeff1[np.where(cocoeff1 < -1)] = -1
        
        cs = m.contourf(x,y,cocoeff1[i,:,:],np.arange(-1.,1.1,.1))
        cs.set_cmap('RdYlBu')
    
        ax.text(0.16,0.015,'Member %i' % (member[i]+1),size='8',horizontalalignment= 'center',
                backgroundcolor='white',verticalalignment= 'center',
                bbox=dict(facecolor='white',edgecolor='black',alpha=0.9),
                transform=ax.transAxes)    
    plt.tight_layout()
    fig.subplots_adjust(bottom=0.098)
    cbar_ax = fig.add_axes([0.15, 0.08, 0.7, 0.01])
    cbar = fig.colorbar(cs, cax=cbar_ax, orientation = 'horizontal',
                        extend='both',extendfrac='auto',ticks=np.arange(-1.,1.2,0.2))
    cbar.set_label('Correlation Coefficient') 
    figure_title = 'LENS 1920-2005, First Leaf Index PDO Correlations'
    fig.text(0.5, .97, figure_title,
     horizontalalignment='center',
     fontsize=14)
    plt.savefig('/Users/zlabe/documents/CESMspring/Fig12.eps',dpi=400,format='eps')
    
PlotCorr(cocoeff1,lat,lon)
    
def ModesTrend(mode):
    """
    Calculates trends in climate modes for historical/future LENS
    
    Parameters
    ----------
    modeline : mean trend line
    modetrend : mean trend
    modeslopeq : trend slope
    modelineq : trend line
    
    """
    
    modemean = sts.nanmean(mode)
    timemean = np.asarray(list(xrange(modemean.shape[0])))
    modeslopeq, modeinterceptq, r_valueq, p_valueq, std_errq = sts.linregress(timemean,modemean)
    modelineq = modeslopeq*timemean+modeinterceptq
    
    modetrend = np.empty((mode.shape[0]))
    for i in xrange(mode.shape[0]):
        timeq = np.asarray(list(xrange(mode.shape[1])))
        modeslope, modeintercept, r_value, p_value, std_err = sts.linregress(timeq,mode[i,:])
        modeline = modeslope*timeq+modeintercept
        modetrend[i] = modeslope
    return modeline, modetrend, modeslopeq, modelineq

#PDOline,PDOtrend, PDOtrendmean, PDOslopemean = ModesTrend(PDOave)
#NINOline,NINOtrend, NINOtrendmean, NINOslopemean = ModesTrend(NINOave)
#NAOline,NAOtrend, NAOtrendmean, NAOslopemean = ModesTrend(NAOave)

def PlotTrends(PDOline,NINOline,NAOline,PDOave,NINOave,NAOave):
    """
    Plots time series for climate modes in historical/future LENS
    
    
    Parameters
    ----------
    PDOline : all PDO values averaged for each year
    PDOave : average PDO for each ensemble
    NAOline : all NAO values averaged for each year
    NAOave : average NAO for each ensemble
    NINOline : all NINO values averaged for each year
    NINOave : average NINO for each ensemble
    """
    fig = plt.figure()
    fig.suptitle('2013-2080 LENS Climate Mode Trends',fontsize=20)
    plt.title('Ensemble Mean',fontsize=13)
    plt.plot(PDOline,label='PDO',color='darkred',linewidth=2)
    plt.plot(sts.nanmean(PDOave),color='k',linewidth=1,linestyle='solid',alpha=0.5)
    plt.plot(sts.nanmean(NINOave),color='k',linewidth=1,linestyle='solid',alpha=0.5)
    plt.plot(sts.nanmean(NAOave),color='k',linewidth=1,linestyle='solid',alpha=0.5)
    plt.plot(NINOline,label='Nino3.4',color='darkorange',linewidth=2)
    plt.plot(NAOline,label='NAO',color='darkblue',linewidth=2)
    plt.xlim([0,len(PDOline)-1])
    plt.ylim([-2,2])
    plt.yticks([-2,-1,0,1,2])
    plt.ylabel('Index')
    plt.legend(loc=4,prop={'size':10},shadow=True)
    x1 = np.arange(0,len(PDOline)+1,17)
    labels = ['2013','2030','2047','2064','2080']
    plt.xticks(x1,labels) 
    plt.grid(True)
    plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/LENS_meanmodetrends_future.png',dpi=300)
    
#PlotTrends(PDOline,NINOline,NAOline,PDOave,NINOave,NAOave)
