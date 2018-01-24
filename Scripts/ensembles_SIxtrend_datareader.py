"""
*Functions read netcdf files for LENS SIx trends
"""

import numpy as np
from netCDF4 import Dataset
#import Ngl 
#import Nio 
#import os

def histtrends():
    """
    Reads in future LENS 1920-2005 SI-x data

    
    Returns
    ----------
    lstfrz_trendh : array of last freeze trends
    damage_trendh : array of damage trends
    leaf_trendh : array of leaf trends 
    lat : array of latitudes
    lon : array of longitudes 
    """
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/cesmclimo/'
    name = 'damagetrends_2005.nc'
    filename = directory + name
    data = Dataset(filename,'r')
    damage_trendh = data.variables['trend'][:]
    lat = data.variables['latitude'][:]
    lon = data.variables['longitude'][:]
    data.close()
    
    name2 = 'leaftrends_2005.nc'
    filename2 = directory + name2
    data2 = Dataset(filename2,'r')
    leaf_trendh = data.variables['trend'][:]
    data2.close()
    
    name3 = 'lstfrztrends_2005.nc'
    filename3 = directory + name3
    data3 = Dataset(filename3,'r')
    lstfrz_trendh = data.variables['trend'][:]
    data3.close()
    
    return lstfrz_trendh,damage_trendh,leaf_trendh, lat, lon
    
def histtrends2():
    """
    Reads in future LENS 1970-2005 SI-x data

    
    Returns
    ----------
    lstfrz_trendh2 : array of last freeze trends
    leaf_trendh2 : array of leaf trends 
    """
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/cesmclimo/'
    name2 = 'Leaftrends_1970.nc'
    filename2 = directory + name2
    data2 = Dataset(filename2,'r')
    leaf_trendh2 = data2.variables['trend'][:]
    data2.close()
    
    name3 = 'lSTFRZtrends_1970.nc'
    filename3 = directory + name3
    data3 = Dataset(filename3,'r')
    lstfrz_trendh2 = data3.variables['trend'][:]
    data3.close()
    
    return lstfrz_trendh2,leaf_trendh2

def futtrends():
    """
    Reads in future LENS 2006-2080 SI-x data

    
    Returns
    ----------
    lstfrz_trendf : array of last freeze trends
    damage_trendf : array of damage trends
    leaf_trendf : array of leaf trends 
    lat : array of latitudes
    lon : array of longitudes 
    """
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/cesmclimo/'
    name = 'damagetrends_0680.nc'
    filename = directory + name
    data = Dataset(filename,'r')
    damage_trendf = data.variables['trend'][:]
    lat = data.variables['latitude'][:]
    lon = data.variables['longitude'][:]
    data.close()
    
    name2 = 'leaftrends_0680.nc'
    filename2 = directory + name2
    data2 = Dataset(filename2,'r')
    leaf_trendf = data.variables['trend'][:]
    data2.close()
    
    name3 = 'lstfrztrends_0680.nc'
    filename3 = directory + name3
    data3 = Dataset(filename3,'r')
    lstfrz_trendf = data.variables['trend'][:]
    data3.close()
    
    return lstfrz_trendf,damage_trendf,leaf_trendf, lat, lon
    
def futtrends2():
    """
    Reads in future LENS 2006-2040 SI-x data

    
    Returns
    ----------
    lstfrz_trendf2 : array of last freeze trends
    leaf_trendf2 : array of leaf trends 
    """
    directory = '/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/cesmclimo/'
    
    name2 = 'Leaf_2040.nc'
    filename2 = directory + name2
    data2 = Dataset(filename2,'r')
    leaf_trendf2 = data2.variables['trend'][:]
    data2.close()
    
    name3 = 'lSTFRZ_2040.nc'
    filename3 = directory + name3
    data3 = Dataset(filename3,'r')
    lstfrz_trendf2 = data3.variables['trend'][:]
    data3.close()
    
    return lstfrz_trendf2,leaf_trendf2
    
#lstfrz_trendh,damage_trendh,leaf_trendh, lat, lon = histtrends()
#lstfrz_trendf,damage_trendf,leaf_trendf, lat, lon = futtrends()


################
### PyNGL Example

#filename = '/users/zml5/desktop/trends'
#wks_type = 'png'
#os.system('rm ' + filename + '.' + wks_type)
#wks = Ngl.open_wks(wks_type,filename)
#
#where = np.isnan(leaf_trendf)
#leaf_trendf[where] = 999.
#
#resources = Ngl.Resources()
#resources.sfXCStartV = float(min(lon))
#resources.sfXCEndV   = float(max(lon))
#resources.sfYCStartV = float(min(lat))
#resources.sfYCEndV   = float(max(lat))
#resources.mpProjection = 'Mercator'
#resources.mpFillOn = True
#resources.mpLimitMode = "LatLon"    # Limit the map view.
#resources.mpMinLonF   = float(min(lon))
#resources.mpMaxLonF   = float(max(lon))
#resources.mpMinLatF   = float(min(lat))
#resources.mpMaxLatF   = float(max(lat))
#resources.mpPerimOn   = True        # Turn on map perimeter.
#resources.cnFillOn              = True     # Turn on contour fill.
#resources.cnLineLabelsOn        = False    # Turn off line labels.
#resources.cnInfoLabelOn         = False    # Turn off info label.
#resources.cnFillPalette        = 'BlWhRe'
#resources.mpFillColors = ['white','gray','transparent','transparent']
#resources.mpMaskAreaSpecifiers  = "USStatesLand"
##resources.pmLabelBarDisplayMode = "Never"
#resources.tiMainString = 'Leaf Trends 2006-2080'
#resources.cnLinesOn      = True
#resources.lbTitleString  = "~F4~days/decade"
#resources.sfMissingValueV = 999.
#
#resources.cnLevelSelectionMode = "ExplicitLevels" # Define own levels.
#resources.cnLevels             = np.arange(-5,5.5,.5)
#
#resources.mpOutlineBoundarySets = "usstates"      # "geophysical"
#resources.mpFillBoundarySets    = "geophysical"
##resources.cnFillDrawOrder       = "Predraw"
#
#maps = Ngl.contour_map(wks,leaf_trendf[0,:,:],resources)
#
#
#Ngl.end()