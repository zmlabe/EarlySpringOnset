"""
*Script plots climate modes for CESM Control in a time series (Selected case year).*
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

### Path may be broken...
directory = '/volumes/data/zml5/cesm/'
files = 'CESM1_CAM5.2_LENS_Control.cvdp_data.401-2200.nc'
filename = directory + files

values = Dataset(filename)
time = values.variables['time'][:] # months since 1/15/401
nino34 = values.variables['nino34'][:]
pdo = values.variables['pdo_timeseries_mon'][:]
nao = values.variables['nao_pc_mon'][:]
values.close()

nino34_shape = np.reshape(nino34,(1800,12))
nao_shape = np.reshape(nao,(1800,12))
pdo_shape = np.reshape(pdo,(1800,12))

nino34_shape = nino34_shape[:599,:]
nao_shape = nao_shape[:599,:]
pdo_shape = pdo_shape[:599,:]
    
######## LF
### Plots for Year 653 [251]
nao653 = nao_shape[251,:]
pdo653 = pdo_shape[251,:]
nino653 = nino34_shape[251,:]

### Plots for Year 418 [16]
nao418 = nao_shape[16,:]
pdo418 = pdo_shape[16,:]
nino418 = nino34_shape[16,:]

### Plots for Year 846 [444]
nao846 = nao_shape[444,:]
pdo846 = pdo_shape[444,:]
nino846 = nino34_shape[444,:]

fig = plt.figure()
ax = plt.plot(nao653,color='darkred',linewidth=3)
ax1 = plt.plot([0]*len(nao653),color='k',linewidth=2)
x = list(xrange(0,len(nao653)))
labels = ['Jan','Feb','March','April','May','June','July','Aug',
          'Sept','Oct','Nov','Dec']
plt.ylabel('Standard Deviation')
plt.xlabel('Months')
plt.xticks(x,labels)
plt.xlim([0,11])
plt.ylim([-3,3])
plt.grid(True)
plt.title('CESM Control Year 653, NAO Index',fontsize=16)
plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/nao_lf846.png',dpi=300)

######## BL
### Plots for Year 406 [4]
nao406 = nao_shape[4,:]
pdo406 = pdo_shape[4,:]
nino406 = nino34_shape[4,:]

### Plots for Year 606 [204]
nao606 = nao_shape[204,:]
pdo606 = pdo_shape[204,:]
nino606 = nino34_shape[204,:]

### Plots for Year 958 [756]
nao958 = nao_shape[556,:]
pdo958 = pdo_shape[556,:]
nino958 = nino34_shape[556,:]

#fig = plt.figure()
#ax = plt.plot(nao653,color='darkred',linewidth=3)
#ax1 = plt.plot([0]*len(nao653),color='k',linewidth=2)
#x = list(xrange(0,len(nao653)))
#labels = ['Jan','Feb','March','April','May','June','July','Aug',
#          'Sept','Oct','Nov','Dec']
#plt.ylabel('Standard Deviation')
#plt.xlabel('Months')
#plt.xticks(x,labels)
#plt.xlim([0,11])
#plt.ylim([-2,2])
#plt.grid(True)
#plt.title('CESM Control Year 653, NAO Index',fontsize=16)
#plt.savefig('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_SpringOnset/Results/nao_bl406.png',dpi=300)