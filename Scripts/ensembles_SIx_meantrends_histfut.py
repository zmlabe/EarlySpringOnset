"""
*Script creates plot for all trends in CESM control and LENS future/historical SIx*
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import nanmean

hist_trendlf,hist_trendbl = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/1920-2005trends.txt',delimiter=',')
fut_trendlf,fut_trendbl = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/2006_2008trends.txt',delimiter=',')
cont_trendlf,cont_trendbl = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/controltrends.txt',delimiter=',')

hist_lf = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/1920-2005lf.txt',delimiter=',')
hist_bl = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/1920-2005bl.txt',delimiter=',')
hist_lstfrz = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/1920-2005lstfrz.txt',delimiter=',')
fut_lf = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/2006_2008lf.txt',delimiter=',')
fut_bl = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/2006_2008bl.txt',delimiter=',')
fut_lstfrz = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/2006_2008lstfrz.txt',delimiter=',')
cont_lf = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/controllf.txt',delimiter=',')
cont_lst = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/controllstfrz.txt',delimiter=',')
cont_bl = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/controlbl.txt',delimiter=',')
best_lf = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/bestlf.txt',delimiter=',')
best_bl = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/bestbl.txt',delimiter=',')
best_lstfrz = np.genfromtxt('/volumes/eas-shared/ault/ecrl/spring-indices/LENS_springonset/data/bestlstfrz.txt',delimiter=',')

print 'DONE!'

cont_lf = cont_lf[-200:]
cont_lst = cont_lst[-200:]

historical_length = len(hist_lf)
future_length = len(fut_lf)
control_length = len(cont_lf)

mean_hislf = nanmean(hist_lf)
mean_hisbl = nanmean(hist_bl)
mean_futlf = nanmean(fut_lf)
mean_futbl = nanmean(fut_bl)

control_dates = list(xrange(402,999))
historical_dates = list(xrange(1920,2006))
future_dates = list(xrange(2006,2081))

c = list(xrange(270,356))
mm = list(xrange(270,364))
d = list(xrange(356,431))
#e = list(xrange(598,691))

zeros = [75]*len(control_dates)
a = np.append(zeros,mean_hislf)
zeros2 = [75]*len(a)
b= np.append(zeros2,mean_futlf)

leaf_line1 = np.append(cont_lf,mean_hislf)
leaf_line = np.append(leaf_line1,mean_futlf)

fig = plt.figure()
ax=fig.add_subplot(111)
fig.suptitle('CESM Control and LENS Members 2-30',fontsize=20)
plt.title('First Leaf and Last Freeze Indices',fontsize=14)

may = [140]*(len(leaf_line)+139)
maxq = [121]*(len(leaf_line)+139)
april = [91]*(len(leaf_line)+139)
base = [70]*(len(leaf_line)+139)
t = np.arange(361+139)

for i in xrange(len(hist_lf)):
    plt.plot(c,hist_lf[i,:],color='gray',linewidth=0.3,alpha=0.3)
    plt.plot(d,fut_lf[i,:],color='gray',linewidth=0.3,alpha=0.3)

plt.plot(april,color='k',linewidth=2,linestyle='--',alpha=0)
plt.plot(cont_lf, color='teal',linewidth=1.1,label='Control Leaf')

plt.plot(cont_lst, color='darkblue',linewidth=1.1,label='Control LSRFRZ')
plt.plot(c,mean_hislf,color='darkorange',linewidth=1.1,label='Historical Leaf')

plt.plot(d,mean_futlf,color='red',linewidth=1.1,label='Future Leaf')
#plt.plot(e,best_lf[41:],color='w',linewidth=1.2,label='Observations LF')
plt.plot(d,fut_lstfrz,color ='maroon',linewidth=1.1,label ='Future LSTFRZ')
#5plt.plot(e,best_lstfrz[41:],color ='gray',linewidth=1.2,label='Observations LSTFRZ')
plt.plot(c,hist_lstfrz,color='saddlebrown',linewidth=1.1,label='Historical LSTFRZ')

labels = ['','','','','','','','','','','','','','','','','','','','','1850',
          '','','','','','','1920','','','','','','','','   2006','','','','','','2062','','','','2100']
x = list(xrange(0,len(leaf_line)+100,10))

plt.fill_between(t,may,maxq,alpha=0.5,color='yellow')
plt.fill_between(t,april,base,alpha=0.5,color='palegoldenrod')
plt.fill_between(t,maxq,april,alpha=0.3,color='salmon')

threshold3 = [81.5]*(len(leaf_line)+100)

plt.axvline(200,linestyle='--',color='k',alpha=0.4)
plt.axvline(270,linestyle='--',color='k',alpha=0.4)
plt.axvline(356,linestyle='--',color='k',alpha=0.4)
line = plt.plot(threshold3,'k',linewidth=3,label='-3$\sigma$ Threshold',
                 linestyle='dashed')
plt.legend(handles=line,loc=1,prop={'size':9},shadow=True)

ax.text(0.22,0.89,'1850 Control',size=11,horizontalalignment='center',
            backgroundcolor='white',verticalalignment= 'center',
            bbox=dict(facecolor='',alpha=0.0),
           transform=ax.transAxes)
ax.text(0.693,0.8,'LENS Historical',size=11,horizontalalignment='center',
            backgroundcolor='orangered',verticalalignment= 'center',
            bbox=dict(facecolor='',alpha=0.0),
           transform=ax.transAxes)
ax.text(0.90,0.631,'LENS Future',size=11,horizontalalignment='center',
            backgroundcolor='white',verticalalignment= 'center',
            bbox=dict(facecolor='',alpha=0.0),
           transform=ax.transAxes)
a=plt.xticks(x,labels)
plt.tick_params(axis='x',which='major',top='off',bottom='off',labelbottom='on',labeltop='off')

plt.ylabel('DOY',fontsize=14)
#plt.grid(True)
plt.vlines(356+55,ymin=0,ymax=81.5,linewidth=3,linestyle='--',color='k',zorder=10)
plt.xlabel('Years',fontsize=14)
plt.xlim([0,452])
plt.ylim([70,130])
ax2=plt.twinx()
ax2.set_ylabel('      March                        April                    May',fontsize=14)
plt.tick_params(axis='y',which='both',right='off',left='off',labelright='off')
plt.savefig('/Users/zlabe/documents/CESMspring/Fig16.eps',dpi=400,format='eps')