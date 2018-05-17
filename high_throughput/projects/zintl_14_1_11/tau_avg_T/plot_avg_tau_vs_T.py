import json
from high_throughput.vasp.plot import plt_set
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


plt_set(plt,scale=12)
fig,(ax0,ax1) = plt.subplots(2,1)
Mn = json.load(open('Yb14MnSb11.json'))
Zn = json.load(open('Yb14ZnSb11.json'))

xmajorLocator   = MultipleLocator(200)
ax0.xaxis.set_major_locator(xmajorLocator)

ax0.plot(Mn['ac']['temperatures'],Mn['ac']['2.9e+20'],label='$ac$')
ax0.plot(Mn['ionic']['temperatures'],Mn['ionic']['2.9e+20'],label='$ii$')
ax0.legend()
ax0.set_xlim([300,1300])
ax0.set_ylim([0,160])
ax0.set_xlabel('Temperature (K)')
ax0.set_ylabel('$ \overline{\\tau} (fs)$')
ax0.set_title('(a) $Yb_{14}MnSb_{11}$',loc='left')
#ax0.annotate('test',xy=(305,Mn['ac']['2.9e+20'][5]+0.5),xytext=(400,140),arrowprops=dict(width=5))



ax1.plot(Zn['ac']['temperatures'],Zn['ac']['3.5e+20'],label='$ac$')
ax1.plot(Zn['ionic']['temperatures'],Zn['ionic']['3.5e+20'],label='$ii$')
ax1.legend()
ax1.set_xlim([300,800])
ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('$ \overline{\\tau} (fs)$')
ax1.set_title('(b) $Yb_{14}ZnSb_{11}$',loc='left')
fig.set_tight_layout({'pad':0.2, 'h_pad':0.1, 'w_pad':0.0})
plt.savefig('tau_avg.svg')


