from high_throughput.vasp.plot import plt_set
from matplotlib import pyplot as plt
plt_set(plt,width=25,height=20,scale=15)
fig,ax=plt.subplots(1,1)
def line(x,q,c):
    return q*x+c

x=[i*0.01 for i in range(-7,87)]
ax.plot(x,[line(i,1,0.02) for i in x],label='q=+1')
ax.scatter(-0.02/1.,0)
ax.text(-0.02/1.,0,'$E_{pin}^p$',fontsize=60,verticalalignment='bottom',horizontalalignment='right')
ax.plot(x,[line(i,-1,0.8) for i in x],label='q=-1')
ax.scatter(0.8/1.,0)
ax.text(0.8/1.,0,'$E_{pin}^n$',fontsize=60,verticalalignment='bottom')

s = 6
ax.axvline(0.1,color='black',linewidth=s,linestyle='dashed')
ax.text(0.1,0.4,'$VBM$',fontsize=60)
ax.axvline(0.7,color='black',linewidth=s,linestyle='dashed')
ax.text(0.7,0.4,'$CBM$',fontsize=60)
ax.axhline(0.0,color='black',linewidth=s)
ax.set_xlim(-0.1,0.9)
ax.set_ylim(-0.1,0.9)
ax.set_xlabel('Fermi level (eV)')
ax.set_ylabel('Defect formation energy (eV)')
ax.legend()
fig.set_tight_layout({'pad':0.2, 'h_pad':0.1, 'w_pad':0.0})
fig.savefig('Epin_model.png')

