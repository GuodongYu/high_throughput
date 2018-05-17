from matplotlib import pyplot as plt
from high_throughput.utils.utils import formate_chemical_formula
from high_throughput.vasp.plot import plt_set 
import glob
import os
import json

plt_set(plt)
fig,ax = plt.subplots(1,1)

plot = 'zt'

files = glob.glob('./%s/*' % plot)
files = glob.glob('./%s/*' % plot)

for f in files:
    comp = os.path.split(f)[-1].replace('.json','')
    compound = formate_chemical_formula(comp)
    pf = json.load(open(f))
    comm = ['units','relaxation_time']
    T = [i for i in pf.keys() if i not in comm][0]
    ax.plot(pf[T]['carriers'],pf[T][plot],label=compound + ' ' + T +' K')
ax.set_xlim([1e+19,1e+21])
ax.set_ylim([0,1])
plt.legend()
plt.savefig(plot+'_vs_carrier.jpg')
    

