import glob
import json
import os
from high_throughput.vasp.plot import plt_set
from high_throughput.utils.utils import formate_chemical_formula

zt_dir = './zt'
pf_dir = './power_factor'

def plot_comparison(prop,xlim=None,ylim=None):
    import matplotlib.pyplot as plt
    plt_set(plt)
    fig,ax = plt.subplots(1,1)
    if prop == 'zt':
        dirs = glob.glob(os.path.join(zt_dir,'*'))
    if prop == 'power_factor':
        dirs = glob.glob(os.path.join(pf_dir,'*'))

    for f in dirs:
        compound = os.path.split(f)[-1].split('.')[0]
        data = json.load(open(f))
        unit = data['units'][prop]
        carrs = data['carriers']
        vals = data[prop+'_max']
        ax.plot(carrs,vals,marker='o',label=formate_chemical_formula(compound))
    if xlim is None:
        xlim = ax.get_xlim()
        avgx = (xlim[0]+xlim[1])/2.
    if ylim is None:
        ylim = ax.get_ylim()
        avgy = (ylim[0]+ylim[1])/2.
    ax.text(avgx/2.,avgy/2., 'Temps: 300 - 1300 K',fontsize=50)
    ax.legend(loc=0)
    ax.set_xlabel('Doping level ($cm^{-3}$)')
    ax.set_ylabel(prop+ '_max %s' % unit)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.savefig(prop+'_max_vs_carr_comparison.jpg')
    plt.clf()

def main(prop):
    plot_comparison(prop)

if __name__ == '__main__':
    main('power_factor')
    main('zt')
