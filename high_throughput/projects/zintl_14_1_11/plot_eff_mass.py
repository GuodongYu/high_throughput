import glob,os
from operator import itemgetter
import json
from high_throughput.vasp.plot import plt_set
from high_throughput.utils.utils import formate_chemical_formula
import numpy as np

root_dir = '/home/gyu/workplace/My_Work/14-1-11/'

def ylabel(mass_type):
    if mass_type == 'cond_eff_mass':
        return 'Conductive effective mass ($m_e$)'
    elif mass_type in ['seeb_eff_mass','seeb_eff_mass_Hall_carr']:
        return 'Seebeck effective mass ($m_e$)'
    elif mass_type == 'comp_fac':
        return 'Fermi surface complexity factor'
    else:
        raise TypeError()

def plot_eff_mass_comparison(mass_type,T,cut=True,xlim=[-0.35,0.0],ylim=[0,2],fmt='png'):
    import matplotlib.pyplot as plt
    plt_set(plt)
    files = glob.glob(root_dir+'/'+mass_type+'/*')
    files.sort()
    legends=[]
    colors = ['black','blue','red','green','purple','orange']
    k=0
    for i in files:
        name = os.path.split(i)[-1].split('.')[0]
        name = formate_chemical_formula(name)
        data = json.load(open(i))
        try:
            mus = data['mu_steps']
        except:
            mus = data['mu_step']
        mass = data[T]
        mass_mus = zip(mus,mass)
        mass_mus = [i for i in mass_mus if data['vbm']-0.4 <= i[0] <= data['vbm']]
        mass = [i[1] for i in mass_mus]
        mus = [i[0] for i in mass_mus]
        mus_labels = data['carr_Ef_pairs'].values()
        if cut:
            max_ = max(mass)
            diff = np.abs(np.array(mass)-max_)
            ind = min(enumerate(diff), key=itemgetter(1))[0]
            mass = mass[ind:]
            mus = mus[ind:]    
        legend,=plt.plot(mus,mass,label=name,color=colors[k])
        plt.axvline(mus_labels[0],color=colors[k],linestyle='dashed', linewidth=3)
        plt.axvline(mus_labels[1],color=colors[k],linestyle='solid', linewidth=3)
        legends.append(legend)
        k+=1
    plt.legend(handles=legends,loc=1)
    plt.xlabel('Fermi level (eV)')
    plt.ylabel(ylabel(mass_type))
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.savefig(mass_type+'_vs_Ef.%s' % fmt,format=fmt)

plot_eff_mass_comparison('seeb_eff_mass_Hall_carr','300')

