import os
import numpy as np
from math import log,exp
from high_throughput.vasp.plot import plt_set

root_dir='/home/gyu/workplace/My_Work/14-1-11/experimental_data/resistivity'
def scattering_type_determine(compound,fmt='png'):
    from matplotlib import pyplot as plt
    plt_set(plt)
    resis_file = os.path.join(root_dir, compound+'.csv')
    with open(resis_file) as f:
        contents = f.read().split('\n')
        Ts = [float(i.split(',')[0]) for i in contents if i]
        T0=round(Ts[0]);T1=round(Ts[-1])
        temps_exp = [log(float(i.split(',')[0])) for i in contents if i]
        if compound == 'Yb14ZnSb11':
            temps_exp = [i for i in temps_exp if i <=log(800)]
        T0 = exp(min(temps_exp));T1=exp(max(temps_exp))
        print temps_exp
        resis_exp = [float(i.split(',')[1]) for i in contents if i]
        resis_exp = [resis_exp[i] for i in range(len(temps_exp))]
        cond = [log(1.0/i) for i in resis_exp]
        x0 = temps_exp[-1];y0=cond[-1]
        x1 = temps_exp[-2];y1=cond[-2]
        k0 = (y1-y0)/(x1-x0)
    legends=[]
    def f(x,k):
        return y0+k*(x-x0)
    for k in [-0.75,-0.7,-0.6,-0.5]:
        legend, = plt.plot(temps_exp, [f(i,k) for i in temps_exp], label='k=%.2f' % k)
        legends.append(legend)
    plt.xlabel('log(T) T:%.f-%.fK' % (T0,T1))
    plt.ylabel('log($\sigma$)')
    plt.title(compound)
    plt.legend(handles=legends, loc=1)
    plt.scatter(temps_exp,cond)
    plt.tight_layout()
    plt.savefig('%s.%s' % (compound,fmt),format=fmt)

scattering_type_determine('Yb14MnSb11')


    
