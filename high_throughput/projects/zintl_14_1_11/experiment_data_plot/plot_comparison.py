import os,glob
from high_throughput.vasp.plot import plt_set

root_dir = '/home/gyu/workplace/My_Work/14-1-11/experimental_data'
names={'Yb14MgSb11':'Yb14MgSb11','Yb14MnSb11':'Yb14MnSb11','CoSb3':'CoSb2.85Te0.15','SiGe':'n_Si80Ge20_nano','PbTe':'Na_doped_PbTe'}

def get_files():
    outs = {}
    outs['seebeck']=[]
    outs['resistivity']=[]
    outs['conductivity']=[]
    outs['thermal_condtivity']=[]
    for prop in ['seebeck','conductivity','thermal_conductivity','resistivity']:
        files= glob.glob(root_dir+'/'+prop+'/*')
        for i in files:
            name = os.path.split(i)[-1][:-4]
            if name in names.values():
                outs[prop].append(i)
    return outs

def unit(prop):
    if prop == 'seebeck':
        return '($\mu VK^{-1}$)'
    elif prop == 'resistivity':
        return '($m\Omega cm$)'
    elif prop == 'thermal_conductivity':
        return '($W/mK$)'


def plot_comparison(prop,fmt='png'):
    import matplotlib.pyplot as plt
    plt_set(plt)
    files_all = get_files()
    files = files_all[prop]
    legends=[]
    Tmin = 1000; Tmax =0
    for i in files:
        with open(i) as f:
            name = os.path.split(i)[-1][:-4]
            contents = f.read().split('\n')
            temps_exp = [float(i.split(',')[0]) for i in contents if i]
            Tmin = Tmin if min(temps_exp) > Tmin else min(temps_exp)
            Tmax = Tmax if max(temps_exp) < Tmax else max(temps_exp)
            quantity_exp = [abs(float(i.split(',')[1])) for i in contents if i]
            legend,=plt.plot(temps_exp, quantity_exp,marker='o',
                linewidth=0,label=name)
            legends.append(legend)
    if prop == 'resistivity':
        for i_add in files_all['conductivity']:
            with open(i_add) as f:
                name = os.path.split(i)[-1][:-4]
                contents = f.read().split('\n')
                temps_exp = [float(i.split(',')[0]) for i in contents if i]
                Tmin = Tmin if min(temps_exp) > Tmin else min(temps_exp)
                Tmax = Tmax if max(temps_exp) < Tmax else max(temps_exp)
                quantity_exp = [1./float(i.split(',')[1]) for i in contents if i]
                legend,=plt.plot(temps_exp, quantity_exp,marker='o',
                    linewidth=0,label=name)
                legends.append(legend)
    plt.legend(handles=legends)
    plt.xlabel('Temperature (K)')
    plt.ylabel(prop + unit(prop))
    #plt.title(compound)
    plt.xlim((Tmin-50,Tmax+50))
    plt.savefig('%s.%s' % (prop,fmt),format=fmt)
    
plot_comparison('seebeck')
