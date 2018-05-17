from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import glob,os,json
from high_throughput.vasp.plot import plt_set
from high_throughput.utils.utils import formate_chemical_formula
seebeck_root_dir='./seebeck'
resistivity_root_dir='./resistivity'
seebeck_dir = '/home/gyu/workplace/My_Work/14-1-11/experimental_data/seebeck'
resistivity_dir = '/home/gyu/workplace/My_Work/14-1-11/experimental_data/resistivity'

def sort_(scatters):
    scatters.sort()
    out = []
    for i in scatters:
        if 'crta' in i:
            ww = i
        else:
            out.append(i)
    out.append(ww)
    return out

def plot_seebeck_comparison(ax,compound,fmt='png'):
    if compound == 'Yb14ZnSb11':
        xmajorLocator   = MultipleLocator(100)
        ymajorLocator   = MultipleLocator(10)
    else:
        xmajorLocator   = MultipleLocator(200)
        ymajorLocator   = MultipleLocator(20)
    scatters = glob.glob(seebeck_root_dir+'/'+compound+'/*')
    scatters = sort_(scatters) 
    for scat_file in scatters:
        scat_type = os.path.split(scat_file)[-1][:-5]
        if scat_type not in ['ac','crta']:
            continue
        if scat_type == 'ac':
            scat_type = 'DP'
        if 'E' in scat_type:
            continue
            scat_type = ''
        if scat_type == 'crta':
            scat_type = 'CRTA'
        data = json.load(open(scat_file))
        temps = data['Temps']
        try:
            seebecks = data['seebeck'][0]
        except:
            seebecks = data['seebeck']['average'][0]
        if scat_type in ['ac','ii']:
            ax.plot(temps,seebecks,label=scat_type,linestyle='dashed')
        else:
            ax.plot(temps,seebecks,label=scat_type)

    exp_file = seebeck_dir+'/'+compound+'.csv'
    with open(exp_file) as f:
        contents = f.read().split('\n')
        temps_exp = [float(i.split(',')[0]) for i in contents if i]
        quantity_exp = [float(i.split(',')[1]) for i in contents if i]
    ax.plot(temps_exp, quantity_exp,
                color='black',marker='o',linewidth=0)
    ax.legend(loc=2)
    ax.set_xticks(xrange(10))
    #ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Seebeck ($\mu V/K$)')
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    #ax.set_title(compound)
    if compound == 'Yb14ZnSb11':
        ax.set_xlim((300,800))
        ax.set_ylim((0,50))
    elif compound == 'Yb14MnSb11':
        ax.set_xlim(400,1300)
        ax.set_ylim((60,200))
    elif compound == 'Yb14MgSb11':
        ax.set_xlim(300,1300)
        ax.set_ylim((40,250))
    return ax
        
def plot_resistivity_comparison(ax,compound,fmt='png'):
    if compound == 'Yb14ZnSb11':
        xmajorLocator   = MultipleLocator(100)
    else:
        xmajorLocator   = MultipleLocator(200)
    ymajorLocator   = MultipleLocator(2)
    scatters = glob.glob(resistivity_root_dir+'/'+compound+'/*')
    scatters = sort_(scatters)
    for scat_file in scatters:
        scat_type = os.path.split(scat_file)[-1][:-5]
        if scat_type not in ['ac','crta']:
            continue
        if scat_type == 'ac':
            scat_type = 'DP'
        if 'E' in scat_type:
            continue
            scat_type = ''
        data = json.load(open(scat_file))
        temps = data['Temps']
        try:
            seebecks = data['resistivity'][0]
        except:
            seebecks = data['resistivity']['average'][0]

        tau0 = data['relaxation_time']
        if scat_type == 'crta':
            scat_type = 'CRTA'
            tau0 = tau0*1e+15
        if scat_type in ['ac','ii']:
            ax.plot(temps,seebecks,label=scat_type,linestyle='dashed')
        else:
            ax.plot(temps,seebecks,label=scat_type)

    exp_file = resistivity_dir+'/'+compound+'.csv'
    with open(exp_file) as f:
        contents = f.read().split('\n')
        temps_exp = [float(i.split(',')[0]) for i in contents if i]
        quantity_exp = [float(i.split(',')[1]) for i in contents if i]
    ax.plot(temps_exp, quantity_exp,
                color='black',marker='o',linewidth=0)
    ax.legend(loc=2)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Resistivity ($m\Omega cm$)')
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    #ax.set_title(compound)
    if compound == 'Yb14ZnSb11':
        ax.set_xlim((300,800))
        ax.set_ylim((0.,10))
    elif compound == 'Yb14MnSb11':
        ax.set_xlim((min(temps_exp)-50,max(temps_exp)+50))
        ax.set_ylim(0,20)
    elif compound == 'Yb14MgSb11':
        ax.set_xlim((min(temps_exp)-50,max(temps_exp)+50))
        ax.set_ylim(0,20)
    return ax
from matplotlib import pyplot as plt
plt_set(plt,2,1,width=120,height=30,scale=35)
title_size = 150
fig,([ax0,ax2,ax4],[ax1,ax3,ax5]) = plt.subplots(2,3)
ax0 = plot_seebeck_comparison(ax0,'Yb14MnSb11')
ax0.set_title('(a)$Yb_{14}MnSb_{11}$',loc='left',fontdict={'fontsize':title_size})
ax1 = plot_resistivity_comparison(ax1,'Yb14MnSb11')
ax2 = plot_seebeck_comparison(ax2,'Yb14ZnSb11')
ax2.set_title('(b)$Yb_{14}ZnSb_{11}$',loc='left',fontdict={'fontsize':title_size})
ax3 = plot_resistivity_comparison(ax3,'Yb14ZnSb11')
ax4 = plot_seebeck_comparison(ax4,'Yb14MgSb11')
ax4.set_title('(c)$Yb_{14}MgSb_{11}$',loc='left',fontdict={'fontsize':title_size})
ax5 = plot_resistivity_comparison(ax5,'Yb14MgSb11')

fig.set_tight_layout({'pad':0.2, 'h_pad':0.1, 'w_pad':0.0})
#plt.show()
plt.savefig('seebeck_resistivity_fit.png')
