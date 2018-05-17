from pymatgen.electronic_structure.bandstructure import Spin
from pymatgen.electronic_structure.boltztrap import BoltztrapAnalyzer
from pymatgen.io.vasp.outputs import Dos,Vasprun
from pymatgen.core.structure import Structure
from scipy import integrate
from scipy.integrate import simps
from scipy.optimize import curve_fit, bisect
from scipy.interpolate  import interp1d
from math import pi,exp,log
import os
import numpy as np
k_b = 8.6173303e-5

class SeebeckDOS(object):
    def __init__(self,bzdir='/home/gyu/workplace/intrinsic_defects/Zintl_Antimonides/Yb14AlSb11/Yb_2/bulk/sta/dos/boltztrap'):
        self.an = BoltztrapAnalyzer.from_files(bzdir)
        self.volume = self.an.vol*1e-24
        Tmin =min(self.an.get_carrier_concentration().keys())
        self.nelec = self.an.get_carrier_concentration()[Tmin][0]*self.volume
        self.cbm,self.vbm = self.an.dos.get_cbm_vbm()
        self.mid=(self.cbm+self.vbm)/2.0
        self.energy_min = min(self.an.dos.energies)
        self.energy_max = max(self.an.dos.energies)
        self.nstep = len(self.an.dos.energies)
        self.step = self.an.dos.energies[1]-self.an.dos.energies[0]
        vbm_diff_list=list(np.abs(self.an.dos.energies - self.vbm))
        cbm_diff_list=list(np.abs(self.an.dos.energies - self.cbm))
        self.vbm_ind = vbm_diff_list.index(min(vbm_diff_list))
        self.cbm_ind = cbm_diff_list.index(min(cbm_diff_list))
        self.mid_ind = (self.vbm_ind+self.cbm_ind)/2
        self.dos_base = self.normalization_base_dos()
        self.Ed = self.vbm-0.16
        self.gamma = 0.04
        self.scale = 0.0
        self.update_distorted_dos()

    def normalization_base_dos(self):
        def nelec_fake():
            xs = [E for E in self.an.dos.energies if E < self.mid]
            ys = [self.an.dos.get_densities()[i] for i in range(self.nstep) if self.an.dos.energies[i] < self.mid]
            return simps(ys,xs)
        nelec_fake = nelec_fake()
        scale = float(self.nelec)/nelec_fake
        dos_base = [i*scale for i in self.an.dos.get_densities()]
        return dos_base

    def update_resonantlevel_parameters(self,Ed,gamma,scale):
        self.Ed = Ed
        self.gamma = gamma
        self.scale = scale
        self.update_distorted_dos()

    def update_distorted_dos(self):
        """
        plus background dos and resonant level dos distortion
        Args:
            scale: the scale of the lorentz_func
        """
        tdos=[]
        for i in range(self.nstep):
            E = self.an.dos.energies[i]
            tot = (self.dos_base[i] +self.scale*self.lorentz_func(E))/self.volume
            tdos.append(tot)
        self.distorted_dos=tdos
        
    def total_nelec(self,Ef,T):
        summ = 0.
        for i in range(self.nstep):
            E = self.an.dos.energies[i]
            summ += self.distorted_dos[i]*self.fermi_dirac_func(E,Ef,T)
        return summ*step*self.volume

    def lorentz_func(self,E):
        """
        Args:
            Ed: resonant level position
            gamma: resonant level width
        """
        return 2.*(1/pi)*(0.5*self.gamma)/((E-self.Ed)**2+(0.5*self.gamma)**2)
    
    def integrate_lorentz(self):
        left = self.Ed-40*self.gamma
        right = self.Ed + 40*self.gamma
        return integrate.quad(self.lorentz_func,left,right)
 
    
    def fermi_dirac_func(self,E,Ef,T):
        return 1./(1.+exp((E-Ef)/(k_b*T)))

    def occupation(self, E, Ef, T, carrier_type='hole'):
        if carrier_type=='hole':
            return 1.-self.fermi_dirac_func(E,Ef,T)
        elif carrier_type == 'electron':
            return self.fermi_dirac_func(E,Ef,T)
   
    def get_doping_level(self,Ef,T):
        summ=0.
        for i in range(self.nstep):
            E = self.an.dos.energies[i]
            summ+=self.distorted_dos[i]*self.fermi_dirac_func( E, Ef, T)
        nelec_now = summ*self.step*self.volume
        return (self.nelec-nelec_now)/self.volume

    def get_Fermi_level(self,carr_con,T):
        """
        Args:
            carr_con: carrier concentration, positive and negative means hole and electrons 
            T: temperature
        """
        def zero_point(Ef):
            return self.get_doping_level(Ef,T) - carr_con
        Ef = bisect(zero_point,self.energy_min,self.energy_max)
        rcarr=self.get_doping_level(Ef,T)
        print 'at %s K, doping_level: %s' % (str(T),str(rcarr))
        return Ef

    def get_Fermi_level_0K(self,carr_con):
        neut = self.nelec/self.volume
        summ = 0.0
        diff = 1e+1000
        for i in range(self.nstep):
            summ += self.distorted_dos[i]
            carr = summ*self.step
            doping = neut-carr
            diff_new = abs(carr_con-doping)
            if diff_new < diff:
                diff = diff_new
                Ef = self.an.dos.energies[i]
                Ef_ind = i
            elif diff_new > diff:
                return Ef, Ef_ind 

    def calcu_Seebeck(self,carr_con,T):
        Ef0,Ef0_ind = self.get_Fermi_level_0K(carr_con)
        Ef = self.get_Fermi_level(carr_con,T)
        Ef_diff = list(np.abs(self.an.dos.energies-Ef))
        Ef_ind = Ef_diff.index(min(Ef_diff))
        if carr_con > 0:
            q=1
            Eedge = self.vbm
            Ef_ = (Eedge - Ef)/(k_b*T)
        else:
            q=-1
            Eedge = self.cbm
            Ef_=(Ef-Eedge)/(k_b*T)
        
        slope1 = (self.distorted_dos[Ef_ind+1] - self.distorted_dos[Ef_ind])/self.step
        slope2 = (self.distorted_dos[Ef_ind] - self.distorted_dos[Ef_ind-1])/self.step
        slope = abs(slope1+slope2)/2.
        dos_term = k_b*T*slope/self.distorted_dos[Ef_ind]
        dos_term = k_b*T*self.distorted_dos[Ef_ind]/abs(carr_con)
        numerator = (pi**2/3.0)*dos_term* (1-0.25*(2*Ef_+Ef_**2)*exp(-Ef_)) + log(2.)*(2+Ef_)*exp(-Ef_)
        denominator = 1 + dos_term*log(2.)*(2+Ef_)*exp(-Ef_)
        S = (k_b/q) * numerator/denominator
        S2 = (k_b/q) * (pi**2/3.0) *dos_term
        return S*1e+6, Ef

class SeebeckSPB(SeebeckDOS):
    def __init__(self):
        SeebeckDOS.__init__(self)

    def effective_mass(self,carr_type='hole'):
        if carr_type == 'hole':
            pass






def plt_set(plt,nrow=1,ncolumn=1):
    fig_width = 20*ncolumn
    fig_height = 16*nrow
    fig_size = [fig_width, fig_height]
    params = {'backend': 'ps',
              'axes.labelsize': 40,
              'font.size': 40,
              'font.weight': 'normal',
              'legend.fontsize': 40,
              'xtick.labelsize': 40,
              'ytick.labelsize': 40,
              'text.usetex': False,
              'figure.figsize': fig_size,
              'lines.markersize': 40,
              'lines.linewidth':10}
    plt.rcParams.update(params)

class ResonantLevelPlotter(object):
    def __init__(self,ResonantLevel):
        self.rl = ResonantLevel

    def plot_backgroupd_dos(self,save_dir='./',fmt='png'):
        from matplotlib import pyplot as plt
        xs = [self.rl.an.dos.energies[i] for i in range(self.rl.nstep)]
        ys = [self.rl.dos_base[i] for i in range(self.rl.nstep)]
        plt.plot(xs,ys)
        plt.savefig(os.path.join(save_dir,'dos.%s' % fmt),format=fmt)
        plt.clf()

    def plot_resonant_level(self, save_dir='./',fmt='png',xlim=None):
        from matplotlib import pyplot as plt
        nstep = 1000
        if xlim is None:
            xlim=(self.rl.vbm-0.5,self.rl.cbm+0.5)
        step = (xlim[1]-xlim[0])/nstep
        xs = [xlim[0]+i*step for i in range(nstep+1)]
        ys = [self.rl.lorentz_func(x) for x in xs]
        plt.plot(xs,ys)
        plt.savefig(os.path.join(save_dir,'Lorentz_function.%s' % fmt),format=fmt)
        plt.clf()

    def plot_distored_dos(self, save_dir='./',fmt='png',xlim=None):
        from matplotlib import pyplot as plt
        nstep = 1000
        if xlim is None:
            xlim=(self.rl.vbm-0.5,self.rl.cbm+0.5)
        step = (xlim[1]-xlim[0])/nstep
        xs = [xlim[0]+i*step for i in range(nstep+1)]
        ys = [self.rl.distorted_dos(x) for x in xs]
        plt.plot(xs,ys)
        plt.savefig(os.path.join(save_dir,'distorted_dos.%s' % fmt),format=fmt)
        plt.clf()

    def plot_dos_fitting(self,save_dir='./',fmt='png',xlim=None):
        from matplotlib import pyplot as plt
        nstep = 1000
        if xlim is None:
            xlim=(self.rl.vbm-0.5,self.rl.cbm+0.5)
        step = (xlim[1]-xlim[0])/nstep
        xs = [xlim[0]+i*step for i in range(nstep+1)]
        ys = [self.rl.dos_background(x) for x in xs]
        yfits = [self.rl.dos_fit(x) for x in xs]
        plt.plot(xs,ys,color='blue')
        plt.plot(xs,yfits,color='red')
        plt.savefig(os.path.join(save_dir,'distorted_dos.%s' % fmt),format=fmt)
        plt.clf()
    
    def plot_seebeck_vs_T(self,carr_cons,save_dir='./',fmt='jpg'):
        from matplotlib import pyplot as plt
        plt_set(plt)
        Ts=range(300,950,50)
        legends=[]
        colors=['red','blue','green','gray','yellow']
        c=0
        for carr_con in carr_cons:
            Ss=[]
            Ef_comms = []
            for T in Ts:
                S,Ef = self.rl.calcu_Seebeck(carr_con,T)
                Ss.append(S)
                if carr_con > 0:
                    Ef_comm = 'v+%.2f' % (Ef-self.rl.vbm) if self.rl.vbm-Ef <0 else 'v%.2f' % (Ef-self.rl.vbm)
                else:
                    Ef_comm = 'c+%.2f' % (Ef-self.rl.cbm) if self.rl.cbm-Ef <0 else 'c%.2f' % (Ef-self.rl.cbm)
                Ef_comms.append(Ef_comm)
            legend,=plt.plot(Ts,Ss, color=colors[c], marker='o', label = 'carr_con:%.4g' % carr_con)
            c+=1
            legends.append(legend)
            for i in range(len(Ts)):
                plt.text(Ts[i],Ss[i],Ef_comms[i],fontsize=15)
        plt.legend(handles=legends, loc=1)
        plt.savefig(os.path.join(save_dir,'Seebeck_vs_T.%s' % fmt),format=fmt)
        plt.clf()

