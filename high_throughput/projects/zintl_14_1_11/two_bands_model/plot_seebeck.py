from scipy.integrate import quad
from scipy.optimize import fsolve, curve_fit
from fdint import fdk
import numpy as np
from math import exp,pi
from comm.unit.unit_convert import *
from comm.unit.unit import kb_hau


def seebeck_spb(eta,r = 0.0):
    """
    using the formulas in j mat chem c 4 209-214 (2016)
    """
    return kb_hau*((2.5 + r) * fdk( 1.5+ r, eta)/((1.5+ r)*fdk(r+0.5, eta))- eta)

def Fermi_Dirac(ene,Ef,T):
    return 1./(1.0 + exp((ene-Ef)/(kb_hau*T)))

def dos_1band(m,ene,edge,type_):
    """
    ene: Hartree
    """
    if type_ == 'vb':
        if ene > edge:
            return 0.
        else:
            return 1./(2*pi**2) * (2*m)**1.5 * (edge-ene)**0.5
    elif type_ == 'cb':
        if ene < edge:
            return 0.
        else:
            return 1./(2*pi**2) * (2*m)**1.5 * (ene-edge)**0.5

def Lorentz_function(e,gamma,Ed):
    return 1./pi*(0.5*gamma/((e-Ed)**2+(0.5*gamma)**2))

def get_carrier(T,dos,Ef,time=10,type_='vb'):
    e_min = Ef - time*kb_hau*T
    e_max = Ef + time*kb_hau*T
    def g(ene):
        if type_ == 'vb':
            return (1-Fermi_Dirac(ene,Ef,T))*dos(ene)
        elif type_ == 'cb':
            return Fermi_Dirac(ene,Ef,T)*dos(ene)
    return quad(g,e_min,e_max)[0]
    
def get_Ef_from_carrier(T,n,dos):
    def carr_diff(Ef):
        return (get_carrier(T,dos,Ef)-n)**2
    out = fsolve(carr_diff, 0.,full_output=True)
    Ef = out[0][0]
    print T/Kelvin_to_hau, 'K' , Ef / eV_to_hau ,get_carrier(T,dos,Ef)/(1.e6*metre_to_hau**(-3))
    return out[0][0]



def main(plot_dos=True):
    from matplotlib import pyplot as plt
    fig,(ax,ax_dos) = plt.subplots(2,1)
    ns = [1.e+20,1.5e+20,3.7e+20,1.1e+21]
    ns = [i * 1.e6 * metre_to_hau**(-3) for i in ns]
    def add_dos(m,edge,n0,gamma,Ed):
        def f(e):
            band_vb = dos_1band(m,e,edge,'vb')
            band_cb = dos_1band(m,e,gap, 'cb')
            band_LF = n0*Lorentz_function(e,gamma,Ed)
            return band_vb+band_LF
        return f


    m = 0.5
    edge = 0.0 * eV_to_hau
    gap = 0.5 * eV_to_hau
    n0 = 1.3e+21 * 1.e6 * metre_to_hau**(-3)
    gamma = 0.05 * eV_to_hau
    Ed = -0.2 * eV_to_hau
    dos = add_dos(m, edge, n0,gamma,Ed)
    es = [i*0.001 for i in range(-50,50)]
    ax_dos.plot([i/eV_to_hau for i in es],[dos(i) for i in es])
    ax_dos.set_xlabel('energy (eV)')
    ax_dos.set_ylabel('density of state')
        
    for n in ns:
        Ss = []
        Ts = range(300,1350,50)
        for T in Ts:
            T = T * Kelvin_to_hau
            Ef = get_Ef_from_carrier(T,n,dos)
            eta = (edge - Ef)/(kb_hau*T)
            S = seebeck_spb(eta,r=0)
            S = S * Kelvin_to_hau/Volt_to_hau * 1e+6
            Ss.append(S)
        n_label = n/(1.e6*metre_to_hau**(-3))
        ax.plot(Ts,Ss,label='%.1e' % n_label)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Seebeck ($\mu V/K$)')
    ax.set_ylim((0,500))
    ax.legend()
    fig.savefig('Seebeck_vs_Temp.svg')

    

if __name__=='__main__':
    main()
