from pymatgen.core.physical_constants import BOLTZMANN_CONST, PLANCK_CONSTANT, ELECTRON_CHARGE,ELECTRON_MASS
from scipy.integrate import quad
from scipy.optimize import fsolve, curve_fit
from fdint import fdk
import numpy as np

hbar = PLANCK_CONSTANT/(2*np.pi)
def seebeck_spb(eta,Lambda=0.5):
    return BOLTZMANN_CONST/ELECTRON_CHARGE * ((2. + Lambda) * fdk( 1.+ Lambda, eta)/ 
                 ((1.+Lambda)*fdk(Lambda, eta))- eta) * 1e+6

def eta_from_seebeck(seeb,Lambda):
    """ takes a value of seebeck and adjusts seebeck until it's equal
    Returns: eta (reduced chemical potential)"""
    out = fsolve(lambda x: (seebeck_spb(x,Lambda) - abs(seeb)) ** 2, 1.,full_output=True)
    return out[0][0]


def seebeck_eff_mass_from_carr(eta, n, T, Lambda):
    """  eta in kB*T units, n in cm-3, T in K
        returns mass in m0 units """
    return (2 * np.pi**2 * abs(n) * 10 ** 6 / (fdk(0.5,eta))) ** (2. / 3)\
            / (2 * ELECTRON_MASS * BOLTZMANN_CONST * T / (PLANCK_CONSTANT/2/np.pi) ** 2)

                
def seebeck_eff_mass_from_seebeck_carr(seeb, n, T, Lambda):
    eta = eta_from_seebeck(seeb,Lambda)
    mass = seebeck_eff_mass_from_carr(eta, n, T, Lambda)
    return mass

def calculate_eff_mass(slope,n):
    return slope*1e-6 * (3.*ELECTRON_CHARGE*PLANCK_CONSTANT**2/(8.*np.pi**2*BOLTZMANN_CONST**2)) * (3.*n*1e+6/np.pi)**(2./3)/ELECTRON_MASS
