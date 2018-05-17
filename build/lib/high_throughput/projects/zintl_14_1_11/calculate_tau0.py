from math import pi
from comm.unit.unit_convert import * 

system = 'Yb14ZnSb11'

if system == 'Yb14MnSb11':
    Cll = 682.6  #kBar
    DP = 2.317 #eV 
    m_hau = 1.23 # me
if system == 'Yb14ZnSb11':
    Cll = 690.3  #kBar
    DP =  8.8 # eV
    m_hau = 0.88 # me


T0 = 300 # K
kb = 8.6173303e-5 # eV/K
kb_hau = kb*eV_to_hau/Kelvin_to_hau

Cll_hau = Cll* 1.e+8 * Pascal_to_hau
DP_hau = DP * eV_to_hau
T0_hau = T0 * Kelvin_to_hau

tau0_hau = pi * Cll_hau/ (2**0.5 * DP_hau**2 * m_hau**(1.5)) * (kb_hau*T0_hau)**(-3./2)

tau0_SI = tau0_hau/second_to_hau

tau0_fs = tau0_SI*1.e+15

print 'tau0 is %f fs' % tau0_fs
