from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.outputs import Vasprun 
from math import ceil,floor,sin,cos,pi
import numpy as np

def formate_chemical_formula(name):
    """
    Converte a chemical formula to be the one for matplotlib plot.
    For example, convert BaBiTe2.9Se0.1 to $BaBiTe_{2.9}Se_{0.1}$, which
        is formated for matplotlib plot
    Args:
        name: the chemical formula you want to convert, e.g. Yb14Al0.25Mn0.75Sb11
              return $Yb_{14}Al_{0.25}Mn_{0.75}Sb_{11}$
    """
    new_name='';new_name=new_name+name[0]
    for i in range(1,len(name)):
        if name[i].isdigit() and not (name[i-1].isdigit() or name[i-1]=='.'):
            new_name=new_name+'$_{'+name[i]
        elif not (name[i].isdigit() or name[i]=='.') and name[i-1].isdigit():
            new_name=new_name+'}$'+name[i]
        else:
            new_name=new_name+name[i]
    if name[-1].isdigit():
        new_name=new_name+'}$'
    return new_name

def digital_deleted_str(name):
    """
    return a string with all digital deleted
    """
    results=''
    for i in name:
        if not i.isdigit():
            results+=i
    return results


def get_rotation_matrix_from_axis_angle(axis_vec,angle,det):
    """
    Args:
        axis: rotation axis 
        angle: rotation angle
    """
    angle = angle*pi/180.
    x,y,z=axis_vec
    m11=cos(angle) + (1-cos(angle))*(x**2)
    m12=(1-cos(angle))*x*y-sin(angle)*z
    m13=(1-cos(angle))*x*z + sin(angle)*y
    m21=(1-cos(angle))*y*x + sin(angle)*z
    m22=cos(angle) + (1-cos(angle))*(y**2)
    m23=(1-cos(angle))*y*z - sin(angle)*x
    m31=(1-cos(angle))*z*x - sin(angle)*y
    m32=(1-cos(angle))*z*y + sin(angle)*x
    m33=cos(angle)+ (1-cos(angle))*(z**2)
    matrix = np.array([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])
    if det == -1:
        matrix = -1*matrix
    return matrix

def check_rotation_two_structs(stru1,stru2):
    n1 = stru1.num_sites; n2 = stru2.num_sites
    ref = stru1 if n1<=n2 else stru2
    comp = stru2 if ref == stru1 else stru1
    ref_latt = ref.lattice
    rotation = False
    for out in comp.sites:
        out_specie =  out.specie
        out_cart = out.coords
        match = False
        for host in ref.sites:
            host_specie = host.specie
            host_cart = host.coords
            if host_specie != out_specie:
                continue
            diff_cart = out_cart - host_cart
            diff_frac = ref_latt.get_fractional_coords(diff_cart)
            vec_match = True
            for i in diff_frac:
                diff = abs(i-round(i,0))
                if diff > 0.01:
                    vec_match = False
            if vec_match:
                match = True
                break
        if not match:
            return True
    return False
            

