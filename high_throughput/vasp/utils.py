from pymatgen.io.vasp.outputs import Vasprun,Outcar
from pymatgen.io.vasp.inputs import Kpoints,Potcar
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.boltztrap import BoltztrapRunner
from pymatgen.io.vasp.sets import MPRelaxSet,MPHSEBSSet,MPStaticSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.groups import PointGroup,SymmOp
from monty.json import jsanitize
from high_throughput.utils.utils import get_rotation_matrix_from_axis_angle
from scipy.optimize import fsolve
from high_throughput.config import *
import numpy as np
import json
import os


class InputsGenerators(object):
    def __init__(self,stru_source):
        try:
            stru = m.get_structure_by_material_id(stru_source)
            self.mpid = stru_source
        except:
            stru = Structure.from_file(stru_source)
            self.mpid = ''
        spa = SpacegroupAnalyzer(stru)
        self.stru = spa.get_primitive_standard_structure()
        self.stru_prim = spa.get_primitive_standard_structure()
        self.stru_conv = spa.get_conventional_standard_structure()
        formula = self.stru.composition.reduced_formula
        self.name = self.mpid+'_'+formula if self.mpid else formula
        try:
            self.bs = m.get_bandstructure_by_material_id(self.mpid)
            self.is_spin_polarized = self.bs.is_spin_polarized
            bs_exist = True
            self.pbevbm = self.bs.get_vbm()
            self.pbecbm = self.bs.get_cbm()
            self.kpbevbm = self.pbevbm['kpoint'].frac_coords
            self.kpbecbm = self.pbecbm['kpoint'].frac_coords 
        except:
            bs_exist = False
            self.kpbevbm = None
            self.kpbecbm = None
            self.bs = None
        if not bs_exist:
            self.is_spin_polarized = False
            for elem in stru.composition.elements:
                if elem.is_transition_metal:
                    self.is_spin_polarized = True
                    break

    def elastic_moduli(self,output_dir='./'):
        stru = self.stru_prim.copy()
        inputs = MPStaticSet(stru).all_input
        path = os.path.join(output_dir,self.name)
        transf = {'history':[{'source':self.mpid}],'defect_type':'elastic_moduli'}
        incar = inputs['INCAR']
        if self.is_spin_polarized:
            incar['ISPIN']=2
        else:
            incar['ISPIN']=1
        incar['IBRION']=6
        incar['ISIF'] = 3
        incar['EDIFF'] = 1e-6
        incar['ISMEAR']=0
        incar['POTIM']=0.015
        del incar['NSW'], incar['LVHAR'], incar['LAECHG']
        os.mkdir(path)
        f=open(path+"/transformations.json",'w')
        f.write(json.dumps(jsanitize(transf)))
        enmax = round(max([i.PSCTR['ENMAX'] for i in inputs['POTCAR']])*1.3)
        incar['ENCUT'] = int(enmax)
        inputs['POTCAR'].write_file(path+"/POTCAR")
        incar.write_file(path+"/INCAR")
        inputs['KPOINTS'].write_file(path+"/KPOINTS")
        inputs['POSCAR'].write_file(path+"/POSCAR")

    def deformation_potential_constant_L(self,out_path='./'):
        stru = self.stru_conv.copy()
        if not set(stru.lattice.angles)==set([90.0]):
            return
        kvbm_cart = self.stru_prim.lattice.reciprocal_lattice.get_cartesian_coords(self.kpbevbm)
        kcbm_cart = self.stru_prim.lattice.reciprocal_lattice.get_cartesian_coords(self.kpbecbm)
        kvbm_conv = stru.lattice.reciprocal_lattice.get_fractional_coords(kvbm_cart)
        kcbm_conv = stru.lattice.reciprocal_lattice.get_fractional_coords(kcbm_cart)
        root = os.path.join(out_path,self.name+'_deformation_potential_constants')
        strains=[-0.01, -0.005, 0.0, 0.005, 0.01]
        inputs = MPRelaxSet(stru).all_input
        kpts = inputs['KPOINTS']
        incar = inputs['INCAR']
        if self.is_spin_polarized:
            incar['ISPIN']=2
        else:
            incar['ISPIN']=1
        enmax = round(max([i.PSCTR['ENMAX'] for i in inputs['POTCAR']])*1.3)
        incar['ENCUT'] = int(enmax)
        incar['ALGO']="Normal"
        incar['ISMEAR']= 0
        incar['EDIFF'] = 1e-5
        incar['EDIFFG'] = -0.005
        incar['IBRION']= 2 
        incar['ISIF'] = 2
        incar['SIGMA'] = 0.05
        axes = {'a':[1,0,0],'b':[0,1,0],'c':[0,0,1]}
        os.mkdir(root)
        for i in ['a','b','c']:
            path1 = os.path.join(root,i)
            os.mkdir(path1)
            for strain in strains:
                if strain < 0:
                    path2 = os.path.join(path1,str(strain)[1:]+'-')
                else:
                    path2 = os.path.join(path1,str(strain))
                os.mkdir(path2)
                stru_tmp=stru.copy()
                stru_tmp.apply_strain(np.array(axes[i])*strain)
                incar.write_file(path2+'/INCAR')
                inputs['POTCAR'].write_file(path2+'/POTCAR')
                inputs['KPOINTS'].write_file(path2+'/KPOINTS')
                stru_tmp.to('poscar',path2+'/POSCAR')

    def deformation_potential_constant_V(self,incar_add={},out_path='./'):
        f = lambda x:3.*x+3.*x**2+x**3
        stru = self.stru_prim.copy()
        root = os.path.join(out_path,self.name+'_DP')
        dVs=[-0.01, -0.005, 0.0, 0.005, 0.01]
        inputs = MPRelaxSet(stru).all_input
        kpts = inputs['KPOINTS']
        incar = inputs['INCAR']
        if self.is_spin_polarized:
            incar['ISPIN']=2
        else:
            incar['ISPIN']=1
        enmax = round(max([i.PSCTR['ENMAX'] for i in inputs['POTCAR']])*1.3)
        incar['ENCUT'] = int(enmax)
        incar['ALGO']="Normal"
        incar['ISMEAR']= 0
        incar['EDIFF'] = 1e-5
        incar['EDIFFG'] = -0.005
        incar['IBRION']= 2 
        incar['ISIF'] = 2
        incar['SIGMA'] = 0.05
        for i in incar_add:
            incar[i] = incar_add[i]
        os.mkdir(root)
        for dV in dVs:
            if dV < 0:
                path = os.path.join(root,str(dV)[1:]+'-')
            else:
                path = os.path.join(root,str(dV))
            g = lambda x:f(x)-dV
            dL = fsolve(g,0.0)
            os.mkdir(path)
            stru_tmp=stru.copy()
            stru_tmp.apply_strain(dL)
            incar.write_file(path+'/INCAR')
            inputs['POTCAR'].write_file(path+'/POTCAR')
            inputs['KPOINTS'].write_file(path+'/KPOINTS')
            stru_tmp.to('poscar',path+'/POSCAR')

    def dielectric(self,uc_type='prim',output_dir='./'):
        if uc_type == 'prim':
            uc_type = 'primitive'
            stru = self.stru_prim.copy()
        elif uc_type == 'conv':
            uc_type = 'conventional'
            stru = self.stru_conv.copy()
        path = os.path.join(output_dir,self.name)
        inputs = MPStaticSet(stru).all_input
        transf = {'history':[{'source':self.mpid,'unit_cell':uc_type}],'defect_type':'dielectric'}
        incar = inputs['INCAR']
        kpoints = Kpoints.automatic_gamma_density(stru,2000)
        if self.is_spin_polarized:
            incar['ISPIN']=2
        else:
            incar['ISPIN']=1
        incar['IBRION']=8
        incar['LEPSILON']=True
        incar['LPEAD']=True
        incar['EDIFF']=0.000001
        incar['LWAVE']=False
        incar['LCHARG']=False
        incar['ISMEAR']=0
        incar['ALGO']="Normal"
        incar['SIGMA']=0.01
        del incar['NSW'], incar['LVHAR'], incar['LAECHG']
        os.mkdir(path)
        f=open(path+"/transformations.json",'w')
        f.write(json.dumps(jsanitize(transf)))
        inputs['POTCAR'].write_file(path+"/POTCAR")
        incar.write_file(path+"/INCAR")
        kpoints.write_file(path+"/KPOINTS")
        inputs['POSCAR'].write_file(path+"/POSCAR")

    def PBEsol_relaxation(self,output_dir='./'):
        path = os.path.join(output_dir,self.name)
        inputs = MPRelaxSet(self.stru).all_input
        incar = inputs['INCAR']
        if self.is_spin_polarized:
            incar['ISPIN']=2
        else:
            incar['ISPIN']=1
        enmax = round(max([i.PSCTR['ENMAX'] for i in inputs['POTCAR']])*1.3)
        incar['ENCUT'] = int(enmax)
        incar["SYSTEM"]=self.name + "_PBEsol_Relax"
        incar["ALGO"] = "Normal"
        incar["EDIFF"] = 0.0001
        incar["EDIFFG"] = 0.001
        incar["ISMEAR"] = 0
        incar["GGA"] = "PS"
        incar["NSW"] = 99
        profile = {'source':self.mpid,'calculation':'PBEsol_relaxation'}
        os.mkdir(path)
        f=open(os.path.join(path,'profile.json'),'w')
        f.write(json.dumps(jsanitize(profile)))
        inputs['KPOINTS'].write_file(path+"/KPOINTS")
        inputs['POSCAR'].write_file(path+"/POSCAR")
        inputs['POTCAR'].write_file(path+"/POTCAR")
        incar.write_file(path+"/INCAR")

    def HSE_gap_with_shot(self,stru=None,output_dir='./',**kws):
        path = os.path.join(output_dir,self.name)
        incar_add={}
        if self.is_spin_polarized:
            incar_add["ISPIN"] = 2
        else:
            incar_add["ISPIN"] = 1
        incar_add["EDIFF"] = 1E-4
        incar_add["SYSTEM"] = self.name+"_HSE_StaticRun"
        incar_add["IBRION"] = -1
        incar_add["HFSCREEN"] = 0.2
        incar_add["PREC"] = "Accurate"
        incar_add["PRECFOCK"] = "Fast"
        incar_add["KPAR"] = 2
        if stru is None:
            dict_params = MPHSEBSSet(self.stru, user_incar_settings=incar_add,
                        reciprocal_density=100).all_input
        else:
            dict_params = MPHSEBSSet(stru, user_incar_settings=incar_add,
                    reciprocal_density=100).all_input
        incar=dict_params["INCAR"]
        kpoints = dict_params["KPOINTS"]
        ks=[]
        if self.kpbevbm is None:
            try:
                self.kpbevbm=np.array([float(i) for i in kws['kvbm']])
                self.kpbecbm=np.array([float(i) for i in kws['kcbm']])
            except:
                print 'Error: PBE bandstruct is a metal, please input kvbm and kcbm by hand!\n \
                       For example, kvbm=[0,0,0],kcbm=[0.5,0.5,0.0]'
                return
        if set(self.kpbevbm==self.kpbecbm)==set([True]):
            ks.append(self.kpbevbm)
        else:
            ks.append(self.kpbevbm);ks.append(self.kpbecbm)
        for k in ks:
            exist = False
            for k_host in kpoints.kpts:
                if (np.array(k)==np.array(k_host)).all():
                    exist=True
                    break
            if exist:
                continue
            kpoints.kpts.append(k)
            kpoints.kpts_weights.append(0)
            kpoints.labels.append('')
        kpoints.num_kpts = len(kpoints.kpts_weights)
        transf = {'history':[{'source':self.mpid,'unit_cell':'primitive'}],'defect_type':'HSE_gap'}
        os.mkdir(path)
        f=open(os.path.join(path,'transformations.json'),'w')
        f.write(json.dumps(jsanitize(transf)))
        incar.write_file(path+"/INCAR")
        dict_params["POTCAR"].write_file(path+"/POTCAR")
        dict_params["POSCAR"].write_file(path+"/POSCAR")
        kpoints.write_file(path+"/KPOINTS")

class BoltztrapRunFromDir(object):
    def __init__(self,root_dir):
        xml = Vasprun(os.path.join(root_dir,'vasprun.xml'))
        self.bs = xml.get_band_structure()
        self.nelec = xml.parameters['NELECT']
    def transport_run(self,doping=[],output_dir='./'):
        BoltztrapRunner(self.bs,self.nelec,doping=doping,symprec=0.1).run(output_dir)
    def band_interpolate(self,output_dir='./'):
        BoltztrapRunner(self.bs,self.nelec,run_type='BANDS').run(output_dir)

class OutcarRead(object):
    def __init__(self,outcar_file):
        self.out = Outcar(outcar_file)
        with open(outcar_file) as f:
            self.OUTCAR = f.read().split('\n')

    def read_elastic_constant(self):
        string = 'TOTAL ELASTIC MODULI'
        for line_num, line in enumerate(self.OUTCAR):
            if string in line:
                start_line = line_num
                break
        block = [i.split()[1:] for i in self.OUTCAR[start_line+3:start_line+9]]
        C = [[float(j) for j in i] for i in block]
        return np.mat(C)

class ElasticProperty(object):
    """
    get the elastic properties of porycrystal from elastic constant of single crystal
    See paper: phys. stat. sol. (a) 99, 423, (1987)
    """
    def __init__(self,elastic_constant):
        self.elastic_constant = elastic_constant
        self.elastic_compliance = self.get_elastic_compliance()
        self.Lame_constants_Voigt = self.get_Lame_constants_Voigt()
        self.Lame_constants_Reuss = self.get_Lame_constants_Reuss()

    def get_elastic_compliance(self):
        """
        the inversed matrix of elastic constant
        """
        S = np.linalg.inv(self.elastic_constant)
        return S

    def get_Lame_constants_Voigt(self):
        C = self.elastic_constant
        I1 = 0.0; I1_ = 0.0
        for i in [0,1,2]:
            for j in [0,1,2]:
                I1 += C[i,j]
        for i in range(6):
            for i in [0,1,2]:
                weight = 1.0
            else:
                weight = 2.0
            I1_ += weight * C[i,i]
        mu = 1./30 * (3.*I1_ - I1)
        lambda_ = 1./15 * (2*I1 - I1_)
        return {'mu':mu, 'lambda': lambda_}

    def get_Lame_constants_Reuss(self):
        S = self.elastic_compliance
        J1 = 0.0; J1_ = 0.0
        for i in [0,1,2]:
            for j in [0,1,2]:
                J1 += S[i,j]
        for i in range(6):
            for i in [0,1,2]:
                weight = 1.0
            else:
                weight = 2.0
            J1_ += weight * S[i,i]
        x = J1_/J1
        mu = (15./J1) * (1. / (6*x - 2))
        lambda_ = 1./J1 -(2./3)*mu
        return {'mu':mu, 'lambda':lambda_}

    def get_bulk_modulus(self,style):
        mu_V = self.Lame_constants_Voigt['mu']
        lambda_V = self.Lame_constants_Voigt['lambda']
        bulk_modulus_V = lambda_V + 2./3 * mu_V
        mu_R = self.Lame_constants_Reuss['mu']
        lambda_R = self.Lame_constants_Reuss['lambda']
        bulk_modulus_R = lambda_R + 2./3 * mu_R
        bulk_modulus_H = (bulk_modulus_V + bulk_modulus_R)/2.0

        if style == 'Voigt':
            return bulk_modulus_V
        elif style == 'Reuss':
            return bulk_modulus_R
        elif style == 'Hill':
            return bulk_modulus_H

    def get_Young_modulus(self,style):
        lambda_V = self.Lame_constants_Voigt['lambda']
        lambda_R = self.Lame_constants_Reuss['lambda']
        bulk_modulus_R = self.get_bulk_modulus('Reuss')
        bulk_modulus_V = self.get_bulk_modulus('Voigt')
        Young_modulus_R = 9.*bulk_modulus_R *(bulk_modulus_R - lambda_R) /(3.*bulk_modulus_R-lambda_R)
        Young_modulus_V = 9.*bulk_modulus_V *(bulk_modulus_V - lambda_V) /(3.*bulk_modulus_V-lambda_V)
        Young_modulus_H = (Young_modulus_R+Young_modulus_V)/2.
        if style == 'Voigt':
            return Young_modulus_V
        elif style == 'Reuss':
            return Young_modulus_R
        elif style == 'Hill':
            return Young_modulus_H

def same_ks(vec1_frac,vec2_frac):
    """
    Are vec1 and vec2 the same just after translation 
    """
    diff = np.array(vec1_frac) - np.array(vec2_frac)
    for i in diff:
        if round(i,2)%1:
            return False
    return True

def get_symops_from_OUTCAR(outfile):
    """
    add s1 operator also if it doesn't exist
    """
    with open(outfile,'r') as f:
        data=f.read().split('\n')
        start_line=0
        end_line=0
        for line_num,line in enumerate(data):
            if 'Space group operators:' in line:
                start_line = line_num
                end_line = 0
            elif start_line and not end_line and not len(line.split()):
                end_line = line_num
    ops=[]
    for line_num in range(start_line+2,end_line):
        cont=data[line_num].split()
        op={}
        op['det']=int(float(cont[1]))
        op['angle']=float(cont[2])
        op['axis']=[float(cont[3]),float(cont[4]),float(cont[5])]
        ops.append(op)
    ops_matrix=[]
    s1=np.array([[-1,0,0],[0,-1,0],[0,0,-1]])
    mats=[]
    for op in ops:
        mat=get_rotation_matrix_from_axis_angle(op['axis'],op['angle'],op['det'])
        mat1=np.round(np.array(mat),3)
        mats.append(mat1)
    exist = False
    for mat in mats:
        if not (mat-s1).any():
            exist = True
            break
    for mat in mats:
        op1=SymmOp.from_rotation_and_translation(mat,[0,0,0])
        ops_matrix.append(op1)
        if exist:
            continue
        mat2=np.dot(s1,mat)
        op2=SymmOp.from_rotation_and_translation(mat2,[0,0,0])
        ops_matrix.append(op2)
    return ops_matrix

def print_rotation_matrix_for_boltztrap(outputfile):
    """
    print rotation matrix for boltztrap run from OUTCAR
    """
    poscar=outputfile.replace('OUTCAR','POSCAR')
    stru=Structure.from_file(poscar)
    lat_mat=stru.lattice.reciprocal_lattice.matrix
    inv=np.linalg.inv(lat_mat)
    ops=get_symops_from_OUTCAR(outputfile)
    rotations=[]
    for op in ops:
        m12=np.dot(lat_mat, op.rotation_matrix)
        m=np.round(np.dot(m12,inv),2)
        int_m = np.int_(m) 
        rotations.append(int_m)
        print int_m[0][0], int_m[0][1], int_m[0][2]
        print int_m[1][0], int_m[1][1], int_m[1][2]
        print int_m[2][0], int_m[2][1], int_m[2][2], '\n'
    return rotations

def get_equivalent_kpoints(stru,kvec_frac,outcarfile):
    """
    Get all equivalent kpts to the k
    """
    ops=get_symops_from_OUTCAR(outcarfile)
    reci_latt = stru.lattice.reciprocal_lattice
    kvec_cart = reci_latt.get_cartesian_coords(kvec_frac)
    kpts=[]
    for op in ops:
        kpt_cart = op.apply_rotation_only(kvec_cart)
        kpt_frac = latt.get_fractional_coords(kpt_cart)
        kpts.append(kpt_frac)
    reduced_kpts=[]
    for k in kpts:
        match = False
        for red_k in reduced_kpts:
            if same_ks(red_k,k):
                match = True
                break
        if not match:
            reduced_kpts.append(k)
    return reduced_kpts
       
def are_equivalent_kpoints(stru,vec1_frac,vec2_farc):
    """
    check whether vec1 and vec2 are the equivalent kpoints.
    """
    equi_vec1 = get_equivalent_kpoints(stru,vec1_frac)
    equivalent = False
    for k1 in equi_vec1:
        if same_ks(k1,vec2_farc):
            equivalent = True
            break
    return equivalent


