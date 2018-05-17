#!/usr/bin/env python

"""
TODO create a VaspInputSet instead?
"""

__author__ = "Geoffroy Hautier"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "November 4, 2012"

from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.inputs import Kpoints,Poscar
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from pymatgen.core.structure import Structure
from monty.json import jsanitize
from high_throughput.defects.defectsmaker import IntrinsicDefectsMaker,DefectRemaker,ExtrinsicDefectsMaker
from high_throughput.defects.utils.util import closest_defect_distance
from high_throughput.config import *
import json
import os
import numpy as np

def IncarSetup(incar,defect_type,hse=False):
    incar['ALGO']="Normal" 
    incar['LWAVE']=False
    incar['LCHARG']=False
    incar['ISMEAR']=0
    if defect_type == 'bulk':
        incar['EDIFF']=0.00001
        incar['LVHAR']=True
    elif defect_type == 'dielectric':
        incar['IBRION']=8
        incar['LEPSILON']=True
        incar['LPEAD']=True
        incar['EDIFF']=0.0001
        incar['SIGMA']=0.01
        del incar['NSW'], incar['LVHAR'], incar['LAECHG']
    else:
        incar['IBRION']=2
        incar['ISIF']=0
        incar['ISPIN']=2
        incar['EDIFF']=0.0001
        #incar['EDIFFG']=0.001
        incar['SIGMA']=0.05
        incar['LVHAR']=True
        incar['LORBIT']=11
    if hse:
       incar['LHFCALC']=True
       incar["ALGO"]="All"
       incar["HFSCREEN"]=0.2
       incar["PRECFOCK"]="Fast"
       incar["KPAR"] = 2
    return incar


def make_vasp_defect_files(dictio, path_base, task_id, compo, hse=False, encut_redo=False):
    """
    simple static method creating VASP files ready for defect computations
    Args:
        dictio:
            the defects data as a dictionnary
        path_base:
            where do we write the files
        task_id:
            some id of the bulk computed data
        compo:
            Composition of the bulk computed data
        hse:
            hse run or not
    """
    count=1
    for site in dictio:
        #### bulk ####
        if site['short_name']=='bulk':
            bulk_unitcell = site['bulk_unitcell'].as_dict()
            uc_type = site['uc_type']
            sc_type = site['sc_type']
            for s in site['supercells']:
                defect_dist = round(closest_defect_distance(s['structure']),2)
                bulk_info = '%s_%s_%s' % (uc_type, sc_type, str(defect_dist))
                dict_transf = {'history':[{'source':task_id,'unitcell_type':site['uc_type']}],
                        'defect_type':'bulk','supercell':s['size']}
                structs = {'bulk_unitcell':bulk_unitcell}
                dict_params = MPStaticSet(s['structure']).all_input
                incar_init=dict_params['INCAR']
                incar=IncarSetup(incar_init,'bulk',hse)
                if encut_redo:
                    enmax = round(max([i.PSCTR['ENMAX'] for i in dict_params['POTCAR']])*1.3)
                    incar['ENCUT'] = int(enmax)
                if hse:
                    kpoint=Kpoints.gamma_automatic()
                else:
                    kpoint=Kpoints.monkhorst_automatic()
                path=path_base+"/"+str(task_id)+'_'+compo.reduced_formula+'/bulk/'+bulk_info
                os.makedirs(path)
                f=open(path+"/transformations.json",'w')
                f.write(json.dumps(jsanitize(dict_transf)))
                g=open(path+"/structures.json",'w')
                g.write(json.dumps(jsanitize(structs)))
                dict_params['POTCAR'].write_file(path+"/POTCAR")
                incar.write_file(path+"/INCAR")
                kpoint.write_file(path+"/KPOINTS")
                dict_params['POSCAR'].write_file(path+"/POSCAR")
            continue

        #### dielectric constants ####
        if site['short_name']=='dielectric':
            dict_transf = {'history':[{'source':task_id,'unit_cell':site['uc_type']}],'defect_type':'dielectric'}
            dict_params = MPStaticSet(site['structure']).all_input
            incar=dict_params['INCAR']
            kpoints=Kpoints.automatic_gamma_density(site['structure'],2000)
            try:
                bs=m.get_bandstructure_by_material_id(task_id)
                if not bs.is_spin_polarized:
                    incar['ISPIN']=1
                else:
                    incar['ISPIN']=2
            except:
                incar['ISPIN']=1
            incar = IncarSetup(incar,'dielectric',hse)
            if encut_redo:
                enmax = round(max([i.PSCTR['ENMAX'] for i in dict_params['POTCAR']])*1.3)
                incar['ENCUT'] = int(enmax)
            path=path_base+"/"+str(task_id)+'_'+compo.reduced_formula+"/"+'dielectric/'+site['uc_type']
            os.makedirs(path)
            f=open(path+"/transformations.json",'w')
            f.write(json.dumps(jsanitize(dict_transf)))
            dict_params['POTCAR'].write_file(path+"/POTCAR")
            incar.write_file(path+"/INCAR")
            kpoints.write_file(path+"/KPOINTS")
            dict_params['POSCAR'].write_file(path+"/POSCAR")
            continue
            
        #### defects ####
        uc_type = site['uc_type']
        sc_type = site['sc_type']
        for charge in site['charges']:
            uc = site['bulk_unitcell'].copy()
            bulk_unitcell = uc.as_dict()
            for s in site['supercells']:
                defect_dist = round(closest_defect_distance(s['structure']),2)
                defect_info = '%s_%s_%s' % (uc_type, sc_type, str(defect_dist))
                uc.make_supercell(s['size'],to_unit_cell=True)
                bulk_supercell = uc.as_dict()
                dict_transf={'history':[{'source':task_id,'unit_cell':site['uc_type']}], 'compo': compo.as_dict(), 
                    'defect_type': site['short_name'], 'defect_site': site['unique_sites'].as_dict(), 
                    'charge': charge, 'supercell': s['size']}
                dict_params=MPRelaxSet(s['structure']).all_input
                try:
                    defect_no_relaxation = s['struct_no_move'].as_dict()
                except:
                    defect_no_relaxation = s['structure'].as_dict()
                structs = {'bulk_unitcell':bulk_unitcell,'bulk_supercell':bulk_supercell, 
                        'defect_no_relaxation':defect_no_relaxation}
                incar=dict_params['INCAR']
                incar=IncarSetup(incar,'defect',hse)
                if encut_redo:
                    enmax = round(max([i.PSCTR['ENMAX'] for i in dict_params['POTCAR']])*1.3)
                    incar['ENCUT'] = int(enmax)
                if hse:
                    kpoint=Kpoints.gamma_automatic()
                else:
                    kpoint=Kpoints.monkhorst_automatic()
                path=path_base+"/"+str(task_id)+'_'+compo.reduced_formula+ \
                    '/'+str(site['short_name'])+"/"+"charge"+str(charge)+'/'+defect_info
                os.makedirs(path)
                f=open(path+"/transformations.json",'w')
                f.write(json.dumps(jsanitize(dict_transf)))
                g=open(path+"/structures.json",'w')
                g.write(json.dumps(jsanitize(structs)))
                comp=s['structure'].composition
                sum_elec=0
                elts=set()
                for p in dict_params['POTCAR']:
                    if p.element not in elts:
                        sum_elec+=comp.as_dict()[p.element]*p.nelectrons
                        elts.add(p.element)
                if charge != 0:
                    incar['NELECT']=sum_elec-charge
                dict_params['POTCAR'].write_file(path+"/POTCAR")
                incar.write_file(path+"/INCAR")
                kpoint.write_file(path+"/KPOINTS")
                dict_params['POSCAR'].write_file(path+"/POSCAR")
                count=count+1
                f.close()
                g.close()

def generate_intrinsic_defect_input_files(source,hse=False,path='./',vac=True,sub=True,inter=False,
         min_dis=13.0,extrapolation=False,encut_redo=False,**kws):
    """
    Args:
        source: material id in Materials Project or structure file e.g. 'mp-149' or POSCAR
        hse: whether hse
        path: path to save outputs
        vac, sub, inter: True or False 
        min_dis: the minimum distance between defects
        extrapolation: True or Faslse, whether generating defects with three different siez supercell
        encut_redo: True will setting the ENCUT to be just 30% larger than the ENMAX in POTCAR
                    False will set ENCUT 520 eV
    """
    try:
        struct=m.get_structure_by_material_id(source)
        mpid=source
    except:
        struct=Structure.from_file(source)
        mpid=kws['label']
    interstitials_extr_elts=[]
    
    make_vasp_defect_files(IntrinsicDefectsMaker(struct,vac,sub,inter,
        min_dis,extrapolation,**kws).defects, path, mpid, struct.composition,hse, encut_redo)

def generate_extrinsic_defect_input_files(source,hse=False,path='./',subst={},inter_ext_elts=[],
        min_dis=13.0,extrapolation=False,encut_redo=False,**kws):
    """
    Args:
        source: material id in Materials Project or structure file e.g. 'mp-149' or POSCAR
        hse: whether hse
        path: path to save outputs
    """
    try:
        struct=m.get_structure_by_material_id(source)
        mpid=source
    except:
        struct=Structure.from_file(source)
        mpid=kws['label']
    inter_ext_elts=[Element(i) for i in inter_ext_elts]
    make_vasp_defect_files(ExtrinsicDefectsMaker(struct,subst,inter_ext_elts,min_dis,extrapolation).defects, path, mpid, struct.composition,hse,encut_redo)

def generate_input_files_with_larger_defect_distance(mpid, bulk_sc, defect,dfc_name, dfc_charge, 
        dfc_dist_add=2.0,hse=False,path='./'):
    defects = DefectRemaker(bulk_sc,defect,dfc_name,dfc_charge,dfc_dist_add).defects
    make_vasp_defect_files(defects, path,mpid,bulk_sc.composition,hse)



