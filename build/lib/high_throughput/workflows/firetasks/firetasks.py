from custodian.custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler, UnconvergedErrorHandler,FrozenJobErrorHandler,NonConvergingErrorHandler,\
                                    MeshSymmetryErrorHandler
from custodian.vasp.jobs import VaspJob
from fireworks.core.firework import FireTaskBase, FWAction
from fireworks.utilities.fw_serializers import FWSerializable
from fireworks.utilities.fw_utilities import explicit_serialize
from monty.json import jsanitize
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import Vasprun,Kpoints,Incar,Poscar,Potcar
#from pymatgen.io.vasp.sets import MPBSHSEVaspInputSet, MPVaspInputSet,MPStaticVaspInputSet, MPNonSCFVaspInputSet
from pymatgen.io.vasp.sets import MPRelaxSet
import numpy as np
import os
import json
import time
import subprocess
from monty.os import cd


@explicit_serialize
class CommandRun(FireTaskBase, FWSerializable):
    def run_task(self,fw_spec):
        command=fw_spec['cmd']
        os.system(command)
        
@explicit_serialize
class MakeVaspDefectFile(FireTaskBase, FWSerializable):
    def run_task(self,fw_spec):
        task_id=fw_spec['task_id']
        initial_struct=fw_spec['initial_struct']
        charge=fw_spec['charge']
        site=fw_spec['site']
        s=fw_spec['s']
        hse=fw_spec['hse']
####----------------initialization over-------------
        dict_transf={'history':[{'source':task_id}], 'compo': initial_struct.composition.as_dict(), 'defect_type': site['short_name'], 'defect_site': site['unique_sites'].as_dict(), 'charge': charge, 'supercell': s['size']}
        dict_params=MPRelaxSet(s['structure']).all_input
        incar=dict_params['INCAR']
        incar['IBRION']=2
        incar['ISIF']=0
        incar['ISPIN']=2
        incar['LWAVE']=False
        incar['LCHARG']=False
        incar['EDIFF']=0.0001
	incar['EDIFFG']=0.001
        incar['ISMEAR']=0
        incar['SIGMA']=0.05
        incar['LVHAR']=True
        incar['LORBIT']=11
        incar['ALGO']="Normal"
        if hse == True:
            incar['LHFCALC']=True
            incar["ALGO"]="All"
            incar["HFSCREEN"]=0.2
            incar["PRECFOCK"]="Fast"
        kpoint=Kpoints.monkhorst_automatic()
        f=open("transformations.json",'w')
        f.write(json.dumps(jsanitize(dict_transf)))
        f.close()
        comp=s['structure'].composition
        sum_elec=0
        elts=set()
        for p in dict_params['POTCAR']:
            if p.element not in elts:
                sum_elec+=comp.as_dict()[p.element]*p.nelectrons
                elts.add(p.element)
            #print p.element
            #print comp.as_dict[p.element]
            #print p.valence
        if charge != 0:
            incar['NELECT']=sum_elec-charge
        dict_params['POTCAR'].write_file("POTCAR")
        incar.write_file("INCAR")
        kpoint.write_file("KPOINTS")
        dict_params['POSCAR'].write_file("POSCAR")
        #print Poscar(s['structure'])
        #Poscar(s['structure']).write_file(path+"POSCAR")
                

class VaspDoneError(Exception):
    def __init__(self):
        pass

class VaspReadyError(Exception):
    def __init__(self):
        pass

@explicit_serialize
class VaspRun(FireTaskBase, FWSerializable):
    def run_task(self,fw_spec):
        #workdir=fw_spec['workdir']
        vasp_cmd=fw_spec['vasp_cmd']
        #with cd(workdir):
	incar=Incar.from_file('INCAR')
	kpoints=Kpoints.from_file('KPOINTS')
	poscar=Poscar.from_file('POSCAR')
        potcar=Potcar.from_file('POTCAR')
        try:
	    out=Outcar(work_dir+'/OUTCAR')
	    if len(out.run_stats) != 7:
		raise  VaspDoneError()
        except:
    	    try:
    	        contcar=Structure.from_file('CONTCAR')
    	        os.rename('CONTCAR','POSCAR')
    	    except:
    	        pass
            job=VaspJob(vasp_cmd)
            handlers=[VaspErrorHandler(),UnconvergedErrorHandler(),FrozenJobErrorHandler(),\
                      NonConvergingErrorHandler(nionic_steps=2, change_algo=True),MeshSymmetryErrorHandler()]
            c=Custodian(handlers,[job],max_errors=10)
            c.run()
        else:
            print 'Vasp job was already done well. No need to rerun!'

@explicit_serialize
class VaspRunCheck(FireTaskBase, FWSerializable):
    def run_task(self,fw_spec):
	workdir=fw_spec['workdir']
	os.chdir(workdir)
	xml=Vasprun('./vasprun.xml')
	if not xml.converged:
            raise VaspDoneError()

@explicit_serialize
class VaspDoubleRelax(FireTaskBase, FWSerializable):
    def run_task(self,fw_spec):
        workdir=fw_spec['workdir']
        vasp_cmd=fw_spec['vasp_cmd']
        os.chdir(workdir)
        jobs=VaspJob.double_relaxation_run(vasp_cmd)
        handlers=[VaspErrorHandler(),UnconvergedErrorHandler(),FrozenJobErrorHandler(),\
                  NonConvergingErrorHandler(nionic_steps=5, change_algo=True),MeshSymmetryErrorHandler()]
        c=Custodian(handlers,jobs,max_errors=10)
        c.run()
      
@explicit_serialize
class Generate_VaspInputFiles_for_Relax_with_PBEsol(FireTaskBase, FWSerializable):
    def run_task(self, fw_spec):
        input=['INCAR','KPOINTS','POSCAR','POTCAR']
        for i in input:
            ready=os.path.isfile(i)
            if ready==False: 
                break
        if ready==False:
            material_id=fw_spec['material_id']
            structure = fw_spec["structure"]
            is_spin_polarized=fw_spec['is_spin_polarized']
            dict_params = MPVaspInputSet().get_all_vasp_input(structure)
            incar = dict_params["INCAR"]
            if is_spin_polarized:
                incar["ISPIN"] = 2
            if not is_spin_polarized:
                incar["ISPIN"] = 1
                del incar["MAGMOM"]
            incar["SYSTEM"]=material_id+"_PBEsol_Relax"
            incar["ALGO"] = "Normal"
            incar["EDIFF"] = 0.0001
            incar["EDIFFG"] = -0.01
            incar["ISMEAR"] = 0
            incar["GGA"] = "PS"
            incar["NSW"] = 99
            incar.write_file("INCAR")
            #kpoint = Kpoints.monkhorst_automatic()
            #kpoint.write_file("KPOINTS")
            dict_params["KPOINTS"].write_file("KPOINTS")
            dict_params["POTCAR"].write_file("POTCAR")
            dict_params["POSCAR"].write_file("POSCAR")
        elif ready == True:
            if os.path.isfile('CONTCAR'):
                if os.stat('CONTCAR').st_size != 0:
                    os.rename('CONTCAR','POSCAR')


@explicit_serialize
class Generate_VaspInputFiles_for_Relax_with_PBEsol_icsd(FireTaskBase, FWSerializable):
    def run_task(self, fw_spec):
	structure= fw_spec['structure']
        icsd_id=fw_spec['icsd_id']
        is_spin_polarized=fw_spec['is_spin_polarized']
        dict_params = MPVaspInputSet().get_all_vasp_input(structure)
        incar = dict_params["INCAR"]
        if is_spin_polarized:
            incar["ISPIN"] = 2
        if not is_spin_polarized:
            incar["ISPIN"] = 1
        incar["SYSTEM"]=icsd_id+"_PBEsol_Relax"
        incar["ALGO"] = "Normal"
        incar["EDIFF"] = 0.0001
        incar["EDIFFG"] = 0.001
        incar["ISMEAR"] = 0
        incar["GGA"] = "PS"
        incar["NSW"] = 99
        incar.write_file("INCAR")
        dict_params["KPOINTS"].write_file("KPOINTS")
        dict_params["POTCAR"].write_file("POTCAR")
        dict_params["POSCAR"].write_file("POSCAR")

@explicit_serialize
class Generate_VaspInputfiles_for_Static_Run(FireTaskBase,FWSerializable):
    def run_task(self,fw_spec):
        job_info_array = fw_spec['_job_info']
        prev_job_info = job_info_array[-1]
        material_id=fw_spec['material_id']
        MPStaticVaspInputSet.from_previous_vasp_run(prev_job_info['launch_dir'])
        incar=Incar.from_file('INCAR')
        incar['EDIFF']=1e-06
        incar['SYSTEM']=material_id+'_scf'
        try:
            del incar['KPOINT_BSE']
        except:
            pass
        incar.write_file('INCAR')
        
     
                
@explicit_serialize
class Generate_VaspInputfiles_for_BS(FireTaskBase,FWSerializable):
    def run_task(self,fw_spec):
        job_info_array = fw_spec['_job_info']
        Material_id=fw_spec['material_id']
        prev_job_info = job_info_array[-1]
        MPNonSCFVaspInputSet.from_previous_vasp_run(prev_job_info['launch_dir'],mode='Line')
        incar=Incar.from_file('INCAR')
        incar['SYSTEM']=Material_id+'_BS'
        try:
            del incar['KPOINT_BSE']
        except:
            pass
        incar.write_file('INCAR')


@explicit_serialize
class Generate_VaspInputFiles_for_HSE_scf_icsd(FireTaskBase, FWSerializable):
    def run_task(self, fw_spec):
        job_info_array = fw_spec['_job_info']
        prev_job_info = job_info_array[-1]
        xml = Vasprun(prev_job_info['launch_dir']+'/vasprun.xml')
        struct = xml.final_structure
	ks=fw_spec['ks_add']
        icsd_id=fw_spec['icsd_id']
	is_spin_polarized=fw_spec['is_spin_polarized']
	incar_add={}
	if is_spin_polarized:
            incar_add["ISPIN"] = 2
        if not is_spin_polarized:
            incar_add["ISPIN"] = 1
	incar_add["EDIFF"] = 0.00001
        incar_add["SYSTEM"] = icsd_id+"_HSE_StaticRun"
        incar_add["IBRION"] = -1
        incar_add["HFSCREEN"] = 0.2
        incar_add["PREC"] = "Accurate"
        incar_add["PRECFOCK"] = "Fast"
        incar_add["KPAR"] = 2
        incar_add["ISYM"] = 2
	dict_params = MPBSHSEVaspInputSet(user_incar_settings=incar_add, added_kpoints=None, mode="Line",
                              kpoints_density=None, kpoints_line_density=30).get_all_vasp_input(struct)
        incar=dict_params["INCAR"]
        kpoints = dict_params["KPOINTS"]
        for i, label in enumerate(kpoints.labels):
            if label != None:
                mark = i
                break
        kpoints.kpts = kpoints.kpts[0:mark]
        kpoints.kpts_weights = kpoints.kpts_weights[0:mark]
        kpoints.labels = kpoints.labels[0:mark]
        for k in ks:
            kpoints.kpts.append(np.array(k))
            kpoints.kpts_weights.append(0)
            kpoints.labels.append('')
        kpoints.num_kpts = len(kpoints.kpts_weights)
        incar.write_file("INCAR")
        dict_params["POTCAR"].write_file("POTCAR")
        dict_params["POSCAR"].write_file("POSCAR")
	kpoints.write_file("KPOINTS")


@explicit_serialize
class Generate_VaspInputFiles_for_HSE_certain_shots(FireTaskBase, FWSerializable):
    def run_task(self, fw_spec):
        job_info_array = fw_spec['_job_info']
        prev_job_info = job_info_array[-1]
        material_id=fw_spec['material_id']
        xml = Vasprun(prev_job_info['launch_dir']+'/vasprun.xml')
        structure = xml.final_structure
        bs=xml.get_band_structure()
        bg= bs.get_band_gap()
        shots=[]
        if bg["direct"]:
            cbm = bs.get_cbm()
            shots.append([cbm['kpoint']._fcoords, cbm['kpoint']._label])
        else:
            cbm = bs.get_cbm()
            vbm = bs.get_vbm()
        shots.append([cbm['kpoint']._fcoords, cbm['kpoint']._label])
        shots.append([vbm['kpoint']._fcoords, vbm['kpoint']._label])
        incar_add = {}
        is_spin_polarized=fw_spec['is_spin_polarized']
        if is_spin_polarized:
            incar_add["ISPIN"] = 2
        if not is_spin_polarized:
            incar_add["ISPIN"] = 1
        #incar_add["ISTART"] = 1
        incar_add["SYSTEM"] = material_id+"_HSE_StaticRun"
        #incar_add["AEXX"] = 0.25
        incar_add["IBRION"] = -1
        incar_add["HFSCREEN"] = 0.2
        incar_add["PREC"] = "Accurate"
        incar_add["PRECFOCK"] = "Fast"
        incar_add["KPAR"] = 2
        incar_add["ISYM"] = 2
        dict_params = MPBSHSEVaspInputSet(user_incar_settings=incar_add, added_kpoints=None, mode="Line",
                    kpoints_density=None, kpoints_line_density=30).get_all_vasp_input(structure)
        incar=dict_params["INCAR"]
        if is_spin_polarized:
            incar["ISPIN"] = 2
        elif not is_spin_polarized:
            try:
                del incar["ISPIN"]
                del incar["MAGMOM"]
            except:
                pass
        try:
            del incar["ICHARG"]
            del incar["ISYM"]
        except:
            pass
        kpoints = dict_params["KPOINTS"]
        label_kpts = enumerate(kpoints.labels)
        for i, label in label_kpts:
            if label != None:
                mark = i
                break
        kpoints.kpts = kpoints.kpts[0:mark]
        kpoints.kpts_weights = kpoints.kpts_weights[0:mark]
        kpoints.labels = kpoints.labels[0:mark]
        for pts in shots:
            kpoints.kpts.append(np.array(pts[0]))
            kpoints.kpts_weights.append(0)
            if pts[1] == None: pts[1] = ''
            kpoints.labels.append(pts[1])
        kpoints.num_kpts = len(kpoints.kpts_weights)
        kpoints.write_file("KPOINTS")
        incar.write_file("INCAR")
        dict_params["POTCAR"].write_file("POTCAR")
        dict_params["POSCAR"].write_file("POSCAR")
        
        

@explicit_serialize
class Insert_Gap_into_monogdb(FireTaskBase, FWSerializable):
    def run_task(self, fw_spec):
        from pymongo import MongoClient
        clt=MongoClient('marilyn.pcpm.ucl.ac.be',27017)
        db=clt.results_GY
        db.authenticate('gyu','pulco')
        hse_gap=db.HSE_gaps
        
        job_info_array = fw_spec['_job_info']
        prev_job_info = job_info_array[-1]
                
        result={}
        
        material_id=fw_spec['material_id']
        xml=Vasprun(prev_job_info['launch_dir']+'/vasprun.xml')
        is_converged=xml.converged
        bs=xml.get_band_structure()
        formula=xml.final_structure.composition.reduced_formula
        material={'mateiral_id':material_id,'formula':formula,'is_spin_polarized':bs.is_spin_polarized}
        
        input={}
        incar=xml.incar
        potcar=xml.potcar_spec
        lattice_parameters={}
        lattice_parameters['a'],lattice_parameters['b'],lattice_parameters['c']=xml.lattice.abc
        lattice_parameters['alpha'],lattice_parameters['beta'],lattice_parameters['gamma']=xml.lattice.angles
        lattice_vectors={}
        lattice_vectors['a_vec'],lattice_vectors['b_vec'],lattice_vectors['c_vec']=(list(xml.lattice.matrix[i]) for i in range(3))
        input['incar']=incar
        input['potcar']=potcar
        input['crystal']={'lattice_parameters':lattice_parameters,'lattice_vectors':lattice_vectors}
        

        gap_info={}
        gap_info['gap']=bs.get_band_gap()['energy']
        
        cbm={}
        cbm['energy']=bs.get_cbm()['energy']
        cbm_kpoint={}
        cbm_kpoint['coordinate']=list(bs.get_cbm()['kpoint'].frac_coords)
        cbm_kpoint['label']=bs.get_cbm()['kpoint'].label
        cbm['kpoint']=cbm_kpoint
        
        vbm={}
        vbm['energy']=bs.get_vbm()['energy']
        vbm_kpoint={}
        vbm_kpoint['coordiante']=list(bs.get_vbm()['kpoint'].frac_coords)
        vbm_kpoint['label']=bs.get_vbm()['kpoint'].label
        vbm['kpoint']=vbm_kpoint

        gap_info['cbm']=cbm
        gap_info['vbm']=vbm
        
        result['material']=material
        result['input']=input
        result['gap_info']=gap_info
        result['is_converged']=is_converged
        
        hse_gap.insert(result)
        
@explicit_serialize
class Insert_Gap_into_monogdb_icsd(FireTaskBase, FWSerializable):
    def run_task(self, fw_spec):
        from pymongo import MongoClient
        clt=MongoClient('marilyn.pcpm.ucl.ac.be',27017)
        db=clt.results_GY
        db.authenticate('gyu','pulco')
        hse_gap=db.HSE_gaps_ICSD
        
	pbe_bs=fw_spec['pbe_bs']
        job_info_array = fw_spec['_job_info']
        prev_job_info = job_info_array[-1]
                
        result={}
        
        icsd_id=fw_spec['icsd_id']
        xml=Vasprun(prev_job_info['launch_dir']+'/vasprun.xml')
	ks=Kpoints.from_file(prev_job_info['launch_dir']+'/KPOINTS')
        is_converged=xml.converged
        bs=xml.get_band_structure()
        formula=xml.final_structure.composition.reduced_formula
        material={'icsd_id':icsd_id,'formula':formula,'is_spin_polarized':bs.is_spin_polarized}
        
        input={}
        incar=xml.incar.as_dict()
        potcar=xml.potcar_spec
	crystal=xml.final_structure.as_dict()
	kpoints=ks.as_dict()
        input['incar']=incar
        input['potcar']=potcar
        input['crystal']=crystal
	input['kpoints']=kpoints
        

        gap_info={}
######### HSE #####################
	gap_info['HSE']={}
        gap_info['HSE']['gap']=bs.get_band_gap()['energy']
        
        cbm={}
        cbm['energy']=bs.get_cbm()['energy']
        cbm_kpoint={}
        cbm_kpoint['coordinate']=list(bs.get_cbm()['kpoint']._fcoords)
        cbm_kpoint['label']=bs.get_cbm()['kpoint'].label
        cbm['kpoint']=cbm_kpoint
        
        vbm={}
        vbm['energy']=bs.get_vbm()['energy']
        vbm_kpoint={}
        vbm_kpoint['coordiante']=list(bs.get_vbm()['kpoint']._fcoords)
        vbm_kpoint['label']=bs.get_vbm()['kpoint'].label
        vbm['kpoint']=vbm_kpoint

        gap_info['HSE']['cbm']=cbm
        gap_info['HSE']['vbm']=vbm
########### PBE from ICSD############
	gap_info['PBE_icsd']={}
	gap_info['PBE_icsd']['gap']=pbe_bs.get_band_gap()['energy']

	pbe_cbm={}
	pbe_cbm['energy']=pbe_bs.get_cbm()['energy']
	pbe_cbm_kpoint={}
	pbe_cbm_kpoint['coordinate']=list(pbe_bs.get_cbm()['kpoint']._fcoords)
	pbe_cbm_kpoint['label']=pbe_bs.get_cbm()['kpoint'].label
	pbe_cbm['kpoint']=pbe_cbm_kpoint
	
        pbe_vbm={}
        pbe_vbm['energy']=pbe_bs.get_vbm()['energy']
        pbe_vbm_kpoint={}
        pbe_vbm_kpoint['coordiante']=list(pbe_bs.get_vbm()['kpoint']._fcoords)
        pbe_vbm_kpoint['label']=pbe_bs.get_vbm()['kpoint'].label
        pbe_vbm['kpoint']=pbe_vbm_kpoint

        gap_info['PBE_icsd']['cbm']=pbe_cbm
	gap_info['PBE_icsd']['vbm']=pbe_vbm
#############################################################################	
        result['material']=material
        result['input']=input
        result['gap_info']=gap_info
        result['is_converged']=is_converged
##################################################################################        
        hse_gap.insert(result)
        
        
@explicit_serialize
class test1(FireTaskBase, FWSerializable):
    def run_task(self, fw_spec):
        f=open('tmp','a+')
        for i in range(50):
            print >>f, i,'test1'
            time.sleep(0.3)
@explicit_serialize
class test2(FireTaskBase, FWSerializable):
    def run_task(self, fw_spec):
        f=open('tmp','a+')
        for i in range(50):
            print>>f, i,'test2'
            time.sleep(0.3)
