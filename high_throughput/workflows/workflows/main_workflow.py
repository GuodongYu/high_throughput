from fireworks.core.firework import Firework, Workflow
from fireworks.core.launchpad import LaunchPad
from defect_tools.intrinsic_defects_maker import IntrinsicDefectsMaker
from pymatgen import MPRester
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from defectswork.firetasks.firetasks import VaspRun, CommandRun
from defectswork.firetasks.firetasks import MakeVaspDefectFile
from defectswork.firetasks.firetasks import Generate_VaspInputFiles_for_Relax_with_PBEsol, \
					    Generate_VaspInputFiles_for_Relax_with_PBEsol_icsd,\
                                            Generate_VaspInputfiles_for_Static_Run,Generate_VaspInputfiles_for_BS
from defectswork.firetasks.firetasks import Generate_VaspInputFiles_for_HSE_certain_shots,\
                                            Generate_VaspInputFiles_for_HSE_scf_icsd
from defectswork.firetasks.firetasks import Insert_Gap_into_monogdb,Insert_Gap_into_monogdb_icsd
from high_throuhgput.config import * 
from pymatgen.io.vasp.sets import Vasprun,Outcar
from pymatgen.io.vasp.outputs import Incar, Potcar
import subprocess
import json
import os
import shutil
from monty.os import cd
import re
import glob
from utils import FWConfig

init = FWConfig()
fworker = init.fworker
launchpad = init.launch_pad
class VaspDoneError(Exception):
    def __init__(self):
        pass

def vasp_jobs_scan_and_run(dir, vasp_cmd, label):
    """
    Args:
	dir: directory need to scan
        vasp_cmd: vasp run command executed by subprocess.Popen, e.g. ['mpirun','vasp_std'] or ['srun','vasp_std']
        label: a label for these jobs
    """
    work_dirs = init.get_directories_NeedVaspRun(dir)
    fws=[]
    njobs=0
    for work_dir in work_dirs:
        queue = init.queue_setup(work_dir)
	fw_name = queue['job_name']
        ftask=VaspRun()
        fw=Firework([ftask],spec={'vasp_cmd':vasp_cmd,'_launch_dir':work_dir,'_queueadapter':queue,'_fworker':fworker},name=fw_name)
        fws.append(fw)
	njobs=njobs+1
    wf=Workflow(fws,name=label)
    launchpad.add_wf(wf)
    return njobs

def wavecar_get_for_localization_check(dir, vasp_cmd, label):
    """
    Args:
        dir: directory need to scan
        vasp_cmd: vasp run command executed by subprocess.Popen, e.g. ['mpirun','vasp_std'] or ['srun','vasp_std']
        label: a label for these jobs
    """
    scan=subprocess.Popen(['find',dir,'-name','POSCAR'],stdout=subprocess.PIPE)
    scan.wait()
    dir_coll=[os.path.split(dir)[0] for dir in scan.stdout.read().split()]
    exclude=re.compile('.*charge0|.*charge_0|.*skip')
    include=re.compile('.*subst|.*sub|.*vac|.*inter|.*as')
    work_dirs=[dir for dir in dir_coll if not exclude.match(dir) and 'charge' in dir.split('/')[-1] and include.match(dir)]
    #work_dirs=[dir for dir in dir_coll 
    #              if ('charge0' not in dir or 'charge_0' not in dir) and 'charge' in dir.split('/')[-1] and ('subst' in dir or 'vac' in dir)]
    fws=[]
    njobs=0
    dir_end='wavecar'
    for work_dir in work_dirs:
	with cd(work_dir):
	    os.mkdir(dir_end)
	    incar=Incar.from_file('INCAR')
	    incar['LVHAR']=True;incar['LWAVE']=True;incar['IBRION']=-1;incar['NSW']=0
	    incar['EDIFF']=0.000001
	    incar.write_file(os.path.join(dir_end,'INCAR'))
	    shutil.copyfile('CONTCAR', os.path.join(dir_end,'POSCAR'))
	    inputs_other=['KPOINTS','POTCAR']
     	    for input in inputs_other:
		shutil.copy(input, dir_end)
	    try:
		shutil.copy('transformations.json',dir_end)
		transf=json.load(open('transformations.json'))
	    except:
		shutil.copy('transformation.json',dir_end)
                transf=json.load(open('transformation.json'))
	    queue={}	    
            queue['job_name']=label+'_'+transf['defect_type']+'_'+str(transf['charge'])+'_wavecar'
            fw_name=transf['defect_type']+'_'+str(transf['charge'])
            ftask=VaspRun()
            fw=Firework([ftask],spec={'vasp_cmd':vasp_cmd,'_launch_dir':os.path.join(work_dir,dir_end),'_queueadapter':queue,'_fworker':fworker},name=fw_name)
            fws.append(fw)
            njobs=njobs+1
    wf=Workflow(fws,name=label)
    launchpad.add_wf(wf)
    return njobs

def neutral_locpot_get_for_potalign(dir, vasp_cmd, label):
    """
    Args:
        dir: directory need to scan
        vasp_cmd: vasp run command executed by subprocess.Popen, e.g. ['mpirun','vasp_std'] or ['srun','vasp_std']
        label: a label for these jobs
    """
    scan=subprocess.Popen(['find',dir,'-name','POSCAR'],stdout=subprocess.PIPE)
    scan.wait()
    dir_coll=[os.path.split(dir)[0] for dir in scan.stdout.read().split()]
    exclude=re.compile('.*charge0|.*charge_0|.*skip')
    include=re.compile('.*subst|.*sub|.*vac|.*inter|.*as')
    work_dirs=[dir for dir in dir_coll if not exclude.match(dir) and 'charge' in dir.split('/')[-1] and include.match(dir)]
    #work_dirs=[dir for dir in dir_coll if 'charge0' not in dir and 'charge' in dir.split('/')[-1] and ('_subst' in dir or '_vac' in dir) ]
    fws=[]
    njobs=0
    dir_end='neutral'
    for work_dir in work_dirs:
	with cd(work_dir):
	    os.mkdir(dir_end)
	    print os.path.join(work_dir,dir_end)
	    incar=Incar.from_file('INCAR')
	    incar['LVHAR']=True;incar['LWAVE']=False;incar['IBRION']=-1;incar['NSW']=0
	    incar['EDIFF']=0.000001;del incar['NELECT']
	    incar.write_file(os.path.join(dir_end,'INCAR'))
	    shutil.copyfile('CONTCAR', os.path.join(dir_end,'POSCAR'))
	    inputs_other=['KPOINTS','POTCAR']
     	    for input in inputs_other:
		shutil.copy(input, dir_end)
	    try:
		shutil.copy('transformations.json',dir_end)
		transf=json.load(open('transformations.json'))
	    except:
		shutil.copy('transformation.json',dir_end)
                transf=json.load(open('transformation.json'))
	    queue={}	    
            queue['job_name']=label+'_'+transf['defect_type']+'_'+str(transf['charge'])+'_neutral'
            fw_name=transf['defect_type']+'_'+str(transf['charge'])
            ftask=VaspRun()
            fw=Firework([ftask],spec={'vasp_cmd':vasp_cmd,'_launch_dir':os.path.join(work_dir,dir_end),'_queueadapter':queue,'_fworker':fworker},name=fw_name)
            fws.append(fw)
            njobs=njobs+1
    wf=Workflow(fws,name=label)
    launchpad.add_wf(wf)
    return njobs
#### not a workflow
def get_iband(dir):
    stru=Structure.from_file(os.path.join(dir,'POSCAR'))
    pots=Potcar.from_file(os.path.join(dir,'POTCAR'))
    elm_nele={}
    for pot in pots:
        elm_nele[pot.element]=pot.nelectrons
    elm_natom=stru.composition.as_dict()    
    
    nele_neu=0.0
    for elm in elm_natom.keys():
        nele_neu = nele_neu + elm_natom[elm]*elm_nele[elm]	
    out=Outcar(os.path.join(dir,'OUTCAR'))
    nelec=int(round(out.nelect))
    charge=nele_neu-nelec
    
    eig_file=open(os.path.join(dir,'EIGENVAL'))
    eigs=eig_file.read().split('\n')
    nk=int(eigs[5].split()[1])
    occu={}
    for i in range(len(eigs))[8:]:
        if len(eigs[i].split())==5:
    	    nth_bnd=int(eigs[i].split()[0])
    	    noccu=float(eigs[i].split()[-1])+float(eigs[i].split()[-2])
    	    if int(eigs[i].split()[0]) in occu:
    	        occu[nth_bnd]=occu[nth_bnd]+noccu
    	    else:
    	        occu[nth_bnd]=noccu
    for key in occu:
        occu[key]=round(occu[key]/nk,2)
    iband=[]
    nbnd=len(occu)
    #print occu
    n_add=0.0
    for ind in range(1,nbnd+1):
        if charge < 0:
    	    ind_band=nbnd+1-ind  
    	    left=abs(charge)-n_add
    	    if abs(left)<0.01:
    	        break
    	    nadd_able=occu[ind_band]
            if left >= nadd_able: 
                n_add=n_add+nadd_able
            elif left < nadd_able:
                n_add=n_add+left
    	    if occu[ind_band]>0.:
    	        iband.append(ind_band)
        elif charge > 0:
    	    ind_band=ind
    	    left= abs(charge)-n_add
    	    if abs(left)<0.01:
    	        break
    	    nadd_able=2.0-occu[ind_band]
    	    if left >= nadd_able:
                n_add=n_add+nadd_able
    	    elif left < nadd_able:
    	        n_add=n_add+left
    	    if nadd_able>0.0:
    	        iband.append(ind_band)
    iband.sort()
    return iband
def parchg_get_for_localization_check(dir, vasp_cmd, label):
    """
    Args:
        dir: directory need to scan
        vasp_cmd: vasp run command executed by os.system, e.g. 'vasp_std'
        label: a label for these jobs
    """
    scan=subprocess.Popen(['find',dir,'-name','wavecar'],stdout=subprocess.PIPE)
    scan.wait()
    dir_coll=scan.stdout.read().split()
    exclude=re.compile('.*charge0|.*charge_0|.*skip')
    include=re.compile('.*subst|.*sub|.*vac|.*inter|.*as')
    work_dirs=[dir for dir in dir_coll if not exclude.match(dir) and include.match(dir)]
    fws=[]
    njobs=0
    dir_end='parchg'
    for work_dir in work_dirs:
	with cd(work_dir):
	    print work_dir
	    if os.path.isdir(dir_end):
		if glob.glob(os.path.join('parchg','PARCHG*')) !=[]:
		    print work_dir, ' was already well done! Skip now!'
		    continue 
		else:
		    delect=subprocess.Popen(['rm','-rf',dir_end])
		    delect.wait()
	    os.mkdir(dir_end)
	    iband=get_iband(work_dir)
	    incar=Incar.from_file('INCAR')
	    incar['ICHARG']=0;incar['ISTART']=1;incar['LPARD']=True
	    incar['IBAND']=iband; incar['KPUSE']=1;incar['LSEPB']=True;incar['LSEPK']=True
	    incar.write_file(os.path.join(dir_end, 'INCAR'))
	    shutil.copyfile('IBZKPT', os.path.join(dir_end,'KPOINTS'))
	    inputs_other=['POTCAR','POSCAR','WAVECAR']
     	    for input in inputs_other:
		shutil.copy(input, dir_end)
	    try:
		transf=json.load(open('transformations.json'))
	    except:
                transf=json.load(open('transformation.json'))
	    queue={}	    
            queue['job_name']=label+'_'+transf['defect_type']+'_'+str(transf['charge'])
	    queue['ntasks']=1;queue['vmem']='40000mb';queue['walltime']='01:00:00'
            fw_name=transf['defect_type']+'_'+str(transf['charge'])
            ftask=CommandRun()
            fw=Firework([ftask],spec={'cmd':vasp_cmd,'_launch_dir':os.path.join(work_dir,dir_end),'_queueadapter':queue,'_fworker':fworker},name=fw_name)
            fws.append(fw)
            njobs=njobs+1
    wf=Workflow(fws,name=label)
    launchpad.add_wf(wf)
    return njobs

def HSE_Gap_with_Vasp_from_icsd_data(file,vasp_cmd):
    fws=[]
    all=json.load(open(file))
    njobs=0
    g=open('/SCRATCH/acad/htbase/yugd/fw_workplace/hse_gaps_icsd/ids')
    ids=g.read().split()
    f=open('wrong_list.dat','a+')
    for task in all:
	icsd_id='icsd-'+str(task['icsd_ids'][0])
	if icsd_id not in ids:
	    continue
	try:
	    struct=Structure.from_dict(task['structure'])
	    pbe_bs=BandStructureSymmLine.from_dict(task['band_structures']['line'])
	    is_spin_polarized=pbe_bs.is_spin_polarized
	    k_vbm=pbe_bs.get_vbm()['kpoint']._fcoords
	    k_cbm=pbe_bs.get_cbm()['kpoint']._fcoords
	except:
	    f.write('all['+str(all.index(task))+']')
	    continue 
##########################################################
	njobs=njobs+1
	queue_hse={}
	queue_pbesol={}
        queue_hse['job_name']=icsd_id+'_HSE'
	queue_pbesol['job_name']=icsd_id+'_PBEsol'
	nsites=struct.num_sites
        if nsites<=5:
	    queue_pbesol['ntasks']=6
	    queue_pbesol['walltime']='6:00:00'
            queue_hse['ntasks']=12
            queue_hse['walltime']='24:00:00'
        elif 5<nsites<=10:
	    queue_pbesol['ntasks']=12
            queue_pbesol['walltime']='6:00:00'
            queue_hse['ntasks']=24
            queue_hse['walltime']='48:00:00'
        elif 10<nsites<=20:
	    queue_pbesol['ntasks']=12
            queue_pbesol['walltime']='24:00:00'
            queue_hse['ntasks']=48
            queue_hse['walltime']='48:00:00'
        elif 20<nsites<=50:
	    queue_pbesol['ntasks']=24
            queue_pbesol['walltime']='48:00:00'
            queue_hse['ntasks']=60
            queue_hse['walltime']='72:00:00'
        else:
	    queue_pbesol['ntasks']=36
            queue_pbesol['walltime']='48:00:00'
            queue_hse['ntasks']=96
            queue_hse['walltime']='72:00:00'

        queue_insert={}
        queue_insert['ntasks']=1
        queue_insert['walltime']='01:00:00'
        queue_insert['vmem'] = '1024mb'
        queue_insert['job_name']=icsd_id+'_Insert'
##########################################################

	if pbe_bs.get_band_gap()['direct']:
	    ks_add=[k_vbm]
	else:
	    ks_add=[k_vbm,k_cbm]
        ftask00 = Generate_VaspInputFiles_for_Relax_with_PBEsol_icsd()
        ftask01 = VaspRun()
	fw0= Firework([ftask00,ftask01],spec={'_pass_job_info':True,'_preserve_fworker':True,'_fworker':fworker,
                                             'vasp_cmd':vasp_cmd,'workdir':'./','structure':struct,'icsd_id':icsd_id,
                                             'is_spin_polarized':is_spin_polarized,'_queueadapter':queue_pbesol},name=icsd_id+'_pbesol')

        ftask10=Generate_VaspInputFiles_for_HSE_scf_icsd()
	ftask11=VaspRun()
	fw1=Firework([ftask10,ftask11],parents=[fw0],spec={'_pass_job_info':True,'_preserve_fworker':True,'_fworker':fworker,
                                             'vasp_cmd':vasp_cmd,'workdir':'./','icsd_id':icsd_id,'ks_add':ks_add,
                                             'is_spin_polarized':is_spin_polarized,'_queueadapter':queue_hse},name=icsd_id+'_HSE')
	finsert=Insert_Gap_into_monogdb_icsd()
	fw2=Firework([finsert],parents=[fw1],spec={'icsd_id':icsd_id,'_queueadapter':queue_insert,'pbe_bs':pbe_bs},name=icsd_id+'_Insert')
	fws.append(fw0);fws.append(fw1);fws.append(fw2)
    wf=Workflow(fws,name='HSE_Gap_for_icsds')
    launchpad.add_wf(wf)
    f.close()
    return njobs
    


def PBEsol_Relax_then_HSE_Gap_with_Vasp(material_id, vasp_cmd, struct, is_spin_polarized):
    """
    add the defect relax calculations(by Vasp) into the launchpad.
    return the total calculation number
    Args:
       material_id: the material under consideration with id from materials project
       vasp_cmd1: the vasp command for PBEsol relaxation, e.g. ['mpirun','vasp']
       vasp_cmd2: the vasp command for HSE static run, e.g. ['mpirun', 'vasp']
       hse: True means considering the hse during the vasp calculations. Default: False
    """

    ftaska1 = Generate_VaspInputFiles_for_Relax_with_PBEsol()
    ftaska2 = VaspRun()

    ftaskb1 = Generate_VaspInputfiles_for_Static_Run()
    ftaskb2 = VaspRun()

    ftaskc1 = Generate_VaspInputfiles_for_BS()
    ftaskc2 = VaspRun()

    ftaskd1 = Generate_VaspInputFiles_for_HSE_certain_shots()
    ftaskd2 = VaspRun()

    ftaskf1 = Insert_Gap_into_monogdb()
    
    queue1={} ### for PBE
    queue2={} ### for HSE
    queue1['job_name']=material_id+'_PBEsol'
    queue2['job_name']=material_id+'_HSE'
    nsites=struct.num_sites
    if nsites<=5:
        queue1['ntasks']=6
        queue1['walltime']='6:00:00'
        queue2['ntasks']=24
        queue2['walltime']='24:00:00'       
    elif 5<nsites<=10:
        queue1['ntasks']=12
        queue1['walltime']='12:00:00'
        queue2['ntasks']=48
        queue2['walltime']='48:00:00'
    elif 10<nsites<=20:
        queue1['ntasks']=24
        queue1['walltime']='24:00:00'
        queue2['ntasks']=60
        queue2['walltime']='48:00:00'
    elif 20<nsites<=50:
        queue1['ntasks']=36
        queue1['walltime']='24:00:00'
        queue2['ntasks']=96
        queue2['walltime']='72:00:00'
    else:
        queue1['ntasks']=48
        queue1['walltime']='24:00:00'
        queue2['ntasks']=96
        queue2['walltime']='72:00:00'        
    
    ##insert info to Mongodb 
    queue3={}
    queue3['ntasks']=1
    queue3['walltime']='01:00:00'
    queue3['vmem'] = '1024mb'
    queue3['job_name']=material_id+'_Insert'  
    
    ## PBEsol relax
    fw1=Firework([ftaska1,ftaska2],spec={'_pass_job_info':True,'_preserve_fworker':True,'_fworker':fworker,
                                         'vasp_cmd':vasp_cmd,'workdir':'./','structure':struct,'material_id':material_id,'is_spin_polarized':is_spin_polarized,
                                         '_queueadapter':queue1},name=material_id+'_PBEsol_Relax')
    ## PBEsol static run
    fw2=Firework([ftaskb1,ftaskb2],parents=[fw1],spec={'_pass_job_info':True,
                                                       'vasp_cmd':vasp_cmd,'workdir':'./','material_id':material_id,
                                                        '_queueadapter':queue1},name=material_id+'_PBEsol_StaticRun')
    ## PBEsol band structure
    fw3=Firework([ftaskc1,ftaskc2],parents=[fw2],spec={'_pass_job_info':True,
                                                       'vasp_cmd':vasp_cmd,'workdir':'./','material_id':material_id,
                                                        '_queueadapter':queue1},name=material_id+'_PBEsol_BS')

    ## HSE static run with PBEsol k_VBM and k_CBM
    fw4=Firework([ftaskd1,ftaskd2],parents=[fw3],spec={'_pass_job_info':True,
                                                       'vasp_cmd':vasp_cmd,'workdir':'./','material_id':material_id,'is_spin_polarized':is_spin_polarized,
                                                        '_queueadapter':queue2},name=material_id+'_HSE_scf')

    ## insert results into mongodb
    fw5=Firework([ftaskf1],parents=[fw4],spec={'material_id':material_id,'_queueadapter':queue3},name=material_id+'_Insert')
    
    wf=Workflow([fw1,fw2,fw3,fw4,fw5],{fw1:[fw2],fw2:[fw3],fw3:[fw4],fw4:[fw5]},name=material_id+'_'+struct.composition.reduced_formula+'_HSE_Gap')
    launchpad.add_wf(wf)

