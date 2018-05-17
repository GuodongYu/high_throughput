from pymatgen.io.vasp.inputs import Poscar, Potcar, Incar, Kpoints
from pymatgen.io.vasp.outputs import Outcar
import subprocess
import getpass
import os
from monty.os import cd
import yaml
from fireworks.fw_config import FWORKER_LOC, QUEUEADAPTER_LOC, LAUNCHPAD_LOC
from fireworks.core.launchpad import LaunchPad
import numpy as np
import json


class VaspDoneError(Exception):
    def __init__(self):
        pass


def get_directories_VaspJobNotDone(root_dir):
    with cd(root_dir): ### avoid the link problems
	root_dir_real = os.getcwd()
    scan=subprocess.Popen(['find',root_dir_real,'-name','POSCAR'],stdout=subprocess.PIPE)
    scan.wait()
    pos_coll=scan.stdout.read().split()
    pos_dirs=[os.path.split(i)[0] for i in pos_coll]
    vaspjob_dirs = []
    for dir in pos_dirs:
	try:
	    pos = Poscar.from_file(os.path.join(dir,'POSCAR'))
	    pot = Potcar.from_file(os.path.join(dir,'POTCAR'))
	    incar = Incar.from_file(os.path.join(dir,'INCAR'))
	    kpt = Kpoints.from_file(os.path.join(dir,'KPOINTS'))
	except:
	    print 'input files are not ready in %s' % dir
	else:
	    try:
		out=Outcar(os.path.join(dir,'OUTCAR'))
		if len(out.run_stats) != 7:
		    vaspjob_dir.append(dir)
	    except:
	        vaspjob_dirs.append(dir)
    return vaspjob_dirs


class FWConfig():
    def __init__(self,**kws):
	self._get_fworker_name()
	self._get_queue_type()
        if self.queue_type == 'SLURM':
            self.mpi_commond = 'srun'
        elif self.queue_type == 'PBS':
            self.mpi_commond = 'mpirun'
	self._get_launchpad()
        self.part = 'default'
        if self.fworker == 'manneback':
            set_dir = os.path.split(QUEUEADAPTER_LOC)[0]
            self.part= kws['part']
            queue_yaml = yaml.load(open(QUEUEADAPTER_LOC))
            if kws and kws['part']=='Zoe':
                queue_yaml['_fw_template_file'] = os.path.join(set_dir,'SLURM_UCL_template_Zoe.txt')
            else:
                queue_yaml['_fw_template_file'] = os.path.join(set_dir,'SLURM_UCL_template.txt')
            dumpfn(queue_yaml, os.path.join(set_dir,'my_qadapter.yaml'), indent=2)

    def _get_fworker_name(self):
        fworker = yaml.load(open(FWORKER_LOC))
        self.fworker = fworker['name']
    
    def _get_queue_type(self):
        queue = yaml.load(open(QUEUEADAPTER_LOC))
        self.queue_type = queue['_fw_q_type']

    def _get_launchpad(self):
	launch = yaml.load(open(LAUNCHPAD_LOC))
	self.launch_pad = LaunchPad(host=launch['host'], port=launch['port'],\
            name=launch['name'], username=launch['username'], password=launch['password'])
    
    def get_directories_in_queue(self):
        user_name=getpass.getuser()
        if self.queue_type == 'PBS':
    	    get_job_ids=['qstat','-u',user_name]
    	    check_job_info = ['qstat','-f']
    	    workdir_str='PBS_O_WORKDIR'
    	    split_str=','
        elif self.queue_type == 'SLURM':
    	    get_job_ids=['squeue','-u',user_name]
    	    check_job_info = ['scontrol','show','job']
    	    workdir_str='WorkDir'
    	    split_str='\n'
        scan=subprocess.Popen(get_job_ids,stdout=subprocess.PIPE)
        scan.wait()
        scan_out=scan.stdout.read().split('\n')
        ids=[i.split()[0] for i in scan_out if user_name in i]
        dirs=[]
        for id in ids:
            scan1=subprocess.Popen(check_job_info+[id],stdout=subprocess.PIPE)
            out=scan1.stdout.read().split(split_str)
            line=[i for i in out if workdir_str in i][0]
            dir=line.replace('\n\t','').split('=')[1]
            with cd(dir):
        	    dir_orig=os.getcwd()
        	    dirs.append(dir_orig)
        return dirs

    def get_directories_NeedVaspRun(self,dir_root):
	"""
	Get the directories of not finished vasp jobs without pending or running in queue
	"""
	dirs_not_finished = get_directories_VaspJobNotDone(dir_root)
        dirs_in_queue =self.get_directories_in_queue()
	dirs_shot = [i for i in dirs_not_finished if i not in dirs_in_queue]
	return dirs_shot
 
    def get_nelec_LHFCALC(self, dir):
        pos = Poscar.from_file(os.path.join(dir,'POSCAR'))
        pot = Potcar.from_file(os.path.join(dir,'POTCAR'))
        incar = Incar.from_file(os.path.join(dir,'INCAR'))
	kpt = Kpoints.from_file(os.path.join(dir,'KPOINTS'))
        if 'LHFCALC' in incar:
    	    LHFCALC = incar['LHFCALC']
        else:
    	    LHFCALC = False
        natoms = pos.natoms
        nelec_atom = [i.nelectrons for i in pot]
        nelec_elt = np.array(natoms)*np.array(nelec_atom)
        nelec_sum = nelec_elt.sum()
        return (nelec_sum , LHFCALC)

    def queue_setup(self,dir):
	nelec, LHFCALC = self.get_nelec_LHFCALC(dir)
	queue = {}
        if self.part == 'Zoe':
            queue['walltime'] = '5-00:00:00'
            queue['queue'] = 'Zoe'

	########### ncores #########
	if not LHFCALC:
	        nelec_1core = 40.
	else:
	        nelec_1core = 10.
		 
	ncore = nelec/nelec_1core
	ncore = int(round(ncore/2))*2 
	queue['ntasks'] = ncore
	############# walltime for zenobe ########
	if self.fworker == 'zenobe':
	    if not LHFCALC:
	    	if ncore < 12:
	    	    time = 24
	    	elif 12 <= ncore < 36:
	    	    time = 60
	    	elif 36 <= ncore < 60:
	    	    time = 90
	    	else:
	       	    time = 150
	    else:
	    	if ncore < 12:
	    	    time = 100
	    	elif 12 <= ncore < 36:
	    	    time = 150
	    	elif 36 <= ncore < 60:
	    	    time = 300
	    	else:
	       	    time = 500	
	    queue['walltime'] = '%s:00:00' % str(time)
	########## job name for defect ###########
	defect = False
	for transf_json in ['transformations.json','transformation.json']:
	    transf_file = os.path.join(dir,transf_json)
	    if os.path.isfile(transf_file):
	    	transf = json.load(open(transf_file))
		defect = True
		break
	if defect:
	    mpid = transf['history'][0]['source']
	    if 'charge' in transf:
		queue['job_name'] = mpid+'_'+transf['defect_type']+'_'+str(transf['charge'])
	    else:
		queue['job_name'] = mpid+'_'+transf['defect_type']
	else:
	    queue['job_name'] = 'No_name'
	return queue




