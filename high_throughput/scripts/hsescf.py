from high_throughput.vasp.input_files_generators import Generator
from pymatgen.io.vasp.outputs import Vasprun
import os
import glob
import json
pbesol_root='/home/gyu/hse_gap/pbesol_relax'
pbesol_dirs=glob.glob(os.path.join(pbesol_root,'*'))
hsescf_root = '/home/gyu/hse_gap/hsescf'

for path in pbesol_dirs:
    try:
        d = json.load(open(path+'/profile.json'))
        mpid = d['source']
        xml=Vasprun(path+'/vasprun.xml')
        relaxed_stru = xml.final_structure
        g = Generator(mpid)
        g.HSE_gap_with_shot(relaxed_stru,output_dir=hsescf_root)
    except:
        pass
