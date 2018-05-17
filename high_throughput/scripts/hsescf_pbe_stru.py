from pymatgen.io.vasp.outputs import Vasprun
from high_throughput.defects.utils.batch_tools import get_mpids
from high_throughput.vasp.utils import  InputsGenerators
import os
import glob
import json

mpids = get_mpids('/home/gyu/scripts/mpid.dat')
output_dir = '/home/gyu/hse_gap'

for mpid in mpids:
    print mpid
    g = InputsGenerators(mpid)
    g.HSE_gap_with_shot(g.stru,output_dir)
