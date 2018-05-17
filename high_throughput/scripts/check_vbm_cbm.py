from high_throughput.defects.database import HseGapOperater
from high_throughput.vasp.utils import are_equivalent_kpoints
from high_throughput.config import *
from pymatgen.core.structure import Structure
import sys
hgo=HseGapOperater()

f=open('/home/gyu/workplace/My_Work/defects_thermoelectrons_HT/k_vbm_cbm_comparison.dat','w')
f.write('\t'.join(['mpid','formula','kpbevbm','khsevbm','equal?','kpbecbm','khsecbm','equal?'])+'\n')
all_data = hgo.collection.find()
for i in all_data:
    try:
        mpid = i['material']['mateiral_id']
        stru = Structure.from_dict(i['input']['crystal'])
        print mpid
        formula = i['material']['formula']
        gap_info =  i['gap_info']
        kcbm = gap_info['cbm']['kpoint']['coordinate']
        kvbm = gap_info['vbm']['kpoint']['coordiante']
        try:
            bs = m.get_bandstructure_by_material_id(mpid)
            kvbm_pbe = bs.get_vbm()['kpoint'].frac_coords
            kcbm_pbe = bs.get_cbm()['kpoint'].frac_coords
            vbm_equal='F';cbm_equal='F'
            if are_equivalent_kpoints(stru,kvbm,kvbm_pbe):
                vbm_equal = 'T'
            if are_equivalent_kpoints(stru,kcbm,kcbm_pbe):
                cbm_equal = 'T'
            record =[mpid,formula,str(kvbm_pbe),str(kvbm),vbm_equal,str(kcbm_pbe),str(kcbm),cbm_equal]
            f.write('\t'.join(record)+'\n')
        except KeyError():
            sys.exit()
    except KeyError():
        sys.exit()
        

