from pymatgen.symmetry.groups import SymmOp
from pymatgen import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.groups import PointGroup
from monty.json import jsanitize
m=MPRester()
def sg_pg_compare(mpid='mp-989535'):
    stru=m.get_structure_by_material_id(mpid)
    spa=SpacegroupAnalyzer(stru)
    pg_name=spa.get_symmetry_dataset()['pointgroup']
    pg_ops = PointGroup(pg_name).symmetry_ops
    sp_ops = []
    sp_mats = []
    pg_mats =[]
    for op in pg_ops:
        rotation = op.rotation_matrix
        pg_mats.append(jsanitize(rotation))
    for op in spa.get_symmetry_operations():
        rotation = op.rotation_matrix
        sp_mats.append(jsanitize(rotation))
        sp_ops.append(SymmOp.from_rotation_and_translation(rotation,(0,0,0)))
    sp_ops=set(sp_ops)
    #pg_mats.sort();sp_mats.sort()
    return pg_mats,sp_mats

