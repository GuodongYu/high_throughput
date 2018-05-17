from pymatgen.core.periodic_table import _pt_data
from pymatgen.core.periodic_table import Element
from high_throughput.config  import *

def common_binary_semiconductors(name):
    """
    name: 2-6, 3-5, or 4-4, the semiconductor name in general
    """
    def check_11(str_):
        ret = ''
        for i in str_:
            if i.isdigit():
                ret+=i
        if not ret:
            return True
        else:
            return False

    gs = [int(i) for i in name.split('-')]
    g0 = gs[0] if gs[0] <= 2 else gs[0]+10
    g1 = gs[1] if gs[1] <= 2 else gs[1]+10
    n0=[]
    n1=[]
    for ele in _pt_data:
        Elm = Element(ele)
        if Elm.group == g0 and not Elm.is_rare_earth_metal:
            n0.append(ele)
        elif Elm.group == g1 and not Elm.is_rare_earth_metal:
            n1.append(ele)
    data = []
    for e0 in n0:
        for e1 in n1:
            datai = m.query(criteria={"elements": {"$all": [e0, e1]},'nelements':2}, properties=["band_gap", "task_id", "pretty_formula"])
            datai = [i for i in datai if i['band_gap']>0.1]
            data+=datai
    return [i for i in data if check_11(i["pretty_formula"])]



