import numpy as np
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.maker import PhaseDiagram
from pymatgen.phasediagram.analyzer import PDAnalyzer
from high_throughput.config import *
from high_throughput.defects.constants import *
from pycdt.corrections.utils import closestsites
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.core.structure import Structure

def get_sc_size(bulk_uc,bulk_sc):
    sc_size=[]
    for vec in bulk_sc.lattice.matrix:
        time = list(np.round(bulk_uc.lattice.get_fractional_coords(vec)))
        sc_size.append(time)
    return sc_size

def are_lattices_equal(st1,st2):
    lat1=np.array(st1.lattice.matrix)
    lat2=np.array(st2.lattice.matrix)
    diff=lat1-lat2
    diff=np.round(diff,3)
    if not diff.any():
        return True
    else:
        return False

def are_structures_equal(st1, st2, tol=0.001):
    """
    Compare whether two structures st1 and st2 are equal.
    Here equal means they must have the same axes, cell shape,
    cell size and species.
    """
    if not are_lattices_equal(st1,st2):
        return False
    equal=True
    for i in range(len(st1)):
        targ=closestsites(st1,st2, st1.sites[i].coords)[1]
        spe1=st1.sites[i].specie
        dis=targ[1]
        spe2=targ[0].specie
        if spe1 != spe2 or dis > tol:
            equal = False
            break
    return equal

def check_rotation_two_structs(stru1,stru2):
    n1 = stru1.num_sites; n2 = stru2.num_sites
    ref = stru1 if n1<=n2 else stru2
    comp = stru2 if ref == stru1 else stru1
    ref_latt = ref.lattice
    rotation = False
    nn=0
    for out in comp.sites:
        out_specie =  out.specie
        out_cart = out.coords
        match = False
        for host in ref.sites:
            host_specie = host.specie
            host_cart = host.coords
            if host_specie != out_specie:
                continue
            diff_cart = out_cart - host_cart
            diff_frac = ref_latt.get_fractional_coords(diff_cart)
            vec_match = True
            for i in diff_frac:
                diff = abs(i-round(i,0))
                if diff > 0.01:
                    vec_match = False
            if vec_match:
                match = True
                break
        nn+=1
        if not match:
            return True
    return False

def strutures_group(strus, tol=0.001):
    sts=[i for i in strus]
    done=[sts[0]]
    left= [i for i in sts if i not in done]
    group=[[strus[0]]]
    for st in left:
        unique = True
        for i in range(len(group)):
            sg = group[i][0]
            if are_structures_equal(st, sg):
                group[i].append(st)
                unique = False
                break
        if unique == True:
            group.append([st])
        done.append(st)
        left = [i for i in sts if i not in done]
    return group

def are_blk_dfc_match(st_blk, st_dfc, tol=0.4):
    if not are_lattices_equal(st_blk,st_dfc):
        return False
    num_blk = st_blk.num_sites
    num_dfc = st_dfc.num_sites
    st_base = st_blk if num_blk <= num_dfc else st_dfc
    st_comp = st_blk if st_blk != st_base else st_dfc
    min_dis=[]
    for i in st_base:
        y=[]
        for j in st_comp:
            y.append(i.distance(j))
        min_dis.append(min(y))
    print np.array(min_dis).mean()
    if np.array(min_dis).mean()<tol:
        return True
    else:
        return False

def closest_defect_distance(stru):
    lat_mat=stru.lattice.matrix
    dis=[]
    for i in [-1,0,1]:
        for j in [-1,0,1]:
            for k in [-1,0,1]:
                if not (i or j or k):
                    continue
                vec=i*lat_mat[0] +j*lat_mat[1] + k*lat_mat[2]
                norm=np.linalg.norm(vec)
                dis.append(norm)
    dis_min=min(dis)
    return dis_min

def get_mu_range(mpid, ext_elts=[]):
    if 'hse' in mpid:
        mpid = mpid.split('_')[0]
    try:
        entry = m.get_entry_by_material_id(mpid)
	in_MP = True
    except:
	from high_throughput.defects.database import TasksOperater
	in_MP = False
	TO = TasksOperater()
	id_ = TO.groups['bulk'][mpid][0]
	rec = TO.collection.find_one({'_id':id_},['output'])
	stru_tmp =Structure.from_dict(rec['output']['crystal'])
	energy = rec['output']['final_energy']
	entry = PDEntry(stru_tmp.composition, energy)
    elts = [i.symbol for i in entry.composition.elements]
    for i in ext_elts:
        elts.append(i)
    entries = m.get_entries_in_chemsys(elts)
    if not in_MP:
	entries.append(entry)
    for entry in entries:
        entry.correction = 0.0
    pd=PhaseDiagram(entries)
    pda=PDAnalyzer(pd)
    chempots={}
    decompositions=pda.get_decomposition(entry.composition).keys()
    decomposition_inds=[pd.qhull_entries.index(entry) for entry in decompositions]
    facets_around=[]
    for facet in pd.facets:
        is_facet_around=True
        for ind in decomposition_inds:
            if ind not in facet:
                is_facet_around=False
        if is_facet_around==True:
            facets_around.append(facet)
    for facet in facets_around:
        s=[]
        for ind in facet:
            s.append(str(pd.qhull_entries[ind].name))
        s.sort()
        chempots['-'.join(s)]=pda.get_facet_chempots(facet)
    chempots={i:{j.symbol:chempots[i][j] for j in chempots[i]} for i in chempots}
    return chempots

def get_pinning_energies(form_en,vbm):
    """
    Args:
        da: defect_analyzer with defects added
    """
    summ = {}
    for i in form_en:
        if i['charge'] not in summ:
            summ[i['charge']]=[]
        summ[i['charge']].append(i['energy'])
    zero_pts = {}
    for i in summ:
        if i != 0:
            zero_pts[i]=vbm-min(summ[i])/i
    pins={}
    pins['positive'] = [-1000.]
    pins['negative'] = [1000.]
    for i in zero_pts:
        if i > 0:
            pins['positive'].append(zero_pts[i])
        if i < 0:
            pins['negative'].append(zero_pts[i])
    pin_p = max(pins['positive'])
    pin_n = min(pins['negative'])

    return {'hole':pin_p, 'electron':pin_n}


def get_oxid_states_composition(comp):
    """
    Get the oxidation states by composition, instead of bond valence theory.
    """
    comp = comp.reduced_composition
    elms = sorted(comp.keys(), key=lambda x:(x.group,-x.Z), reverse=True)
    symbs = [i.symbol for i in elms]
    oxi_range_short = {}
    oxi_range_long ={}
    oxi_range_mid ={}
    for elm in elms:
        if elm.is_alkali:
            oxi_range_short[elm]=[1]
            oxi_range_long[elm]=[1]
            oxi_range_mid[elm]=[1]
        elif elm.is_alkaline:
            oxi_range_long[elm]=[2]
            oxi_range_short[elm]=[2]
            oxi_range_mid[elm]=[2]
        elif elm.symbol == 'F':
            oxi_range_short[elm]=[-1]
            oxi_range_long[elm]=[-1]
            oxi_range_mid[elm]=[-1]
        elif elm.symbol == 'O':
            if 'F' in symbs:
                oxi_range_short[elm]=[-2,-1]
                oxi_range_long[elm]=[-2,-1]
                oxi_range_mid[elm]=[-2,-1]
            else:
                oxi_range_short[elm]=[-2]
                oxi_range_long[elm]=[-2]
                oxi_range_mid[elm]=[-2]
        else:
            if elm==elms[-1]:
                oxi_range_short[elm] = [i for i in elm.common_oxidation_states if i>0]
                oxi_range_mid[elm] = [i for i in elm.oxidation_states if i>0]
                oxi_range_long[elm]=[i for i in elm.oxidation_states if i>0]
            elif elm==elms[0]:
                oxi_range_short[elm] = [i for i in elm.common_oxidation_states if i<0]
                oxi_range_mid[elm] = [i for i in elm.common_oxidation_states if i<0]
                oxi_range_long[elm] = [i for i in elm.oxidation_states if i<0]
            else:
                oxi_range_short[elm] = elm.common_oxidation_states
                oxi_range_long[elm] = elm.oxidation_states
                oxi_range_mid[elm] = elm.oxidation_states
    def product(args, **kwds):
        pools = map(tuple, args) * kwds.get('repeat', 1)
        result = [[]]
        for pool in pools:
            result = [x+[y] for x in result for y in pool]
        for prod in result:
            yield tuple(prod)
    
    for i in product([oxi_range_short[elm] for elm in elms]):
        if not np.sum(np.array(i)*np.array([comp[elm] for elm in elms])):
            oxid_states = dict(zip(elms,i))
            return oxid_states
    for i in product([oxi_range_mid[elm] for elm in elms]):
        if not np.sum(np.array(i)*np.array([comp[elm] for elm in elms])):
            oxid_states = dict(zip(elms,i))
            return oxid_states
    for i in product([oxi_range_long[elm] for elm in elms]):
        if not np.sum(np.array(i)*np.array([comp[elm] for elm in elms])):
            oxid_states = dict(zip(elms,i))
            return oxid_states
    print 'Get oxidation states failed, Please input mannuly!'
    oxid_states={}
    for elm in elms:
        oxid_states[elm]=int(raw_input(elm.symbol+':  '))
    return oxid_states

def get_oxidation_states(stru):
    try:
        valences = BVAnalyzer().get_valences(stru)
        elems = [i for i in stru.species]
        oxids = set(zip(elems,valences))
        reduced_elems = set(elems)
        out={i:[] for i in reduced_elems}
        for i in oxids:
            if i[1] not in out[i[0]]:
                out[i[0]].append(i[1])
        out1 = {i:np.sign(out[i][0])*np.max(np.abs(out[i])) for i in out}
        return out1
    except:
        oxi_states = {}
        oxi_range = {}
        reduc_comp = stru.composition.reduced_composition
        oxid_states = get_oxid_states_composition(reduc_comp)
        return oxid_states

def get_defect_charge_list(stru, oxi_range_from=None, fill=False):
    """
    stru: pymatgen Structure object
    oxi_range_from: pmg or others
    fill: bool value 
    """
    elems = stru.composition.elements
    qlist={}
    oxid_states = get_oxidation_states(stru)
    anion = [i for i in oxid_states if oxid_states[i]>0]
    cation = [i for i in oxid_states if oxid_states[i]<0]
    max_min_oxid={}
    for x in stru.composition.elements:
        max_min_oxid[x]=[-100,100]
        if oxi_range_from == 'pmg':
            max_min_oxid[x][0] = min(x.oxidation_states)
            max_min_oxid[x][1] = max(x.oxidation_states)
        else:
            max_min_oxid[x][0] = oxi_state_range[x][0]
            max_min_oxid[x][1] = oxi_state_range[x][1]

        if oxid_states[x] > max_min_oxid[x][1]:
            max_min_oxid[x][1] = oxid_states[x]
        elif oxid_states[x] < max_min_oxid[x][0]:
            max_min_oxid[x][0] = oxid_states[x]

        if oxid_states[x]>0:
            if max_min_oxid[x][0] < 0:
                max_min_oxid[x][0] = 0
        else:
            if max_min_oxid[x][1] > 0:
                max_min_oxid[x][1] = 0
       
    ############################################################
    allowed_subst={}
    if len(cation)>1:
        for i in xrange(len(cation)-1):
            allowed_subst[cation[i]]=cation[:i]+cation[i+1:]
        allowed_subst[cation[len(cation)-1]]=cation[:-1]
    if len(anion)>1:
        for i in xrange(len(anion)-1):
            allowed_subst[anion[i]]=anion[:i]+anion[i+1:]
        allowed_subst[anion[len(anion)-1]]=anion[:-1]

    for elem in elems:
        oxid_stat = oxid_states[elem]
        name_vac = '%s_vac' % elem.symbol
        name_inter = '%s_inter' % elem.symbol
        qs_vac = [-q for q in range(max_min_oxid[elem][0],max_min_oxid[elem][-1]+1)]
        qs_inter = [q for q in range(max_min_oxid[elem][0],max_min_oxid[elem][-1]+1)]
        if fill:
            if min(qs_vac)*max(qs_vac) > 0:
                qs_vac = [np.sign(qs_vac[0])*i for i in range(np.max(np.abs(qs_vac))+1)]
                qs_inter = [np.sign(qs_inter[0])*i for i in range(np.max(np.abs(qs_inter))+1)]
        qlist[name_vac] = qs_vac
        qlist[name_inter] = qs_inter
    for host in allowed_subst:
        for subst in allowed_subst[host]:
            name_subst = '%s_subst_%s' % (host.symbol, subst.symbol)
            qhost = oxid_states[host]
            qsubsts = range(max_min_oxid[subst][0],max_min_oxid[subst][-1]+1)
            qs_subst = [i-qhost for i in qsubsts]
            if fill:
                if max(qs_subst)*min(qs_subst) > 0:
                    sign=np.sign(qs_subst[0])
                    qs_subst = [sign*i for i in range(np.max(np.abs(qs_subst))+1)]
            qlist[name_subst]=qs_subst
    return qlist
