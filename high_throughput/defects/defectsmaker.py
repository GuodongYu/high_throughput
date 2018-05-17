#!/usr/bin/env python

__author__ = "Geoffroy Hautier"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "November 4, 2012"

from pymatgen import Molecule
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer,PointGroupAnalyzer
from pymatgen.core.periodic_table import DummySpecie, Element
from pymatgen.core.structure import PeriodicSite,Structure
from pymatgen.analysis.bond_valence import BVAnalyzer
import numpy as np
from sys import exit
from scipy.spatial import Voronoi,ConvexHull 
from scipy.spatial.qhull import QhullError
import random
from pycdt.corrections.kumagai_correction import find_defect_pos,closestsites
from high_throughput.defects.utils.util import closest_defect_distance
from high_throughput.defects.constants import *
from high_throughput.defects.utils.util import *
from ase.build import find_optimal_cell_shape



def get_nearest_neighbours(index,struct):
    dists=list(struct.distance_matrix[index])
    dists.remove(dists[index])
    min_dist=np.array(dists).min()
    nearest_neighbour_infos=struct.get_neighbors(struct[index],min_dist+0.2,include_index=True)
    nearest_neighbours=[i[2] for i in nearest_neighbour_infos]
    return nearest_neighbours

def get_nearest_neighbours_from_point(vec_cartesian,struct):
    site = PeriodicSite(DummySpecie(),vec_cartesian,struct.lattice,coords_are_cartesian=True)
    dists = [i[1] for i in struct.get_neighbors(site,10,include_index=True) if i[1]>0.1]
    min_dist = min(dists)
    nearest_neighbour_infos=struct.get_neighbors(site,min_dist+0.2,include_index=True)
    nearest_neighbours=[i[2] for i in nearest_neighbour_infos if i[1]>0.1]
    return nearest_neighbours

def get_struct_with_nearest_neighbours_movement(vec_cart,stru):
    structure=stru.copy()
    nearest_neighbours=get_nearest_neighbours_from_point(vec_cart,structure)
    smin=0.1;smax=0.3
    la,lb,lc = structure.lattice.abc
    for i in nearest_neighbours:
        atom=structure[i]
        da=random.choice([random.uniform(smin/la, smax/la),random.uniform(-smax/la, -smin/la)])
        db=random.choice([random.uniform(smin/lb, smax/lb),random.uniform(-smax/lb, -smin/lb)])
        dc=random.choice([random.uniform(smin/lc, smax/lc),random.uniform(-smax/lc, -smin/lc)])
        structure[i]=atom.specie,[atom.a+da,atom.b+db,atom.c+dc]
    return structure

def get_site_in_UC(stru_dfc,stru_bulk_sc,stru_bulk_uc):
    spa = SpacegroupAnalyzer(stru_bulk_uc)
    stru_bulk_uc_symm = spa.get_symmetrized_structure()
    def get_neighbors(site,stru):
        neighs = stru.get_neighbors(site,20.)
        neighs = sorted(neighs,key=lambda x:x[1])
        neighs = [(i[0].specie.symbol,round(i[1],2)) for i in neighs]
        dists = list(set([i[1] for i in neighs]))
        dists.sort()
        cut = (dists[1]+dists[2])/2.0
        neighs = [i for i in neighs if i[1]<=cut]
        return set(neighs)
    pos_blk, pos_dfc = find_defect_pos(stru_bulk_sc, stru_dfc)
    if pos_blk is None and pos_dfc is None:
        raise KeyError('Can not determine defect type')
    elif pos_blk is None:
        defect_type = 'inter'
    elif pos_dfc is None:
        defect_type = 'vac'
    else:
        defect_type = 'subst'
    site, dist, index = closestsites(stru_bulk_sc, stru_dfc,pos_blk)[0]
    if defect_type == 'subst':
        subst = closestsites(stru_bulk_sc, stru_dfc,pos_blk)[1][0].specie
    else:
        subst = None
    neighs_sc =  get_neighbors(site, stru_bulk_sc)
    for s in stru_bulk_uc_symm.equivalent_sites:
        if s[0].specie != site.specie:
            continue
        neighs_uc = get_neighbors(s[0], stru_bulk_uc_symm)
        if neighs_uc == neighs_sc:
            return s[0]

def inds_after_sym(p_inds,stru):
    spa = SpacegroupAnalyzer(stru,symprec=0.5)
    se = spa.get_symmetrized_structure()    
    equ_inds = se.equivalent_indices
    ind_dict ={}
    for i in range(len(equ_inds)):
        ind_dict[i]=equ_inds[i]
    out = []
    for i in p_inds:
        for j in ind_dict:
            if i in ind_dict[j]:
                out.append(j)
                break
    return out

def polyhedral_largest(frac_cord, struct):
    cart = struct.lattice.get_cartesian_coords(frac_cord)
    neigs_mess = struct.get_neighbors_in_shell(cart,0,10,include_index=True)
    neigs = sorted(neigs_mess,key=lambda i:(i[1],i[0].specie.symbol))
    k=0
    OK = False
    while not OK:
        sites = neigs[0:3+k]
        k+=1
        inds = [i[2] for i in sites]
        dis = [round(i[1],1) for i in sites]
        species = [i[0].specie.symbol for i in sites]
        dis_ = [dis[0]]
        for i in dis[1:]:
            equ = False
            for j in dis_:
                if abs(i-j)<=0.2:
                    dis_.append(j)
                    equ = True
                    break
            if not equ:
                dis_.append(i)
        pts = [i[0].coords for i in sites]
        try:
            poly = ConvexHull(pts)
        except QhullError:
            d3 = False
            continue
        else:
            if len(sites) == len(poly.vertices):
                poly_out = poly
                inds_out = inds
                dis_out = dis_
                species_out = species
                d3 = True
            elif len(sites) != len(poly.vertices) and d3:
                OK = True
                if d3:
                    poly={}
                    molecule = Molecule(species_out,poly_out.points)
                    pga = PointGroupAnalyzer(molecule)
                    nsymmop = len(pga.symmops)
                    #ind_dis = zip(inds_after_sym(inds_out,struct),dis_out)
                    #ind_dis = sorted(ind_dis, key=lambda i:(i[1],i[0]))
                    poly['frac_coord']=frac_cord
                    poly['poly']=poly_out
                    poly['species'] = species_out
                    poly['site_ind']=  inds_after_sym(inds_out,struct)
                    poly['dist'] = dis_out
                    poly['nsymmop'] = nsymmop
                    poly['struct'] = struct
                    return poly
                else:
                    return 

def polyhedral_cut(poly):
    out={}
    out['frac_coord']=poly['frac_coord']
    dist_n = {}
    dis = [i for i in poly['dist']]
    for i in dis:
        if i not in dist_n:
            dist_n[i]=0
        dist_n[i]+=1
    dis = list(set(dis))
    dis.sort()
    n=0
    for i in range(len(dis)):
        if abs(dis[i]-dis[0])<=0.5:
            n+=dist_n[dis[i]]
    pts = poly['poly'].points[0:n]
    try:
        convex = ConvexHull(pts)
        verts = list(convex.vertices)
    except QhullError:
        return
    else:
        card=poly['struct'].lattice.get_cartesian_coords(poly['frac_coord'])
        pts_add = np.array(list(pts)+[card])
        conv_ = ConvexHull(pts_add)
        if conv_.volume > convex.volume:
            return
        else:
            out['poly']=convex
            out['species'] = poly['species'][0:n]
            out['dist'] = poly['dist'][0:n]
            out['site_ind']=poly['site_ind'][0:n]
            mole = Molecule(out['species'],pts)
            pga = PointGroupAnalyzer(mole)
            nsymmop = len(pga.symmops)
            out['nsymmop']=nsymmop
            return out

def distance(frac1,frac2,stru):
    dis = stru.lattice.get_distance_and_image(frac1, frac2)
    return dis[0]

def screen_polys_dist(polys,stru,dis_bar=2.0):
    ref=[polys[0]]
    for i in polys[1:]:
        match = False
        for j in range(len(ref)):
            if distance(i['frac_coord'],ref[j]['frac_coord'],stru)<=dis_bar:
                match = True
                #if i[1].volume > ref[j][1].volume:
                if i['nsymmop'] >1 and ref[j]['nsymmop']==1:
                    ref[j]=i
                elif i['poly'].volume > ref[j]['poly'].volume:
                    ref[j]=i
                break
        if not match:
            ref.append(i)
    ref.sort()
    return ref

def screen_polys(polys,num,screen_criteria='symmetry'):
    """
    num: the total number of interstial sites 
    screen_criteria: symmetry or distance, which means the site with high symmetry(distance to neighbors) has the higher priority
                                          
    """
    ref = [polys[0]]
    for i in polys[1:]:
        match = False
        for j in ref:
                check = abs(i['poly'].volume-j['poly'].volume) < 0.1 and abs(i['poly'].area-j['poly'].area)< 0.1
                if check:
                    if i['poly'].volume > j['poly'].volume:
                        ref[ref.index(j)] = i 
                    match =True
                    break
        if not match:
            ref.append(i)  
    for i in ref:
        print i['poly'].volume, i['poly'].area, i['site_ind'],i['species'], i['dist'], i['nsymmop']
    out_symm = sorted(ref,key=lambda i:(i['nsymmop'],min(i['dist'])),reverse=True)
    out_dist = sorted(ref,key=lambda i:min(i['dist']),reverse=True)
    
    if num == 'all' or num >=len(out_symm):
        return out_symm
    elif num <= 2:
        if out_symm[0] != out_dist[0]:
            return [out_symm[0],out_dist[0]]
        else:
            return out_symm[0:num]
    else:
        out = [out_symm[0],out_dist[0]]
        for i in range(1,len(out_symm)):
            if len(out) == num:
                return out
            if out_symm[i] not in out:
                out.append(out_symm[i])

def get_interstitial_sites(structure, dis_bar=2.0, num=2):
    """
    get interstitial positions(IntrPoints) for a specific structure. 
    method:
    1st, collect all atoms in a 3x3x3 supercell.  
    2nd, get all Voronoi cell cornors by Voronoi analysis 
    3rd, take all the atoms and Voronoi cell cornors (DummySpecie atoms) in the centeral unit cell as a new unit cell to make a crystal.
    4th, keep nonequivalent DummySpecie atoms by space group symmetry analysis.
    5th, if the distances among some DummySpecie atoms are less than 'criter', just leave one of them.
    6th, screen the final DummySpecie atom to make each one have different neighbours. 
            
    Args:
        struct: the structure under consideration for obtaining the interstitial positions.
        standadized: True or False. must be the same as the standidized parameter in DefectsMaker function.  
        save_inters: True or False. If True, the unitcell structure with all interstitial positions occupied by 'Lr' or 'No' atoms are output with cif format.
    """
    criter=2.0
    struct = structure.copy()
    struct_orig = structure.copy()
    lat=struct.lattice.matrix
    def struct_with_voronoi_verts(structure):
        struct = structure.copy()
        sites_cart=struct.cart_coords
        s=[]
        for x in sites_cart:
            for i in range(-2,4):
                for j in range(-2,4):
                    for k in range(-2,4):
                        s.append(x+i*lat[0]+j*lat[1]+k*lat[2])
        vor=Voronoi(s)
        verts=vor.vertices
        inters=[]
        for i in verts:
            fracs= struct.lattice.get_fractional_coords(i)
            inbox=True
            for j in fracs:
                if j<0 or j>=1.0:
                    inbox = False
                    break
            if inbox:
                inters.append(list(fracs)) 
                struct.append(DummySpecie(),fracs)
        return struct, inters
    def get_inter_sites_after_symmetry(struct_with_voronoi_vects):
        spa=SpacegroupAnalyzer(struct_with_voronoi_vects,symprec=0.001, angle_tolerance=1)
        stru=spa.get_symmetrized_structure()
        eqs=stru.equivalent_sites
        re_inters=[]
        for i in eqs:
            if i[0].specie==DummySpecie():
                re_inters.append(list(i[0].frac_coords))
        re_inters.sort()
        return re_inters
    struct_withvv, inters =  struct_with_voronoi_verts(struct)
    try:
        inters_sym = get_inter_sites_after_symmetry(struct_withvv,)
    except:
        inters_sym = inters    
    polys_all = [polyhedral_largest(i,struct_orig) for i in inters_sym]
    polys_cut = [polyhedral_cut(i) for i in polys_all if i]
    polys = [i for i in polys_cut if i]
    OK = False
    while not OK:
        polys_dis = screen_polys_dist(polys,struct_orig,dis_bar)
        if polys == polys_dis:
            OK = True
        else:
            polys = polys_dis
    inter_sites =  [i['frac_coord'] for i in screen_polys(polys_dis,num)]
    for i in inter_sites:
        struct_orig.append(DummySpecie(),i)
    struct_orig.to('poscar','POSCAR')
    return inter_sites

class DefectStructMaker(object):
    """
    a class to generate defective structures in supercells in view of their computations with a PW code
    The standard defects such as antisites, vacancies are generated
    and appropriate supercells are made
    TODO: develop a better way to find interstitials
    """
    def __init__(self, struct_bulk, min_dis=13.0, extrapolation=False):
        """
        Args:
            struct_bulk_source: the bulk structure, pymatgen Structure object
            struct_defect: unconverged defect structure, the pymatgen Structure object, the one without relaxation is better. 
            defect_type: the name of the defect on which the potential alignment isn't converged
            charge: the defect charge
            dfc_dist_add: the added distance between defects based on the origin one
        """
        self.extrapolation = extrapolation
        self.struct_backup = struct_bulk
        finder = SpacegroupAnalyzer(struct_bulk,symprec=0.5)
        self.st_prim=finder.get_primitive_standard_structure()
        self.st_conv=finder.get_conventional_standard_structure()
        self.min_dis = min_dis
        self.find_optimized_supercell()
        best_supercell = self.optimized_supercell
        finder = SpacegroupAnalyzer(self.struct_orig,symprec=0.5) 
        self.struct = finder.get_symmetrized_structure()
        self.inter_sites = get_interstitial_sites(self.struct_orig)
        
    def find_optimized_supercell(self):
        def defect_dis(lat):
            diss=[]
            for i in [-1,0,1]:
                for j in [-1,0,1]:
                    for k in [-1,0,1]:
                        if not (i or j or k):
                            continue
                        vec=lat[0]*i+lat[1]*j+lat[2]*k
                        dis=np.linalg.norm(vec)
                        diss.append(dis)
            return min(diss)
        def lat_contants_diff(lat):
            a=np.linalg.norm(lat[0])
            b=np.linalg.norm(lat[1])
            c=np.linalg.norm(lat[2])
            dab=abs(a-b)
            dac=abs(a-c)
            dbc=abs(b-c)
            diffs=[dab,dac,dbc]
            return max(diffs)
        trans=[]
        strus={'primitive':self.st_prim,'conventional':self.st_conv}
        for uc_type in strus:
            stru = strus[uc_type]
            lat_uc = stru.lattice.matrix
            for sc_type in ['sc','fcc']:
                gp={}
                gp['uc_type']=uc_type
                gp['sc_type']=sc_type
                dis=0.0
                i=1
                while dis < self.min_dis:
                    p=find_optimal_cell_shape(lat_uc, i, sc_type)
                    lat_sc=np.dot(p,lat_uc)
                    dis=defect_dis(lat_sc)
                    lat_diff=lat_contants_diff(lat_sc)
                    i+=1
                gp['dis']=dis
                gp['dis-']=-1.0*dis
                gp['lat_diff']=lat_diff
                gp['num']=(i-1)*stru.num_sites
                gp['supercell']=p
                trans.append(gp)
        trans.sort(key=lambda x:(x['num'],x['dis-'],x['lat_diff']))
        choice=trans[0] 
        self.struct_orig = strus[choice['uc_type']]
        self.optimized_supercell=[choice['supercell']]
        self.uc_type=choice['uc_type']
        self.sc_type=choice['sc_type']
        dis=choice['dis']
        num=choice['num']
        if self.extrapolation:
            j=2
            lat_uc = self.struct_orig.lattice.matrix
            for dis_new in [self.min_dis-3., self.min_dis-6.]:
                dis=0.0
                i=1
                while dis < dis_new:
                    p=find_optimal_cell_shape(lat_uc, i, self.sc_type)
                    lat_sc=np.dot(p,lat_uc)
                    dis=defect_dis(lat_sc)
                    lat_diff=lat_contants_diff(lat_sc)
                    i+=1
                self.optimized_supercell.append(p)
                j+=1

    def make_defect_cell_vacancy(self, target_site, supercell, move_neighbours=True):
        """
        Create a supercell for a vacancy
        TODO: This method needs to be changed to something smarter
        """
        structure=self.struct_orig.copy()
        for i in range(structure.num_sites):
            s = structure._sites[i]
            if s.distance_from_point(target_site.coords)<0.001:
                index_uc=i
        defect_site=[False]*structure.num_sites
        defect_site[index_uc]=True
        structure.add_site_property('defect_site',defect_site)
        structure.make_supercell(supercell,to_unit_cell=True)
        for i in range(structure.num_sites):
            if structure[i].properties['defect_site']:
                index=i
                break
        la,lb,lc=structure.lattice.abc
        if move_neighbours==True:
            nearest_neighbours=get_nearest_neighbours(index,structure)
            smin=0.1;smax=0.3
            for i in nearest_neighbours:
                atom=structure[i]
                da=random.choice([random.uniform(smin/la, smax/la),random.uniform(-smax/la, -smin/la)])
                db=random.choice([random.uniform(smin/lb, smax/lb),random.uniform(-smax/lb, -smin/lb)])
                dc=random.choice([random.uniform(smin/lc, smax/lc),random.uniform(-smax/lc, -smin/lc)])
                structure[i]=atom.specie,[atom.a+da,atom.b+db,atom.c+dc]

        structure.replace(index, DummySpecie())
        structure.remove_species([DummySpecie()])
        return structure
    
    def make_defect_interstitial(self, target_site, supercell):
        structure=self.struct_orig.copy()
        structure.append(DummySpecie(), target_site.frac_coords)
        structure.make_supercell(supercell,to_unit_cell=True)
        for i in range(structure.num_sites):
            if structure[i].specie == DummySpecie():
                structure.replace(i,target_site.specie)
                break
        structure.remove_species([DummySpecie()])
        return structure
            
    def make_defect_cell_intrinsic_subst_defects(self, target_site, subst, supercell, move_neighbours=True):
        """
        subst defines elements that you allow on this site
        charge is the original oxidation state of the element on the target site
        """
        struct=self.struct_orig.copy()
        for i in range(struct.num_sites):
            s = struct._sites[i]
            if s.distance_from_point(target_site.coords)<0.001:
                index_uc=i
        defect_site=[False]*struct.num_sites
        defect_site[index_uc]=True
        struct.add_site_property('defect_site',defect_site)
        struct.make_supercell(supercell,to_unit_cell=True)
        for i in range(struct.num_sites):
            if struct[i].properties['defect_site']:
                index=i
                break
        la,lb,lc=struct.lattice.abc
        if move_neighbours==True:
            nearest_neighbours=get_nearest_neighbours(index,struct)
            smin=0.1;smax=0.3
            for i in nearest_neighbours:
                atom=struct[i]
                da=random.choice([random.uniform(smin/la, smax/la),random.uniform(-smax/la, -smin/la)])
                db=random.choice([random.uniform(smin/lb, smax/lb),random.uniform(-smax/lb, -smin/lb)])
                dc=random.choice([random.uniform(smin/lc, smax/lc),random.uniform(-smax/lc, -smin/lc)])
                struct[i]=atom.specie,[atom.a+da,atom.b+db,atom.c+dc]
        struct.replace(index, subst)
        return struct

class DefectsMaker(DefectStructMaker):
    """
    a class to generate defective structures in supercells in view of their computations with a PW code
    The standard defects such as antisites, vacancies are generated
    and appropriate supercells are made
    TODO: develop a better way to find interstitials
    """
    def __init__(self, structure, max_min_oxid, allowed_subst, oxid_states, min_dis=13.,
                 inter=False, interstitials_extr_elts=[], vac=True, extrapolation=True):
        """
        Args:
            structure:
                the bulk structure
            max_min_oxid:
                the minimal and maximum oxidation state of each element as a dict. For instance
                {Element("O"):(-2,0)}
            allowed_subst:
                the allowed substitutions between elements. Including intrinsic (e.g., anti-sites) and
                extrinsic. Example: {Element("Co"):[Element("Zn"),Element("Mn")]} means Co sites can be substituted
                by Mn or Zn.
            oxid_states:
                the oxidation state of the elements in the compound e.g. {Element("Fe"):2,Element("O"):-2}
            size_limit:
                maximum atom numbers allowed in the supercell
            interstitials_sites:
                a list of PeriodicSites in the bulk structure on which we put an interstitial
            standardized:
                True means we use a standardized primitive cell
        """
        DefectStructMaker.__init__(self, structure, min_dis, extrapolation)
        self.defects = []
        nb_per_elts = {e:0 for e in structure.composition.elements}
        
        ### bulk###
        def get_bulk_sc(stru,sc_size):
            s=stru.copy()
            s.make_supercell(sc_size,to_unit_cell=True)
            return s
        self.defects.append({'short_name':'bulk','uc_type':self.uc_type,'sc_type':self.sc_type,'bulk_unitcell':self.struct_orig,'supercells':
               [{'size': s_size,'structure':get_bulk_sc(self.struct_orig, s_size)} for s_size in self.optimized_supercell]})
        self.defects.append({'short_name':'dielectric','structure':self.struct_orig,'uc_type':self.uc_type})
        charges_limit=range(-100,100)
        for s in self.struct.equivalent_sites:
            nb_per_elts[s[0].specie] = nb_per_elts[s[0].specie]+1
            #vac
            if vac==True:
                list_charges=[]
                for c in range(max_min_oxid[s[0].specie][0], max_min_oxid[s[0].specie][1]+1):
                    list_charges.append(-c)
                if 0 not in list_charges:
                    list_charges.append(0)
                list_charges=range(min(list_charges),max(list_charges)+1)
                list_charges=list(set(charges_limit)&set(list_charges))
                self.defects.append({'short_name': s[0].specie.symbol+str(nb_per_elts[s[0].specie])+"_vac",
                                     'unique_sites': s[0],
                                     'supercells':[{'size': s_size,'structure':self.make_defect_cell_vacancy(s[0], s_size), 
                                     'struct_no_move':self.make_defect_cell_vacancy(s[0], s_size,move_neighbours=False)} 
                                     for s_size in self.optimized_supercell],
                                     'charges':list_charges,
                                     'uc_type':self.uc_type,'sc_type':self.sc_type,
                                     'bulk_unitcell':self.struct_orig})
            #sub
            if s[0].specie in allowed_subst:
                for subst in allowed_subst[s[0].specie]:
                    list_charges_subst=[c-oxid_states[s[0].specie] for c in range(max_min_oxid[subst][0], max_min_oxid[subst][1]+1)]
                    if 0 not in list_charges_subst:
                        list_charges_subst.append(0)
                    list_charges_subst=range(min(list_charges_subst),max(list_charges_subst)+1)
                    list_charges_subst=list(set(charges_limit)&set(list_charges_subst))
                    self.defects.append({'short_name':s[0].specie.symbol+str(nb_per_elts[s[0].specie])+"_subst_"
                                +subst.symbol, 'unique_sites':s[0],'supercells':[{'size':s_size,'structure':
                                self.make_defect_cell_intrinsic_subst_defects(s[0], subst, s_size),
                                'struct_no_move':self.make_defect_cell_intrinsic_subst_defects(
                                s[0], subst, s_size,move_neighbours=False)} for s_size in self.optimized_supercell],
                                'charges':list_charges_subst,'uc_type':self.uc_type,'sc_type':self.sc_type,'bulk_unitcell':self.struct_orig})
            
        #interstitials
        if inter:
            for elt in self.struct.composition.elements:
                count = 1
                list_charges_inter=[c for c in range(max_min_oxid[elt][0],max_min_oxid[elt][1]+1)]
                if 0 not in list_charges_inter:
                    list_charges_inter.append(0)
                list_charges_inter=range(min(list_charges_inter),max(list_charges_inter)+1)
                list_charges_inter=list(set(charges_limit)&set(list_charges_inter))
                for frac_coord in self.inter_sites:
                    self.defects.append({'short_name':elt.symbol+str(count)+"_inter",
                                         'unique_sites':PeriodicSite(elt, frac_coord, structure.lattice),
                                         'supercells':[{'size':s_size,
                                                        'structure':self.make_defect_interstitial(
                                                            PeriodicSite(elt, frac_coord, self.struct_orig.lattice), s_size)}
                                                       for s_size in self.optimized_supercell],
                                        'charges':list_charges_inter,
                                        'uc_type':self.uc_type,'sc_type':self.sc_type,
                                        'bulk_unitcell':self.struct_orig})
                    count = count+1
        for elt in interstitials_extr_elts:
            count = 1
            list_charges_subst=[c for c in range(max_min_oxid[elt][0],max_min_oxid[elt][1]+1)]
            if 0 not in list_charges_subst:
                list_charges_subst.append(0)
            list_charges_subst=range(min(list_charges_subst),max(list_charges_subst)+1)
            list_charges_subst=list(set(charges_limit)&set(list_charges_subst))
            for frac_coord in self.inter_sites:
                self.defects.append({'short_name':elt.symbol+str(count)+"_inter",
                                     'unique_sites':PeriodicSite(elt, frac_coord, structure.lattice),
                                     'supercells':[{'size':s_size,
                                                    'structure':self.make_defect_interstitial(
                                                        PeriodicSite(elt, frac_coord, self.struct_orig.lattice), s_size)}
                                                   for s_size in self.optimized_supercell],
                                     'charges':list_charges_subst,'uc_type':self.uc_type,'sc_type':self.sc_type,
                                    'bulk_unitcell':self.struct_orig})
                count = count+1

def IntrinsicDefectsMaker(struct, vac=True, sub=True, inter=False,\
                             min_dis=13., extrapolation=False,**kwargs):
    """
    generate intrinsic defects for a given structure including vacancies, substitutions and interstitials 
    Args:
        struct: the structure for generating the intrinstic defects
        standardized: whether converting the structrue into the priv
        vac: whether considering intrinsic vacancies 
        sub: whether considering intrinsic substitutions
        inter: whether considering intrinsic interstitials
        save_inters: whether output all interstitial sites into a cif file. The tag only makes sense when inter=True
        kwargs: if Bond velence theory can not give the oxidation state, give them manuelly by kwargs. 
    """
    oxid_states = get_oxidation_states(struct)
    anion = [i for i in oxid_states if oxid_states[i]>0]
    cation = [i for i in oxid_states if oxid_states[i]<0]
    max_min_oxid={}
    oxi_range_from = 'pmg'
    for x in struct.composition.elements:
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

########## Define dict allowed_subst ###############
    allowed_subst={}
    if sub:
        if len(cation)>1:
            for i in xrange(len(cation)-1):
                allowed_subst[cation[i]]=cation[:i]+cation[i+1:]
            allowed_subst[cation[len(cation)-1]]=cation[:-1]
        if len(anion)>1:
            for i in xrange(len(anion)-1):
                allowed_subst[anion[i]]=anion[:i]+anion[i+1:]
            allowed_subst[anion[len(anion)-1]]=anion[:-1]
    interstitials_extr_elts=[]
    return DefectsMaker(struct,max_min_oxid,allowed_subst,oxid_states, min_dis, inter, interstitials_extr_elts, vac, extrapolation)

def ExtrinsicDefectsMaker(struct, subst={}, inter_ext_elts=[], min_dis=13., extrapolation=False):
    """
    generate intrinsic defects for a given structure including vacancies, substitutions and interstitials 
    Args:
        struct: the structure for generating the intrinstic defects
        subst: whether considering intrinsic substitutions
        inter: whether considering intrinsic interstitials
    """
    subst = {Element(i):[Element(j) for j in subst[i]] for i in subst}
    ext_elts = [i for i in inter_ext_elts]
    for i in subst:
        ext_elts+=subst[i]
    ext_elts = list(set(ext_elts))
    oxid_states = get_oxidation_states(struct)
    for i in ext_elts:
        oxid_states[i] = i.common_oxidation_states[0]
    max_min_oxid={}
    oxi_range_from = 'pmg'
    for x in struct.composition.elements+ext_elts:
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
    vac = False
    inter = False
    return DefectsMaker(struct,max_min_oxid,subst,oxid_states, min_dis, inter, inter_ext_elts, vac, extrapolation)

class DefectRemaker(DefectStructMaker):
    """
    a class to generate a new defect structure based on the existed defects tructure when increasing
    the defect distance. 
    """
    def __init__(self, struct_bulk_sc, struct_defect, defect_name, charge, dfc_dist_add=2.0, dist_min=13.0,**kws):
        """
        Args:
            struct_bulk_sc: the supercell bulk structure with the same cell size as struct_defect, pymatgen Structure object
            struct_defect: defect structure, the pymatgen Structure object, the one without relaxation is better. 
            defect_type: the name of the defect
            charge: the defect charge
            dfc_dist_add: the added distance between defects based on the origin one
        """
        if kws:
            self.old_bulk_uc=kws['stru_bulk_uc']
        else:
            self.old_bulk_uc=None
        dfc_dist_previous = closest_defect_distance(struct_defect)-0.01
        dfc_dist = dfc_dist_previous + dfc_dist_add
        dfc_dist = max(dfc_dist, dist_min)
        DefectStructMaker.__init__(self,struct_bulk_sc,dfc_dist,extrapolation=False)
        self.old_bulk_sc = struct_bulk_sc
        self.old_dfc_stru = struct_defect
        self.defect_name = defect_name
        self.charge = charge
        self.defects=[]
        self.defect_backup=[]
        self.site_inUC, self.dfc_type, self.subst = self.get_defect_site_inUC()
        if 'subst' in defect_name:
            if self.subst != Element(defect_name.split('_')[-1]):
                raise KeyError("substition elements don't match between defect_name and defect structure")
        if self.dfc_type not in self.defect_name:
            raise KeyError("defect_name doesn't match defect structure!!")
        def get_bulk_sc(stru,sc_size):
            s=stru.copy()
            s.make_supercell(sc_size,to_unit_cell=True)
            return s
        bulk_sc = {'short_name':'bulk','uc_type':self.uc_type,'sc_type':self.sc_type,'bulk_unitcell':self.struct_orig,'supercells':
                 [{'size': s_size,'structure':get_bulk_sc(self.struct_orig, s_size)} for s_size in self.optimized_supercell]}
        self.defects.append(bulk_sc)
        self.defects.append({'short_name':'dielectric','structure':self.struct_orig,'uc_type':self.uc_type})
        if self.dfc_type == 'vac':
            defect = {'short_name': self.defect_name,'unique_sites': self.site_inUC,'sc_type':self.sc_type,
                             'supercells':[{'size': s_size,'structure':self.make_defect_cell_vacancy(self.site_inUC, s_size), 
                             'struct_no_move':self.make_defect_cell_vacancy(self.site_inUC, s_size,move_neighbours=False)} 
                                            for s_size in self.optimized_supercell],
                             'charges':[self.charge], 'uc_type':self.uc_type,'bulk_unitcell':self.struct_orig}
            self.defects.append(defect)
        if self.dfc_type == 'subst':
            defect = {'short_name':self.defect_name, 'unique_sites':self.site_inUC,'sc_type':self.sc_type,
                            'supercells':[{'size':s_size,
                            'structure':self.make_defect_cell_intrinsic_subst_defects(self.site_inUC, self.subst, s_size),
                            'struct_no_move':self.make_defect_cell_intrinsic_subst_defects(
                            self.site_inUC, self.subst, s_size,move_neighbours=False)}
                            for s_size in self.optimized_supercell],
                            'charges':[self.charge],'uc_type':self.uc_type,'bulk_unitcell':self.struct_orig}
            self.defects.append(defect)
        if self.old_bulk_uc:
            sc_size=get_sc_size(self.old_bulk_uc,self.old_bulk_sc)
            pos_blk, pos_dfc = find_defect_pos(self.old_bulk_sc, self.old_dfc_stru)
            if pos_dfc is None:
                pos_dfc = pos_blk
            defect_backup = {'short_name':self.defect_name, 
                            'unique_sites':get_site_in_UC(self.old_dfc_stru,self.old_bulk_sc,self.old_bulk_uc),
                            'sc_type':'unknown',
                            'supercells':[{'size':sc_size, 'structure': get_struct_with_nearest_neighbours_movement(
                                    pos_dfc,self.old_dfc_stru),'struct_no_move':self.old_dfc_stru}],
                            'charges':[self.charge],'uc_type':'unknown','bulk_unitcell':self.old_bulk_uc}
            self.defect_backup.append(defect_backup)


    def get_defect_site_inUC(self):    
        def get_neighbors(site,stru):
            neighs = stru.get_neighbors(site,20.)
            neighs = sorted(neighs,key=lambda x:x[1])
            neighs = [(i[0].specie.symbol,round(i[1],2)) for i in neighs]
            dists = list(set([i[1] for i in neighs]))
            dists.sort()
            cut = (dists[1]+dists[2])/2.0
            neighs = [i for i in neighs if i[1]<=cut]
            return set(neighs)
        pos_blk, pos_dfc = find_defect_pos(self.old_bulk_sc, self.old_dfc_stru)
        if pos_blk is None and pos_dfc is None:
            raise KeyError('Can not determine defect type')
        elif pos_blk is None:
            defect_type = 'inter'
        elif pos_dfc is None:
            defect_type = 'vac'
        else:
            defect_type = 'subst'
        site, dist, index = closestsites(self.old_bulk_sc, self.old_dfc_stru,pos_blk)[0]
        if defect_type == 'subst':
            subst = closestsites(self.old_bulk_sc, self.old_dfc_stru,pos_blk)[1][0].specie
        else:
            subst = None
        neighs_sc =  get_neighbors(site, self.old_bulk_sc)
        for s in self.struct.equivalent_sites:
            if s[0].specie != site.specie:
                continue
            neighs_uc = get_neighbors(s[0], self.struct)
            if neighs_uc == neighs_sc:
                return (s[0], defect_type, subst)

