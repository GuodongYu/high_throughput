from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.sites import PeriodicSite
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pycdt.core.defects_analyzer import ComputedDefect
from high_throughput.defects.utils.plotter import DefectPlotter
from pycdt.core.defects_analyzer import DefectsAnalyzer
from high_throughput.config import *
from high_throughput.defects.utils.util import *
import json
import subprocess
from os.path import *
import os
from pymongo.collection import Collection,ObjectId
from pycdt.corrections.utils import closestsites
from high_throughput.defects.kumagai_correction import KumagaiCorrection, KumagaiBulkInit, read_ES_avg
from high_throughput.defects.db_config import *
from high_throughput.config import *
from monty.json import jsanitize, MontyEncoder, MontyDecoder
from monty.serialization import dumpfn
import numpy as np
import itertools
from monty.os import cd
import logging
logging.basicConfig()


class TasksIdsClassify(object):
    def __init__(self):
        self.collection_name = tasks_collection_name
        self.collection = Collection(db,self.collection_name)
        self._classify()
        self._num_unsuccessuful()

    def _classify(self):
        groups={}
        for i in ['defect','bulk','dielectric','band']:
            groups[i]={}
        for i in ['unsuccessful','unknown']:
            groups[i]=[]
        import re
        patt=re.compile('(.*vac|.*subst|.*inter|.*sub)')
        for i in self.collection.find({},['transformations','state','nsites']):
            if i['state'] != 'successful':
                groups['unsuccessful'].append(i['_id'])
                continue
            try:
                mpid=i['transformations']['history'][0]['source']
                defect_type=i['transformations']['defect_type']
                if defect_type=='bulk':
                    if mpid not in groups['bulk']:
                        groups['bulk'][mpid]=[]
                    groups['bulk'][mpid].append((i['_id'],i['nsites']))
		elif defect_type == 'band':
		    if mpid not in groups['band']:
			groups['band'][mpid]=[]
		    groups['band'][mpid].append(i['_id'])
                elif defect_type=='dielectric':
                    if mpid not in groups['dielectric']:
                        groups['dielectric'][mpid]=[]
                    groups['dielectric'][mpid].append(i['_id'])
                elif patt.match(defect_type):
                    if mpid not in groups['defect']:
                        groups['defect'][mpid]={}
                    charge = i['transformations']['charge']
                    name='%s:%s' % (defect_type,str(charge))
                    if name not in groups['defect'][mpid]:
                        groups['defect'][mpid][name]=[]
                    groups['defect'][mpid][name].append((i['_id'],i['nsites']))
                else:
                    groups['unknown'].append(i['_id'])
            except:
                groups['unknown'].append(i['_id'])
        self.groups_detail = groups.copy()
        groups_short={}
        groups_short['dielectric']=groups['dielectric']
        groups_short['band']=groups['band']
        groups_short['defect']={}
        groups_short['bulk']={}
        for mpid in groups['bulk']:
            groups_short['bulk'][mpid]=[]
            for i in groups['bulk'][mpid]:
                groups_short['bulk'][mpid].append(i[0])
        for mpid in groups['defect']:
            groups_short['defect'][mpid]=[]
            for name in groups['defect'][mpid]:
                id_nsite = groups['defect'][mpid][name]
                id_nsite.sort(key=lambda x:x[1],reverse=True)
                groups_short['defect'][mpid].append(id_nsite[0][0])
        self.groups=groups_short

    def _num_unsuccessuful(self):
        num=0
        for i in self.collection.find({},['state']):
            if i['state'] != 'successful':
                num=num+1
        self.num_unsuccessuful=num


    def _ids_strus_group(self, ids):
        id_sts={}
        for id in ids:
            stru_dict=self.collection.find_one({'_id':id})['output']['crystal']
            stru=Structure.from_dict(stru_dict)
            id_sts[id]=stru
        done=[]
        left= [id for id in id_sts if id not in done]
        group=[]
        for id in left:
            st = id_sts[id]
            unique = True
            for i in range(len(group)):
                sg = id_sts[group[i][0]]
                if are_structures_equal(st, sg):
                    group[i].append(id)
                    unique = False
                    break
            if unique == True:
                group.append([id])
            done.append(id)
            left = [id for id in id_sts if id not in done]
        return [i[0] for i in group]

    def _group_blk_diel(self, ids_blk, ids_diel):
        group=[]
        ids_blk_st_sc={}
        for id in ids_blk:
            record=self.collection.find_one({'_id':id})
            stru_dict=record['output']['crystal']
            stru=Structure.from_dict(stru_dict)
            sc=record['transformations']['supercell']
            if len(sc)==1:
                sc=sc[0]
            ids_blk_st_sc[id]=[stru,sc]

        ids_diel_st={}
        for id in ids_diel:
            record=self.collection.find_one({'_id':id})
            stru_dict=stru_dict=record['output']['crystal']
            stru=stru=Structure.from_dict(stru_dict)
            ids_diel_st[id]=stru
            
        for id_blk in ids_blk_st_sc:
            st_blk = ids_blk_st_sc[id_blk][0]
            sc = ids_blk_st_sc[id_blk][1]
            for id_diel in ids_diel:
                st_diel = ids_diel_st[id_diel].copy()
                st_diel.make_supercell(sc)
                if not check_rotation_two_structs(st_blk, st_diel) or are_lattices_equal(st_blk, st_diel):
                    group.append({'bulk':id_blk,'dielectric':id_diel})
                    break
        return group

    def _dfc_vs_blkdiel_group(self, ids_dfc, group_blkdiel):
        out = []
        ids_left = [i for i in ids_dfc]
        group = [i for i in group_blkdiel]
        for nth in range(len(group)):
            if not ids_left:
                break
            pair = group[nth]
            group[nth]['defect']=[]
            id_blk=pair['bulk']
            record = self.collection.find_one({'_id':id_blk})
            stru_blk_dict = record['output']['crystal']
            stru_blk = Structure.from_dict(stru_blk_dict)
            sc_blk= record['transformations']['supercell']
            if len(sc_blk)==1:
                sc_blk=sc_blk[0]
            for id_dfc in ids_dfc:
                record_dfc = self.collection.find_one({'_id':id_dfc})
                #stru_dfc_dict = record_dfc['structures']['defect_no_relaxation']
                stru_dfc_dict = record_dfc['input']['crystal']
                stru_dfc = Structure.from_dict(stru_dfc_dict)
                sc_dfc = record_dfc['transformations']['supercell']
                if are_blk_dfc_match(stru_blk, stru_dfc):
                    group[nth]['defect'].append(id_dfc)
                    ids_left.remove(id_dfc)
            if group[nth]['defect']:
                out.append(group[nth])
        if ids_left:
            print 'Warnning: number of rest defects: %s' % str(len(ids_left))
        return out 
            
    def get_group_on_mpid(self, mpid):
        """
        Classify of all calculations in databse.
        """
        if mpid not in self.groups['bulk']:
            raise KeyError("No %s bulk in database" % mpid)
        if mpid not in self.groups['dielectric']:
            raise KeyError("No %s dielectric in database" % mpid)
        if mpid not in self.groups['defect']:
            raise KeyError("No %s defects in database" % mpid)
        ids_blk= self.groups['bulk'][mpid]
        ids_blk_g = self._ids_strus_group(ids_blk)
        ids_diel=self.groups['dielectric'][mpid]
        ids_diel_g = self._ids_strus_group(ids_diel)
        blk_diel_g= self._group_blk_diel(ids_blk_g,ids_diel_g)
        ids_dfc=self.groups['defect'][mpid]
        group= self._dfc_vs_blkdiel_group(ids_dfc, blk_diel_g)
        return group

class TasksOperater(TasksIdsClassify):
    def __init__(self):
        """
        The operaters for tasks collection, which is the original collection 
        by pymatgen-db inseration, including the defects, bulk, and dielectric.
        """
        TasksIdsClassify.__init__(self)

    def get_structure(self,mpid):
        try:
            blk_id = self.groups['bulk'][mpid][0]
            rec = self.collection.find_one({'_id':blk_id},['input.crystal'])
            stru = Structure.from_dict(rec['input']['crystal'])
        except:
            stru = m.get_structure_by_material_id(mpid)
        spa = SpacegroupAnalyzer(stru)
        prim_stru = spa.get_primitive_standard_structure()
        return prim_stru

    def del_part(self, mpid, part, *kw):
        """
        Args:
            mpid: material id in Materials Project
            part: 'bulk', 'dielectric' or 'defect'
            kw: case 1) defect_type and charge when part is 'defect'
                        for example, 'Cu1_vac', -1 or 'all' 
                case 2) all for delecting all defect for mpid
        """
        if part=='dielectric':
            ids = self.groups_detail[part][mpid]
            for idd in ids:
                self.collection.delete_one({'_id':idd})
        elif part == 'bulk':
            ids = [i[0] for i in self.groups_detail['bulk'][mpid]]
            for idd in ids:
                self.collection.delete_one({'_id':idd})
        if part == 'defect':
            dfcs = self.groups_detail['defect'][mpid]
            id_num = []
            for i in dfcs:
                id_num+=dfcs[i]
            ids = [i[0] for i in id_num]
            for idd in ids:
                if kw[0] == 'all':
                    self.collection.delete_one({'_id':idd})
                else:
                    record=self.collection.find_one({'_id':idd})
                    defect_type = record['transformations']['defect_type']
                    charge = record['transformations']['charge']
                    if kw[0] == defect_type and kw[1]==charge:
                        self.collection.find_one({'_id':idd})

    def del_records_no_mpid(self):
        for i in self.collection.find({},['transformations.history.source','state']):
            try:
                mpid=i['transformations']['history'][0]['source']
		if not mpid.replace(' ', ''):
		    raise KeyError()
            except:
                self.collection.delete_one({'_id':i['_id']})

    def del_rawdata(self,mpid):
        try:
            self.del_part(mpid,'bulk')
        except:
            pass
        try:
            self.del_part(mpid,'dielectric')
        except:
            pass
        try:
            self.del_part(mpid,'defect','all')
        except:
            pass

    def del_all_unsuccessful(self):
        """
        delect all unsuccessful jobs
        """
        for i in self.collection.find({},['state']):
            if i['state'] !='successful':
                self.collection.delete_one({'_id':i['_id']})

    def del_unsuccessful(self,mpid):
        """
        delect all unsuccessful jobs about mpid
        """
        for i in self.collection.find({},['transformations.history.source','state']):
            if i['transformations']['history'][0]['source']==mpid:
                if i['state'] !='successful':
                    self.collection.delete_one({'_id':i['_id']})

    def calculate_defect_correction(self, mpid, db_update=False, path_save_potalign=None):
        """
        mpid: material id in Materials Project
        plot: whether plot potential aligenment pictures
        kws: path = 'path' to save the plots
        """
        groups=self.get_group_on_mpid(mpid)
        if not len(groups):
            raise IOError('Can not match among dielectric bulk and defects')
        elif len(groups)==1:
            print 'There is 1 group for %s' % mpid
        else:
            print 'There are %s groups for %s' % (str(len(groups)),mpid)
        for g in groups:
            if not ( g['bulk'] and g['dielectric'] and g['defect']):
                print 'Data is not completed for this group'
                continue
            id_blk=g['bulk']
            id_diel=g['dielectric']
            ids_dfc=g['defect']
            diel_record = self.collection.find_one({'_id':g['dielectric']})
            epsilon_static = diel_record['calculations'][0]['output']['epsilon_static']
            epsilon_ionic = diel_record['calculations'][0]['output']['epsilon_ionic']
            epsilon =(np.array(epsilon_static)+np.array(epsilon_ionic)).tolist()
            blk_record = self.collection.find_one({'_id':g['bulk']})
            dim_blk = blk_record['output']['site_potentials']['from_outdar']['ngxf_dims']
            stru_blk_dict = blk_record['output']['crystal']
            stru_blk = Structure.from_dict(stru_blk_dict)
            puredat = blk_record['output']['site_potentials']['from_outdar']
            dim_blk=puredat['ngxf_dims']
            bulkinit = KumagaiBulkInit(stru_blk, dim_blk, epsilon)
            gamma = bulkinit.gamma
            g_sum = bulkinit.g_sum
            print 'Defect_distance: %f' % closest_defect_distance(stru_blk)
            for id_dfc in g['defect']:
                dfc_record=self.collection.find_one({'_id':id_dfc})
                charge = dfc_record['transformations']['charge']
                dfc_type = dfc_record['transformations']['defect_type']
                #stru_dfc_dict = dfc_record['output']['crystal'] ## structure with relaxation
                stru_dfc_dict = dfc_record['calculations'][0]['input']['crystal']
                #stru_dfc_dict = dfc_record['calculations'][0]['output']['crystal']
                defdat = dfc_record['output']['site_potentials']['from_outdar']
                dim_dfc = defdat['ngxf_dims']
                stru_dfc = Structure.from_dict(stru_dfc_dict)
                kumagaicorr = KumagaiCorrection(dim_blk, epsilon, charge, gamma, g_sum, stru_blk, stru_dfc, puredat, defdat)
                if path_save_potalign:
                    title='%s_%s%s' % (mpid,dfc_type,str(charge))
                    new_path = join(path_save_potalign,mpid)
                    if not isdir(new_path):
                        os.makedirs(new_path)
                    with cd(new_path):
                        corr = kumagaicorr.correction(title, partflag='AllSplit')
                else:
                    corr = kumagaicorr.correction(partflag='AllSplit')
                pc = corr[0]; potalign=corr[1]; tot=corr[2]; std_devia=corr[3]
                if std_devia>0.06:
                    text='>0.06'
                else:
                    text=''
                if db_update:
                    self.collection.update({'_id':id_dfc}, 
                            {'$set':{'kumagai_correction':{'pc':pc,'potalign':potalign, \
                                    'std_devia':std_devia, 'total':tot}}})
                print '%-10s %-16s pc: %10s potalign: %11s  std_devia: %10s %6s total: %11s' % \
                        (mpid, dfc_type+str(charge), str(pc)+' eV', str(potalign)+' eV', 
                                str(std_devia)+' eV',text, str(tot)+' eV')

    def _gather_computed_defects(self, mpid, corr=True):
        collection=[]
        grps = self.get_group_on_mpid(mpid)
        for gp in grps:
            colli={}
            ids_dfc = gp['defect']
            id_blk = gp['bulk']
            record_blk = self.collection.find_one({'_id':id_blk},['output'])
            stru_blk_dict = record_blk['output']['crystal']
            stru_blk = Structure.from_dict(stru_blk_dict)
            energy_blk = record_blk['output']['final_energy']
            blk_entry = ComputedStructureEntry(stru_blk, energy_blk)
            blk_entry_as_dict = blk_entry.as_dict()
            colli['bulk_entry']=blk_entry_as_dict
            colli['computed_defects']=[]
            for id_dfc in ids_dfc:
                record_dfc = self.collection.find_one({'_id':id_dfc})
                stru_dfc_dict = record_dfc['output']['crystal']
                stru_dfc = Structure.from_dict(stru_dfc_dict)
                dfc_site_dict = record_dfc['transformations']['defect_site']
                dfc_site = PeriodicSite.from_dict(dfc_site_dict)
                energy_dfc = record_dfc['output']['final_energy']
                sc = record_dfc['transformations']['supercell']
                dfc_type = record_dfc['transformations']['defect_type']
                if corr:
                    correction = record_dfc['kumagai_correction']['total']
                else:
                    correction = 0.0
                charge = record_dfc['transformations']['charge']
                data={'potalign':{'std_devia':record_dfc['kumagai_correction']['std_devia']}}
                dfc_entry = ComputedStructureEntry(stru_dfc, energy_dfc,data=data)
                defect = ComputedDefect(dfc_entry, dfc_site, supercell_size=sc, 
                        charge=charge, charge_correction=correction, name=dfc_type)
                defect_dict = defect.as_dict()
                colli['computed_defects'].append(defect_dict)
            collection.append(colli)
        return collection

class DefectOperater(object):
    def __init__(self):
        """
        The operaters on defect collection, which stores the the bulk entries,
        computed defects and defect results: pinning and formations.
        """
        self.collection_name = defect_collection_name
        self.collection = Collection(db, self.collection_name)
        #self._gather_defect_names()
        self._gather_all_mpids()
        self.to = TasksOperater() 
        self.hgo = HseGapOperater() 
    def _gather_all_mpids(self):
        mpids=[]
        for i in self.collection.find({},['mpid']):
            mpids.append(i['mpid'])
        self.all_mpids=mpids

    def has_mpid(self,mpid):
        if mpid in self.all_mpids:
            return True
        else:
            self._gather_all_mpids()
            if mpid in self.all_mpids:
                return True
            else:
                return False

    def get_defects(self,mpid):
        existed = self.has_mpid(mpid)
        if not existed:
            raise IOError("%s doesn't exist in database" % mpid)
        for i in self.collection.find({'mpid':mpid},['defects']):
            dfc={}
            for gp in i['defects']:
                for j in gp['computed_defects']:
                    if j['name'] not in dfc:
                        dfc[j['name']]=[]
                    dfc[j['name']].append(j['charge'])
                    #defects[i['mpid']].append(j['full_name'])
            return dfc

    def _gather_defect_names(self):
        defects={}
        for i in self.collection.find({},['mpid','defects']):
            dfc={}
            for gp in i['defects']:
                for j in gp['computed_defects']:
                    if j['name'] not in dfc:
                        dfc[j['name']]=[]
                    dfc[j['name']].append(j['charge'])
                    #defects[i['mpid']].append(j['full_name'])
            defects[i['mpid']]=dfc
        self.all_defects=defects

    def _get_vbm_gap(self, mpid, label):
        if label not in ['dft','hse']:
            raise KeyError('label must be hse or dft')
        if label=='hse':
            hsegapinfo = self.hgo.get_hsegap(mpid)
            vbm = hsegapinfo['vbm']
            gap = hsegapinfo['gap']
        elif label=='dft':
            try:
                bs = self.collection.find_one({'mpid':mpid},['results.dft'])
                vbm = bs['results']['dft']['vbm']
                gap = bs['results']['dft']['gap']
            except:
		try:
                    bs = m.get_bandstructure_by_material_id(mpid)
                    vbm = bs.get_vbm()['energy']
                    gap = bs.get_band_gap()['energy']
		except:
		    id_ = self.to.groups['band'][mpid][0]
		    recd = self.to.collection.find_one({'_id':id_},['calculations.output.vbm','calculations.output.cbm'])
		    vbm = recd['calculations'][0]['output']['vbm']
		    cbm = recd['calculations'][0]['output']['cbm']
		    gap = cbm - vbm
        return ( vbm, gap )
    
    def compute_formation_energies(self, mpid, label, corr=True, ext_elts=[],**kws):
        if label not in ['hse','dft']:
            raise KeyError("label must be 'hse' or 'dft'!!")
        try:
            defects = kws['defects']
        except:
            if corr:
                defects=self.to._gather_computed_defects(mpid)
            else:
                defects=self.to._gather_computed_defects(mpid, corr=False)
        vbm,gap =self._get_vbm_gap(mpid, label)
        #### results will be inserted to database if db_update==True ####
        results = {}
        results_tmp = []
        results['vbm']=vbm
        results['gap']=gap
        mu_range = get_mu_range(mpid,ext_elts)
        elts = [i.symbol for i in self.to.get_structure(mpid).composition.elements]
        elts = elts + ext_elts
        for gp in defects:
            regionsi={}
            bulk_entry_dict = gp['bulk_entry']
            bulk_entry = ComputedStructureEntry.from_dict(bulk_entry_dict)
            for region, mu in mu_range.items():
                da = DefectsAnalyzer(bulk_entry, vbm, mu, gap)
                for dfc_dict in gp['computed_defects']:
                    defect = ComputedDefect.from_dict(dfc_dict)
                    dfc_elts = [i.symbol for i in defect.entry.composition.elements]
                    include = True
                    for el in dfc_elts:
                        if el not in elts:
                            include = False
                            break
                    if include:
                        da.add_computed_defect(defect)
                regionsi[region]={}
                form_ens = da.get_formation_energies()
                regionsi[region]['formation_energies'] = form_ens
                forms_q0 = [i['energy'] for i in form_ens if not i['charge']]
                forms_q0.append(1000)
                if min(forms_q0)<=0:
                    regionsi[region]['stable']=False
                else:
                    regionsi[region]['stable']=True
                regionsi[region]['chem_pot'] = da._mu_elts
            results_tmp.append(regionsi)
        def combine():
            out={}
            regions=results_tmp[0].keys()
            for region in regions:
                out[region]={}
                out[region]['chem_pot'] = results_tmp[0][region]['chem_pot']
                stable = set([i[region]['stable'] for i in results_tmp])
                if stable == set([True]):
                    out[region]['stable']=True
                else:
                    out[region]['stable']=False
                energies=[]
                for i in results_tmp:
                    energies+=i[region]['formation_energies']
                out[region]['formation_energies']= energies
            return out
        results['regions'] = combine()
        return results

    def _update_defects(self, mpid, hse=True, dft=True):
        defects=self.to._gather_computed_defects(mpid, corr=True)
        kws={'defects':defects}
        record={}
        record['mpid']=mpid
        record['defects']=defects
        record['results']={}
        if hse:
            results_hse = self.compute_formation_energies(mpid, 'hse', corr=True,**kws) 
            record['results']['hse'] = results_hse
        if dft:
            results_dft = self.compute_formation_energies(mpid, 'dft', corr=True,**kws)
            record['results']['dft'] = results_dft
        if mpid in self.all_mpids:
            self.collection.delete_one({'mpid':mpid})
        self.collection.insert(record)

    def _get_potalign_std_devia(self,mpid):
        record=self.collection.find_one({'mpid':mpid},['defects.computed_defects.entry.data',
            'defects.computed_defects.full_name'])
        dfc_devia={}
        for gp in record['defects']:
            for i in gp['computed_defects']:
                name = i['full_name']
                std_devia = i['entry']['data']['potalign']['std_devia']
                dfc_devia[name]=std_devia
        return dfc_devia
    
    def get_lost_defects(self, mpid):
        """
        get the defects without calculation 
        """
        defects_lost = {}
        stru = self.to.get_structure(mpid)
        qlist = get_defect_charge_list(stru)
        defects = self.get_defects(mpid)
        def del_num(s):
            out=''
            for i in s:
                if i.isdigit():
                    continue
                out+=i
            return out
        qlist_in = defects
        #for dfc in defects:
        #    name = '_'.join(dfc.split('_')[:-1])
        #    charge = int(dfc.split('_')[-1])
        #    if name not in qlist_in:
        #        qlist_in[name]=[]
        #    qlist_in[name].append(charge)
        names = []
        for dfc in qlist_in:
            name=del_num(dfc)
            if name not in qlist:
                continue
            name_list = [i for i in name.split('_') if i]
            if 'vac' in name_list:
                name_list.remove('vac')
                name = name_list[0]+'_vac'
            if 'as' in name_list:
                name_list.remove('as')
                name_list.remove('on')
                name = name_list[1]+'_subst_'+name_list[0]
            names.append(name)
            for q in qlist[name]:
                if q not in qlist_in[dfc]:
                    if dfc not in defects_lost:
                        defects_lost[dfc]=[]
                    defects_lost[dfc].append(q)
        for name in qlist:
            if name not in names and 'inter' not in name:
                defects_lost[name]=qlist[name]
        return defects_lost

    def get_potalign_unconverged_defects(self,mpid,potalign_tol):
        out={}
        std_devias = self._get_potalign_std_devia(mpid)
        for dfc in std_devias:
            if std_devias[dfc]>potalign_tol:
                name = '_'.join(dfc.split('_')[:-1])
                charge = int(dfc.split('_')[-1])
                if name not in out:
                    out[name]=[]
                out[name].append(charge)
        return out
   
    def shrink_defect_charge_range(self,mpid,vbmgap_from,corr=True,potalign_tol=0.06):
        """
        Get the smallest safe charge range based on the transition levels, 
        the defect with charge outside this range can be discarded safely.
        """
        transition_levels = self.get_transition_levels(mpid,vbmgap_from,corr=True,potalign_tol=0.06)
        gap=self._get_vbm_gap(mpid,vbmgap_from)[1]
        output={}
        for dfc_type in transition_levels:
            output[dfc_type]=[]
            qminus = -100
            qplus = 100
            for pair in transition_levels[dfc_type]:
                tran_lev = transition_levels[dfc_type][pair]
                qs =pair.split('|')
                q1 = int(qs[0]);q2 = int(qs[1])
                if q1<0 or q2<0:
                    qmin = min([q1,q2])
                    if tran_lev > gap:
                        if qmin>qminus:
                            qminus = qmin
                if q1>0 or q2>0:
                    qmax = max([q1,q2])
                    if tran_lev < 0:
                        if qmax<qplus:
                            qplus = qmax
            output[dfc_type]=[qminus,qplus]
        return output

    def report_dopability(self, mpid, label, corr=True, potalign_tol=0.06, in_detail=False):
        """
        label: 'hse' or 'dft'
        """
        if label not in ['dft','hse']:
            raise KeyError("label must be 'hse' or 'dft'!!")
        if corr:
            results = self.collection.find_one({'mpid':mpid})['results'][label]
        else:
            results = self.compute_formation_energies(mpid, label, corr=False)
        dfc_devia=self._get_potalign_std_devia(mpid)
        vbm=results['vbm'];gap=results['gap']
        report_in_detail={}
        for region in results['regions']:
            report_in_detail[region]={}
            formation_energies = results['regions'][region]['formation_energies']
            formation_energies_screen=[]
            for i in formation_energies:
                name=i['name'];charge=str(i['charge'])
                full_name=name+'_'+charge
                if dfc_devia[full_name]<=potalign_tol:
                    formation_energies_screen.append(i)
            pinning_energies = get_pinning_energies(formation_energies_screen,vbm)
            diff_pinp_vbm = pinning_energies['hole'] -vbm
            diff_pinn_cbm = pinning_energies['electron'] - vbm - gap
            stable = results['regions'][region]['stable']
            if diff_pinp_vbm < 0.0:
                report_in_detail[region]['dopability_p_type'] = 'Good'
            else:
                report_in_detail[region]['dopability_p_type'] = 'Bad'
            if diff_pinn_cbm > 0.0:
                report_in_detail[region]['dopability_n_type'] = 'Good'
            else:
                report_in_detail[region]['dopability_n_type'] = 'Bad'
            report_in_detail[region]['stable']='Good' if stable else 'Bad'
        report_in_short={}
        dopability_p_all=[report_in_detail[i]['dopability_p_type'] for i in report_in_detail]
        dopability_n_all=[report_in_detail[i]['dopability_n_type'] for i in report_in_detail]
        stable_all=[report_in_detail[i]['stable'] for i in report_in_detail]
        if set(dopability_p_all)==set(['Good']):
            report_in_short['p_type_dopability'] = 'Robust'
        elif set(dopability_p_all)==set(['Bad']):
            report_in_short['p_type_dopability'] = 'Bad'
        else:
            report_in_short['p_type_dopability'] = 'Good'
        if set(dopability_n_all)==set(['Good']):
            report_in_short['n_type_dopability'] = 'Robust'
        elif set(dopability_n_all)==set(['Bad']):
            report_in_short['n_type_dopability'] = 'Bad'
        else:
            report_in_short['n_type_dopability'] = 'Good'
        if set(stable_all)==set(['Good']):
            report_in_short['stable'] = 'Robust'
        elif set(stable_all)==set(['Bad']):
            report_in_short['stable'] = 'Bad'
        else:
            report_in_short['stable'] = 'Good'
        if in_detail:
            return report_in_detail
        else:
            return report_in_short
            
    def _pick_defect_types_with_lowest_energy(self, ens_input):
        def del_num(name):
            new=''
            for i in name:
                if not i.isdigit():
                    new+=i
            return new
        def get_min_lines(keys):
            keys_min=[i for i in keys]
            for i in keys:
        	for j in keys:
        	    if j==i:
        		continue
        	    yi=np.array(ens[i]);yj=np.array(ens[j])
                    if set(yi>=yj)==set([True]):
                        keys_min.remove(i)
        		break 
            keys_min_absolute=[i for i in keys_min]

            all=set(keys_min_absolute)
            done=set([])
            left= all-done
            for i in all:
        	done.add(i)
        	for j in all-done:
        	    yi=np.array(ens[i])+0.1;yj=np.array(ens[j])
        	    if set(yi>=yj)==set([True]):
        		keys_min.remove(i)
        		break
        	    else:
        		done=done-set([i])
            return keys_min
        ens = {i:ens_input[i] for i in ens_input}
        keys = ens.keys()
        main_keys=set([del_num(i) for i in keys])
        key_group={}
        for i in main_keys:
            key_group[i]=[]
        for key in keys:
            main_key=del_num(key)
            key_group[main_key].append(key)
        left_keys=[]
        for main_key in key_group:
            left=get_min_lines(key_group[main_key])
            left_keys = left_keys + left
        return left_keys

    def get_transition_levels(self,mpid,vbmgap_from,corr=True,potalign_tol=0.06):
        """
        transition_level with respect to vbm
        """
        dfc_devia=self._get_potalign_std_devia(mpid)
        if corr:
            results = self.collection.find_one({'mpid':mpid})['results'][vbmgap_from]
        else:
            results = self.compute_formation_energies(mpid, vbmgap_from, corr=False)
        vbm=results['vbm'];gap=results['gap']
        region=results['regions'].keys()[0]
        form_ens_one_region=results['regions'][region]['formation_energies']
        form_ens={}
        for i in form_ens_one_region:
            name=i['name']
            charge=i['charge']
            energy=i['energy']
            full_name= name+'_'+str(charge)
            if dfc_devia[full_name]>potalign_tol:
                continue
            if name not in form_ens:
                form_ens[name]={}
            form_ens[name][charge]=energy
        transition_levels={}
        for name in form_ens:
            transition_levels[name]={}
            for pair in itertools.combinations(form_ens[name].keys(),2):
                q1 = pair[0]; q2 = pair[1]
                e1 = form_ens[name][q1]; e2 = form_ens[name][q2]
                pair_name=str(q1)+'|'+str(q2)
                tran_lev = (e1-e2)/(q2-q1)
                transition_levels[name][pair_name]=tran_lev
        return transition_levels



    def plot_formation_energies(self, mpid, label, corr=True, ext_elts=[], potalign_tol=0.06, xlim=None, ylim=(-1,3), path='.', fmt='png'):
        stru=self.to.get_structure(mpid)
        reduced_formula = stru.composition.reduced_formula
        if label not in ['hse','dft']:
            raise KeyError("label must be 'hse' or 'dft'!!")
        new_path = path+'/'+mpid+'_'+reduced_formula+'/'+label
        if not isdir(new_path):
            os.makedirs(new_path)
        ########## comments for plots ######
        comments={}
        oxida_stats = get_oxidation_states(stru)
        comments['oxidation_states'] = {i.symbol:oxida_stats[i] for i in oxida_stats}
        comments['potalign_convergence_criterion'] = potalign_tol
        dfc_devia=self._get_potalign_std_devia(mpid)
        dfc_names = set(['_'.join(i.split('_')[:-1]) for i in dfc_devia])
        std_devias ={i:{} for i in dfc_names}
        defects_in_plots = {}
        defects_not_in_plots = {}
        for i in dfc_devia:
            name = '_'.join(i.split('_')[:-1])
            charge = int(i.split('_')[-1])
            std_devia = dfc_devia[i]
            std_devias[name][charge]=std_devia
            if corr:
                if std_devia <= potalign_tol:
                    if name not in defects_in_plots:
                        defects_in_plots[name]=[]
                    defects_in_plots[name].append(charge)
                else:
                    if name not in defects_not_in_plots:
                        defects_not_in_plots[name]=[]
                    defects_not_in_plots[name].append(charge)
            else:
                if name not in defects_in_plots:
                    defects_in_plots[name]=[]
                defects_in_plots[name].append(charge)
        comments['potalign_standard_deviation'] = std_devias
        comments['defects_in_plots'] = defects_in_plots
        comments['defects_not_in_plots'] = defects_not_in_plots
        file_name=new_path+'/comments.json'
        #with open(new_path+'/comments.json','w') as comm:
        #    comm.write(json.dumps(jsanitize(comments)))
        dumpfn(comments, file_name, cls=MontyEncoder, indent=2,sort_keys=True)
        ####################################
        if ext_elts:
            results = self.compute_formation_energies(mpid, label, corr, ext_elts)
        else:
            if corr:
                results = self.collection.find_one({'mpid':mpid})['results'][label]
            else:
                results = self.compute_formation_energies(mpid, label, corr=False)
        vbm=results['vbm'];gap=results['gap']
        if not xlim:
            xlim=(-0.2, gap+0.5)
        x_min = vbm + xlim[0]
        x_max = vbm + xlim[1]
        x_step = 0.0001
        nstep = int((x_max - x_min)/x_step)
        x_pts = [x_min+i*x_step for i in xrange(nstep+1)]
        picked_defect_type = False
        for region in results['regions']:
            form_ens_split=results['regions'][region]['formation_energies']
            form_ens={}
            for i in form_ens_split:
                name=i['name']
                charge=i['charge']
                energy=i['energy']
                if name in defects_not_in_plots:
                    if charge in defects_not_in_plots[name]:
                        continue
                ens = [ energy + charge * ( x - vbm ) for x in x_pts]
                if name not in form_ens:
                    form_ens[name]={}
                form_ens[name][charge]=ens
            #### form_ens save plots for every detailed defect e.g. Cu1_vac0, Cu1_vac1 ...
            ###########################################################################
            defect_types = form_ens.keys()
            ens_opt = {}
            qlabels_opt = {}
            inds_qlabel_opt = {}
            for t in defect_types:
                charges=[]; enss=[]
                for charge, ens in form_ens[t].items():
                    charges.append(charge)
                    enss.append(ens)
                ens_opt_t = np.array(enss).min(axis=0)
                ens_opt[t]=ens_opt_t ### formation energies for plot 
                
                q_indexs = np.argmin(enss, axis=0)
                charges_opt = [charges[i] for i in q_indexs]
                def get_inds_for_q_labels(charges):
                    qs_backup=list(set(charges)); qs_backup.sort()
                    qs=[i for i in qs_backup]
                    inds=[]
                    for q in qs_backup:
                        inds_first = np.where(np.array(charges)==q)[0]
                        inds_second = [i for i in inds_first if ylim[0]<=ens_opt_t[i]<=ylim[1]]
                        if inds_second:
                            ind = int(np.array(inds_second).mean())
                            inds.append(ind)
                        else:
                            qs.remove(q)
                    return (inds,qs)
                inds_qs = get_inds_for_q_labels(charges_opt)
                inds_qlabel_t = inds_qs[0]
                qlabels_t = inds_qs[1]
                inds_qlabel_opt[t] = inds_qlabel_t
                qlabels_opt[t] = qlabels_t
            ###### ens_opt save the Eform plot for each defect_type, use this for plot in detail ############
            ###### qlabels_opt save the charges appearing in plot for each defect_type #####
            ###### inds_qlabel_opt save the index of plot of each q label #######
            if not picked_defect_type:
                left_keys = self._pick_defect_types_with_lowest_energy(ens_opt)
                picked_defect_type = True
            for key in defect_types:
                if key in left_keys:
                    continue
                del ens_opt[key], qlabels_opt[key], inds_qlabel_opt[key]

            defect_types_after_pick = ens_opt.keys()
            #############################################################################
            import matplotlib
            from matplotlib import pyplot as plt
            from math import ceil
            #matplotlib.style.use('classic')
            ncolumn = 3
            nline_in_subplot=20
            adict=locals()
            nsubplot=int(ceil(len(ens_opt)/float(nline_in_subplot)))  #### subplot numbers
            nrow=int(ceil(nsubplot/float(ncolumn)))
            fig_width = 8*ncolumn
            fig_height = 6*nrow
            if nsubplot <= ncolumn:
                fig_width = 8*nsubplot
                fig_height = 6
            fig_size = [fig_width, fig_height]
            params = {'backend': 'ps',
                    'axes.labelsize': 25,
                    'font.size': 25,
                    'font.weight': 'normal',
                    'legend.fontsize': 15,
                    'xtick.labelsize': 25,
                    'ytick.labelsize': 25,
                    'text.usetex': False,
                    'figure.figsize': fig_size,
                    'lines.markersize': 5}
            plt.rcParams.update(params)
            if nsubplot > 1:
                fig,axes=plt.subplots(nrow,ncolumn,sharey=True)
                n=0
                for i in range(nrow):
                    for j in range(ncolumn):
                        if nrow>1:
                            adict['ax%s' % n]=axes[i][j]
                        elif nrow==1:
                            adict['ax%s' % n]=axes[j]
                        adict['legend%s' % n]=[]
                        n=n+1
                if nrow*ncolumn>nsubplot:
                    for i in range(nsubplot,nrow*ncolumn):
                        fig.delaxes(adict['ax%s' % str(i)])
            elif nsubplot==1:
                fig,adict['ax0']=plt.subplots(1,1)
                adict['legend0']=[]
            import matplotlib.cm as cm

            def get_legends(types):
                legends = []
                for name in types:
                    legends.append(name)
                return legends

	    #defect_types.sort(key=lambda x:(x.split('_')[-1][:-1],x.split('_')[0]))
	    defect_types_after_pick.sort()
	    colors=['blue','green','red','cyan','pink','gray']
	    linetypes=['-','--','-.',':']
	    i=0
            for c,cnt in zip(defect_types_after_pick,range(len(defect_types_after_pick))):
	        j=i%len(colors)
	        k=i/len(colors)
                ind=int(ceil((cnt+1)/float(nline_in_subplot)))-1
                adict['legend%s' % str(ind)].append(c)
                adict['ax%s' % str(ind)].plot(np.array(x_pts)-vbm, ens_opt[c], linewidth=2,
                          color=colors[j],linestyle=linetypes[k])#, color=colors[cnt])
 	        for ll in range(len(qlabels_opt[c])):
	            adict['ax%s' % str(ind)].text(x_pts[inds_qlabel_opt[c][ll]]-vbm,ens_opt[c][inds_qlabel_opt[c][ll]],
                                str(qlabels_opt[c][ll]),color=colors[j],horizontalalignment='center')
	        if i==nline_in_subplot-1:
	    	    i=0
	    	    continue
	        i=i+1
                #plt.plot(x, y[c], linewidth=3, color=colors[cnt])
            for i in range(nsubplot):
                adict['ax%s' % i].plot([x_min-vbm, x_max-vbm], [0, 0], 'k-')
                adict['ax%s' % i].axvline(x=0.0, linestyle='--', color='k', linewidth=2)
                adict['ax%s' % i].axvline(x=gap, linestyle='--', color='k', linewidth=2)
                adict['ax%s' % i].set_ylim(ylim)
                adict['ax%s' % i].set_xlim(xlim)
                adict['ax%s' % i].legend(get_legends(adict['legend%s' % i]),loc='upper right')
                #adict['ax%s' % i].set_xlabel("Fermi energy (eV)"s
                #adict['ax%s' % i].set_ylabel("Formation Energy (eV)")
            adict['ax%s' % 0].set_ylabel("Formation energy (eV)")
            adict['ax%s' % str(nsubplot-1)].set_xlabel('Fermi energy (eV)')
            fig.set_tight_layout({'pad':0.2, 'h_pad':0.1, 'w_pad':0.0})
            fig.savefig(join(new_path,'%s.%s' % (region,fmt)))
            plt.clf()

    def _plot_formation_energies_old(self, mpid, label, xlim=None, ylim=(-1,3), path='', fmt='jpg',  **kws):
        try:
            vbm = kws['vbm']
            gap = kws['gap']
            if not (vbm and gap):
                raise KeyError()
        except:
            if label=='hse':
                hsegapinfo = self.hgo.hse_gaps[mpid]
                vbm = hsegapinfo['vbm']
                gap = hsegapinfo['gap']
            elif label=='dft':
                bs = m.get_bandstructure_by_material_id(mpid)
                vbm = bs.get_vbm()['energy']
                gap = bs.get_band_gap()['energy']
            else:
                raise KeyError('label tag only can be hse or dft')
        record = self.collection.find_one({'mpid':mpid})
        bulk_entry_dict = record['bulk_entries'][0]
        bulk_entry = ComputedStructureEntry.from_dict(bulk_entry_dict)
        mu_range = get_mu_range(mpid)
        new_path = path+'/'+mpid 
        if not isdir(new_path):
            os.makedirs(new_path)
        for region, mu in mu_range.items():
            da = DefectsAnalyzer(bulk_entry, vbm, mu, gap)
            for dfc_dict in record['computed_defects']:
                defect = ComputedDefect.from_dict(dfc_dict) 
                da.add_computed_defect(defect)
            plotter = DefectPlotter(da)
            form_en_plot = plotter.get_plot_form_energy(xlim=xlim, ylim=ylim)
            form_en_plot.savefig(join(new_path, region+'.'+fmt),dpi=200)

class HseGapOperater(object):
    def __init__(self):
        self.collection_name = hsegap_collection_name
        self.collection = Collection(db,self.collection_name)
        #self._gather_hse_gaps()
        self._gather_all_mpids()
    
    def _gather_all_mpids(self):
        mpids=[]
        for i in self.collection.find({},['mpid']):
            mpids.append(i['mpid'])
        self.all_mpids=mpids

    def has_mpid(self,mpid):
        if mpid in self.all_mpids:
            return True
        else:
            self._gather_all_mpids()
            if mpid in self.all_mpids:
                return True
            else:
                return False

    def get_hsegap(self,mpid):
        rec= self.collection.find_one({'mpid':mpid},['gap_info'])
        if 'gap_info' in rec:
            if 'vbm' in rec['gap_info']:
                return{'vbm':rec['gap_info']['vbm']['energy'],'gap':rec['gap_info']['gap']}
        
    
    def _gather_hse_gaps(self):
        all_gaps={}
        for i in self.collection.find({},['material.mateiral_id','gap_info']):
            if 'gap_info' in i:
                if 'vbm' in i['gap_info']:
                    all_gaps[i['material']['mateiral_id']]={'vbm':i['gap_info']['vbm']['energy'],'gap':i['gap_info']['gap']}
        self.hse_gaps=all_gaps

    def insert_new_record(self, calc_dir):
        transf=json.load(open(os.path.join(calc_dir,'transformations.json')))
        mpid=transf['history'][0]['source']
        existed = self.has_mpid(mpid)
        if existed:
            print '%s is already in database' % mpid
            return 
        try:
            xml=Vasprun(os.path.join(calc_dir,'vasprun.xml'))
        except:
            print '%s: calculation is not finished!' % mpid
            return 
        bs=xml.get_band_structure()
        formula=xml.final_structure.composition.reduced_formula
        material={'mateiral_id':mpid,'formula':formula,'is_spin_polarized':bs.is_spin_polarized}
         
        record=xml.as_dict()
        record['mpid']=mpid
        record['band_structure'] = bs.as_dict()
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
        
        record['material']=material
        record['gap_info']=gap_info
        record = jsanitize(record) 
        self.collection.insert(record)
        
